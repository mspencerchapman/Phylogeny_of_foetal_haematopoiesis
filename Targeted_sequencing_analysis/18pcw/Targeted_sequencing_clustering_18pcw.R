#Cluster the posterior probability vectors from the different tissues
library(data.table)
library(viridis)
library(dplyr)
library(tidyr)

my_working_directory="~/R Work/Fetal HSPCs/Phylogeny_of_foetal_haematopoiesis/"
setwd(my_working_directory)

#Set the colour scheme for the heat maps
viridis <- viridis::viridis(11)
colour.scale <- c("lightgrey", brewer.pal(9, name = "YlOrRd"))

#Cluster the cell fraction results for each aggregate tissue for each mutation
samples=tissues_comb

#DEFINE SIMILARITY USING SOFT COSINE SIMILARITY SCORE
###Choose the kind of data to input for the analysis
#NB. clearest outcomes displaying relationship of PB to mesenchyme/ placental CD34+s is by using
# post.prob with tissues_comb/tissues and censoring values <0.5 (i.e. more likely absent than present)

#(1) Generate the vaf_present/clean.post.prob matrix - then keep only those with the samples desired
cell_frac=calculate_cell_frac(NV=NV,NR=NR)
cell_frac_present<-cell_frac[details_targ$mut_ref,gsub("_comb","",samples)]

##Include mutations present with highest confidence. Set any with a prob < threshold of 0.9 (in the combined prob vectors) to 0.
prob_threshold_to_include<-0.9
cell_frac_present[clean.post.prob[details_targ$mut_ref,samples]<prob_threshold_to_include]<-0
#Log transform and scale the data to be between 0.01 and 1
log_cell_frac_present=cell_frac_present
log_cell_frac_present[cell_frac_present != 0] <- log(log_cell_frac_present[cell_frac_present != 0])
log_cell_frac_present_scaled=log_cell_frac_present
scale_range=c(0.05,1) #The scale range impacts on results, as it sets the difference between "present at low level" and "absent"
log_cell_frac_present_scaled[cell_frac_present!=0] = plotrix::rescale(log_cell_frac_present[cell_frac_present!=0],newrange = scale_range) #scale these figures between 0.01 and 1

#Transform the data for the rest of the analysis
data=t(log_cell_frac_present_scaled[rowSums(clean.post.prob[,germ_layers_comb]>prob_threshold_to_include)>0,])
sum(rowSums(clean.post.prob[,germ_layers_comb]>prob_threshold_to_include)>0) #Print the number of muts included in analysis

#(2) Get the names of the mutations that are not all 0 (i.e. that "differentiate" the tissues)
diff_muts=colnames(data) #get a vector of mut_refs for these differing mutations

#(3) Calculate the phylogenetic distances for these mutations
mut_pairs=as.data.frame(data.table::CJ(a=diff_muts,b=diff_muts,unique=FALSE,sorted=FALSE)) #create a df of all possible pairs
mut_pairs$uid=apply(mut_pairs[,1:2],1,paste,collapse="-") #create a unique id for each pair

dim(mut_pairs) #How many are included - gives an idea of how long it will take

new_tree=TRUE
if(new_tree) {
  #Generate data frame of all possible node pairs
  nodes=unique(tree.multi$edge[,2])
  node_combs=as.data.frame(data.table::CJ(a=nodes,b=nodes,unique=TRUE))
  #Calculate genetic distances for all of these node pairs - save as a list
  nodeheights=nodeHeights(tree.multi) #The lapply function needs the nodeheights object
  node_comb_dist=lapply(1:nrow(node_combs),function(i) {
    node_a=node_combs[i,1]
    node_b=node_combs[i,2]
    node_height_a=nodeheights[tree.multi$edge[,2]==node_a,2]
    node_height_b=nodeheights[tree.multi$edge[,2]==node_b,2]
    if(node_a==node_b) {stop(return(0))}
    ancestral_nodes_a=get_ancestral_nodes(node_a,tree.multi$edge) #get all the ancestral nodes of each node & include the node itself
    ancestral_nodes_b=get_ancestral_nodes(node_b,tree.multi$edge)
    common_ancestors=dplyr::intersect(c(node_a,ancestral_nodes_a),c(node_b,ancestral_nodes_b)) #pull out those that are common ancestors
    if(length(common_ancestors)==0) {
      genetic_dist <-(node_height_a+node_height_b) #If there are no common ancestors listed, then the closest common ancestor is the tree root & the genetic distance is the sum of the node heights
    } else {
      common_ancestors_heights=sapply(common_ancestors, function(x) {
        nodeheights[tree.multi$edge[,2]==x,2] #otherwise find the heights of all common ancestors
      })
      genetic_dist = node_height_a+node_height_b-2*max(common_ancestors_heights) #the common ancestor with the maximum nodeheight is the most recent
    }
    return(genetic_dist)
  })
  #Name the list with the node pair in the format "node_1-node_2"
  node_comb_names=apply(node_combs,1,paste,collapse="-")
  names(node_comb_dist) <- node_comb_names
}

#Now using this list as a reference, get the distance for all individual mutation pairs
shortest_dist = function(mut_a,mut_b, tree, details) { #the function to look up the distance
  node_a=details$node[details$mut_ref==mut_a]
  node_b=details$node[details$mut_ref==mut_b]
  genetic_dist<-node_comb_dist[[paste(node_a,node_b,sep="-")]]
  return(genetic_dist)
}

mut_pairs$dist<-NA
mut_pairs$dist[mut_pairs$uid %in% rownames(ref_dist)]=ref_dist[mut_pairs$uid[mut_pairs$uid %in% rownames(ref_dist)],"dist"]
empties=which(is.na(mut_pairs$dist))
for(i in empties) {
  mut_pairs$dist[i] <-shortest_dist(mut_a=mut_pairs[i,1],mut_b=mut_pairs[i,2],tree=tree.multi,details=details_targ)
  if(i%%1000==0) {print(i)}
}

save_dist=TRUE
if(save_dist) {
  ref_dist=mut_pairs[,c("uid","dist")]
  rownames(ref_dist)<-ref_dist$uid
}

##(4) Convert the "distance" into a "similarity", by taking the inverse.
mut_pairs$sim=1/(mut_pairs$dist+1) #Add one to denominator to avoid dividing by 0, and make max of 1.

##(5) Define functions to calculate the cosine similarity for tissue pairs
sum_tissue_scores = function(tissue_1,tissue_2,mut_pairs_sim) {
  tissue_pairs_df=as.data.frame(data.table::CJ(a=data[tissue_1,diff_muts],b=data[tissue_2,diff_muts],unique=FALSE,sorted = FALSE))
  tissue_pairs_df$sim=mut_pairs_sim
  tissue_pairs_df$numerator=apply(tissue_pairs_df[,1:3],1,prod)
  return(sum(tissue_pairs_df$numerator)/2)
}

get_soft_cosim=function(tissue_1,tissue_2,mut_pairs_sim) {
  num=sum_tissue_scores(tissue_1,tissue_2,mut_pairs_sim = mut_pairs_sim)
  denom=sum_tissue_scores(tissue_1,tissue_1,mut_pairs_sim = mut_pairs_sim)^0.5 * sum_tissue_scores(tissue_2,tissue_2,mut_pairs_sim = mut_pairs_sim)^0.5
  soft_cosim=num/denom
  return(soft_cosim)
}

##Apply this across all possible tissue pairs
tissue_pairs=as.data.frame(data.table::CJ(tissues,tissues))
colnames(tissue_pairs) <- c("tissue_1","tissue_2")
tissue_pairs$soft_cosim = NA
for(i in 1:nrow(tissue_pairs)) {
  if(i%%10==0) {print(i)}
  tissue_pairs$soft_cosim[i] <-get_soft_cosim(tissue_1 = tissue_pairs[i,1],tissue_2=tissue_pairs[i,2],mut_pairs_sim = mut_pairs$sim)
}

#(6) Convert to wide format for the heatmap
sim_mat<-tissue_pairs%>%
  pivot_wider(names_from = tissue_2,values_from = soft_cosim) %>%
  dplyr::select(-1) %>%
  as.matrix()
rownames(sim_mat) <- colnames(sim_mat)

pdf("Figures/18pcw/Tissue_clustering_heatmap_18pcw.pdf")
heatmap(sim_mat,
        scale = "none",
        col=viridis,
        distfun = function(x) dist(x,method="euclidean"),
        hclustfun = function(x) hclust(x, method="complete"),
        labRow = gsub("_"," ",gsub("_comb","",tissues_comb)),
        labCol = gsub("_"," ",gsub("_comb","",tissues_comb)),
        margins=c(10,10),
        cexRow = 0.8,
        cexCol = 0.8)
dev.off()

#Make a dendrogram of the dist_mat
dist_mat=dist(sim_mat,method="euclidean")
tissues.hclust=hclust(dist_mat,method="complete")
dd=as.dendrogram(tissues.hclust)
plot(dd)

pdf("Figures/18pcw/Tissue_mutation_clustering_heatmap_18pcw.pdf")
#Use this dendrogram to apply to the original
heatmap(data,
        main="Heatmap of mutation cell fractions in different tissues",
        Rowv = dd,
        scale="none",
        labCol=NA,
        margins=c(5,5),
        labRow = gsub("_"," ",gsub("_comb","",tissues_comb)),
        #ylab="Tissue",
        xlab="Mutation",
        cexRow = 0.8,
        cexCol = 0.8,
        col=colour.scale)
dev.off()


##PLOT THE SCALE BARS FOR THESE ANALYSES##
#This is to plot a colour scale bar to visualize how the colour scales on the log_vaf tree corresponds to the vaf

#This bit rescales the scale bar using the same transformations applied to the vaf data
scale_bar=rev(1/exp(seq(0,8))) #These are the tick points that will go on the scale
log_scale_bar=log(scale_bar)

#Get the vaf_present & log_vaf_present objects (performed as within the generate_targ_seq_plots function)
cell_frac=calculate_cell_frac(NV=NV,NR=NR)
cell_frac_present<-cell_frac[details_targ$mut_ref,gsub("_comb","",samples)]

##Include mutations present with highest confidence. Set any with a prob < threshold of 0.9 (in the combined prob vectors) to 0.
cell_frac_present[clean.post.prob[details_targ$mut_ref,samples]<prob_threshold_to_include]<-0 #taken from earlier in script
log_cell_frac_present=cell_frac_present
log_cell_frac_present[cell_frac_present != 0] <- log(log_cell_frac_present[cell_frac_present != 0])
log_cell_frac_present_scaled=log_cell_frac_present

#Use these this to perform the identical transformations on the scale bar
xrange=range(log_cell_frac_present[cell_frac_present != 0])
mfac=(scale_range[2]-scale_range[1])/(xrange[2]-xrange[1]) #scale range is taken from earlier in script
log_scale_bar_scaled = scale_range[1] + (log_scale_bar - xrange[1])*mfac

#Map these
tick_positions=round(log_scale_bar_scaled,digits=2)
tick_values=signif(scale_bar,digits=2)

#Remove those for which the positions are outside the scale bar
tick_positions_trimmed<-tick_positions[tick_positions>0 & tick_positions<=1]
tick_values_trimmed<-tick_values[tick_positions>0 & tick_positions<=1]

#Use the same colour scale, this process mirrors the add_var_col function in its processing
#of the colour scale
colfunc = colorRampPalette(colour.scale)
lut=colfunc(101)
scale = (length(lut)-1)
pdf("Figures/18pcw/Scale_bar_tissuemut_clustering_18pcw")
plot(c(0,10), c(0,1), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main="")
axis(side=2,pos=1,las=1,at=tick_positions_trimmed,labels=as.character(tick_values_trimmed))
for (i in 1:(length(lut)-1)) {
  y = (i-1)/scale
  rect(1,y,2,y+1/scale, col=lut[i], border=NA)
}
dev.off()