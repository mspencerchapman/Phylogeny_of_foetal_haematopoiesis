library(stringr)
library(ape)
library(seqinr)
library(ggtree)
library(tidyr)
library(dplyr)
library(ggplot2)
library(plotrix)
library(phangorn)
library(RColorBrewer)

#Functions required in script
aggregate_cols=function(NV,NR,lcm_summary_table,aggregate_feature,exclude_class=NULL) {
  classes=names(table(lcm_summary_table[,aggregate_feature]))
  if(!is.null(exclude_class)){classes<-classes[!classes%in%exclude_class]}
  NV_list=lapply(classes, function(class) {
    return(apply(NV,1,function(x) {sum(x[colnames(NV) %in% lcm_summary_table$Sample_ID[lcm_summary_table[,aggregate_feature] == class]])}))
  })
  NR_list=lapply(classes, function(class) {
    return(apply(NR,1,function(x) {sum(x[colnames(NR) %in% lcm_summary_table$Sample_ID[lcm_summary_table[,aggregate_feature] == class]])}))
  })
  NV_agg=Reduce(cbind,NV_list)
  NR_agg=Reduce(cbind,NR_list)
  dimnames(NV_agg)=dimnames(NR_agg)=list(rownames(NV),classes)
  return(list(NR=NR_agg,NV=NV_agg))
}

combine_post_prob=function(post.prob,lcm_summary_table,aggregate_feature,exclude_class=NULL) {
  classes=names(table(lcm_summary_table[,aggregate_feature]))
  if(!is.null(exclude_class)){classes<-classes[!classes%in%exclude_class]}
  comb_list=lapply(classes,function(class) {
    sample_cols=c(lcm_summary_table$Sample_ID[lcm_summary_table[,aggregate_feature]==class])
    return(apply(post.prob[,sample_cols,drop=FALSE], MARGIN=1,FUN=function(x) {1-prod(1-x)}))
  })
  post.prob_comb=Reduce(cbind,comb_list)
  rownames(post.prob_comb)<-rownames(post.prob)
  colnames(post.prob_comb)<-paste0(classes,"_comb")
  return(post.prob_comb)
}

my_working_directory="~/R Work/Fetal HSPCs/Phylogeny_of_foetal_haematopoiesis/"
treemut_dir="~/R Work/R_scripts/treemut/"
setwd(my_working_directory)

#Define the file paths for the data files
tree_file_path="Data/18pcw/Tree_18pcw.tree"
file_annot="Data/18pcw/Filtered_mut_set_annotated_18pcw"
cgpvaf_LCM_SNV_file ="Data/18pcw/PD41768.LCM.snp.tsv"
cgpvaf_colonies_SNV_file ="Data/18pcw/PD41768.colonies.snp.tsv"
lcm_smry_path = "Data/18pcw/targeted_seq_ids_18wks.csv"
clean.post.prob_path="Data/18pcw/clean.post.prob_18wks.tsv" #This is generated in this script, but is saved as it takes some time

R_function_files = list.files("R_functions",pattern=".R",full.names=TRUE)
sapply(R_function_files,source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)

#--------------------------------------------------------------------------------------------------
##IMPORTING THE CGPVAF TABLES (SNVs and INDELs), AND AGGREGATING COUNTS ACROSS TISSUES/ GERM LAYERS
#--------------------------------------------------------------------------------------------------
#Import the cgpvaf files and combine the LCM and single-cell colony (SCC) matrices into single matrix
LCM_mats<-import_cgpvaf_SNV_and_INDEL(SNV_output_file=cgpvaf_LCM_SNV_file,INDEL_output_file=NULL)
SCC_mats<-import_cgpvaf_SNV_and_INDEL(SNV_output_file=cgpvaf_colonies_SNV_file,INDEL_output_file=NULL)

#Import the LCM summary table
lcm_smry = read.csv(lcm_smry_path,stringsAsFactors = FALSE)

#Create aggregated read counts for each tissue/ germ-layer/ pre- or post-gastrulation category (as defined in the lcm_smry table)
Tissue_agg=aggregate_cols(LCM_mats$NV,LCM_mats$NR,lcm_summary_table = lcm_smry,aggregate_feature = "Tissue")
Germ_layer_agg=aggregate_cols(LCM_mats$NV,LCM_mats$NR,lcm_summary_table = lcm_smry,aggregate_feature = "Germ_layer")

#Combine the individual LCM sample counts with the aggregated counts
NV<-cbind(LCM_mats$NV,Tissue_agg$NV,Germ_layer_agg$NV)
NR<-cbind(LCM_mats$NR,Tissue_agg$NR,Germ_layer_agg$NR)

#Import the details matrix and tree (needed to interpret the targeted sequencing data)
load(file_annot); tree.18wk <- read.tree(tree_file_path); tree.18wk$tip.label=gsub("_hum","",tree.18wk$tip.label)
details <- filtered_muts$COMB_mats.tree.build$mat

#In details file, create column for whether it is included in the targeted sequencing baits
details$targeted_tree = ifelse(details$mut_ref %in% rownames(NV),"YES","NO")

#Now subset the tree and the details matrix to only include those mutations that are in the bait set
tree_targ = get_subset_tree(tree = tree.18wk, details = details, v.field = "targeted_tree",value = "YES")
details_targ = details[details$targeted_tree=="YES",]

NV=NV[details_targ$mut_ref,]
NR=NR[details_targ$mut_ref,]

#--------------------------------------------------------------------------------------------------
##ESTIMATE THE POSTERIOR PROBABILITY OF EACH MUTATION BEING PRESENT IN A GIVEN TISSUE##
#--------------------------------------------------------------------------------------------------
depth.cols=SCC_mats$NR[details_targ$mut_ref,]; mtr.cols=SCC_mats$NV[details_targ$mut_ref,] #The matrices for calculating background error (single cell colonies)

#Only use single cell colonies that are included in the tree
depth.cols <- depth.cols[,names(depth.cols) %in% tree.18wk$tip.label]
mtr.cols <- mtr.cols[,names(mtr.cols) %in% tree.18wk$tip.label]

# Check matrix structure is identical
print(all(row.names(NR) == row.names(NV)))
print(all(names(NR) == names(NV)))
print(all(row.names(NR) == details_targ$mut_ref))
print(all(row.names(depth.cols) == row.names(NR)))
print(all(row.names(mtr.cols) == row.names(NR)))

# Set parameters
alpha_mut <- rep(1.5, nrow(NR))
beta_mut <- rep(10, nrow(NR))
min_theta <- 0.001
max_theta <- 0.1

# Estimate alpha and beta parameters for the error distribution
# Uses empirically weighted method of moments, as described by Keinman, JASA 1973
# Loop through each mutation -> find colonies not carrying that variant -> estimate alpha / beta from those colonies

alpha_error <- beta_error <- rep(0,nrow(NR))
descendants_per_mut <- Descendants(x = tree.18wk, node = details_targ$node, type = "tips")
for (i in 1:nrow(NR)) {
  desc_tips <- tree.18wk$tip.label[descendants_per_mut[[i]]]
  
  curr.n <- unlist(depth.cols[i, !(names(depth.cols) %in% desc_tips)])
  curr.x <- unlist(mtr.cols[i, !(names(mtr.cols) %in% desc_tips)])
  curr.x <- curr.x[curr.n > 0]
  curr.n <- curr.n[curr.n > 0]
  
  min_error <- max((sum(curr.x) + 1) / (sum(curr.n) + 2), 0.001)
  moment.ests.1 <- moment.est(curr.x, curr.n, curr.n)
  if (moment.ests.1[1] < min_error) {
    p_hat <- min_error
    gamma_hat <- 0
  } else {
    new.wts <- curr.n / (1 + moment.ests.1[2] * (curr.n-1))
    moment.ests.2 <- moment.est(curr.x, curr.n, new.wts)
    p_hat <- moment.ests.2[1]
    gamma_hat <- moment.ests.2[2]
  }
  theta_hat <- min(max(gamma_hat / (1-gamma_hat), min_theta), max_theta)
  
  alpha_error[i] <- p_hat / theta_hat
  beta_error[i] <- (1 - p_hat) / theta_hat
}

# Define prior prob as fraction of colonies in the haematopoietic tree with the mutation
prior_mut <- sapply(1:nrow(NR), function(i) {length(Descendants(tree.18wk, details_targ$node[i], "tips")[[1]]) / length(tree.18wk$tip.label)})

# Calculate posterior probabilities
post.prob <- cbind(sapply(1:ncol(NR), function(i) {bbprob.calculator(x_i = NV[,i], n_i = NR[,i], alpha_error = alpha_error, beta_error = beta_error, alpha_mut = alpha_mut, beta_mut = beta_mut, prior_mut = prior_mut)}))
colnames(post.prob) <- names(NR); row.names(post.prob) <- row.names(NR)
post.prob[post.prob<0.05]<-0

#########################################################################
#In a given tissue, individual cuts may be clonal/ oligoclonal, and therefore in particular lineages may be
#represented at a higher vaf than i the aggregated read counts for the tissue.  Therefore, reassessing aggregated
#counts may dilute the ability to call mutations. Therefore call mutations based on 1-prod(1-x) where x is the individual
#sample probability i.e. 1 - the probability that the mutation is absent in all individual samples.

post.prob<-cbind(post.prob,
                 combine_post_prob(post.prob,lcm_smry,aggregate_feature="Tissue"),
                 combine_post_prob(post.prob,lcm_smry,aggregate_feature="Germ_layer"))

#########################################################################
#CREATE FUNCTION TO ALTER CERTAIN POST.PROBS BASED ON PHYLOGENY i.e. "clean it up" to make sense with the phylogeny
#This function does two main alterations:
#(1) Removes positive posterior probabiities when there are no other mutations called on same branch or ancestral branch
#(2) Boosts the posterior probabilities when there are high posterior probabilities in descendant mutations & other mutations on the same branch & ancestral mutations
re_run=FALSE
if(re_run) {
  clean.post.prob=apply(post.prob,2,clean_up_post,details_targ,tree.18wk)
  write.table(clean.post.prob,file = clean.post.prob_path,quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
} else {
  clean.post.prob=as.matrix(read.delim(clean.post.prob_path,sep="\t"))
  colnames(clean.post.prob)<-gsub("\\."," ",colnames(clean.post.prob)) #Saving replaces the spaces with \. - change these back
}

#########################################################################
## CREATE THE GENERATION TREES##
#Create versions of tree where branch lengths are proportional to the estimated number of cell divisions by this point
tree.gen=tree_targ
for(node in unique(tree.gen$edge[,1])) {
  multi=sum(tree.gen$edge[,1]==node)
  tree.gen$edge.length[tree.gen$edge[,1]==node] <- log2(multi) #set branch lengths according to generations according to the formula log2(number of daughters)
}

#Manually change nodes that are >1 division from targeted sequencing (see the "lineage_loss_boot.R" script for calculation)
tree.gen$edge.length[which(tree.gen$edge[,2]==238)] <-2
tree.gen$edge.length[which(tree.gen$edge[,2]==383)] <-1+tree.gen$edge.length[which(tree.gen$edge[,2]==383)]
tree.gen$edge.length[which(tree.gen$edge[,2]==368)] <-1+tree.gen$edge.length[which(tree.gen$edge[,2]==368)]
tree.gen$edge.length[which(tree.gen$edge[,2]==367)] <-2+tree.gen$edge.length[which(tree.gen$edge[,2]==367)]

#Generate two further derivatives - one made ultra-metric, and then "cut" at a given number of generations
tree.gen.ultra<-force.ultrametric(tree.gen,method = "extend")
tree.gen.ultra.squash=squash_tree(tree.gen.ultra,cut_off = 10.9)

####NOW READY TO PLOT THE TREES IN A VARIETY OF WAYS ####
#For this function to run, need several objects in the environment:
#(1) post.prob & clean.post.prob matrices
#(2) the vaf matrix (this gets transformed to the log_vaf_present matrix & scaled within the function)
#(3) the vectors: tissues, germ_layers, tissues_comb and germ_layers_comb
#(4) note that the post.prob & clean.post.prob matrices are expected to include the "_comb" versions for tissues & germ_layers
#(6) A colour.scale object containing the colour scale for the trees (and the scale bar)

#Create vectors of the two "comb" data types
tissues_comb=paste0(tissues,"_comb")
germ_layers_comb=paste0(germ_layers,"_comb")

#Define the colour scale
colour.scale <- c("lightgrey", brewer.pal(9, name = "YlOrRd"))
colour.scale<-colour.scale[c(3:9,9)]

#Do some plots using the "generate_targ_seq_plots" function
pdf("Figures/18pcw/Ectoderm_samples.18wks.pdf",width=13,height=7)
samples=c("HAIR_FOLLICLE_comb","PERIPHERAL_NERVE_comb","EPIDERMIS_comb","ECTODERM_comb")
generate_targ_seq_plots(samples,
                        tree=tree.gen.ultra.squash,
                        details_targ=details_targ,
                        matrices=list(NV=NV[details_targ$mut_ref,],NR=NR[details_targ$mut_ref,]),
                        post_prob_type="clean", #other option is "raw"
                        info_type="log_vaf", #other option is "post.prob" or "vaf"
                        prob_threshold_to_include=0.5,
                        plot_cell_frac=FALSE,
                        plot_donut=T, #plot a donut for positive nodes
                        donut_info="cell_frac", #options are "cell_frac" or "lineages_lost"
                        CI=0.95,
                        radius=4,
                        scale_muts_to_branch = T)

dev.off()

pdf("Figures/18pcw/Germ_layers_cellfracpie.18wks.pdf",width=13,height=7)
samples=germ_layers_comb
generate_targ_seq_plots(samples,
                        tree=tree.gen.ultra.squash,
                        details_targ=details_targ,
                        matrices=list(NV=NV[details_targ$mut_ref,],NR=NR[details_targ$mut_ref,]),
                        post_prob_type="clean", #other option is "raw"
                        info_type="log_cell_frac", #options are "post.prob", "cell_frac", "log_cell_frac"
                        prob_threshold_to_include=0.5,
                        plot_cell_frac=FALSE,
                        plot_donut=TRUE, #plot a donut for positive nodes
                        donut_info="cell_frac", #options are "cell_frac" or "lineages_lost"
                        CI=0.95,
                        radius=4,
                        scale_muts_to_branch = TRUE)

dev.off()

#------------------------------------------------------------------------------------------------
#MUTATION TIMING TREE
#------------------------------------------------------------------------------------------------
#Each divergence event is assigned a number from 1 (earliest diverging) to 0 (latest diverging).
#These numbers are arbitrary, but facilitate plotting on the tree & tracking the divergence.
mutation_timing_values=c(1,0.85,0.6,0.5,0.4,0)
names(mutation_timing_values)=c("pre_ICM","pre_epiblast","pre_mesoderm","pre_lateral_plate","pre_haem","post_haem")
raw_colours=c("black","purple","magenta","blue","green","orange","red")
colfunc = colorRampPalette(raw_colours)
col.scale = colfunc(101)
my_pal=col.scale[mutation_timing_values*100]

#Define the timing of each mutation based on presence/ absence in different tissues going from the earliest diverging (trophoblast) to latest (heart & limb)
clean.post.prob <- clean.post.prob[details_targ$mut_ref,]
mutation_timing = sapply(1:nrow(clean.post.prob), function(i) {
  if(clean.post.prob[i,"ECTODERM_comb"] >0.5) {
    return(mutation_timing_values["pre_mesoderm"])
  } else if(clean.post.prob[i,"KIDNEY_comb"] >0.5) {
    return(mutation_timing_values["pre_lateral_plate"])
  } else if(clean.post.prob[i,"HEART_comb"] > 0.5 |clean.post.prob[i,"MESODERM_comb"] > 0.5|clean.post.prob[i,"ALL_LCM"] > 0.5) {
    return(mutation_timing_values["pre_haem"])
  } else {
    return(mutation_timing_values["post_haem"])
  }
})
details_targ$mutation_timing=mutation_timing

pdf("Figures/18pcw/mutation_timing_tree_18pcw.pdf",width=10,height=6)
tree.multi.squash=plot_tree(tree.multi.squash, cex.label = 0)
add_annotation(tree=tree.multi.squash,
               details=details_targ,
               list(mtr=NV[details_targ$mut_ref,],dep=NR[details_targ$mut_ref,]),
               annot_function=function(tree,details,matrices,node) {
                 add_var_col(tree,
                             details,
                             matrices,
                             node,
                             var_field = "mutation_timing",
                             pval_based=FALSE,
                             colours=raw_colours,
                             scale_muts_to_branch = FALSE,
                             lwd=3)
               }
)
draw.circle(1,10,radius=3,nv=100,col=my_pal[3])
text(3,10,"= Pre mesoderm commitment",pos=4)
draw.circle(1,8,radius=3,nv=100,col=my_pal[4])
text(3,8,"= Pre lateral plate mesoderm commitment",pos=4)
draw.circle(1,6,radius=3,nv=100,col=my_pal[5])
text(3,6,"= Pre haematopoietic commitment",pos=4)
dev.off()

#------------------------------------------------------------------------------------------------
#ANALYSIS TO CALCULATE MUTATION BURDEN DISTRIBUTIONS & NUMBERS OF LINEAGES BY DEVELOPMENTAL STAGE
#------------------------------------------------------------------------------------------------
#Get subset tree of pre-gastrulation mutations
#Need to redefine the "get_subset_tree" function to filter based on those mutations with a timing > a given value
get_subset_tree=function(tree, details, v.field = "Mut_type", value = "SNV") {
  get_new_edge_length = function(node, tree, details,v.field,value) {
    sum(details$node == node & details[v.field]>=value)
  }
  tree_subset = tree
  tree_subset$edge.length = sapply(tree$edge[,2], get_new_edge_length, tree = tree, details = details,v.field = v.field,value=value)
  return(tree_subset)
}

#Function to "prune" the tree of lineages that don't exist at the given developmental stage
#If a tip has zero branch length & so does its direct antecedent branch -> remove that tip, as these lineages did not exist at this time point
prune_tree_of_zero_tips = function(tree) {
  current_tree=NULL
  for(i in 1:20) {
    if(is.null(current_tree)){current_tree<-tree} else {current_tree<-tree_pruned}
    current_tree$tip.label<-1:length(current_tree$tip.label)
    tips_to_remove=NULL
    for(i in 1:length(current_tree$tip.label)) { #iterate through all the terminal branches
      private_branch_length=current_tree$edge.length[current_tree$edge[,2]==i] #get the terminal branch length
      if(private_branch_length>0) { #if this has mutations, this branch must be retained
        next
      } else {
        ancestor=current_tree$edge[current_tree$edge[,2]==i,1] #get the direct ancestor node
        if(ancestor==current_tree$edge[1,1]) { 
          ancestral_branch_length<-1 #i.e. if the ancestor node is the root, make the branch length positive
        } else {
          ancestral_branch_length=current_tree$edge.length[current_tree$edge[,2]==ancestor] #otherwise get the ancestral branch length
          }
        all_daughters=current_tree$edge[current_tree$edge[,1]==ancestor,2] #get all the daughter branches from this direct ancestor node
        daughter_branch_lengths<-current_tree$edge.length[current_tree$edge[,2]%in%all_daughters] #get the daughter branch lengths
        if(ancestral_branch_length==0|all(daughter_branch_lengths==0)) { #if the ancestral branch length is 0 or all the daughter branches are zero, "prune" this tip
          tips_to_remove=c(tips_to_remove,i)
        }
      }
    }
    if(is.null(tips_to_remove)) {stop(return(current_tree))} #if no branches are removed in an iteration, this is the final tree to be returned
    tree_pruned=drop.tip(current_tree,trim.internal=FALSE,current_tree$tip.label[tips_to_remove])
  }
}

get_mut_burden_of_ancestral_cells = function(tree,details, mutation_timing_value) {
  #Shorten all branches so that only mutations with a given "mutation timing value" or above are retained
  tree.ancestral = get_subset_tree(tree, details, v.field = "mutation_timing",value = mutation_timing_value)
  #Now drop the tips that have been shortened to zero using the "pruning" function
  tree_pruned=prune_tree_of_zero_tips(tree.ancestral)
  #However, the pruned tree can still have multiple length 0 terminal nodes at polytomies.
  #We do not know if each of these was a separate lineage at the time of commitment, or if only one gave rise to the others. Therefore, remove all but one.
  check_nodes<-unique(tree_pruned$edge[tree_pruned$edge[,2]%in%1:length(tree_pruned$tip.label),1]) #Get a list of all the unique "pre-terminal" nodes i.e. all those nodes that give rise to terminal nodes.
  drop_tips=NULL #Create a empty vector of tips that will be dropped
  for(i in check_nodes) { #Go through each "pre-terminal node"
    all_daughters=tree_pruned$edge[tree_pruned$edge[,1]==i,2] #Get all the daughter branches
    daughter_branch_lengths<-tree_pruned$edge.length[tree_pruned$edge[,2]%in%all_daughters] #Get the branch lengths of these
    matched_tips=tree_pruned$edge[tree_pruned$edge[,2]%in%all_daughters,2] #Get the terminal node numbers in the same order
    zero_tips=matched_tips[daughter_branch_lengths==0] #Select the terminal nodes that have a branch length of 0
    drop_tips=c(drop_tips,zero_tips[-1]) #If there is more than 1 with a branch length of 0, keep only one, drop the others by adding to the "drop_tips" vector
  }
  tree_pruned=drop.tip(tree_pruned,trim.internal=FALSE,tree_pruned$tip.label[drop_tips])
  mut_burden=get_mut_burden(tree_pruned)
  return(mut_burden)
}

#Run this function for each "developmental time point"
ancestral_mut_burdens=lapply(mutation_timing_values[mutation_timing_values%in%details_targ$mutation_timing],function(i) {
  get_mut_burden_of_ancestral_cells(tree=tree.multi,details=details_targ,mutation_timing_value = i)
})

#Combine these into a single dataframe
df=Reduce(rbind,mapply(FUN=function(x,timing) return(data.frame(timing=timing,mutation_burden=x)),x=ancestral_mut_burdens,timing=names(ancestral_mut_burdens),SIMPLIFY=FALSE))
df$timing=factor(df$timing,levels=names(mutation_timing_values))

p1<-df%>%
  filter(timing!="post_haem")%>%
  ggplot(aes(x=mutation_burden,fill=timing)) +
  geom_density()+
  labs(fill="Developmental stage")+
  facet_grid(rows=vars(timing),scales="free_y") +
  theme_classic() +
  theme(legend.position = "none")+
  scale_fill_manual(values=my_pal[mutation_timing_values%in%details_targ$mutation_timing]) +
  ggtitle("Mutation burden distribution by developmental stage")

p2<-df%>%
  filter(timing!="post_haem")%>%
  group_by(timing)%>%
  dplyr::summarise(mean_mutation_burden=mean(mutation_burden),median=median(mutation_burden),n=n()) %>%
  ggplot(aes(y=mean_mutation_burden,x=timing,fill=timing)) +
  geom_bar(stat="identity",col="black") +
  labs(fill="Developmental stage")+
  theme_classic() +
  theme(legend.position = "none")+
  scale_fill_manual(values=my_pal[mutation_timing_values%in%details_targ$mutation_timing]) +
  ggtitle("Mean mutation burden by developmental stage")

p3<-df%>%
  filter(timing!="post_haem")%>%
  group_by(timing)%>%
  dplyr::summarise(mean_mutation_burden=mean(mutation_burden),median=median(mutation_burden),lineages_present=n()) %>%
  ggplot(aes(y=lineages_present,x=timing,fill=timing)) +
  geom_bar(stat="identity",col="black") +
  labs(fill="Developmental stage")+
  theme_classic() +
  scale_fill_manual(values=my_pal[mutation_timing_values%in%details_targ$mutation_timing]) +
  ggtitle("Number of lineages captured by developmental stage")

#Combine and save this plot
comb_plot<-arrangeGrob(p1,p2,p3,ncol=3)
ggsave(comb_plot,filename = "Figures/18pcw/Ancestral_mutations_and_lineages.pdf",device= "pdf",width = 15,height=6)

#---------------------------------------------------------------------------
#SUMMING CELL LINEAGES CAPTURED FOR A GIVEN TISSUE ACROSS ENTIRE GENERATIONS
#---------------------------------------------------------------------------

#Try to pick apart branches that represent multiple generations
#All combined LCM counts give the best power to pick these out - therefore add this column
NV<-cbind(NV,matrix(apply(NV[,colnames(NV)%in%lcm_smry$Sample_ID],1,sum),ncol=1,dimnames=list(NULL,"ALL_LCM")))
NR<-cbind(NR,matrix(apply(NR[,colnames(NR)%in%lcm_smry$Sample_ID],1,sum),ncol=1,dimnames=list(NULL,"ALL_LCM")))

#Work out branches/ mutations for each generation
#Gens 1 & 2 are straight-foward
gen_branches_1=236 #Don't include the minor branch as this is multiple generations & does not contribute to embryo.
gen_branches_2=tree_targ$edge[tree_targ$edge[,1]%in%236,2] #237 and 305

#Gen 3 includes: early muts of 238, 284, early muts on node 305 daughters (306/341/361/367/368/383)
#STRATEGY IS: find all branches with any mutations contributing & then exclude those that are not that generation
gen_branches_3=tree_targ$edge[tree_targ$edge[,1]%in%gen_branches_2,2]
#exclude the late muts from 238
late_238=find_early_muts_from_branch("ALL_LCM",node=238,tree=tree_targ,details=details_targ,matrices = list(NV=NV,NR=NR),return_late_muts = T)
#exclude the late muts from the node 305 daughters
daughters<-Descendants(tree_targ,305,type = "children")
single_dist=sapply(daughters,function(x) check_branch_distribution("ALL_LCM",node=x,tree=tree_targ,details=details_targ,matrices=list(NV=NV,NR=NR)))
names(single_dist)<-daughters
mixed_dist_branches=single_dist[single_dist<0.05 & !is.na(single_dist)] #3 branches (367,368 & 383) have mixed distributions
#binom_mix function fails for branch 383 - therefore only apply to 367 & 368
late_367_368=Reduce(c,lapply(c(367,368),function(x) {find_early_muts_from_branch("ALL_LCM",node=x,tree=tree_targ,details=details_targ,matrices = list(NV=NV,NR=NR),return_late_muts = T)}))
#Manual execution for branch 383
node_df=check_branch_distribution("ALL_LCM",383,tree_targ,details_targ,list(NV=NV,NR=NR),return_counts = TRUE)
node_df$vaf=node_df$NV/node_df$NR #highest vaf is idx 1010, therefore exclude 608 & 2287
late_383=c(608,2287)
not_gen_3=c(late_238,late_367_368,late_383)

#Gen 4 includes: late muts of 238, early muts from node 284 daughters (285/293/91/298/301/106), node 305 daughters (306/341/361/367/368/383) - later muts only if split
#STRATEGY IS: exclude early muts of 238,late muts from 284 daughters, early muts from 305 daughters (NB, some muts from 305 daughters will be included in gen 3 & gen 4)
gen_branches_4=c(238,Descendants(tree_targ,284,"children"),Descendants(tree_targ,305,"children"))
early_238=find_early_muts_from_branch("ALL_LCM",node=238,tree=tree_targ,details=details_targ,matrices = list(NV=NV,NR=NR),return_late_muts = F)
#exclude the late muts from the node 284 daughters
daughters<-Descendants(tree_targ,284,type = "children")
single_dist=sapply(daughters,function(x) check_branch_distribution("ALL_LCM",node=x,tree=tree_targ,details=details_targ,matrices=list(NV=NV,NR=NR)))
names(single_dist)<-daughters
mixed_dist_branches=single_dist[single_dist<0.05 & !is.na(single_dist)] #3 branches (367,368 & 383) have mixed distributions
late_284_daughters=Reduce(c,lapply(names(mixed_dist_branches),function(x) {find_early_muts_from_branch("ALL_LCM",node=x,tree=tree_targ,details=details_targ,matrices = list(NV=NV,NR=NR),return_late_muts = T)}))
#exclude early muts from mixed 305 daughters (this is 367,368 and 383)
node_df=check_branch_distribution("ALL_LCM",367,tree_targ,details_targ,list(NV=NV,NR=NR),return_counts = TRUE)
not_gen_4_367=c(274,367,2474)
early_368=find_early_muts_from_branch("ALL_LCM",368,tree_targ,details_targ,list(NV=NV,NR=NR),return_late_muts = F)
node_df=check_branch_distribution("ALL_LCM",383,tree_targ,details_targ,list(NV=NV,NR=NR),return_counts = TRUE)
early_383=1010
not_gen_4=c(early_238,late_284_daughters, not_gen_4_367,early_368,early_383)

#Create summary df to get total generation cell 
summary_df=data.frame(tissue=NA,generation=NA,median_cell_frac=NA,CI_95_lower=NA,CI_95_upper=NA)

tissues=names(table(lcm_smry$Tissue))
germ_layers=names(table(lcm_smry$Germ_layer))

for(sample in c(tissues,germ_layers)) {
  #SUM ACROSS THE TREE FOR A GIVEN GENERATION
  print(sample)
  counts=list()
  
  #Get counts df for each generation
  counts$gen_branches_1=Reduce(rbind,lapply(gen_branches_1,get_node_read_counts,sample=sample,tree=tree_targ,details=details_targ,matrices=list(NV=NV,NR=NR)))
  counts$gen_branches_2=Reduce(rbind,lapply(gen_branches_2,get_node_read_counts,sample=sample,tree=tree_targ,details=details_targ,matrices=list(NV=NV,NR=NR)))
  counts$gen_branches_3=Reduce(rbind,lapply(gen_branches_3,get_node_read_counts,sample=sample,tree=tree_targ,details=details_targ,matrices=list(NV=NV,NR=NR),exclude_mut_indexes=not_gen_3))
  counts$gen_branches_4=Reduce(rbind,lapply(gen_branches_4,get_node_read_counts,sample=sample,tree=tree_targ,details=details_targ,matrices=list(NV=NV,NR=NR),exclude_mut_indexes=not_gen_4))
  
  #For each tissue have a "Generation 0" entry with a cell fraction of 1
  summary_df=rbind(summary_df,c(sample,0,1,1,1))
  
  for(generation in 1:4) {
    print(generation)
    #Subset the counts df for that generation
    gen_df=counts[[generation]]
    #Bootstrap the counts to get cell fracs for each lineage
    cell_fracs=Reduce(rbind,lapply(1:nrow(gen_df), function(i) {bootstrap_counts(gen_df[i,])}))
    #Sum across bootstraps to get bootstrapped sums
    if(nrow(gen_df)==1) {total_cell_fracs<-cell_fracs} else {total_cell_fracs=colSums(cell_fracs)}
    #Total fractions > 1 don't make sense, therefore cooerce all values >1 to 1
    total_cell_fracs[total_cell_fracs>1]<-1
    #Add
    summary_df=rbind(summary_df,data.frame(tissue=sample,generation=generation,median_cell_frac=median(total_cell_fracs),CI_95_lower=quantile(total_cell_fracs,0.025),CI_95_upper=quantile(total_cell_fracs,0.975)))
  }
}

summary_df[,2:5]<-apply(summary_df[,2:5],2,as.numeric)
str(summary_df)

#Visualization options
to_visualize=tissues
to_visualize=germ_layers
to_visualize=c("ECTODERM","KIDNEY","HEART")
to_visualize=c("HAIR FOLLICLE","EPIDERMIS")

#Use ggplot2 to create "lineage loss over generations"
summary_df%>%
  filter(tissue%in%to_visualize) %>%
  ggplot(aes(col=tissue,fill=tissue,x=generation,y=median_cell_frac)) +
  geom_point() +
  geom_line(lwd=1) +
  geom_ribbon(aes(col=NULL,fill=tissue,ymin=CI_95_lower,ymax=CI_95_upper),alpha=0.2) +
  scale_y_continuous(limits=c(0,1)) +
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2") +
  theme_classic()+
  ggtitle("Lineages captured in different tissues through cell generations")

ggsave(p1,file=filename = "Figures/18pcw//Lineages_captured_through_cell_gens_18pcw.pdf",width = 7,height=5)

#------------------------------------------------------------------------------------------------------------------------------------------------
##CALCULATE THE PROPORTION OF BLOOD ANTECEDENTS FROM DIFFERENT GENERATIONS FROM THE ZYGOTE WITH DETECTABLE PROGENY IN NON-HAEMATOPOIETIC TISSUES
#------------------------------------------------------------------------------------------------------------------------------------------------

#(1) Take the ultrametric generation time tree
#"Remove" the minor branch by setting all edge lenths to 0 - this is not represented in any of the targeted sequencing tissues and confuses things
minor_branch_nodes = c(389,get_all_node_children(389,tree.gen.ultra))
tree.gen.ultra$edge.length[tree.gen.ultra$edge[,2] %in% minor_branch_nodes] <- 0
tree.gen.ultra$coords<-NULL #reset the coords

#(2)
nodeheights <- nodeHeights(tree.gen.ultra) #get node heights of each edge
generation_cut_offs = seq(from=0,to=7.8,by=0.2) #decide what "generations" to assess. Beyond 8 becomes a bit meaningless as true generations is impossible to know.

#(3) The main function to capture the "lineages represented" statistic
#This iterates through the generations defined above by "cutting" horizontally across the generations tree at each heigh, and finds which branches are present
#Returns the proportion of branches that have any detectable mutations in each tissue at that point
captured_by_gen=lapply(generation_cut_offs, function (x) {
  #Define the branches to include at this cut off
  branches_to_include=tree.gen.ultra$edge[which(nodeheights[,1] <= x & !nodeheights[,2] <= x),2] #get the branches crossing each cut_off
  #Capture the total number of branches - the denominator
  total_lineages=length(branches_to_include)
  #Find whether any of the mutations in these branches are likely present (post.prob>0.5) in each tissue
  #If any are, this whole branch is counted as "represented" in the tissue
  lineages_represented=rowSums(sapply(branches_to_include,function(node) {
    apply(clean.post.prob[details_targ$node==node,,drop=FALSE],2,function(x) {any(x>prob_threshold_to_include)})
  }))
  fraction_represented=lineages_represented/total_lineages
  return(fraction_represented)
})

#Clean up output into dataframe
captured_by_gen=as.data.frame(Reduce(rbind,captured_by_gen))

#Make the colnames the corresponding tissue
colnames(captured_by_gen) <- colnames(clean.post.prob)

#Add a generation column
captured_by_gen$generation<-generation_cut_offs

#Gather into tidy df & visualize with ggplot
p1<-gather(captured_by_gen,key = sample,value = lineages_represented,-generation) %>%
  filter(sample%in%germ_layers_comb & sample!="UNKNOWN_comb")%>%
  ggplot(aes(x=generation,y=lineages_represented,col=sample)) +
  geom_line(lwd=1.5) +
  scale_x_continuous(breaks = seq(0,12,1)) +
  scale_y_continuous(breaks=seq(0,1,0.1),limits = c(0,1)) +
  theme_classic() +
  ggtitle("Proportion of blood lineages found in tissue through generations")

ggsave(p1,filename = "Figures/18pcw/Fraction_of_lineages_represented_18pcw.pdf",device = "pdf",width=7,height=5)

#---------------------------------------------------------------------------
### PLOT A SCALE BAR ###
#---------------------------------------------------------------------------
#This is to plot a colour scale bar to visualize how the colour scales on the log_vaf tree corresponds to the vaf

#This bit rescales the scale bar using the same transformations applied to the vaf data
scale_bar=rev(1/exp(seq(0,8))) #These are the tick points that will go on the scale
log_scale_bar=log(scale_bar)
prob_threshold_to_include=0.5 #This should be same as for the within the "generate_targeted_seq_plots" function

#This bit takes the elements of the plotrix rescale function
scale_range=c(0.01,1)

#Get the vaf_present & log_vaf_present objects (performed as within the generate_targ_seq_plots function)
cell_frac=calculate_cell_frac(NV,NR)
cell_frac_present<-cbind(cell_frac,cell_frac[,gsub("_comb","",colnames(post.prob.mat)[grep("_comb",colnames(post.prob.mat))])]) #To double up the "_comb" results, so that matched the post.prob.mat
cell_frac_present=cell_frac_present[details_targ$mut_ref,]
colnames(cell_frac_present)<- colnames(post.prob)
cell_frac_present[post.prob.mat<prob_threshold_to_include]<-0
log_cell_frac_present=cell_frac_present
log_cell_frac_present[cell_frac_present != 0] <- log(log_cell_frac_present[cell_frac_present != 0])#change all the non 0 vafs to the log of their vaf

#Use these this to perform the identical transformations on the scale bar
xrange=range(log_cell_frac_present[cell_frac_present != 0])
mfac=(scale_range[2]-scale_range[1])/(xrange[2]-xrange[1])
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
pdf("Figures/18pcw/Trees_colour_scale_bar_18pcw.pdf")
plot(c(0,10), c(0,1), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main="")
axis(side=2,pos=1,las=1,at=tick_positions_trimmed,labels=as.character(tick_values_trimmed))
for (i in 1:(length(lut)-1)) {
  y = (i-1)/scale
  rect(1,y,2,y+1/scale, col=lut[i], border=NA)
}
dev.off()
