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

my_working_directory="~/R_work/Phylogeny_of_foetal_haematopoiesis/"
treemut_dir="~/R_work/treemut" #Path for the cloned "treemut" directory
setwd(my_working_directory)

#Define the file paths for the data files
tree_file_path="Data/8pcw/Tree_8pcw.tree"
file_annot="Data/8pcw/Filtered_mut_set_annotated_8pcw"
cgpvaf_LCM_SNV_file ="Data/8pcw/PD43947.LCM.snp.tsv"
cgpvaf_LCM_INDEL_file ="Data/8pcw/PD43947.LCM.indel.tsv"
cgpvaf_colonies_SNV_file ="Data/8pcw/PD43947.colonies.snp.tsv"
cgpvaf_colonies_INDEL_file ="Data/8pcw/PD43947.colonies.indel.tsv"
lcm_smry_path = "Data/8pcw/targeted_seq_ids_8wks.csv"

R_function_files = list.files("R_functions",pattern=".R",full.names=TRUE)
sapply(R_function_files,source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)

#--------------------------------------------------------------------------------------------------
##IMPORTING THE CGPVAF TABLES (SNVs and INDELs), AND AGGREGATING COUNTS ACROSS TISSUES/ GERM LAYERS
#--------------------------------------------------------------------------------------------------
#Import the cgpvaf files and combine the LCM and single-cell colony (SCC) matrices into single matrix
LCM_mats<-import_cgpvaf_SNV_and_INDEL(SNV_output_file=cgpvaf_LCM_SNV_file,INDEL_output_file=cgpvaf_LCM_INDEL_file)
SCC_mats<-import_cgpvaf_SNV_and_INDEL(SNV_output_file=cgpvaf_colonies_SNV_file,INDEL_output_file=cgpvaf_colonies_INDEL_file)

#Import the LCM summary table
lcm_smry = read.csv(lcm_smry_path,stringsAsFactors = FALSE)

#Create aggregated read counts for each tissue/ germ-layer/ pre- or post-gastrulation category (as defined in the lcm_smry table)
Tissue_agg=aggregate_cols(LCM_mats$NV,LCM_mats$NR,lcm_summary_table = lcm_smry,aggregate_feature = "Tissue")
Germ_layer_agg=aggregate_cols(LCM_mats$NV,LCM_mats$NR,lcm_summary_table = lcm_smry,aggregate_feature = "Germ_layer")
Gastrulation_agg=aggregate_cols(LCM_mats$NV,LCM_mats$NR,lcm_summary_table = lcm_smry,aggregate_feature = "Gastrulation",exclude_class = "UNKNOWN")

#Combine the individual LCM sample counts with the aggregated counts
NV<-cbind(LCM_mats$NV,Tissue_agg$NV,Germ_layer_agg$NV,Gastrulation_agg$NV)
NR<-cbind(LCM_mats$NR,Tissue_agg$NR,Germ_layer_agg$NR,Gastrulation_agg$NR)

#Import the details matrix and tree (needed to interpret the targeted sequencing data)
load(file_annot); tree.8wk <- read.tree(tree_file_path); tree.8wk$tip.label=gsub("_hum","",tree.8wk$tip.label)
details <- filtered_muts$COMB_mats.tree.build$mat

#In details file, create column for whether it is included in the targeted sequencing baits
details$targeted_tree = ifelse(details$mut_ref %in% rownames(NV),"YES","NO")

#Now subset the tree and the details matrix to only include those mutations that are in the bait set
tree_targ = get_subset_tree(tree = tree.8wk, details = details, v.field = "targeted_tree",value = "YES")
details_targ = details[details$targeted_tree=="YES",]

NV=NV[details_targ$mut_ref,]
NR=NR[details_targ$mut_ref,]

#--------------------------------------------------------------------------------------------------
##ESTIMATE THE POSTERIOR PROBABILITY OF EACH MUTATION BEING PRESENT IN A GIVEN TISSUE##
#--------------------------------------------------------------------------------------------------
depth.cols=SCC_mats$NR[details_targ$mut_ref,]; mtr.cols=SCC_mats$NV[details_targ$mut_ref,] #The matrices for calculating background error (single cell colonies)

#Only use single cell colonies that are included in the tree
depth.cols <- depth.cols[,names(depth.cols) %in% tree.8wk$tip.label]
mtr.cols <- mtr.cols[,names(mtr.cols) %in% tree.8wk$tip.label]

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
descendants_per_mut <- Descendants(x = tree.8wk, node = details_targ$node, type = "tips")
for (i in 1:nrow(NR)) {
  desc_tips <- tree.8wk$tip.label[descendants_per_mut[[i]]]
  
  # Deal with case of lesion segregation!
  if (row.names(NR)[i] == "22-37053571-G-A") {
    desc_tips <- c(desc_tips, tree.8wk$tip.label[Descendants(tree.8wk, node=527, type="tips")[[1]]])
  }
  
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
prior_mut <- sapply(1:nrow(NR), function(i) {length(Descendants(tree.8wk, details_targ$node[i], "tips")[[1]]) / length(tree.8wk$tip.label)})

# Calculate posterior probabilities
post.prob <- cbind(sapply(1:ncol(NR), function(i) {bbprob.calculator(x_i = NV[,i], n_i = NR[,i], alpha_error = alpha_error, beta_error = beta_error, alpha_mut = alpha_mut, beta_mut = beta_mut, prior_mut = prior_mut)}))
colnames(post.prob) <- names(NR); row.names(post.prob) <- row.names(NR)
post.prob[post.prob<0.05]<-0

#########################################################################
#In a given tissue, individual cuts may be clonal/ oligoclonal, and therefore in particular lineages may be
#represented at a higher vaf than the aggregated read counts for the tissue.  Therefore, reassessing aggregated
#counts may dilute the ability to call mutations. Therefore call mutations based on 1-prod(1-x) where x is the individual
#sample probability i.e. 1 - the probability that the mutation is absent in all individual samples.

post.prob<-cbind(post.prob,
                 combine_post_prob(post.prob,lcm_smry,aggregate_feature="Tissue"),
                 combine_post_prob(post.prob,lcm_smry,aggregate_feature="Germ_layer"),
                 combine_post_prob(post.prob,lcm_smry,aggregate_feature="Gastrulation",exclude_class = "UNKNOWN"))

#########################################################################
#CREATE FUNCTION TO ADJUST POST.PROBS BASED ON PHYLOGENY i.e. "clean it up" to make sense with the phylogeny
#This function does two main alterations:
#(1) Removes positive posterior probabilities when there are no other mutations called on same branch or ancestral branch
#(2) Boosts the posterior probabilities when there are high posterior probabilities in descendant mutations & other mutations on the same branch & ancestral mutations
re_run=F
if(re_run) {
  clean.post.prob=apply(post.prob,2,clean_up_post.multi,details_targ,tree.8wk)
  write.table(clean.post.prob,file = "Data/8pcw/clean.post.prob_8wks.tsv",quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
} else { 
  clean.post.prob=as.matrix(read.delim("Data/8pcw/clean.post.prob_8wks.tsv",sep="\t")) #Re-import it if you have previously run the analysis
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
tree.gen$edge.length[tree.gen$edge[,2]==280] <-2 #Rescale node 280 - we have good evidence that this is two generations from the trophoblast targeted data
tree.squash=squash_tree(tree_targ,cut_off=22) #Create a "squashed" tree for displaying this way

#Extend the branches of the generations tree for clarity, and cut at 10 generations
tree.gen.ultra=force.ultrametric(tree.gen,method="extend")
tree.gen.ultra.squash=squash_tree(tree.gen.ultra,cut_off=10)

#########################################################################
####NOW READY TO PLOT THE TREES IN A VARIETY OF WAYS ####

#For this function to run, need several objects in the environment:
#(1) post.prob & clean.post.prob matrices
#(2) the vaf matrix (this gets transformed to the log_vaf_present matrix & scaled within the function)
#(3) the vectors: tissues, germ_layers, tissues_comb and germ_layers_comb
#(4) note that the post.prob & clean.post.prob matrices are expected to include the "_comb" versions for tissues & germ_layers
#(6) A colour.scale object containing the colour scale for the trees (and the scale bar)

#Create vectors of the two "comb" data types
tissues_comb=paste0(names(table(lcm_smry$Tissue)),"_comb")
germ_layers_comb=paste0(names(table(lcm_smry$Germ_layer)),"_comb")

#Do some plots using the "generate_targ_seq_plots" function
pdf("Figures/8pcw/Gut_samples_cellfrac.8wks.pdf",width=13,height=7)
colour.scale <- c("lightgrey", brewer.pal(9, name = "YlOrRd")); colour.scale<-colour.scale[c(3:9,9)]
samples=lcm_smry$Sample_ID[lcm_smry$Tissue=="GUT"]
generate_targ_seq_plots(samples,
                        tree=tree.gen,
                        details_targ=details_targ,
                        matrices=list(NV=NV[details_targ$mut_ref,],NR=NR[details_targ$mut_ref,]),
                        post_prob_type="clean", #other option is "raw"
                        info_type="log_cell_frac", #options are "post.prob", "cell_frac", "log_cell_frac"
                        prob_threshold_to_include=0.5,
                        plot_cell_frac=F,
                        plot_donut=F, #plot a donut for positive nodes
                        donut_info="cell_frac", #options are "cell_frac" or "lineages_lost"
                        CI=0.95,
                        radius=4,
                        scale_muts_to_branch = TRUE)

dev.off()

pdf("Figures/8pcw/Aggregated_tissues_cellfrac.8wks.pdf",width=13,height=7)
colour.scale <- c("lightgrey", brewer.pal(9, name = "YlOrRd")); colour.scale<-colour.scale[c(3:9,9)]
samples=tissues_comb
generate_targ_seq_plots(samples,
                        tree=tree.squash,
                        details_targ=details_targ,
                        matrices=list(NV=NV[details_targ$mut_ref,],NR=NR[details_targ$mut_ref,]),
                        post_prob_type="clean", #other option is "raw"
                        info_type="log_cell_frac", #options are "post.prob", "cell_frac", "log_cell_frac"
                        prob_threshold_to_include=0.5,
                        plot_cell_frac=FALSE,
                        plot_donut=FALSE, #plot a donut for positive nodes
                        donut_info="cell_frac", #options are "cell_frac" or "lineages_lost"
                        CI=0.95,
                        radius=4,
                        scale_muts_to_branch = TRUE)

dev.off()


pdf("Figures/8pcw/Germ_layers_cellfracpie.8wks.pdf",width=13,height=7)
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

#MUTATION TIMING TREE
#Each 
mutation_timing_values=c(1,0.85,0.6,0.5,0.4,0)
names(mutation_timing_values)=c("pre_ICM","pre_epiblast","pre_mesoderm","pre_lateral_plate","pre_haem","post_haem")

colfunc = colorRampPalette(c("black","purple","magenta","blue","green","orange","red"))
col.scale = colfunc(101)
my_pal=col.scale[mutation_timing_values*100]

#Define the timing of each mutation based on presence/ absence in different tissues going from the earliest diverging (trophoblast) to latest (heart & limb)
clean.post.prob <- clean.post.prob[details_targ$mut_ref,]
mutation_timing = sapply(1:nrow(clean.post.prob), function(i) {
  if(clean.post.prob[i,"TROPHECTODERM"] >0.5) {
    return(mutation_timing_values["pre_ICM"])
  } else if(clean.post.prob[i,"EXTRA_EMBRYONIC_MESODERM_comb"] >0.5) {
    return(mutation_timing_values["pre_epiblast"])
  } else if(clean.post.prob[i,"ECTODERM_comb"] >0.5) {
    return(mutation_timing_values["pre_mesoderm"])
  } else if(clean.post.prob[i,"ENDODERM_comb"] >0.5) {
    return(mutation_timing_values["pre_mesoderm"]) 
  } else if(clean.post.prob[i,"KIDNEY_comb"] >0.5|clean.post.prob[i,"INTER_VERTEBRAL_DISC_comb"] >0.5) {
    return(mutation_timing_values["pre_lateral_plate"])
  } else if(clean.post.prob[i,"HEART_comb"] >0.5|clean.post.prob[i,"LIMB_comb"]>0.5) {
    return(mutation_timing_values["pre_haem"])
  } else {
    return(mutation_timing_values["post_haem"])
  }
})
details_targ$mutation_timing=mutation_timing

pdf("Figures/8pcw/mutation_timing_tree.pdf",width=10,height=6)
tree.squash=plot_tree(tree.squash, cex.label = 0)
add_annotation(tree=tree.squash,
               details=details_targ,
               list(mtr=mtr,dep=dep),
               annot_function=function(tree,details,matrices,node) {
                 add_var_col(tree,
                             details,
                             matrices,
                             node,
                             var_field = "mutation_timing",
                             pval_based=FALSE,
                             colours=c("black","purple","magenta","blue","green","orange","red"),
                             scale_muts_to_branch = FALSE,
                             lwd=3)
               }
)
draw.circle(1,14,radius=3,nv=100,col=my_pal[1])
text(3,14,"= Pre ICM commitment",pos=4)
draw.circle(1,12,radius=3,nv=100,col=my_pal[2])
text(3,12,"= Pre epiblast commitment",pos=4)
draw.circle(1,10,radius=3,nv=100,col=my_pal[3])
text(3,10,"= Pre mesoderm commitment",pos=4)
draw.circle(1,8,radius=3,nv=100,col=my_pal[4])
text(3,8,"= Pre lateral plate mesoderm commitment",pos=4)
draw.circle(1,6,radius=3,nv=100,col=my_pal[5])
text(3,6,"= Pre haematopoietic commitment",pos=4)
dev.off()

#ANALYSIS TO CALCULATE MUTATION BURDEN DISTRIBUTIONS BY DEVELOPMENTAL STAGE
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
    for(i in 1:length(current_tree$tip.label)) {
      private_branch_length=current_tree$edge.length[current_tree$edge[,2]==i]
      if(private_branch_length>0) {
        next
      } else {
        ancestor=current_tree$edge[current_tree$edge[,2]==i,1]
        ancestral_branch_length=current_tree$edge.length[current_tree$edge[,2]==ancestor]
        all_daughters=current_tree$edge[current_tree$edge[,1]==ancestor,2]
        daughter_branch_lengths<-current_tree$edge.length[current_tree$edge[,2]%in%all_daughters]
        if(ancestral_branch_length==0|all(daughter_branch_lengths==0)) {
          tips_to_remove=c(tips_to_remove,i)
        }
      }
    }
    if(is.null(tips_to_remove)) {stop(return(current_tree))}
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
ancestral_mut_burdens=lapply(mutation_timing_values,function(i) {
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
  scale_fill_manual(values=my_pal) +
  ggtitle("Mutation burden distribution by developmental stage")

p2<-df%>%
  filter(timing!="post_haem")%>%
  group_by(timing)%>%
  dplyr::summarise(mean_mutation_burden=mean(mutation_burden),median=median(mutation_burden),n=n()) %>%
  ggplot(aes(y=mean_mutation_burden,x=timing,fill=timing)) +
  geom_bar(stat="identity",col="black") +
  labs(fill="Developmental stage")+
  theme_classic() +
  scale_fill_manual(values=my_pal) +
  ggtitle("Mean mutation burden by developmental stage")

p3<-df%>%
  filter(timing!="post_haem")%>%
  group_by(timing)%>%
  dplyr::summarise(mean_mutation_burden=mean(mutation_burden),median=median(mutation_burden),lineages_present=n()) %>%
  ggplot(aes(y=lineages_present,x=timing,fill=timing)) +
  geom_bar(stat="identity",col="black") +
  labs(fill="Developmental stage")+
  theme_classic() +
  scale_fill_manual(values=my_pal) +
  ggtitle("Number of lineages capture by developmental stage")

#Combine and save this plot
comb_plot<-arrangeGrob(p1,p2,p3,ncol=3)
ggsave(comb_plot,filename = "Figures/8pcw/Ancestral_mutations_and_lineages.pdf",device= "pdf",width = 15,height=6)

#---------------------------------------------------------------------------
#SUMMING CELL LINEAGES CAPTURED FOR A GIVEN TISSUE ACROSS ENTIRE GENERATIONS
#---------------------------------------------------------------------------

#Work out branches for each generation
gen_branches_1=tree.multi$edge[tree.multi$edge[,1]==278,2]
gen_branches_2=tree.multi$edge[tree.multi$edge[,1]%in%gen_branches_1,2]
gen_branches_3=c(tree.multi$edge[tree.multi$edge[,1]%in%gen_branches_2[-1],2],280) #Get daughter branches from all except node 280 (2nd half of 280 is 3rd generation)
gen_branches_4=tree.multi$edge[tree.multi$edge[,1]%in%gen_branches_3,2]

#Generations 2 & 3 are complicated -as branch 280 is half 2nd & half 3rd generation
#Work out which muts are which by which are present/ absent in trophoblast
gen_3_muts=names(post.prob[details_targ$node==280,"TROPH"][post.prob[details_targ$node==280,"TROPH"]==0])
gen_2_muts=names(post.prob[details_targ$node==280,"TROPH"][post.prob[details_targ$node==280,"TROPH"]==1])
gen_3_muts_idx=which(details_targ$mut_ref%in%gen_3_muts)
gen_2_muts_idx=which(details_targ$mut_ref%in%gen_2_muts)

#Many mutations on the gen_4_branches are from later generations (i.e. the branches have mutations from multiple generations)
#These lines of code find which mutations are from later generations and therefore not genuinely from gen_4_branches
#The ALL_LCM counts (aggregated counts across all LCM biopsies) give the best power to extract this information
single_dist=sapply(gen_branches_4,function(x) check_branch_distribution("ALL_LCM",node=x,tree=tree.multi,details=details_targ,matrices=list(NV=NV,NR=NR)))
names(single_dist)=gen_branches_4
mixed_dist_branches=single_dist[single_dist<0.05 & !is.na(single_dist)]
not_gen_4=Reduce(c,lapply(names(mixed_dist_branches), function(x) find_early_muts_from_branch("ALL_LCM",node=x,tree=tree.multi,details=details_targ,matrices = list(NV=NV,NR=NR),return_late_muts = TRUE)))

#Create summary df to get total generation cell 
summary_df=data.frame(tissue=NA,generation=NA,median_cell_frac=NA,CI_95_lower=NA,CI_95_upper=NA)

for(sample in c(tissues,germ_layers,same_germ_layer)) {
  #SUM ACROSS THE TREE FOR A GIVEN GENERATION
  print(sample)
  counts=list()
  
  #Get counts df for each generation
  counts$gen_branches_1=Reduce(rbind,lapply(gen_branches_1,get_node_read_counts,sample=sample,tree=tree.multi,details=details_targ,matrices=list(NV=NV,NR=NR)))
  counts$gen_branches_2=Reduce(rbind,lapply(gen_branches_2,get_node_read_counts,sample=sample,tree=tree.multi,details=details_targ,matrices=list(NV=NV,NR=NR),exclude_mut_indexes=gen_3_muts_idx))
  counts$gen_branches_3=Reduce(rbind,lapply(gen_branches_3,get_node_read_counts,sample=sample,tree=tree.multi,details=details_targ,matrices=list(NV=NV,NR=NR),exclude_mut_indexes=gen_2_muts_idx))
  counts$gen_branches_4=Reduce(rbind,lapply(gen_branches_4,get_node_read_counts,sample=sample,tree=tree.multi,details=details_targ,matrices=list(NV=NV,NR=NR),exclude_mut_indexes=not_gen_4))
  
  #For each tissue have a "Generation 0" entry with a cell fraction of 1
  summary_df=rbind(summary_df,c(sample,0,1,1,1))
  
  for(generation in 1:4) {
    print(generation)
    #Subset the counts df for that generation
    gen_df=counts[[generation]]
    #Bootstrap the counts to get cell fracs for each lineage
    cell_fracs=Reduce(rbind,lapply(1:nrow(gen_df), function(i) {bootstrap_counts(gen_df[i,])}))
    #Sum across bootstraps to get bootstrapped sums
    total_cell_fracs=colSums(cell_fracs)
    #Total fractions > 1 don't make sense, therefore cooerce all values >1 to 1
    total_cell_fracs[total_cell_fracs>1]<-1
    #Add
    summary_df=rbind(summary_df,data.frame(tissue=sample,generation=generation,median_cell_frac=median(total_cell_fracs),CI_95_lower=quantile(total_cell_fracs,0.025),CI_95_upper=quantile(total_cell_fracs,0.975)))
  }
}

#Make the numeric columns numeric
summary_df[,2:5]<-apply(summary_df[,2:5],2,as.numeric)

#Visualization options - several options, but the figure in the paper is just germ layers
to_visualize=c("TROPH","MESENCHYME","PERIPHERAL_BLOOD","POST_GASTRULATION")
# extra_embryonic=c("PERIPHERAL_BLOOD","TROPH","MESENCHYME")
# to_visualize=lcm_smry$Sample_ID[lcm_smry$Tissue=="PERIPHERAL_BLOOD"]
# to_visualize=tissues
# to_visualize=c("MESODERM","ENDODERM","ECTODERM")
# to_visualize=germ_layers[!germ_layers=="UNKNOWN"]

#Use ggplot2 to create "lineage loss over generations"
color.cat=c("#d7529f", "#68374f", "#ae7196", "#8c0250", "#609111", "#d6061a", "#798872", "#743502", "#b5753e", "#1f5113")

p2<-summary_df%>%
  filter(tissue%in%to_visualize) %>%
  ggplot(aes(col=tissue,fill=tissue,x=generation,y=median_cell_frac)) +
  geom_point() +
  geom_line(lwd=1) +
  geom_ribbon(aes(col=NULL,fill=tissue,ymin=CI_95_lower,ymax=CI_95_upper),alpha=0.2) +
  scale_y_continuous(limits=c(0,1)) +
  scale_discrete_manual(color.cat) +
  theme_classic()+
  labs(col="Tissue",fill="Tissue") +
  ggtitle("Lineages captured in different tissues through cell generations")

ggsave(p2,file="Figures/8pcw/Lineages_captured_through_cell_gens_8pcw.pdf",width = 7,height=5)

#------------------------------------------------------------------------------------------------------------------------------------------------
##CALCULATE THE PROPORTION OF BLOOD ANTECEDENTS FROM DIFFERENT GENERATIONS FROM THE ZYGOTE WITH DETECTABLE PROGENY IN NON-HAEMATOPOIETIC TISSUES
#------------------------------------------------------------------------------------------------------------------------------------------------

#(1) Use the tree.gen.ultra tree previously generated
nodeheights <- nodeHeights(tree.gen.ultra) #get node heights of each edge
generation_cut_offs = seq(from=0,to=7.8,by=0.2) #decide what "generations" to assess. Beyond 8 becomes a bit meaningless as true generations is impossible to know.


#(2) The main function to capture the statistic
captured_by_gen=lapply(generation_cut_offs, function (x) {
  #Define the branches to include at this cut off
  branches_to_include=tree.gen.ultra$edge[which(nodeheights[,1] <= x & !nodeheights[,2] <= x),2] #get the branches crossing each cut_off
  #Capture the total number of branches - the denominator
  total_lineages=length(branches_to_include)
  #Find whether any of the mutations in these branches are likely present in each tissue (as previous)
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

#Need to manually clean-up that the second half of branch 280 is absent from troph (the above function cannot split branches)
captured_by_gen[captured_by_gen$generation>2,"TROPH"] <-0
captured_by_gen[captured_by_gen$generation>2,"TROPHECTODERM_comb"] <-0
captured_by_gen[captured_by_gen$generation>2,"TROPHECTODERM"] <-0

#Gather into tidy df & visualize with ggplot
p1<-gather(captured_by_gen,key = sample,value = lineages_represented,-generation) %>%
  filter(sample%in%germ_layers_comb & sample!="UNKNOWN_comb")%>%
  ggplot(aes(x=generation,y=lineages_represented,col=sample)) +
  geom_line(lwd=1.5) +
  scale_x_continuous(breaks = seq(0,14,1)) +
  theme_classic() +
  ggtitle("Proportion of blood lineages found in tissue through generations")

ggsave(p1,filename = "Figures/8pcw/Fraction_of_lineages_represented_8pcw.pdf",device = "pdf",width=7,height=5)

#---------------------------------
### PLOT A SCALE BAR FOR THE TREES###
#---------------------------------
#This is to plot a colour scale bar to visualize how the colour scales on the log_cell_frac tree corresponds to the actual cell fraction

#This bit rescales the scale bar using the same transformations applied to the vaf data
scale_bar=rev(1/exp(seq(0,8))) #These are the tick points that will go on the scale
log_scale_bar=log(scale_bar)
prob_threshold_to_include=0.9 #This should be same as for the within the "generate_targeted_seq_plots" function

#This bit takes the elements of the plotrix rescale function
scale_range=c(0.01,1)

#Get the vaf_present & log_vaf_present objects (performed as within the generate_targ_seq_plots function)
cell_frac=calculate_cell_frac(matrices$NV,matrices$NR)
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
pdf("Figures/8pcw/Trees_colour_scale_bar_8pcw.pdf")
plot(c(0,10), c(0,1), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main="")
axis(side=2,pos=1,las=1,at=tick_positions_trimmed,labels=as.character(tick_values_trimmed))
for (i in 1:(length(lut)-1)) {
  y = (i-1)/scale
  rect(1,y,2,y+1/scale, col=lut[i], border=NA)
}
dev.off()
