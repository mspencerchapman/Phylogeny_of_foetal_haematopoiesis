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

my_working_directory="~/R Work/Fetal HSPCs/Foetal_phylogeny/"
treemut_dir="~/R Work/R_scripts/treemut/"
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

##IMPORTING THE LCM CGPVAF TABLES (snps and indels)
#Import the cgpvaf files and combine the LCM and single-cell colony (SCC) matrices into single matrix
LCM_mats<-import_cgpvaf_SNV_and_INDEL(SNV_output_file=cgpvaf_LCM_SNV_file,INDEL_output_file=cgpvaf_LCM_INDEL_file)
SCC_mats<-import_cgpvaf_SNV_and_INDEL(SNV_output_file=cgpvaf_colonies_SNV_file,INDEL_output_file=cgpvaf_colonies_INDEL_file)
NV = cbind(LCM_mats$NV,SCC_mats$NV); NR = cbind(LCM_mats$NR,SCC_mats$NR)

#Import the LCM sample info
lcm_smry = read.csv(lcm_smry_path,stringsAsFactors = FALSE)

#Add aggregate columns for each tissue type (as defined in the lcm_smry table)
tissues = names(table(lcm_smry$Tissue))
for(i in tissues) {
  NV[,i] <- apply(NV,1,function(x) {sum(x[colnames(NV) %in% lcm_smry$Sample_ID[lcm_smry$Tissue == i]])})
  NR[,i] <- apply(NR,1,function(x) {sum(x[colnames(NR) %in% lcm_smry$Sample_ID[lcm_smry$Tissue == i]])})
}

#Add aggregate columns for germ layers (as defined in the lcm_smry table)
germ_layers=names(table(lcm_smry$Germ_layer))
for(i in germ_layers) {
  NV[,i] <- apply(NV,1,function(x) {sum(x[colnames(NV) %in% lcm_smry$Sample_ID[lcm_smry$Germ_layer == i]])})
  NR[,i] <- apply(NR,1,function(x) {sum(x[colnames(NR) %in% lcm_smry$Sample_ID[lcm_smry$Germ_layer == i]])})
}

#Add aggregate columns for mesoderm vs non-mesoderm (as defined in the lcm_smry table)
same_germ_layer=names(table(lcm_smry$Gastrulation))
for(i in same_germ_layer[-3]) { #Don't do for the "unknown" again
  NV[,i] <- apply(NV,1,function(x) {sum(x[colnames(NV) %in% lcm_smry$Sample_ID[lcm_smry$Gastrulation == i]])})
  NR[,i] <- apply(NR,1,function(x) {sum(x[colnames(NR) %in% lcm_smry$Sample_ID[lcm_smry$Gastrulation == i]])})
}

#Add a column to aggregate all LCM samples
NV[,"ALL_LCM"] <- apply(NV,1,function(x) {sum(x[colnames(NV) %in% lcm_smry$Sample_ID])})
NR[,"ALL_LCM"] <- apply(NR,1,function(x) {sum(x[colnames(NR) %in% lcm_smry$Sample_ID])})

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


##SCRIPT INCLUDING FUNCTIONS TO CALCULATE LINEAGE LOSS AT A GIVEN NODE IN A GIVEN TISSUE
##AND TO WORK OUT LINEAGES CAPTURED THROUGH THE FIRST FOUR ROUNDS OF CELL DIVISION

#Start with working from root on generations tree
NV=NV[details_targ$mut_ref,]
NR=NR[details_targ$mut_ref,]

par(mfrow=c(1,1))
tree.multi=plot_tree(tree.multi,cex.label=0)
add_annotation(tree=tree.multi,
               details=details_targ,
               matrices=list(mtr=NV,dep=NR),
               annot_function=plot_node_number)

#---------------------------------------------------
#PLOTTING LINEAGE LOSS DONUT CHARTS AT THE KEY NODES
#---------------------------------------------------

#Plot the tree for the tissue of interest
pdf("~/Desktop/8wk_germ_layers_tree.clean.vaf.pdf",width=10,height=7)
samples=germ_layers_comb
generate_targ_seq_plots(samples,
                        tree=tree.multi,
                        details_targ=details_targ,
                        matrices=list(NV=NV,NR=NR),
                        post_prob_type="clean", #other option is "raw"
                        info_type="log_vaf", #other option is "post.prob" or "vaf"
                        prob_threshold_to_include=0.5,
                        plot_cell_frac=FALSE,
                        plot_donut=FALSE, #plot a donut for positive nodes
                        donut_info="cell_frac", #other option is "lineages_lost"
                        CI=0.8,
                        radius=4)
dev.off()


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

#Also, lots of mutations on the gen_4_branches are from later generations
#These lines of code find which ones are from later branches and therefore not genuinely from gen_4_branches
#The ALL_LCM counts give the best power to pick these out in one tissue
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

summary_df[,2:5]<-apply(summary_df[,2:5],2,as.numeric)
str(summary_df)

#Visualization options
extra_embryonic=c("PERIPHERAL_BLOOD","TROPH","MESENCHYME")
to_visualize=lcm_smry$Sample_ID[lcm_smry$Tissue=="PERIPHERAL_BLOOD"]
to_visualize=tissues
to_visualize=c("MESODERM","ENDODERM","ECTODERM")
to_visualize=germ_layers[!germ_layers=="UNKNOWN"]
to_visualize=c("TROPH","MESENCHYME","PERIPHERAL_BLOOD","POST_GASTRULATION")

#Use ggplot2 to create "lineage loss over generations"
summary_df$tissue=factor(summary_df$tissue,levels=c("POST_GASTRULATION","MESODERM",""))
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


p_comb<-arrangeGrob(p1,p2,ncol=2)
plot(p_comb)
