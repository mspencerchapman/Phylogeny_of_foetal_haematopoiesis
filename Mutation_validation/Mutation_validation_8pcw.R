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

my_working_directory="~/R Work/Fetal HSPCs/Phylogeny_of_foetal_haematopoiesis/"
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
SCC_mats<-import_cgpvaf_SNV_and_INDEL(SNV_output_file=cgpvaf_colonies_SNV_file,INDEL_output_file=cgpvaf_colonies_INDEL_file)

#Import the details matrix and tree (needed to interpret the targeted sequencing data)
load(file_annot); tree<- read.tree(tree_file_path); tree$tip.label=gsub("_hum","",tree$tip.label)
details <- filtered_muts$COMB_mats.tree.build$mat

#In details file, create column for whether it is included in the targeted sequencing baits
details$targeted_tree = ifelse(details$mut_ref %in% rownames(SCC_mats$NV),"YES","NO")

#Now subset the tree and the details matrix to only (1) SNVs that are (2) included in the bait set
tree_targ = get_subset_tree(tree = tree, details = details, v.field = "targeted_tree",value = "YES")
details_targ = details[details$targeted_tree=="YES",]

#REVIEW MUTATION SETS FOR THE BINOMIAL MIXTURE MODEL
#1.Define vectors of mutation references ("mut_refs") for private vs shared, SNVs vs INDELs, autosomal vs xy, as these need to be analysed separately

private_SNVs=details_targ$mut_ref[details_targ$Mut_type=="SNV" & details_targ$node %in% 1:length(tree_targ$tip.label)]
auto_private_SNVs = private_SNVs[!grepl("X",private_SNVs)&!grepl("Y",private_SNVs)]
xy_private_SNVs = private_SNVs[grepl("X",private_SNVs)|grepl("Y",private_SNVs)]

shared_SNVs=details_targ$mut_ref[details_targ$Mut_type=="SNV" & !details_targ$node %in% 1:length(tree_targ$tip.label)]
auto_shared_SNVs = shared_SNVs[!grepl("X",shared_SNVs)&!grepl("Y",shared_SNVs)]
xy_shared_SNVs = shared_SNVs[grepl("X",shared_SNVs)|grepl("Y",shared_SNVs)]

#2.Extract the validation result counts, and reformat into a dataframe of NV, NR and mut_ref
#Using the "counts_only" argument, the validate_mutation function just gives a total read/ depth count for a given mutation aggregate across all colonies expected to have that mutation
validation_results_counts=lapply(details_targ$mut_ref,validate_mutation,tree=tree_targ,details=details_targ,NV=SCC_mats$NV,NR=SCC_mats$NR,counts_only=TRUE)
names(validation_results_counts)=details_targ$mut_ref
mut_counts_df=Reduce(rbind,mapply(x=validation_results_counts,names=names(validation_results_counts),FUN = function(x,names) {data.frame(NV=x[1],NR=x[2],mut_ref=names)},SIMPLIFY = FALSE))

#Define a function to plot binomial mixture model output filtering the mutation set in various ways
binomial_mix_plot=function(mut_counts_df,mut_ref_set,depth_cutoff,nrange=1:4,title=NULL,return_res=FALSE,...) {
  NV=mut_counts_df$NV[mut_counts_df$mut_ref%in%mut_ref_set & mut_counts_df$NR>=depth_cutoff]
  NR=mut_counts_df$NR[mut_counts_df$mut_ref%in%mut_ref_set & mut_counts_df$NR>=depth_cutoff]
  cols=c("blue","red","green","orange","magenta")
  hist(NV/NR,breaks=seq(0,1,0.05),xlim=c(0,1),col='gray',freq=F,xlab="VAF",main=title,cex.main=0.7,cex.lab=0.7,cex.axis=0.7,...)
  lines(density(NV/NR),lwd=2,lty='dashed')
  res = binom_mix(NV,NR,nrange=nrange)
  for (i in 1:res$n){
    meancov = round(mean(NR))
    lines(x=(0:meancov)/meancov,
          y=meancov*res$prop[i]*dbinom(0:meancov,meancov,prob=res$p[i]),
          type="l",col=cols[i],lwd=2)
  }
  res$mut_ref=as.character(mut_counts_df$mut_ref[mut_counts_df$mut_ref%in%mut_ref_set & mut_counts_df$NR>=depth_cutoff])
  df=data.frame(binom.peak = res$p[order(res$p)], sample_fraction=res$prop[order(res$p)])
  print(df)
  if(return_res){return(res)}
}

#Apply to the autosomal mutations, using either 8x or 40x cutoffs, and reviewing private or shared
pdf("Figures/8pcw/Mutation_validation_binomial_mix_model_8pcw.pdf",width=7,height=7)
par(mfrow=c(2,2))
binomial_mix_plot(mut_counts_df=mut_counts_df,mut_ref_set = auto_private_SNVs,depth_cutoff = 8,title = "8pcw: Private autosomal SNVs (depth>=8)") #
binomial_mix_plot(mut_counts_df=mut_counts_df,mut_ref_set = auto_private_SNVs,depth_cutoff = 40,title = "8pcw: Private autosomal SNVs (depth>=40)")
binomial_mix_plot(mut_counts_df=mut_counts_df,mut_ref_set = auto_shared_SNVs,depth_cutoff = 8,title = "8pcw: Shared autosomal SNVs (depth>=8)",ylim=c(0,7))
binomial_mix_plot(mut_counts_df=mut_counts_df,mut_ref_set = auto_shared_SNVs,depth_cutoff = 40, title = "8pcw: Shared autosomal SNVs (depth>=40)",ylim=c(0,11))
dev.off()

#For the autosomal private mutations, pull out the subclonal and clonal clusters of mutations for mutational signature analysis
res=binomial_mix_plot(mut_counts_df=mut_counts_df,mut_ref_set = auto_private_SNVs,depth_cutoff = 8,return_res = TRUE)
clonal_clust=which.min(abs(res$p-0.5)) #This will pull out the peak closest to 0.5
subclonal_muts=res$mut_ref[res$Which_cluster!=clonal_clust]
clonal_muts=res$mut_ref[res$Which_cluster==clonal_clust]

#Create and save vcf files of these mutations
subclonal_vcf_file = create_vcf_files(filtered_muts$COMB_mats.tree.build$mat, select_vector = filtered_muts$COMB_mats.tree.build$mat$mut_ref %in% subclonal_muts)
clonal_vcf_file = create_vcf_files(filtered_muts$COMB_mats.tree.build$mat, select_vector = filtered_muts$COMB_mats.tree.build$mat$mut_ref %in% clonal_muts)

write.table(subclonal_vcf_file, sep = "\t", quote = FALSE, file = "Data/8pcw/subclonal_private_muts_8pcw.vcf", row.names = FALSE)
write.table(clonal_vcf_file, sep = "\t", quote = FALSE, file = "Data/8pcw/clonal_private_muts_8pcw.vcf", row.names = FALSE)

#It is surprising that some shared mutations may be called as "subclonal",as the fact that it has been called
#in more than one colony (and that the sharing fits the phylogeny) makes it highly unlikely to be in vitro acquired
#Therefore we can review these manually
res=binomial_mix_plot(mut_counts_df=mut_counts_df,mut_ref_set = auto_shared_SNVs,depth_cutoff = 8,title = "8pcw: Shared autosomal SNVs (depth>=8)",ylim=c(0,10),return_res = TRUE)
clonal_clust=which.min(abs(res$p-0.5)) #This will pull out the peak closest to 0.5
subclonal_muts=res$mut_ref[res$Which_cluster!=clonal_clust]

#This visualization allows you to vizualize the overall vafs and read counts for individual branches
pdf("Figures/8pcw/Shared_mutations_called_as_subclonal_8pcw.pdf",width=15,height=5)
par(mfrow=c(1,1))
for(i in 1:length(subclonal_muts)) {
  tree_targ=plot_tree(tree=tree_targ,cex.label=0.2)
  add_annotation(tree,
                 details=details_targ,
                 matrices=list(NV=SCC_mats$NV, NR=SCC_mats$NR),
                 annot_function = plot_mut_vaf_by_branch,
                 mut=subclonal_muts[i],
                 cex=0.4)
  add_annotation(tree,
                 details=details,
                 matrices=list(NV=NV, NR=NR),
                 annot_function = highlight_nodes, #Highlights the allocated node for the mutation with dashed red
                 nodes=details$node[details$mut_ref==subclonal_muts[i]],
                 lty=2)
}
dev.off()

#-----------------------------------------------------------------------
#INDIVIDUAL COLONY TARGETED SEQUENCING RESULTS OVERLAID ON THE PHYLOGENY
#As a visual check on the accuracy of the phylogeny from the WGS, we overlaid the results of the targeted
#re-sequencing on the original phylogeny. This would quickly highlight any inaccuracies.

#Use a raw vaf colour scale
details_targ_full <- cbind(details_targ,calculate_vaf(SCC_mats$NV[details_targ$mut_ref,],SCC_mats$NR[details_targ$mut_ref,]))

#Define samples that were previously placed on phylogeny ("colonies_to_validate") and those that were not ("new_samps")
colonies = colnames(SCC_mats$NV)
new_samps = colonies[!colonies %in% tree$tip.label]
colonies_to_validate = tree$tip.label[tree$tip.label%in%colonies]

#Save validation plots for all colonies that were including in the original WGS phylogeny
pdf("Figures/8pcw/Colony_validation_trees_8pcw.pdf")
samples=colonies_to_validate
sapply(samples, function(SAMPLE) {
  tree_targ=plot_tree(tree_targ, cex.label = 0)
  text(50,1.5,SAMPLE)
  text(50,0,paste0("Mean depth is ",round(mean(NR[,SAMPLE]),digits = 2)))
  add_annotation(tree=tree_targ,
                 details=details_targ_full,
                 list(mtr=SCC_mats$NV,dep=SCC_mats$NR),
                 annot_function=function(tree,details,matrices,node) {
                   add_var_col(tree,details,matrices,node,var_field = SAMPLE,pval_based=FALSE,lwd = 2.5)
                 }
  )
  plot_sample_tip_point(tree=tree_targ,details=details_targ_full,sample=SAMPLE)
}
)
dev.off()

#Save validation plots for poly-clonal colonies not previously on phylogeny
pdf("Figures/8pcw/New_polyclonal_colonies_phylogeny_placement_8pcw.pdf")
samples=new_samps
sapply(samples, function(SAMPLE) {
  tree_targ=plot_tree(tree_targ, cex.label = 0)
  text(50,1.5,SAMPLE)
  text(50,0,paste0("Mean depth is ",round(mean(SCC_mats$NR[,SAMPLE]),digits = 2)))
  add_annotation(tree=tree_targ,
                 details=details_targ_full,
                 list(mtr=SCC_mats$NV,dep=SCC_mats$NR),
                 annot_function=function(tree,details,matrices,node) {
                   add_var_col(tree,details,matrices,node,var_field = SAMPLE,pval_based=FALSE,lwd = 2.5)
                 }
  )
})
dev.off()

#validate mutation by mutation, aggregating reads across samples expected to have the mutation
validation_results=sapply(details_targ_full$mut_ref,validate_mutation,tree=tree,details=details_targ,NV=SCC_mats$NV,NR=SCC_mats$NR,pval_cutoff=0.045,vaf_cutoff=0.3,depth_cutoff=15)
details_targ_full$validation_results=validation_results
table(validation_results)

#Validation tree plot - annotate the tree with mutations that have been "successfully validated"
pdf("Figures/8pcw/validated_muts_tree_8pcw.pdf",width=20,height=6)
tree_targ=plot_tree(tree_targ, cex.label = 0)
add_annotation(tree=tree_targ,
               details=details_targ_full,
               list(mtr=SCC_mats$NV,dep=SCC_mats$NR),
               annot_function=function(tree,details,matrices,node) {
                 add_categorical_col(tree,
                                     details,
                                     matrices,
                                     node,
                                     var_field = "validation_results",
                                     annot=list(PASS="green",FAIL="red",inadequate_depth="grey"),
                                     lwd = 2.5)
               }
)
dev.off()


