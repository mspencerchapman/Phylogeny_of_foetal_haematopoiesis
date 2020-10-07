#This script is designed to be interactive & quick to allow rapid exploring of filtering parameters
library(stringr)
library(ape)
library(seqinr)
library(ggtree)
library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)

# my_working_directory="~/R Work/Fetal HSPCs/Phylogeny_of_foetal_haematopoiesis/"
# treemut_dir="~/R Work/R_scripts/treemut/"
# 
my_working_directory="/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/Phylogeny_of_foetal_haematopoiesis"
treemut_dir="/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut"

setwd(my_working_directory)

R_function_files = list.files("R_functions",pattern=".R",full.names=TRUE)
sapply(R_function_files,source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)

#Set the file paths for saved files based on these IDs.
mats_and_params_file ="/lustre/scratch119/realdata/mdt1/team154/ms56/fetal_HSC/filtering_runs/mats_and_params/mats_and_params_a1_MS2_both_reduced" #This file is too large to go on github, and must be re-derived from the raw data if needed
filtered_muts_file = "Data/18pcw/Filtered_mut_set_18pcw"
dna_string_file = "Data/18pcw/DNA_string_file_18pcw.fa"
mpboot_tree_file = paste0(dna_string_file,".treefile")
tree_file_path = "Data/18pcw/Tree_18pcw.tree"
vcf_header_path = "Data/vcfHeader.txt"
vcf_path = "Data/18pcw/Filtered_mut_set_vcf_18pcw.vcf"
shared_vcf_path="Data/18pcw/Mutations_shared_18pcw.vcf"
vagrent_input_path = "Data/18pcw/Filtered_mut_set_vcf_with_header_18pcw.vcf"
vagrent_output_path = paste0(vagrent_input_path,".annot")
file_annot="Data/18pcw/Filtered_mut_set_annotated_18pcw"
sensitivity_analysis_path <- "Data/18pcw/sensitivity_analysis_18pcw"  

#If have previous run analysis - load up files
previously_run = T
if(previously_run) {
  load(file_annot); tree <- read.tree(tree_file_path)
} else {
  #Otherwise, load up the "list"mats_and_params" file and the file of "early somatic mutations"
  load(mats_and_params_file)
  load("Data/18pcw/early_somatic_muts")
  
  #Can specify coverage cut-off here, to exclude low coverage samples from analysis - samples too difficult to even place on tree
  min_sample_mean_cov = 4
  other_samples_to_remove = "PD41768id_hum" #must be the sample ID without the "_MTR" or "_DEP" suffix
  COMB_mats$gender="male"
  
  #Remove the low coverage samples and their private mutations
  if(min_sample_mean_cov > 0) {
    output = remove_low_coverage_samples(COMB_mats = COMB_mats,
                                         filter_params = filter_params,
                                         min_sample_mean_cov = min_sample_mean_cov,
                                         other_samples_to_remove = other_samples_to_remove,
                                         min_variant_reads_auto = 3, #these parameters are to remove mutations from the matrix that are no longer positive in any samples
                                         min_variant_reads_xy = 2)
    COMB_mats= output$COMB_mats
    filter_params = output$filter_params
  }
  
  
  
  #REVIEW MEAN DEPTH HISTOGRAMS TO DECIDE MANUAL CUT-OFFS FOR EXCLUDING OUTLIERS
  #hist(filter_params$mean_depth, breaks = 100, xlim = c(0,60))
  #hist(filter_params$mean_depth[COMB_mats$mat$Chrom %in% c("X","Y")], breaks = 200, xlim = c(0,30))
  #hist(filter_params$mean_depth[!COMB_mats$mat$Chrom %in% c("X","Y")], breaks = 200, xlim = c(0,60))
  XY_low_depth_cutoff = 4; XY_high_depth_cutoff = 9; AUTO_low_depth_cutoff = 8; AUTO_high_depth_cutoff = 18
  
  #Get the filtered mutation set - returns object with the filtered mat/ NV/NR matrices, as well as a full genotype matrix, and matrix of shared mutations only
  filtered_muts = get_filtered_mut_set(input_set_ID = "18pcw",  #the Run_ID of the unfiltered mutation set used as input - though won't account for any removal of samples from set
                                       COMB_mats = COMB_mats,  #the main full mutation matrix as outputed by the "HSC_filtering_treebuild_table.R" script
                                       filter_params = filter_params,  #the filter_params matrix as outputed by the "HSC_filtering_treebuild_table.R" script
                                       gender = COMB_mats$gender, #patient's gender
                                       
                                       #These parameters decide whether a mutation is retained in the "true somatic mutation" set
                                       retain_muts = early_somatic_muts,  #any mutations that should be manually retained, despite not meeting filtering criteria, NULL by default
                                       germline_pval = -10,  #the log10 p-value cutoff for likelihood that mutations came from germline distribution. Mutations with value < cutoff are retained.
                                       rho = 0.1,  #rho cutoff for the beta-binomial filter, a measure of how "over-dispersed" the counts are compared to a binomial distribution. Mutations with value > cutoff are retained.
                                       mean_depth = c(AUTO_low_depth_cutoff,AUTO_high_depth_cutoff, XY_low_depth_cutoff, XY_high_depth_cutoff),   #Numeric vector of length 4 defining mean depth at mutation site cut-offs. This is in the order (1) lower threshold for autosomes, (2) upper threshold for autosomes, (3) lower threshold for XY, (4) upper threshold for XY. This removes mis-mapping/ low reliability loci.
                                       pval_dp2=NA,  #the p-value cut-off if using the "pval within pos" filter, with positive samples defined as having >= 2 reads. Mutations with value > cutoff are retained.
                                       pval_dp3=0.001,   #the p-value cut-off if using the "pval within pos" filter, with positive samples defined as having >= 3 reads (allows for more index hopping). Mutations with value > cutoff are retained.
                                       min_depth = c(6,4), #Numeric vector of length 2 defining minimum depths that at least one positive sample must have for mutation to be retained (AUTO and XY)
                                       min_pval_for_true_somatic = 0.1,   #Default: 0.1. the minimum p-value that at least one sample must have for the variant:normal read distribution coming from that expected for a true somatic. Mutations with value > cutoff are retained.
                                       min_vaf = NA, #Numeric vector of length 2 defining minimum vaf in at least one sample for mutation to be retained (AUTO and XY). Mutations with value > cutoff are retained.
                                       
                                       #These parameters decide the genotype for each sample for each "true somatic mutation".  These may be less stringent than the initial parameters.
                                       min_variant_reads_SHARED = 2,  #the minimum number of reads for samples to be assigned a positive genotype
                                       min_pval_for_true_somatic_SHARED = 0.05,  #the p-value for coming from "true somatic mutation" read distribution to be assigned a positive genotype
                                       min_vaf_SHARED = NA) #Numeric vector of length 2, defining minimum vaf to be assigned a positive genotype
  
  
  #Decide an ID for this filtered set, depending on approach taken, and save
  save(filtered_muts, file = filtered_muts_file)
  
  #Write a fasta file of the dummy DNA strings
  write.fasta(filtered_muts$dna_strings, names=names(filtered_muts$dna_strings), dna_string_file)
  
  #BUILD TREE with MPBoot
  system(paste0("/lustre/scratch117/casm/team154/tc16/Programs/mpboot-sse-1.1.0-Linux/bin/mpboot -s ", dna_string_file," -bb 1000"))
  
  #Import the tree into R using ape
  require(ape)
  tree_file = paste0(dna_string_file,".treefile")
  tree <- read.tree(tree_file)
  tree <- drop.tip(tree,"Ancestral")
  tree$edge.length = rep(1, nrow(tree$edge)) #Initially need to assign edge lengths of 1 for the tree_muts package to work
  
  #Assign mutations back to the tree
  df = reconstruct_genotype_summary(tree) #Define df (data frame) for treeshape
  
  #Get matrices in order, and run the main assignment functions
  mtr = as.matrix(filtered_muts$COMB_mats.tree.build$NV)
  depth = as.matrix(filtered_muts$COMB_mats.tree.build$NR)
  p.error = c(rep(0.01, ncol(filtered_muts$COMB_mats.tree.build$NR)))
  res = assign_to_tree(mtr[,df$samples], depth[,df$samples], df, error_rate = p.error) #Get res (results!) object
  
  #Assign edge lengths from the res object
  tree$edge.length <- res$df$df$edge_length
  filtered_muts$COMB_mats.tree.build$mat$node <- tree$edge[res$summary$edge_ml,2] #Add node information to the filtered_muts object
  filtered_muts$COMB_mats.tree.build$mat$pval <- res$summary$pval #Add pval information to the filtered_muts object
  
  treefit_pval_cutoff = 1e-3
  poor_fit = res$summary$pval < treefit_pval_cutoff; print(sum(poor_fit))  #See how many mutations don't have read counts that fit the tree very well
  
  #Remove poor_fit muts if they clearly do not fit the phylogeny
  filter_poor_fit = "yes"
  if(filter_poor_fit == "yes") {
    filtered_muts$COMB_mats.tree.build <- list_subset(filtered_muts$COMB_mats.tree.build,select_vector = !poor_fit)
    tree$edge.length<-sapply(tree$edge[,2], function(node) sum(filtered_muts$COMB_mats.tree.build$mat$node==node))
  }
  
  #Now make multi-furcating version of the tree - used for most downstream analysis
  tree<-di2multi(tree)
  #Assign mutations back to the tree
  df = reconstruct_genotype_summary(tree) #Define df (data frame) for multi treeshape
  mtr = as.matrix(filtered_muts$COMB_mats.tree.build$NV); depth = as.matrix(filtered_muts$COMB_mats.tree.build$NR)
  p.error = c(rep(0.01, ncol(filtered_muts$COMB_mats.tree.build$NR)))
  res = assign_to_tree(mtr[,df$samples], depth[,df$samples], df, error_rate = p.error) #Get res (results!) object
  
  #Assign edge lengths from the res object
  tree$edge.length <- res$df$df$edge_length
  filtered_muts$COMB_mats.tree.build$mat$node <- tree$edge[res$summary$edge_ml,2] #Add node information to the filtered_muts object
  filtered_muts$COMB_mats.tree.build$mat$pval <- res$summary$pval #Add pval information to the filtered_muts object
  
  #Save the multi-furcating tree file
  write.tree(tree, file = tree_file_path)
  
  #ANNOTATING THE MUTATIONS
  #-----------------------------------------------
  #Write vcf files for VariantCaller analysis - separate out "All mutations" & Shared mutations"
  vcf_file = create_vcf_files(filtered_muts$COMB_mats.tree.build$mat)
  write.table(vcf_file, sep = "\t", quote = FALSE, file = vcf_path, row.names = FALSE)
  
  #Also save the shared mutations for mutational signature analysis
  shared_vcf_file = create_vcf_files(filtered_muts$COMB_mats.tree.build$mat,select_vector = which(!filtered_muts$COMB_mats.tree.build$mat$node%in%1:length(tree$tip.label)))
  write.table(shared_vcf_file, sep = "\t", quote = FALSE, file = shared_vcf_path, row.names = FALSE)
  
  #1. paste vcf file to a dummy header file
  system(paste0("cat ",vcf_header_path," ",vcf_path," > ", vagrent_input_path))
  
  #2. commands to run vagrent
  system(paste0("AnnotateVcf.pl -i ",vagrent_input_path," -o ",vagrent_output_path," -sp Human -as NCBI37 -c /lustre/scratch117/casm/team78pipelines/reference/human/GRCh37d5/vagrent/e75/vagrent.cache.gz"))
  
  #3. import vagrent output
  vagrent_output = fread(vagrent_output_path,skip = "#CHROM")
  annot_info = as.data.frame(str_split(vagrent_output$INFO, pattern = ";",simplify = TRUE), stringsAsFactors = FALSE)
  colnames(annot_info) <- c("VT","VD","VC","VW")
  
  annot_info$VC <- gsub(x=annot_info$VC, pattern = "VC=", replacement = "")
  annot_info$VT <- gsub(x=annot_info$VT, pattern = "VT=", replacement = "")
  annot_info$VW <- gsub(x=annot_info$VW, pattern = "VW=", replacement = "")
  annot_info$VD <- gsub(x=annot_info$VD, pattern = "VD=", replacement = "")
  
  filtered_muts$COMB_mats.tree.build$mat <- cbind(filtered_muts$COMB_mats.tree.build$mat,split_vagrent_output(df = annot_info,split_col = "VD"))
  filtered_muts$COMB_mats.tree.build$mat$variant_ID <- paste(filtered_muts$COMB_mats.tree.build$mat$Gene, filtered_muts$COMB_mats.tree.build$mat$Protein, sep = " ") #Add a "variant_ID" that includes that gene & protein coding change
  filtered_muts$COMB_mats.tree.build$mat$Type[filtered_muts$COMB_mats.tree.build$mat$Type == ""] <- "no_annotation" #Label blank field as "no_annotation"
  
  
  #-----------------------------------------------
  #Save the annotated filtered_muts files (post tree filtering)
  save(filtered_muts, file = file_annot)
}

#-----------------------------------------------
#CORRECT THE EDGE LENGTHS BASED ON SAMPLE SENSITIVITY
#NB - for function to work correctly, the sensitivity_df must be set up in exactly the correct format
#If there are no SNVs in the tree, and/or sensitivity analysis for these, set "include_indels" to FALSE

#Create trees of only SNVs or only INDELs (this function does not do any correction)
tree_SNV = get_subset_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, v.field="Mut_type", value = "SNV")
tree_INDEL = get_subset_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, v.field="Mut_type", value = "INDEL")

#Store the uncorrected SNV mutation burden
SNV_burden_u=get_mut_burden(tree_SNV)

##PERFORM REGRESSION TO GET COVERAGE-SPECIFIC IN VITRO MUTS CORRECTION
#Create sensitivity-corrected SNV tree, this is necessary because for the 18pcw, there are a number
#of samples with low coverage (<10X). At this level it becomes harder to distinguish between true
#mutations and in vitro mutations using the read counts, so more in vitro mutations are retained (as seen by
#the negative correlation between private mutation burden and read depth). This code removes this
#correlation while maintaining the same overall proportional reduction in private mutations corresponding
#to the binomial mixture model.

#Import the sensitivity analysis file
sensitivity_df <- read.delim(sensitivity_analysis_path, header = TRUE, stringsAsFactors = FALSE, sep = " ")

#Correct the tree using this df
tree_SNV_c = get_corrected_tree(tree = tree_SNV, details = filtered_muts$COMB_mats.tree.build$mat, include_indels = FALSE, sensitivity_df = sensitivity_df,get_edge_from_tree=TRUE)

##LINEAR REGRESSION MODEL FOR CORRECTION
#This needs the table of sequencing summary statistics used for the regression
smry_seq_18pcw<-read.csv("Data/18pcw/smry_seq_18pcw.csv")
#(1) Get the corrected private branch lengths for each sample
private_mut_burden=tree_SNV_c$edge.length[tree_SNV_c$edge[,2]%in%1:length(tree_SNV_c$tip.label)]
#(2) Get a matched vector of sequencing coverage
sample_coverage=sapply(gsub("_hum","",tree_SNV_c$tip.label),function(x) smry_seq_18pcw$coverage[smry_seq_18pcw$Donor_ID==x])
#(3) Visualize the correlation using the log of the private mutation burden
plot(log(private_mut_burden)~sample_coverage)
#(4) Linear model of this correlation and add to the plot
lin_mod <-lm(log(private_mut_burden)~sample_coverage)
abline(lin_mod)
#(5) Use the correlation coefficient with coverage to correct the log(private_mut_burden)
#The 19.75 is chosen such that the total proportion of in vitro mutations is as per the binomial mixture model
cut_off=19.75
log_corrected_muts=log(private_mut_burden)+lin_mod$coefficients["sample_coverage"]*(cut_off-sample_coverage)
plot(log_corrected_muts~sample_coverage)
abline(lm(log_corrected_muts~sample_coverage))
#(6) Apply the correction such that samples above the threshold remain at the same, they don't gain true muts from the correction
corrected_muts=sapply(1:length(sample_coverage), function(i) {if(sample_coverage[i]<cut_off) {return(exp(log_corrected_muts[i]))} else {return(private_mut_burden[i])}})
#(7) Review that correction has worked as expected by looking for any correlation with corrected muts
plot(corrected_muts~sample_coverage); abline(lm(corrected_muts~sample_coverage))
#(8) This is the check that the proportion of mutations retained matches that from the binomial mixture model
sum(corrected_muts)/sum(private_mut_burden)

#Create corrected SNV tree
tree_SNV_c$edge.length[tree_SNV_c$edge[,2]%in%1:length(tree$tip.label)] <- corrected_muts

tree_c = get_corrected_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, include_indels = TRUE, sensitivity_df = sensitivity_df)
tree_INDEL_c = get_corrected_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, include_SNVs = FALSE, sensitivity_df = sensitivity_df)

#Create a "half-corrected" hybrid tree: Corrected SNVs, uncorrected INDELs (as difficult to do accurately & doesn't make a big difference)
tree_hybrid=tree_SNV_c #Start from corrected SNV tree
tree_hybrid$edge.length <- tree_SNV_c$edge.length + tree_INDEL$edge.length #Combine the edge lengths of the corrected SNV tree & the uncorrected indel tree
pdf("Figures/18pcw/Haematopoietic_phylogeny_18pcw.pdf",width=15,height=7)
tree_hybrid=plot_tree(tree_hybrid,cex.label = 0,lwd=1,default_edge_color = "black")
dev.off()


#------------------------------------------------------------------
#PLOT THE CELL MIGRATION TREE TO SHOW LACK OF CLUSTERING BY LOCATION
#------------------------------------------------------------------
details=filtered_muts$COMB_mats.tree.build$mat

#Label tree by the tissue that they were isolated from (liver, femur 1 or femur 2)
tree_SNV_c$tip.label=gsub("_hum","",tree_SNV_c$tip.label)
tree_SNV_c$tip.label=as.character(sapply(tree_SNV_c$tip.label, function(sample) smry_seq_18pcw$Tissue[which(smry_seq_18pcw$Donor_ID==sample)]))

#Make the tree ultra-metric for visualization
tree_SNV_c_ultra <- make.ultrametric.tree(tree_SNV_c)
tree_SNV_c_ultra$edge.length=42*tree_SNV_c_ultra$edge.length #Multiply by the average mutation burden for visualization

#To simply say if branch is found in colonies from only one or multiple organs.. (therefore found in each femur is also black)
#1.relabel tips to only be the organ in which they were found
pdf("Figures/18pcw/Haematopoietic_phylogeny_tissue_origin_18pcw.pdf",width=7,height=4)
tree_SNV_c_ultra=plot_tree(tree_SNV_c_ultra, cex.label = 0,plot_axis = FALSE)
add_annotation(tree_SNV_c_ultra,
               details,
               matrices,
               annot_function=plot_sharing_multiple,
               sharing_cols=c("black","#17698E","#17A258","brown"),
               lwd=0.8)
dev.off()

#---------------------------------------------------------------------------------------------------
#PLOT WITH ANNOTATION OF MUTATIONS CAUSING CODING CHANGES IN SHARED BRANCHES (the most robust calls)
#---------------------------------------------------------------------------------------------------
#Define the coding changes in shared branches
details$shared_coding_change <- ifelse(!details$node %in% 1:length(tree$tip.label) & details$Type %in% c("protein_coding:exon:CDS:substitution:codon_variant:non_synonymous_codon",
                                                                                                               "protein_coding:exon:CDS:insertion:frameshift_variant",
                                                                                                               "protein_coding:exon:CDS:deletion:frameshift_variant",
                                                                                                               "protein_coding:exon:CDS:substitution:codon_variant:stop_gained",
                                                                                                               "protein_coding:exon:CDS:substitution:codon_variant:initiator_codon_change",
                                                                                                               "protein_coding:exon:CDS:deletion:inframe_variant:inframe_codon_loss"),
                                             "Shared coding change", "no")

pdf("Figures/18pcw/Haematopoietic_phylogeny_shared_CDS_muts_18pcw.pdf",width=15,height=7)
tree_hybrid=plot_tree(tree_hybrid, cex.label = 0)
plot_tree_labels(tree_hybrid,
                 details = details,
                 type = "line",
                 query.field = "shared_coding_change",
                 data.frame(value="Shared coding change",col="red",pch = 17,stringsAsFactors = FALSE),
                 label.field = "variant_ID",
                 cex.label = 0.8)
dev.off()