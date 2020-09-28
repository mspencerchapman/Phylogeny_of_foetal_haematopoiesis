#This script is designed to be interactive & quick to allow rapid exploring of filtering parameters
library(stringr)
library(ape)
library(seqinr)
library(tidyr)
library(dplyr)
library(ggplot2)

my_working_directory="~/R Work/Fetal HSPCs/Phylogeny_of_foetal_haematopoiesis/"
treemut_dir="~/R Work/R_scripts/treemut/"
setwd(my_working_directory)

R_function_files = list.files("R_functions",pattern=".R",full.names=TRUE)
sapply(R_function_files,source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)

#Set the file paths for saved files based on these IDs.
mats_and_params_file = "Data/8pcw/Mutation_matrices_and_parameters_8pcw_reduced"
filtered_muts_file = "Data/8pcw/Filtered_mut_set_8pcw"
dna_string_file = "Data/8pcw/DNA_string_file_8pcw.fa"
mpboot_tree_file = paste0(dna_string_file,".treefile")
tree_file_path = "Data/8pcw/Tree_8pcw.tree"
vcf_header_path = "Data/VCF_header_for_VaGrent.txt"
vcf_path = "Data/8pcw/Filtered_mut_set_vcf_8pcw.vcf"
vagrent_input_path = "Data/8pcw/Filtered_mut_set_vcf_with_header_8pcw.vcf"
vagrent_output_path = paste0(vagrent_input_path,".annot")
file_annot="Data/8pcw/Filtered_mut_set_annotated_8pcw"
sensitivity_analysis_path <- "Data/8pcw/sensitivity_analysis_8wks"  

#If have previous run analysis - load up files, or can rerun the analysis
previously_run = TRUE
if(previously_run) {
  load(file_annot); tree <- read.tree(tree_file_path)
} else {
  #Otherwise start here
  load(mats_and_params_file )
  COMB_mats$nsamp = ncol(COMB_mats$NV)
  COMB_mats$gender="male"
  
  min_sample_mean_cov = 0
  
  #Remove the low coverage samples and their private mutations
  if(min_sample_mean_cov > 0) {
    output = remove_low_coverage_samples(COMB_mats = COMB_mats,
                                         filter_params = filter_params,
                                         min_sample_mean_cov = min_sample_mean_cov,
                                         other_samples_to_remove = NULL,
                                         min_variant_reads_auto = 3, #these parameters are to remove mutations from the matrix that are no longer positive in any samples
                                         min_variant_reads_xy = 2)
    COMB_mats= output$COMB_mats
    filter_params = output$filter_params
  }
  
  
  #REVIEW MEAN DEPTH HISTOGRAMS TO DECIDE MANUAL CUT-OFFS FOR EXCLUDING OUTLIERS
  #hist(filter_params$mean_depth, breaks = 100, xlim = c(0,60))
  XY_low_depth_cutoff = 5; XY_high_depth_cutoff = 17; AUTO_low_depth_cutoff = 10; AUTO_high_depth_cutoff = 40 #for real analysis
  
  #This is the main function - applies the set cut-offs to the mutation set, filtering the mutations, assigning genotypes for each to all the samples, and building the dummy dna strings for tree building.
  filtered_muts = get_filtered_mut_set(input_set_ID = "8pcw",  #the Run_ID of the unfiltered mutation set used as input - though won't account for any removal of samples from set
                                       COMB_mats = COMB_mats,  #the main full mutation matrix as outputed by the "HSC_filtering_treebuild_table.R" script
                                       filter_params = filter_params,  #the filter_params matrix as outputed by the "HSC_filtering_treebuild_table.R" script
                                       gender = COMB_mats$gender, #patient's gender
                                       
                                       #These parameters decide whether a mutation is retained in the "true somatic mutation" set
                                       retain_muts = early_somatic_muts,  #any mutations that should be manually retained, despite not meeting filtering criteria, NULL by default
                                       germline_pval = -10,  #the log10 p-value cutoff for mutations coming from an expected germline distribution
                                       rho = 0.1,  #rho cutoff for the beta-binomial filter, a measure of how "over-dispersed" the counts are compared to a binomial distribution
                                       mean_depth = c(AUTO_low_depth_cutoff,AUTO_high_depth_cutoff, XY_low_depth_cutoff, XY_high_depth_cutoff),   #Numeric vector of length 4 defining mean depth at mutation site cut-offs. This is in the order 1. lower threshold for autosomes, 2. upper threshold for autosomes, 3. lower threshold for XY, 4. upper threshold for XY. This removes mis-mapping/ low reliability loci.
                                       pval_dp2=NA,  #the p-value cut-off if using the "pval within pos" filter, with positive samples defined as having >= 2 reads
                                       pval_dp3=0.01,   #the p-value cut-off if using the "pval within pos" filter, with positive samples defined as having >= 3 reads (allows for more index hopping)
                                       min_depth = c(6,4), #Numeric vector of length 2 defining minimum depths that at least one positive sample must have for mutation to be retained (AUTO and XY)
                                       min_pval_for_true_somatic = 0.1,   #Default: 0.1. the minimum p-value that at least one sample must have for the variant:normal read distribution coming from that expected for a true somatic
                                       min_vaf = NA, #Numeric vector of length 2 defining minimum vaf in at least one sample for mutation to be retained (AUTO and XY)
                                       
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
  tree <- read.tree(mpboot_tree_file)
  tree <- drop.tip(tree,"Ancestral")
  tree <- multi2di(tree)
  tree$edge.length = rep(1, nrow(tree$edge)) #Initially need to assign edge lengths of 1 for the tree_muts package to work
  
  #ASSIGN MUTATIONS TO THE TREE USING THE TREE_MUT PACKAGE
  df = reconstruct_genotype_summary(tree) #Define df (data frame) for treeshape
  
  #Get matrices in order, and run the main assignment functions
  mtr = filtered_muts$COMB_mats.tree.build$NV; mtr = as.matrix(mtr)
  depth = filtered_muts$COMB_mats.tree.build$NR; depth = as.matrix(depth)
  p.error = rep(0.01, ncol(filtered_muts$COMB_mats.tree.build$NR))
  res = assign_to_tree(mtr[,df$samples], depth[,df$samples], df, error_rate = p.error) #Get res (results!) object
  treefit_pval_cutoff = 1e-3
  poor_fit = res$summary$pval < treefit_pval_cutoff  #See how many mutations don't have read counts that fit the tree very well
  sum(poor_fit)
  hist(log(res$summary$pval[poor_fit]), breaks = 100)
  
  #Good point to review the poor fit ones.  Why poor fit? Then define if any "poor fit" mutations should be kept
  retain_poor_fit <- rep(FALSE, sum(poor_fit))
  retain_poor_fit[1] <- TRUE  #can edit this to any that you would like to keep
  poor_fit[which(poor_fit == TRUE)[retain_poor_fit]] <- FALSE
  
  #Remove poor_fit muts if they look dodgy
  filter_poor_fit = "yes"
  if(filter_poor_fit == "yes") {
    filtered_muts$COMB_mats.tree.build <- list_subset(filtered_muts$COMB_mats.tree.build,select_vector = !poor_fit)
    mtr = as.matrix(filtered_muts$COMB_mats.tree.build$NV)
    depth = as.matrix(filtered_muts$COMB_mats.tree.build$NR)
    p.error = c(rep(0.01, ncol(filtered_muts$COMB_mats.tree.build$NR)))
    res = assign_to_tree(mtr[,df$samples], depth[,df$samples], df, error_rate = p.error) #Get res (results!) object
  }
  
  #Assign edge lengths from the res object
  tree$edge.length <- res$df$df$edge_length
  
  #Add node and pval information to the filtered_muts object
  filtered_muts$COMB_mats.tree.build$mat$node <- tree$edge[res$summary$edge_ml,2]
  filtered_muts$COMB_mats.tree.build$mat$pval <- res$summary$pval
  
  #Save the res file and the tree file
  write.tree(tree, file = tree_file_path)
  
  #ANNOTATING THE MUTATIONS
  #-----------------------------------------------
  #Write vcf files for VariantCaller analysis - separate out "All mutations" & Shared mutations"
  vcf_file = create_vcf_files(filtered_muts$COMB_mats.tree.build$mat)
  write.table(vcf_file, sep = "\t", quote = FALSE, file = "Data/8pcw/Mutations_all_8pcw.vcf", row.names = FALSE)
  
  #Also save the shared mutations for mutational signature analysis
  shared_vcf_file = create_vcf_files(filtered_muts$COMB_mats.tree.build$mat,select_vector = which(!filtered_muts$COMB_mats.tree.build$mat$node%in%1:length(tree$tip.label)))
  write.table(shared_vcf_file, sep = "\t", quote = FALSE, file = "Data/8pcw/Mutations_shared_8pcw.vcf", row.names = FALSE)
  
  #-----------------------------------------------
  #Save the annotated filtered_muts files once Vagrent info is incorporated (post tree filtering)
  save(filtered_muts, file = file_annot)
  }

#-----------------------------------------------
#CORRECT THE EDGE LENGTHS BASED ON SAMPLE SENSITIVITY & INVITRO MUTATION PROPORTION
#NB - for function to work correctly, the sensitivity_df must be set up in exactly the correct format
#If there are no SNVs in the tree, and/or sensitivity analysis for these, set "include_indels" to FALSE

#Create trees of only SNVs or only INDELs (this function does not do any correction)
tree_SNV = get_subset_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, v.field="Mut_type", value = "SNV")
tree_INDEL = get_subset_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, v.field="Mut_type", value = "INDEL")

#Before correcting for sensitivity, correct for false positives by removing the likely proportion of in-vitro mutations
invitro_prop=0.27
tree_SNV$edge.length[tree_SNV$edge[,2]%in%1:length(tree_SNV$tip.label)] <-  (1-invitro_prop)*tree_SNV$edge.length[tree_SNV$edge[,2]%in%1:length(tree_SNV$tip.label)]

#Import the sensitivity analysis file
sensitivity_df <- read.delim(sensitivity_analysis_path, header = TRUE, stringsAsFactors = FALSE, sep = " ")

#Create corrected SNV tree
tree_SNV_c = get_corrected_tree(tree = tree_SNV, details = filtered_muts$COMB_mats.tree.build$mat, include_indels = FALSE, sensitivity_df = sensitivity_df,get_edge_from_tree=TRUE)

#Create other trees: (1) tree corrected for SNVs & indels (2) corrected indel tree
tree_c = get_corrected_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, include_indels = TRUE, sensitivity_df = sensitivity_df)
tree_INDEL_c = get_corrected_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, include_SNVs = FALSE, sensitivity_df = sensitivity_df)

#Create the "half-corrected" hybrid tree: Corrected SNVs, uncorrected INDELs (as difficult to do accurately & doesn't make a big difference)
tree_hybrid=tree_SNV_c #Start from corrected SNV tree
tree_hybrid$edge.length <- tree_SNV_c$edge.length + tree_INDEL$edge.length #Combine the edge lengths of the corrected SNV tree & the uncorrected indel tree
tree_hybrid.multi<-di2multi(tree_hybrid)

pdf("Figures/8pcw/Haematopoietic_phylogeny_8pcw.pdf",width=15,height=7)
tree_hybrid.multi=plot_tree(tree_hybrid.multi,cex.label = 0,lwd=1,default_edge_color = "black")
dev.off()

#---------------------------------------------------------------------------------------------------
#PLOT WITH ANNOTATION OF MUTATIONS CAUSING CODING CHANGES IN SHARED BRANCHES (the most robust calls)
#---------------------------------------------------------------------------------------------------

details.multi<-filtered_muts$COMB_mats.tree.build$mat
NV <- filtered_muts$COMB_mats.tree.build$NV
NR <- filtered_muts$COMB_mats.tree.build$NR

details.multi[,1:13]<-apply(details.multi[,1:13],2,as.character)
details.multi$variant_ID <- paste(details.multi$Gene, details.multi$Protein, sep = " ") #Add a "variant_ID" that includes that gene & protein coding change
details.multi$Type[details.multi$Type == ""] <- "no_annotation" #Label blank field as "no_annotation"

##Use treemut to reassign mutations to the multi-furcating tree
tree_hybrid.multi$tip.label<-gsub("_hum","",tree_hybrid.multi$tip.label)
df = reconstruct_genotype_summary(tree_hybrid.multi) #Define df (data frame) for treeshape

#Get matrices in order, and run the main assignment functions
mtr = as.matrix(NV);depth=as.matrix(NR)
colnames(depth)=colnames(mtr)=gsub("_hum","",colnames(mtr))
p.error = rep(0.01, ncol(NR))
res = assign_to_tree(mtr[,df$samples], depth[,df$samples], df, error_rate = p.error) #Get res (results!) object
details.multi$node<- tree_hybrid.multi$edge[res$summary$edge_ml,2]

#Define the coding changes in shared branches
details.multi$shared_coding_change <- ifelse(!details.multi$node %in% 1:length(tree$tip.label) & details.multi$Type %in% c("protein_coding:exon:CDS:substitution:codon_variant:non_synonymous_codon",
                                                                                                         "protein_coding:exon:CDS:insertion:frameshift_variant",
                                                                                                         "protein_coding:exon:CDS:deletion:frameshift_variant",
                                                                                                         "protein_coding:exon:CDS:substitution:codon_variant:stop_gained",
                                                                                                         "protein_coding:exon:CDS:substitution:codon_variant:initiator_codon_change",
                                                                                                         "protein_coding:exon:CDS:deletion:inframe_variant:inframe_codon_loss"),
                                       "Shared coding change", "no")

pdf("Figures/8pcw/Haematopoietic_phylogeny_shared_CDS_muts_8pcw.pdf",width=15,height=7)
tree_hybrid.multi=plot_tree(tree_hybrid.multi, cex.label = 0)
plot_tree_labels(tree_hybrid.multi,
                 details = details.multi,
                 type = "line",
                 query.field = "shared_coding_change",
                 data.frame(value="Shared coding change",col="red",pch = 17,stringsAsFactors = FALSE),
                 label.field = "variant_ID",
                 cex.label = 0.8)
dev.off()
