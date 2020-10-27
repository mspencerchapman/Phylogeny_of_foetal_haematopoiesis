

#Locally
SNV_output_file="~/Mounts/Lustre/fetal_HSC/fetal_8wks/LCM_WGS/cgpVAF_fetal_8wk_WGS_LCM_merged_SNVs.tsv"
colony_SNV_bed_file="~/Mounts/Lustre/fetal_HSC/fetal_8wks/caveman_raw/caveman_pileup/"
final_SNV_mutset="~/R_work/Phylogeny_of_foetal_haematopoiesis/Data/8pcw/Filtered_mut_set_annotated_8pcw"
R_function_files=list.files("~/R_work/Phylogeny_of_foetal_haematopoiesis/R_functions/",pattern=".R",full.names = T)

#On lustre - fetal_8wks folder
SNV_output_file="LCM_WGS/cgpVAF_fetal_8wk_WGS_LCM_merged_SNVs.tsv"
final_SNV_mutset="../Phylogeny_of_foetal_haematopoiesis/Data/8pcw/Filtered_mut_set_annotated_8pcw"
colony_SNV_bed_file="caveman_raw/caveman_pileup/fetal_8wks_all_samps.bed"
MS_filtered_bed_file="LCM_WGS/MS_filters/LCM_WGS_MSFiltered.bed"
R_function_files=list.files("../Phylogeny_of_foetal_haematopoiesis/R_functions/",pattern=".R",full.names = T)

sapply(R_function_files,source)
library(data.table)
library(pheatmap)

#Import
LCM_WGS_mats=import_cgpvaf_SNV_and_INDEL(SNV_output_file = SNV_output_file)
LCM_WGS_mats$gender="male"

colony_SNV_bed_mat=fread(colony_SNV_bed_file)
names(colony_SNV_bed_mat)=c("Chrom","Pos","Ref","Alt")
colony_SNV_bed_mat$Pos <- sapply(colony_SNV_bed_mat$Pos, function(x) gsub("^\\s+|\\s+$","",x)) #strip the white space
colony_SNV_bed_mat$mut_ref=apply(colony_SNV_bed_mat[,1:4],1,paste,collapse="-")

MS_filtered_SNV_bed_mat=fread(MS_filtered_bed_file,header=F)
names(MS_filtered_SNV_bed_mat)=c("Chrom","Pos","Ref","Alt")
MS_filtered_SNV_bed_mat$Pos <- sapply(MS_filtered_SNV_bed_mat$Pos, function(x) gsub("^\\s+|\\s+$","",x)) #strip the white space
MS_filtered_SNV_bed_mat$mut_ref=apply(MS_filtered_SNV_bed_mat[,1:4],1,paste,collapse="-")

sum(LCM_WGS_mats$mat$mut_ref%in%colony_SNV_bed_mat$mut_ref)
sum(LCM_WGS_mats$mat$mut_ref%in%filtered_muts$COMB_mats.tree.build$mat$mut_ref)
load(final_SNV_mutset)

LCM_WGS_new=list_subset(LCM_WGS_mats,select_vector = !LCM_WGS_mats$mat$mut_ref%in%colony_SNV_bed_mat$mut_ref & LCM_WGS_mats$mat$mut_ref%in%MS_filtered_SNV_bed_mat$mut_ref)

NR=LCM_WGS_new$NR
NV=LCM_WGS_new$NV


#Remove low depth mutations
hist(rowMeans(LCM_WGS_new$NR[grepl("X",LCM_WGS_new$mat$mut_ref)|grepl("Y",LCM_WGS_new$mat$mut_ref),]),breaks=50)
hist(rowMeans(LCM_WGS_new$NR[!grepl("X",LCM_WGS_new$mat$mut_ref)&!grepl("Y",LCM_WGS_new$mat$mut_ref),]),xlim=c(0,50),breaks=200);abline(v=45,col="red")

#Remove the new germline mutations & depth outlier mutations
filter_params=data.frame(mean_depth=rowMeans(LCM_WGS_new$NR))
filter_pass=data.frame(
  germline_pval=germline.binomial.filter(LCM_WGS_new) < 0.0001,
  #rho_val=beta.binom.filter(LCM_WGS_new)>0.01,
  mean_depth=sapply(1:nrow(LCM_WGS_new$NV),assess_mean_depth,COMB_mats=LCM_WGS_new,14,45,10,25)
  )
retain=apply(filter_pass,1,function(x) all(x==1))
LCM_WGS_filt=list_subset(LCM_WGS_new,retain)
cf=calculate_cell_frac(LCM_WGS_filt$NV,LCM_WGS_filt$NR)
pheatmap(cf,show_rownames=F)

lcm_wgs_smry=read.csv()
