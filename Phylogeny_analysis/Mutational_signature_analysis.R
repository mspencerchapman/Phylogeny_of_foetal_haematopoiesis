#Script to plot the 96-profile signatures
#The vcf files for this analysis are derived from the "filter_from_table" scripts & the "Mutation_validation" script
library(BSgenome)
library(ggplot2)
library(MutationalPatterns)
library(gridExtra)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome,character.only=TRUE)

my_working_directory="~/R Work/Fetal HSPCs/Phylogeny_of_foetal_haematopoiesis/"
setwd(my_working_directory)

#Import the vcf files
vcf_files <- c(list.files("Data/8pcw",pattern = ".vcf", full.names = TRUE),
               list.files("Data/18pcw",pattern = ".vcf", full.names = TRUE))
sample_names<-c("8 pcw: Clonal","8 pcw: All","8pcw: Shared","8pcw: Subclonal","18 pcw: Clonal","18 pcw: All","18pcw: Shared","18pcw: Subclonal")
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
clonality<-rep(c("Clonal","All","Shared","Subclonal"),2) #Define a vector of the vcf type

#VISUALIZE THE 96_profile PLOTS
#1. Generate the 96_profile counts
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
head(mut_mat)

#2.Plot them using "plot_96_profile" function on this matrix
p1<-plot_96_profile(mut_mat[,which(clonality=="Shared")])
p2<-plot_96_profile(mut_mat[,which(clonality=="Clonal")])
p3<-plot_96_profile(mut_mat[,which(clonality=="Subclonal")])
combined_plot=arrangeGrob(p1,p2,p3,nrow=3)
plot(combined_plot)
ggsave(combined_plot,filename = "Figures/Combined/Mutation_signatures_96profile.pdf",device="pdf",height=11.7,width=8.3)
