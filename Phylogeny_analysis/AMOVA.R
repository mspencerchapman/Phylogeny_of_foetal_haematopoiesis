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
setwd(my_working_directory)
treemut_dir="~/R Work/R_scripts/treemut/"
setwd(my_working_directory)

#Define the file paths for the data files
tree_file_path="Data/18pcw/Tree_18pcw.tree"
file_annot="Data/18pcw/Filtered_mut_set_annotated_18pcw"
sensitivity_analysis_path="Data/18pcw/sensitivity_analysis_18pcw"

R_function_files = list.files("R_functions",pattern=".R",full.names=TRUE)
sapply(R_function_files,source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)

load(file_annot); tree<- read.tree(tree_file_path); tree$tip.label=gsub("_hum","",tree$tip.label)
details <- filtered_muts$COMB_mats.tree.build$mat

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
sensitivity_df$Sample=gsub("_hum","",sensitivity_df$Sample)

#Correct the tree using this df
tree_SNV_c = get_corrected_tree(tree = tree_SNV, details = filtered_muts$COMB_mats.tree.build$mat, include_indels = FALSE, sensitivity_df = sensitivity_df,get_edge_from_tree=TRUE)

##LINEAR REGRESSION MODEL FOR CORRECTION
#This needs the table of sequencing summary statistics used for the regression
smry_seq_18pcw=read.csv("Data/18pcw/smry_seq.csv")
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
tree_hybrid$edge.length <- tree_SNV_c$edge.length + tree_INDEL$edge.length 

#Makes sense for this analysis to make the tree ultra-metric. Use Peter's function (need to source the code to make this work)
utree = make.ultrametric.tree(tree_hybrid)
#utree=force.ultrametric(tree_hybrid,method="extend")

cellkey <- data.frame(Sample = utree$tip.label, stringsAsFactors = FALSE)
cellkey$Cell_type <- sapply(utree$tip.label, function(x) smry_seq$Tissue[which(smry_seq$Donor_ID == x)])

#Build a distance matrix
myvcv <- vcv(utree)
mymax <- max(nodeHeights(utree))
distmat <- matrix(sapply(myvcv, function(cell) ((mymax-cell)*2)), nrow=nrow(myvcv))
colnames(distmat) <- colnames(myvcv)
rownames(distmat) <- rownames(myvcv)
diag(distmat) <- 0

#Setup functions for amova
# function to perform amova.
amova.fn <- function(distmat, groupnames, cell_key) {
  groupnums <- length(groupnames)
  dw <- c()
  cellnums <- c()
  for (i in 1:groupnums) {
    tgroup <- groupnames[i]
    tcells <- which(rownames(distmat) %in% cell_key$Sample[cell_key$Cell_type==tgroup])
    cellnums <- c(cellnums, length(tcells))
    dw <- c(dw, sum(distmat[tcells, tcells]))
  }
  tdist <- distmat[colnames(distmat) %in% cellkey$Sample[cell_key$Cell_type %in% groupnames], rownames(distmat) %in% cell_key$Sample[cell_key$Cell_type %in% groupnames]]
  dap <- (sum(tdist) - sum(dw))/2
  
  dfAP <- length(groupnames) - 1 
  dfWP <- sum(cellnums-1)
  N <- sum(cellnums)
  SSwp <- sum(dw/(cellnums*2))
  SSap <- sum(((dw + dap)/(2*N)) - (dw/(2*cellnums)))
  
  MSwp <- SSwp/dfWP
  MSap <- SSap/dfAP
  nc <- (N - (sum(cellnums^2)/N))/dfAP
  varwp <- MSwp
  varap <- (MSap - MSwp)/nc
  obsphi <- varap/(varwp + varap)
  return(obsphi)
}

# function to randomise sample labels and repeat
randamova.fn <- function(distmat, groupnames, cell_key) {
  
  # change added 2018.01.22: when randomizing, only include the part of the distance matrix that involves the cell types being considered.
  tcells <- which(rownames(distmat) %in% cell_key$Sample[cell_key$Cell_type %in% groupnames])
  #
  
  randmat <- distmat[tcells, tcells]
  colnames(randmat) <- sample(colnames(randmat))
  rownames(randmat) <- colnames(randmat)
  randphi <- amova.fn(distmat=randmat, groupnames=groupnames, cell_key=cell_key)
  return(randphi)
}

# function tying it all in together
amovapval.fn <- function(distmat, groupnames, cell_key, iterations, plottitle) {
  # calculate observed
  obsphi <- amova.fn(distmat=distmat, groupnames=groupnames, cell_key=cell_key)
  # calculate null
  randphis <- sapply(1:iterations, function(cell) randamova.fn(distmat = distmat, groupnames=groupnames, cell_key=cell_key))
  # calculate pval
  pval <- length(which(randphis>obsphi))/length(randphis) 
  
  hist(randphis, col="grey", 100, main=plottitle, xlab="Phi statistic")
  abline(v=obsphi, col="red", lwd=2)
  legend("topright", legend=paste0("Observed\n p = ", signif(pval,digits=2)), lwd=2, col="red", bty="n")
}

#RUN FOR THE 18 PCW SAMPLE LOCATIONS
pdf("Figures/18pcw/AMOVA_cell_location_18pcw.pdf",width=5,height=5)
amovapval.fn(distmat = distmat, groupnames = c("F1", "F2", "L"), cell_key = cellkey, iterations = 30000, plottitle = "Phylogenetic clustering of HSPC locations")
dev.off()


