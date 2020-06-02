library(ape)
library(ggplot2)

#Set working directories and import custom functions
my_working_directory="~/R Work/Fetal HSPCs/Phylogeny_of_foetal_haematopoiesis/"
treemut_dir="~/R Work/R_scripts/treemut/"
setwd(my_working_directory)

R_function_files = list.files("R_functions",pattern=".R",full.names=TRUE)
sapply(R_function_files,source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)

#Functions to define sample sensitivity for indels
#Take SNV branch lengths as surrogate for branch "time" and therefore take overall sample sensitivity as a mean of the
#sensitivity of each contributing branch scaled by the SNV branch length
get_branch_indel_sensitivity=function(node,tree,sensitivity_df) {
  daughters <- getTips(tree = tree, node = node)
  all_sens_INDELs <- sensitivity_df[sensitivity_df$Sample %in% daughters,"INDEL_sensitivity"]
  branch_INDEL_sens = 1 - prod(1-all_sens_INDELs) #branch sensitivity is calculated assuming coverage is randomly distributed over genome
  return(branch_INDEL_sens)
}
get_sample_indel_sensitivity=function(sample,tree,sensitivity_df){
  sample_nodes=get_ancestral_nodes(which(tree$tip.label==sample),tree$edge) #get the nodes contributing to that sample
  node_sensitivities=sapply(sample_nodes,get_branch_indel_sensitivity,tree,sensitivity_df) #get the overall indel sensitivity for each of these branches
  sample_edge_lengths = sapply(sample_nodes, function(node) tree$edge.length[tree$edge[,2]==node]) #get the corrected SNV edge lengths for each branch
  sample_indel_sensitivity=sum(sample_edge_lengths*node_sensitivities)/sum(sample_edge_lengths) #get a weighted mean sensitivity by summing the product of each edge length by the edge indel sensitivity, and dividing by the total edge length
  return(sample_indel_sensitivity)
}

#------------ANALYSIS FOR 8PCW FOETUS---------------------#
#Import files and set up the 8pcw corrected SNV tree
sensitivity_analysis_path <- "Data/8pcw/sensitivity_analysis_8wks"
tree_file_path = "Data/8pcw/Tree_8pcw.tree"
file_annot="Data/8pcw/Filtered_mut_set_annotated_8pcw"
tree=read.tree(tree_file_path); load(file_annot); sensitivity_df <- read.delim(sensitivity_analysis_path, header = TRUE, stringsAsFactors = FALSE, sep = " ")
tree_SNV = get_subset_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, v.field="Mut_type", value = "SNV")
invitro_prop=0.27
tree_SNV$edge.length[tree_SNV$edge[,2]%in%1:length(tree_SNV$tip.label)] <-  (1-invitro_prop)*tree_SNV$edge.length[tree_SNV$edge[,2]%in%1:length(tree_SNV$tip.label)]
tree_SNV_c = get_corrected_tree(tree = tree_SNV, details = filtered_muts$COMB_mats.tree.build$mat, include_indels = FALSE, sensitivity_df = sensitivity_df,get_edge_from_tree=TRUE)

make_ultramatric=TRUE
if(make_ultrametric) {
  #Make this tree ultrametric so that branch lengths best reflect time
  #First make 0 terminal branches 0.1 to avoid dividing by 0 in the algorithm
  tree_SNV_c$edge.length[tree_SNV_c$edge.length==0 & tree_SNV_c$edge[,2]%in%1:length(tree_SNV_c$tip.label)]<-0.1
  tree_SNV_c=make.ultrametric.tree(tree_SNV_c)
}

#Extract the raw mean indels per sample from the data
tree_INDEL = get_subset_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, v.field="Mut_type", value = "INDEL")
mean_indels.data=mean(get_mut_burden(tree_INDEL))

#NB- may actually be better to do this weighting with the ultrametric tree
#Using these function we can get a vector of the calculated sample sensitivities
sample_indel_sensitivities=sapply(tree_SNV_c$tip.label,get_sample_indel_sensitivity,tree=tree_SNV_c,sensitivity_df)

#Now, assuming a poisson distribution, we can infer the underlying lambda using a bayesian approach
#Run simulations with an underlying "true lambda" and infer how many indels will be detected for each sample and get the overall mean indel burden
#Use the mean indel burden as the "summary statistic" to gate the simulations and find the most likely range of true underlying lambdas
nsim=100000
lambda_vec=runif(nsim,min=1,max=10) #Prior probability is uniform distribution of 1 - 10 indels/ cell
results_df=data.frame(lambda=lambda_vec,mean_indels=NA)
for(i in 1:nsim) {
  lambda=lambda_vec[i]
  true_indel_no=rpois(n=length(tree$tip.label),lambda=lambda)
  names(true_indel_no)<- tree$tip.label
  indel_vec<-NULL
  for(sample in tree$tip.label) {
    indel_vec = c(indel_vec,rbinom(n=1,size=true_indel_no[sample],prob = sample_indel_sensitivities[sample]))
  }
  results_df$mean_indels[i] <- round(mean(indel_vec),digits = 2)
  if(i%%5000==0) {print(i)}
}

#Gate the results on those that give the same mean indel number as the data
lambda_results=results_df$lambda[results_df$mean_indels==round(mean_indels.data,digits = 2)]

#Get the 95% confidence interval on the true underlying mean indels
results_8pcw=quantile(lambda_results,probs = c(0.025,0.5,0.975))
data.frame(indel_mutation_rate=lambda_results)%>%
  ggplot(aes(x=indel_mutation_rate)) +
  geom_histogram(fill="dark blue",col="black",alpha=0.5) +
  theme_classic() +
  labs(title = "8 pcw foetus: posterior distribution of indel mutation rate")

 
#------------ANALYSIS FOR 18PCW FOETUS---------------------#
#Import files and set up the 8pcw corrected SNV tree
sensitivity_analysis_path <- "Data/18pcw/sensitivity_analysis_18pcw"
tree_file_path = "Data/18pcw/Tree_18pcw.tree"
file_annot="Data/18pcw/Filtered_mut_set_annotated_18pcw"
tree=read.tree(tree_file_path); load(file_annot); sensitivity_df <- read.delim(sensitivity_analysis_path, header = TRUE, stringsAsFactors = FALSE, sep = " ")

#Create trees of only SNVs or only INDELs (this function does not do any correction)
tree_SNV = get_subset_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, v.field="Mut_type", value = "SNV")

#Store the uncorrected SNV mutation burden
SNV_burden_u=get_mut_burden(tree_SNV)

#Correct the tree using the sensitivity df
tree_SNV_c = get_corrected_tree(tree = tree_SNV, details = filtered_muts$COMB_mats.tree.build$mat, include_indels = FALSE, sensitivity_df = sensitivity_df,get_edge_from_tree=TRUE)

##LINEAR REGRESSION MODEL FOR SNV CORRECTION (see "filtering_from_table_18pcw.R" script for more detail)
#This needs the table of sequencing summary statistics used for the regression
smry_seq_18pcw=read.csv("Data/18pcw/smry_seq_18pcw.csv")
private_mut_burden=tree_SNV_c$edge.length[tree_SNV_c$edge[,2]%in%1:length(tree_SNV_c$tip.label)]
sample_coverage=sapply(gsub("_hum","",tree_SNV_c$tip.label),function(x) smry_seq_18pcw$coverage[smry_seq_18pcw$Donor_ID==x])
lin_mod <-lm(log(private_mut_burden)~sample_coverage)
cut_off=19.75
log_corrected_muts=log(private_mut_burden)+lin_mod$coefficients["sample_coverage"]*(cut_off-sample_coverage)
corrected_muts=sapply(1:length(sample_coverage), function(i) {if(sample_coverage[i]<cut_off) {return(exp(log_corrected_muts[i]))} else {return(private_mut_burden[i])}})

#Create corrected SNV tree
tree_SNV_c$edge.length[tree_SNV_c$edge[,2]%in%1:length(tree$tip.label)] <- corrected_muts

make_ultrametric=TRUE
if(make_ultrametric) {
  #Make this tree ultrametric so that branch lengths best reflect time
  #First make 0 terminal branches 0.1 to avoid dividing by 0 in the algorithm
  tree_SNV_c$edge.length[tree_SNV_c$edge.length==0 & tree_SNV_c$edge[,2]%in%1:length(tree_SNV_c$tip.label)]<-0.1
  tree_SNV_c=make.ultrametric.tree(tree_SNV_c)
}

#Extract the raw mean indels per sample from the data
tree_INDEL = get_subset_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, v.field="Mut_type", value = "INDEL")
mean_indels.data=mean(get_mut_burden(tree_INDEL))

#NB- may actually be better to do this weighting with the ultrametric tree
#Using these function we can get a vector of the calculated sample sensitivities
sample_indel_sensitivities=sapply(tree_SNV_c$tip.label,get_sample_indel_sensitivity,tree=tree_SNV_c,sensitivity_df)

#Now, assuming a poisson distribution, we can infer the underlying lambda using a bayesian approach
#Run simulations with an underlying "true lambda" and infer how many indels will be detected for each sample and get the overall mean indel burden
#Use the mean indel burden as the "summary statistic" to gate the simulations and find the most likely range of true underlying lambdas
nsim=100000
lambda_vec=runif(nsim,min=1,max=10) #Prior probability is uniform distribution of 1 - 10 indels/ cell
results_df=data.frame(lambda=lambda_vec,mean_indels=NA)
for(i in 1:nsim) {
  lambda=lambda_vec[i]
  true_indel_no=rpois(n=length(tree$tip.label),lambda=lambda)
  names(true_indel_no)<- tree$tip.label
  indel_vec<-NULL
  for(sample in tree$tip.label) {
    indel_vec = c(indel_vec,rbinom(n=1,size=true_indel_no[sample],prob = sample_indel_sensitivities[sample]))
  }
  results_df$mean_indels[i] <- round(mean(indel_vec),digits = 2)
  if(i%%5000==0) {print(i)}
}

#Gate the results on those that give the same mean indel number as the data
lambda_results=results_df$lambda[results_df$mean_indels==round(mean_indels.data,digits = 2)]

#Get the 95% confidence interval on the true underlying mean indels
results_18pcw=quantile(lambda_results,probs = c(0.025,0.5,0.975))
data.frame(indel_mutation_rate=lambda_results)%>%
  ggplot(aes(x=indel_mutation_rate)) +
  geom_histogram(fill="dark blue",col="black",alpha=0.5) +
  theme_classic() +
  labs(title = "18 pcw foetus: posterior distribution of indel mutation rate")

#Print the results
print(rbind(results_8pcw,results_18pcw))
