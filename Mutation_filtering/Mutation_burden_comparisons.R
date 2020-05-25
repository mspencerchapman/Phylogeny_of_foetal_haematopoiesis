my_working_directory = "~/R Work/Fetal HSPCs/Foetal_phylogeny/"
setwd(my_working_directory)

library(stringr)
library(ape)
library(seqinr)
library(ggtree)
library(tidyr)
library(dplyr)
library(ggplot2)
library(forcats)

R_function_files = list.files("R_functions",pattern=".R",full.names=TRUE)
sapply(R_function_files,source)
setwd("R_functions/treemut/"); source("treemut.R"); setwd(my_working_directory)

##----Import and correct the 8pcw mutation burden info-----
tree_file_path="Data/8pcw/Tree_8pcw.tree"
file_annot="Data/8pcw/Filtered_mut_set_annotated_8pcw"
sensitivity_analysis_path <- "Data/8pcw/sensitivity_analysis_8wks"

load(file_annot); tree <- read.tree(tree_file_path)
sensitivity_df <- read.delim(sensitivity_analysis_path, header = TRUE, stringsAsFactors = FALSE, sep = " ")

#Create trees of only SNVs or only INDELs (this function does not do any correction)
tree_SNV = get_subset_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, v.field="Mut_type", value = "SNV")
tree_INDEL = get_subset_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, v.field="Mut_type", value = "INDEL")

#Store the uncorrected mutation burden
SNV_burden_u =get_mut_burden(tree_SNV)

#Before correcting for sensitivity, correct for false positives by removing the propotion of invitro mutations in the private branches
invitro_prop=0.27
tree_SNV$edge.length[tree_SNV$edge[,2]%in%1:length(tree_SNV$tip.label)] <-  (1-invitro_prop)*tree_SNV$edge.length[tree_SNV$edge[,2]%in%1:length(tree_SNV$tip.label)]

#Create corrected SNV tree
tree_SNV_c = get_corrected_tree(tree = tree_SNV, details = filtered_muts$COMB_mats.tree.build$mat, include_indels = FALSE, sensitivity_df = sensitivity_df,get_edge_from_tree=TRUE)
tree_c = get_corrected_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, include_indels = TRUE, sensitivity_df = sensitivity_df)
tree_INDEL_c = get_corrected_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, include_SNVs = FALSE, sensitivity_df = sensitivity_df)

#Get mutation burdens
SNV_burden_c = get_mut_burden(tree_SNV_c)
INDEL_burden = get_mut_burden(tree_INDEL)
TOTAL_burden = SNV_burden_c+INDEL_burden
PROP_shared = nodeHeights(tree_SNV_c)[tree_SNV_c$edge[,2] %in% 1:length(tree_SNV_c$tip.label),1]/nodeHeights(tree_SNV_c)[tree_SNV_c$edge[,2] %in% 1:length(tree_SNV_c$tip.label),2]

df = data.frame(total_muts = TOTAL_burden, SNV_muts_c= SNV_burden_c, SNV_muts_u=SNV_burden_u,INDEL_muts = INDEL_burden,PROP_shared=PROP_shared)
df$samples = tree$tip.label
df.tidy = gather(df, key = "mutation_burden",value = "number_of_mutations", -samples)
df.tidy.8pcw <- df.tidy
df.tidy.8pcw$foetus <- "8 pcw"

##----Now import and correct the 18pcw mutation burden info-----

#Set the file paths for saved files based on these IDs.
tree_file_path = "Data/18pcw/Tree_18pcw.tree"
file_annot = "Data/18pcw/Filtered_mut_set_annotated_18pcw"
sensitivity_analysis_path <- "Data/18pcw/sensitivity_analysis_18pcw"

load(file_annot); tree <- read.tree(tree_file_path)
sensitivity_df <- read.delim(sensitivity_analysis_path, header = TRUE, stringsAsFactors = FALSE, sep = " ")
smry_seq_18wks<-read.csv("Data/18pcw/smry_seq.csv")

#Create trees of only SNVs or only INDELs (this function does not do any correction)
tree_SNV = get_subset_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, v.field="Mut_type", value = "SNV")
tree_INDEL = get_subset_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, v.field="Mut_type", value = "INDEL")

#Store the uncorrected mutation burden
SNV_burden_u =get_mut_burden(tree_SNV)

##PERFORM REGRESSION TO GET COVERAGE-SPECIFIC IN VITRO MUTS CORRECTION
#Create sensitivity-corrected SNV tree
tree_SNV_c = get_corrected_tree(tree = tree_SNV, details = filtered_muts$COMB_mats.tree.build$mat, include_indels = FALSE, sensitivity_df = sensitivity_df,get_edge_from_tree=TRUE)

##LINEAR REGRESSION MODEL FOR CORRECTION
#(1) Get the corrected private branch lengths for each sample
private_mut_burden=tree_SNV_c$edge.length[tree_SNV_c$edge[,2]%in%1:length(tree_SNV_c$tip.label)]
#(2) Get a matched vector of sequencing coverage
sample_coverage=sapply(gsub("_hum","",tree_SNV_c$tip.label),function(x) smry_seq_18wks$coverage[smry_seq_18wks$Donor_ID==x])
#(3) Visualize the correlation using the log of the private mutation burden
plot(log(private_mut_burden)~sample_coverage)
#(4) Linear model of this correlation and add to the plot
lin_mod <-lm(log(private_mut_burden)~sample_coverage); abline(lin_mod)
#(5) Use the correlation coefficient with coverage to correct the log(private_mut_burden)
#The 19.75 is chosen such that the total proportion of in vitro mutations is as per the binomial mixture model
cut_off=19.75
log_corrected_muts=log(private_mut_burden)+lin_mod$coefficients["sample_coverage"]*(cut_off-sample_coverage)
plot(log_corrected_muts~sample_coverage); abline(lm(log_corrected_muts~sample_coverage))
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

#Get mutation burdens
SNV_burden_c = get_mut_burden(tree_SNV_c)
INDEL_burden = get_mut_burden(tree_INDEL)
TOTAL_burden = SNV_burden_c+INDEL_burden
PROP_shared = nodeHeights(tree_SNV_c)[tree_SNV_c$edge[,2] %in% 1:length(tree_SNV_c$tip.label),1]/nodeHeights(tree_SNV_c)[tree_SNV_c$edge[,2] %in% 1:length(tree_SNV_c$tip.label),2]

df = data.frame(total_muts = TOTAL_burden, SNV_muts_c= SNV_burden_c, SNV_muts_u=SNV_burden_u,INDEL_muts = INDEL_burden,PROP_shared=PROP_shared)
df$samples = tree$tip.label
df.tidy = gather(df, key = "mutation_burden",value = "number_of_mutations", -samples)
df.tidy.18pcw <- df.tidy
df.tidy.18pcw$foetus <- "18 pcw"

#Combine the 18pcw and 8pcw dataframes
df.tidy=rbind(df.tidy.8pcw,df.tidy.18pcw)

#Plot SNV histogram
palette=c("#E30613","#1D71B8")
my_theme = theme_classic(base_size=7)

p2.1 <- df.tidy %>%
  mutate(foetus=fct_rev(foetus)) %>%
  filter(mutation_burden == "SNV_muts_c") %>%
  ggplot(aes(x = number_of_mutations,col=foetus,fill=foetus)) + 
  geom_histogram(alpha=0.8) +
  xlim(c(0,70)) +
  my_theme +
  theme(legend.position = "none") +
  facet_grid(rows=vars(foetus)) +
  labs(x = "Single nucleotide variants (corrected)", y = "Colonies") + 
  scale_color_manual(values=c("black","black")) +
  scale_fill_manual(values=palette)

#Plot INDEL histogram
p2.2 <- df.tidy %>%
  mutate(foetus=fct_rev(foetus)) %>%
  filter(mutation_burden == "INDEL_muts") %>%
  ggplot(aes(x = number_of_mutations, col=foetus,fill=foetus)) + 
  geom_histogram(binwidth=1) +
  scale_x_continuous(limits=c(-0.5,6.5),breaks = 0:6)+
  my_theme +
  theme(legend.position="none") +
  facet_grid(rows=vars(foetus)) +
  labs(x = "Indels", y = "Colonies") + 
  scale_color_manual(values=c("black","black")) +
  scale_fill_manual(values=palette)

#Plot the mean mutation burden for each foetus, separating into private & shared mutations
mean_SNV_burden <- df.tidy %>%
  filter(mutation_burden=="SNV_muts_c") %>%
  group_by(foetus) %>%
  summarise(mean=mean(number_of_mutations))

mean_shared<- df.tidy %>%
  filter(mutation_burden=="PROP_shared") %>%
  group_by(foetus) %>%
  summarise(mean=mean(number_of_mutations))

df2=data.frame(Foetus=c("18 pcw","8 pcw"),Shared_mutations=mean_shared$mean*mean_SNV_burden$mean,Private_mutations=(1-mean_shared$mean)*mean_SNV_burden$mean)

p2.3<- df2 %>%
  gather(key = "Type",value="Mutation_burden",-Foetus)%>%
  mutate(Foetus=fct_rev(Foetus)) %>%
  ggplot(aes(x=Foetus,y=Mutation_burden,col=Type,fill=Type)) +
  geom_bar(stat="identity")+
  labs(y="Mean single nucleotide variant burden")+
  scale_fill_manual(values=c("#C6C6C6","#878787"))+
  scale_color_manual(values=c("black","black")) +
  my_theme

#Arrange the plots and save
my_grid_2<-arrangeGrob(p2.1,p2.2,p2.3,ncol=3,widths=c(2,1.5,1.5))
plot(my_grid_2)
ggsave(file="Figures/Combined/Mutation_burden_comparisons.pdf",my_grid_2,device="pdf",width = 8,height=2.5)
