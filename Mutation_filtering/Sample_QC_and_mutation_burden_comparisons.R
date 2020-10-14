my_working_directory="~/R_work/Phylogeny_of_foetal_haematopoiesis/"
treemut_dir="~/R_work/treemut/"
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
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)

#Import the 18pcw smry seq object and add a "foetus" column
smry_seq_18pcw<-read.csv("Data/18pcw/smry_seq_18pcw.csv")
smry_seq_18pcw$Foetus="18 pcw"

#Add the 8pcw smry seq object and create placeholder columns for the features that are absent from the 8wk
smry_seq_8pcw<-read.csv("Data/8pcw/smry_seq_8pcw.csv")
smry_seq_8pcw$Tissue="L"
smry_seq_8pcw$Cell_type="Progenitor"
smry_seq_8pcw$Foetus="8 pcw"

#Combine the two data frames
smry_seq=rbind(smry_seq_18pcw,smry_seq_8pcw[,colnames(smry_seq_18pcw)])

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
sample_coverage=sapply(gsub("_hum","",tree_SNV_c$tip.label),function(x) smry_seq$coverage[smry_seq$Donor_ID==x])
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
df.tidy$samples<-gsub("_hum","",df.tidy$samples)

#NOW READY TO CREATE THE PLOTS
palette=c("#E30613","#1D71B8")
my_theme = theme_classic(base_family="Arial")+theme(text=element_text(size=7),
                                                    axis.title = element_text(size=7),
                                                    axis.text = element_text(size = 7),
                                                    legend.text = element_text(size=7),
                                                    strip.text = element_text(size=7))

p1.1<-smry_seq %>%
  mutate(Foetus=fct_rev(Foetus)) %>%
  ggplot(aes(y=Percentage,x=Foetus,col=Foetus)) +
  geom_boxplot(outlier.size=0.5) +
  labs(x="Foetus",y="Human read percentage prior to filtering (%)") +
  scale_color_manual(values=palette) +
  #theme(legend.position = "none") +
  my_theme


p1.2<-smry_seq %>%
  mutate(Foetus=fct_rev(Foetus)) %>%
  #filter(Foetus=="18pcw") %>%
  ggplot(aes(x=Percentage,y=coverage,col=Foetus)) +
  labs(x="Human read percentage prior to filtering (%)",y="Coverage achieved post-filtering (x)") +
  scale_color_manual(values=palette) +
  geom_smooth(method="lm",size=0.5)+
  geom_point(size=0.5,alpha=0.5) +
  my_theme

p1.3<-smry_seq %>%
  mutate(Foetus=fct_rev(Foetus)) %>%
  ggplot(aes(y=coverage,x=Foetus,col=Foetus)) +
  geom_boxplot(outlier.size = 0.5,outlier.alpha = 0.75) +
  labs(y="Coverage achieved post-filtering (x)") +
  scale_color_manual(values=palette) +
  #theme(legend.position = "none") +
  my_theme

#Now do correlations with mutation burden, corrected and uncorrected
#Add information from df.tidy into the smry_seq df
smry_seq$SNV_burden_u=sapply(smry_seq$Donor_ID,function(sample) ifelse(sample%in%df.tidy$samples,df.tidy$number_of_mutations[df.tidy$samples==sample & df.tidy$mutation_burden=="SNV_muts_u"],NA))
smry_seq$SNV_burden_c=sapply(smry_seq$Donor_ID,function(sample) ifelse(sample%in%df.tidy$samples,df.tidy$number_of_mutations[df.tidy$samples==sample & df.tidy$mutation_burden=="SNV_muts_c"],NA))
smry_seq$INDEL_burden=sapply(smry_seq$Donor_ID,function(sample) ifelse(sample%in%df.tidy$samples,df.tidy$number_of_mutations[df.tidy$samples==sample & df.tidy$mutation_burden=="INDEL_muts"],NA))

#SNV BURDEN (uncorrected) vs COVERAGE plots
p1.4<-smry_seq %>%
  mutate(Foetus=fct_rev(Foetus)) %>%
  #filter(coverage>10)%>%
  ggplot(aes(x=coverage,y=SNV_burden_u,col=Foetus)) +
  geom_jitter(alpha=0.5,width = 0.1,height = 0.4,size=0.5) +
  geom_smooth(method="lm",fullrange=FALSE,size=0.5) +
  scale_color_manual(values=palette) +
  ylim(c(0,70)) +
  labs(y="Uncorrected SNV burden",x="Mean sequencing coverage") +
  my_theme


#SNV BURDEN (corrected) vs COVERAGE plots
p1.5<-smry_seq %>%
  mutate(Foetus=fct_rev(Foetus)) %>%
  #filter(coverage>10)%>%
  ggplot(aes(x=coverage,y=SNV_burden_c,col=Foetus)) +
  geom_point(alpha=0.5,size=0.5) +
  geom_smooth(method="lm",fullrange=FALSE,size=0.5) +
  scale_color_manual(values=palette) +
  ylim(c(0,75)) +
  labs(y="Corrected SNV burden",x="Mean sequencing coverage") +
  my_theme

#INDEL BURDEN vs COVERAGE plots
p1.6<-smry_seq %>%
  mutate(Foetus=fct_rev(Foetus)) %>%
  #filter(coverage>10)%>%
  ggplot(aes(x=coverage,y=INDEL_burden,col=Foetus)) +
  geom_jitter(alpha=0.5,width = 0,height = 0.25,size=0.5) +
  geom_smooth(method="lm",fullrange=FALSE,size=0.5) +
  scale_color_manual(values=palette) +
  labs(y="INDEL burden",x="Mean sequencing coverage") +
  my_theme

#Arrange the plots in a grid - this is the final figure
my_grid<-arrangeGrob(p1.1,p1.2,p1.3,p1.4,p1.5,p1.6,ncol=3,nrow=2)
plot(my_grid)
ggsave(file="Figures/Combined/Sample_contamination_and_coverage.pdf",my_grid,device=cairo_pdf,width=10,height=8)


#Plot SNV histogram
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
  summarise(mean=mean(number_of_mutations),n=n(),sd=sd(number_of_mutations))

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
ggsave(file="Figures/Combined/Mutation_burden_comparisons.pdf",my_grid_2,device=cairo_pdf,width = 8,height=2)
