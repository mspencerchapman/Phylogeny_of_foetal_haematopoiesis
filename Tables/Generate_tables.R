#Generate the tables summarising the number of biopsies and mean depth for the LCM targeted sequencing
library(dplyr)

setwd("~/R Work/Fetal HSPCs/Foetal_phylogeny/")

#For the 8pcw foetus
lcm_smry<-read.csv("Data/8pcw/targeted_seq_ids_8wks.csv")
gpvaf_LCM_SNV_file ="Data/8pcw/PD43947.LCM.snp.tsv"
cgpvaf_LCM_INDEL_file ="Data/8pcw/PD43947.LCM.indel.tsv"
LCM_mats<-import_cgpvaf_SNV_and_INDEL(SNV_output_file=cgpvaf_LCM_SNV_file,INDEL_output_file=cgpvaf_LCM_INDEL_file)

lcm_smry$mean_coverage=sapply(lcm_smry$Sample_ID,function(Sample_ID) {mean(LCM_mats$NR[,Sample_ID])})
lcm_smry%>%group_by(Tissue)%>%summarise(n_biopsies=n(),aggregate_depth=sum(mean_coverage)) %>% write.csv(file = "Tables/8pcw/lcm_tissues_summary_8pcw.csv")
lcm_smry%>%arrange(Tissue)%>%dplyr::select(Sample_ID,Tissue,mean_coverage)%>% write.csv(file = "Tables/8pcw/lcm_tissues_summary_8pcw.csv")

#For the 18pcw foetus
lcm_smry<-read.csv("Data/18pcw/targeted_seq_ids_18wks.csv")
cgpvaf_LCM_SNV_file ="Data/18pcw/PD41768.LCM.snp.tsv"
LCM_mats<-import_cgpvaf_SNV_and_INDEL(SNV_output_file=cgpvaf_LCM_SNV_file)

lcm_smry$mean_coverage=sapply(lcm_smry$Sample_ID,function(Sample_ID) {mean(LCM_mats$NR[,Sample_ID])})
lcm_smry%>%group_by(Tissue)%>%summarise(n_biopsies=n(),aggregate_depth=sum(mean_coverage)) %>% write.csv(file = "Tables/18pcw/lcm_tissues_summary_18pcw.csv")
lcm_smry%>%arrange(Tissue)%>%dplyr::select(Sample_ID,Tissue,mean_coverage)%>% write.csv(file = "Tables/18pcw/lcm_individual_biopsy_summary_18pcw.csv")
