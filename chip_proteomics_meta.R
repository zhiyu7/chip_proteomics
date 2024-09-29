gc()
rm(list=ls())

library(dplyr)
library(data.table)
library(rmeta)
library(tableone)
library(readxl)   
library(tidyr)
library(stringr)

####meta analysis using meta.summaries() in rmeta
#########JHS, MESA, CHS, FhS
map=fread("/medpop/esp2/zyu/chip_protemoics/data/JHS_protein_ID_key.csv")
#map=fread("/medpop/esp2/zyu/proteomics_pathway/data/chs/5K_Processed_ANML_Normalized/CHS_Markers_5K_ANML-Norm_2022_04.csv")
#map=map%>%select(SeqId, SomaId, TargetFullName, Target, UniProt, EntrezGeneSymbol)
jhs_univariate_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_JHS_proteomics_CHIP_noadjust.txt")
jhs_univariate_chip_proteomics=jhs_univariate_chip_proteomics%>%select(jhs_estimate=beta, jhs_se=se, jhs_p=p, jhs_n=N, Protein=SomaId, Chip=chip)
mesa_univariate_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_MESA_proteomics_CHIP_noadjust.txt")
mesa_univariate_chip_proteomics=mesa_univariate_chip_proteomics%>%select(mesa_estimate=beta, mesa_se=se, mesa_p=p, mesa_n=N, Protein=SomaId, Chip=chip)
chs_univariate_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_CHS_proteomics_CHIP_noadjust.txt")
chs_univariate_chip_proteomics=chs_univariate_chip_proteomics%>%select(chs_estimate=beta, chs_se=se, chs_p=p, chs_n=N, Protein=SomaId, Chip=chip)

jhs_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_JHS_proteomics_CHIP_fully_adjust.txt")
jhs_chip_proteomics=jhs_chip_proteomics%>%select(jhs_estimate=beta, jhs_se=se, jhs_p=p, jhs_n=N, Protein=SomaId, Chip=chip)
mesa_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_MESA_proteomics_CHIP_fully_adjust.txt")
mesa_chip_proteomics=mesa_chip_proteomics%>%select(mesa_estimate=beta, mesa_se=se, mesa_p=p, mesa_n=N, Protein=SomaId, Chip=chip)
chs_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_CHS_proteomics_CHIP_fully_adjust.txt")
chs_chip_proteomics=chs_chip_proteomics%>%select(chs_estimate=beta, chs_se=se, chs_p=p, chs_n=N, Protein=SomaId, Chip=chip)

jhs_nopeer_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_JHS_proteomics_CHIP_nopeer.txt", head=T)
jhs_nopeer_chip_proteomics=jhs_nopeer_chip_proteomics%>%select(jhs_estimate=beta, jhs_se=se, jhs_p=p, jhs_n=N, Protein=SomaId, Chip=chip)
mesa_nopeer_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_MESA_proteomics_CHIP_nopeer.txt", head=T)
mesa_nopeer_chip_proteomics=mesa_nopeer_chip_proteomics%>%select(mesa_estimate=beta, mesa_se=se, mesa_p=p, mesa_n=N, Protein=SomaId, Chip=chip)
chs_nopeer_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_CHS_proteomics_CHIP_nopeer.txt", head=T)
chs_nopeer_chip_proteomics=chs_nopeer_chip_proteomics%>%select(chs_estimate=beta, chs_se=se, chs_p=p, chs_n=N, Protein=SomaId, Chip=chip)

jhs_female_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_JHS_proteomics_CHIP_female_fully_adjust.txt", head=T)
jhs_female_chip_proteomics=jhs_female_chip_proteomics%>%select(jhs_estimate=beta, jhs_se=se, jhs_p=p, jhs_n=N, Protein=SomaId, Chip=chip)
mesa_female_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_MESA_proteomics_CHIP_female_fully_adjust.txt", head=T)
mesa_female_chip_proteomics=mesa_female_chip_proteomics%>%select(mesa_estimate=beta, mesa_se=se, mesa_p=p, mesa_n=N, Protein=SomaId, Chip=chip)
chs_female_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_CHS_proteomics_CHIP_female_fully_adjust.txt", head=T)
chs_female_chip_proteomics=chs_female_chip_proteomics%>%select(chs_estimate=beta, chs_se=se, chs_p=p, chs_n=N, Protein=SomaId, Chip=chip)

jhs_male_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_JHS_proteomics_CHIP_male_fully_adjust.txt", head=T)
jhs_male_chip_proteomics=jhs_male_chip_proteomics%>%select(jhs_estimate=beta, jhs_se=se, jhs_p=p, jhs_n=N, Protein=SomaId, Chip=chip)
mesa_male_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_MESA_proteomics_CHIP_male_fully_adjust.txt", head=T)
mesa_male_chip_proteomics=mesa_male_chip_proteomics%>%select(mesa_estimate=beta, mesa_se=se, mesa_p=p, mesa_n=N, Protein=SomaId, Chip=chip)
chs_male_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_CHS_proteomics_CHIP_male_fully_adjust.txt", head=T)
chs_male_chip_proteomics=chs_male_chip_proteomics%>%select(chs_estimate=beta, chs_se=se, chs_p=p, chs_n=N, Protein=SomaId, Chip=chip)

mesa_black_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Finalfinal0227_MESA_proteomics_CHIP_black_fully_adjust.txt", head=T)
mesa_black_chip_proteomics=mesa_black_chip_proteomics%>%select(mesa_estimate=beta, mesa_se=se, mesa_p=p, mesa_n=N, Protein=SomaId, Chip=chip)
chs_black_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_CHS_proteomics_CHIP_black_fully_adjust.txt", head=T)
chs_black_chip_proteomics=chs_black_chip_proteomics%>%select(chs_estimate=beta, chs_se=se, chs_p=p, chs_n=N, Protein=SomaId, Chip=chip)

mesa_white_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Finalfinal0227_MESA_proteomics_CHIP_white_fully_adjust.txt", head=T)
mesa_white_chip_proteomics=mesa_white_chip_proteomics%>%select(mesa_estimate=beta, mesa_se=se, mesa_p=p, mesa_n=N, Protein=SomaId, Chip=chip)
chs_white_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_CHS_proteomics_CHIP_white_fully_adjust.txt", head=T)
chs_white_chip_proteomics=chs_white_chip_proteomics%>%select(chs_estimate=beta, chs_se=se, chs_p=p, chs_n=N, Protein=SomaId, Chip=chip)

jhs_interaction_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_JHS_proteomics_interaction_CHIP_fully_adjust.txt")
jhs_interaction_chip_proteomics=jhs_interaction_chip_proteomics%>%select(jhs_estimate=beta, jhs_se=se, jhs_p=p, jhs_n=N, Protein=SomaId, Chip=chip)
mesa_interaction_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_MESA_proteomics_interaction_CHIP_fully_adjust.txt")
mesa_interaction_chip_proteomics=mesa_interaction_chip_proteomics%>%select(mesa_estimate=beta, mesa_se=se, mesa_p=p, mesa_n=N, Protein=SomaId, Chip=chip)
chs_interaction_chip_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_CHS_proteomics_interaction_CHIP_fully_adjust.txt")
chs_interaction_chip_proteomics=chs_interaction_chip_proteomics%>%select(chs_estimate=beta, chs_se=se, chs_p=p, chs_n=N, Protein=SomaId, Chip=chip)

#jhs_cad_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/JHS_proteomics_CADt_CAD_all_agesexracepc_withbatch.txt")
#jhs_cad_proteomics=jhs_cad_proteomics%>%select(jhs_estimate=beta,jhs_se=se,jhs_p=p,jhs_fdr=p_FDR,jhs_n=N,Protein=SomaId)
#mesa_cad_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/MESA_proteomics_CADt_CAD_all_agesexracepc_adjust.txt")
#mesa_cad_proteomics=mesa_cad_proteomics%>%select(mesa_estimate=beta,mesa_se=se,mesa_p=p,mesa_fdr=p_FDR,mesa_n=N,Protein=SomaId)
#chs_cad_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/CHS_proteomics_CADt_CAD_all_agesexracepc_adjust.txt")
#chs_cad_proteomics=chs_cad_proteomics%>%select(chs_estimate=beta,chs_se=se,chs_p=p,chs_fdr=p_FDR,chs_n=N,Protein=SomaId)

jhs_nopeer_cad_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/JHS_proteomics_CAD_t_CAD_all_fully_adjust.txt")
jhs_nopeer_cad_proteomics=jhs_nopeer_cad_proteomics%>%select(jhs_estimate=beta,jhs_se=se,jhs_p=p,jhs_n=N,Protein=SomaId)
mesa_nopeer_cad_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/MESA_proteomics_CADt_CAD_all_fully_adjust.txt")
mesa_nopeer_cad_proteomics=mesa_nopeer_cad_proteomics%>%select(mesa_estimate=beta,mesa_se=se,mesa_p=p,mesa_n=N,Protein=SomaId)
chs_nopeer_cad_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/CHS_proteomics_CADt_CAD_all_fully_adjust.txt")
chs_nopeer_cad_proteomics=chs_nopeer_cad_proteomics%>%select(chs_estimate=beta,chs_se=se,chs_p=p,chs_n=N,Protein=SomaId)

#######ARIC
aric_map=fread("/medpop/esp2/zyu/chip_protemoics/data/ARIC_protein_ID_key.csv")
 
read_excel_allsheets <- function(filename, tibble = FALSE) {
	sheets <- readxl::excel_sheets(filename)
	x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
	if(!tibble) x <- lapply(x, as.data.frame)
	names(x) <- sheets
	x
}


aric_univariate_ea_chip_proteomics <- read_excel_allsheets("/medpop/esp2/zyu/chip_protemoics/output/aric/univariate_results/aric_chip_proteom_univariate_results_ea_2023-01-26.xlsx")
aric_univariate_ea_chip_proteomics=bind_rows(aric_univariate_ea_chip_proteomics, .id = "CHIP")
aric_univariate_ea_chip_proteomics=merge(aric_univariate_ea_chip_proteomics, aric_map, by.x="seqid_in_sample", by.y="SeqId")
aric_univariate_ea_chip_proteomics=aric_univariate_ea_chip_proteomics%>%mutate(aric_ea_fdr=NA)%>%select(aric_ea_estimate=beta, aric_ea_se=se, aric_ea_p=p, aric_ea_fdr, aric_ea_n=N, Protein=SomaId, Chip=CHIP)
aric_univariate_ea_chip_proteomics$Chip <- str_replace(aric_univariate_ea_chip_proteomics$Chip, "large_ASXL1", "ASXL1_large")
aric_univariate_ea_chip_proteomics$Chip <- str_replace(aric_univariate_ea_chip_proteomics$Chip, "^((?!_large).)*ASXL1$", "ASXL1_all")
aric_univariate_ea_chip_proteomics$Chip <- str_replace(aric_univariate_ea_chip_proteomics$Chip, "large_DNMT3A", "DNMT3A_large")
aric_univariate_ea_chip_proteomics$Chip <- str_replace(aric_univariate_ea_chip_proteomics$Chip, "^((?!_large).)*DNMT3A$", "DNMT3A_all")
aric_univariate_ea_chip_proteomics$Chip <- str_replace(aric_univariate_ea_chip_proteomics$Chip, "large_TET2", "TET2_large")
aric_univariate_ea_chip_proteomics$Chip <- str_replace(aric_univariate_ea_chip_proteomics$Chip, "^((?!_large).)*TET2$", "TET2_all")
aric_univariate_ea_chip_proteomics$Chip <- str_replace(aric_univariate_ea_chip_proteomics$Chip, "large_JAK2", "JAK2_large")
aric_univariate_ea_chip_proteomics$Chip <- str_replace(aric_univariate_ea_chip_proteomics$Chip, "^((?!_large).)*JAK2$", "JAK2_all")
aric_univariate_ea_chip_proteomics$Chip <- str_replace(aric_univariate_ea_chip_proteomics$Chip, "large_CHIP", "vf_01")
aric_univariate_ea_chip_proteomics$Chip <- str_replace(aric_univariate_ea_chip_proteomics$Chip, "^((?!_large).)*CHIP$", "chip")


aric_univariate_aa_chip_proteomics <- read_excel_allsheets("/medpop/esp2/zyu/chip_protemoics/output/aric/univariate_results/aric_chip_proteom_univariate_results_aa_2023-01-26.xlsx")
aric_univariate_aa_chip_proteomics=bind_rows(aric_univariate_aa_chip_proteomics, .id = "CHIP")
aric_univariate_aa_chip_proteomics=merge(aric_univariate_aa_chip_proteomics, aric_map, by.x="seqid_in_sample", by.y="SeqId")
aric_univariate_aa_chip_proteomics=aric_univariate_aa_chip_proteomics%>%mutate(aric_aa_fdr=NA)%>%select(aric_aa_estimate=beta, aric_aa_se=se, aric_aa_p=p, aric_aa_fdr, aric_aa_n=N, Protein=SomaId, Chip=CHIP)
aric_univariate_aa_chip_proteomics$Chip <- str_replace(aric_univariate_aa_chip_proteomics$Chip, "large_ASXL1", "ASXL1_large")
aric_univariate_aa_chip_proteomics$Chip <- str_replace(aric_univariate_aa_chip_proteomics$Chip, "^((?!_large).)*ASXL1$", "ASXL1_all")
aric_univariate_aa_chip_proteomics$Chip <- str_replace(aric_univariate_aa_chip_proteomics$Chip, "large_DNMT3A", "DNMT3A_large")
aric_univariate_aa_chip_proteomics$Chip <- str_replace(aric_univariate_aa_chip_proteomics$Chip, "^((?!_large).)*DNMT3A$", "DNMT3A_all")
aric_univariate_aa_chip_proteomics$Chip <- str_replace(aric_univariate_aa_chip_proteomics$Chip, "large_TET2", "TET2_large")
aric_univariate_aa_chip_proteomics$Chip <- str_replace(aric_univariate_aa_chip_proteomics$Chip, "^((?!_large).)*TET2$", "TET2_all")
aric_univariate_aa_chip_proteomics$Chip <- str_replace(aric_univariate_aa_chip_proteomics$Chip, "large_JAK2", "JAK2_large")
aric_univariate_aa_chip_proteomics$Chip <- str_replace(aric_univariate_aa_chip_proteomics$Chip, "^((?!_large).)*JAK2$", "JAK2_all")
aric_univariate_aa_chip_proteomics$Chip <- str_replace(aric_univariate_aa_chip_proteomics$Chip, "large_CHIP", "vf_01")
aric_univariate_aa_chip_proteomics$Chip <- str_replace(aric_univariate_aa_chip_proteomics$Chip, "^((?!_large).)*CHIP$", "chip")


aric_ea_chip_proteomics <- read_excel_allsheets("/medpop/esp2/zyu/chip_protemoics/output/aric/covariate_and_PEERS_adjusted_results/aric_chip_proteom_cov_PEER_results_ea_2023-01-26.xlsx")
aric_ea_chip_proteomics=bind_rows(aric_ea_chip_proteomics, .id = "CHIP")
aric_ea_chip_proteomics=merge(aric_ea_chip_proteomics, aric_map, by.x="seqid_in_sample", by.y="SeqId")
aric_ea_chip_proteomics=aric_ea_chip_proteomics%>%mutate(aric_ea_fdr=NA)%>%select(aric_ea_estimate=beta, aric_ea_se=se, aric_ea_p=p, aric_ea_fdr, aric_ea_n=N, Protein=SomaId, Chip=CHIP)
aric_ea_chip_proteomics$Chip <- str_replace(aric_ea_chip_proteomics$Chip, "large_ASXL1", "ASXL1_large")
aric_ea_chip_proteomics$Chip <- str_replace(aric_ea_chip_proteomics$Chip, "^((?!_large).)*ASXL1$", "ASXL1_all")
aric_ea_chip_proteomics$Chip <- str_replace(aric_ea_chip_proteomics$Chip, "large_DNMT3A", "DNMT3A_large")
aric_ea_chip_proteomics$Chip <- str_replace(aric_ea_chip_proteomics$Chip, "^((?!_large).)*DNMT3A$", "DNMT3A_all")
aric_ea_chip_proteomics$Chip <- str_replace(aric_ea_chip_proteomics$Chip, "large_TET2", "TET2_large")
aric_ea_chip_proteomics$Chip <- str_replace(aric_ea_chip_proteomics$Chip, "^((?!_large).)*TET2$", "TET2_all")
aric_ea_chip_proteomics$Chip <- str_replace(aric_ea_chip_proteomics$Chip, "large_JAK2", "JAK2_large")
aric_ea_chip_proteomics$Chip <- str_replace(aric_ea_chip_proteomics$Chip, "^((?!_large).)*JAK2$", "JAK2_all")
aric_ea_chip_proteomics$Chip <- str_replace(aric_ea_chip_proteomics$Chip, "large_CHIP", "vf_01")
aric_ea_chip_proteomics$Chip <- str_replace(aric_ea_chip_proteomics$Chip, "^((?!_large).)*CHIP$", "chip")


aric_aa_chip_proteomics <- read_excel_allsheets("/medpop/esp2/zyu/chip_protemoics/output/aric/covariate_and_PEERS_adjusted_results/aric_chip_proteom_cov_PEER_results_aa_2023-01-26.xlsx")
aric_aa_chip_proteomics=bind_rows(aric_aa_chip_proteomics, .id = "CHIP")
aric_aa_chip_proteomics=merge(aric_aa_chip_proteomics, aric_map, by.x="seqid_in_sample", by.y="SeqId")
aric_aa_chip_proteomics=aric_aa_chip_proteomics%>%mutate(aric_aa_fdr=NA)%>%select(aric_aa_estimate=beta, aric_aa_se=se, aric_aa_p=p, aric_aa_fdr, aric_aa_n=N, Protein=SomaId, Chip=CHIP)
aric_aa_chip_proteomics$Chip <- str_replace(aric_aa_chip_proteomics$Chip, "large_ASXL1", "ASXL1_large")
aric_aa_chip_proteomics$Chip <- str_replace(aric_aa_chip_proteomics$Chip, "^((?!_large).)*ASXL1$", "ASXL1_all")
aric_aa_chip_proteomics$Chip <- str_replace(aric_aa_chip_proteomics$Chip, "large_DNMT3A", "DNMT3A_large")
aric_aa_chip_proteomics$Chip <- str_replace(aric_aa_chip_proteomics$Chip, "^((?!_large).)*DNMT3A$", "DNMT3A_all")
aric_aa_chip_proteomics$Chip <- str_replace(aric_aa_chip_proteomics$Chip, "large_TET2", "TET2_large")
aric_aa_chip_proteomics$Chip <- str_replace(aric_aa_chip_proteomics$Chip, "^((?!_large).)*TET2$", "TET2_all")
aric_aa_chip_proteomics$Chip <- str_replace(aric_aa_chip_proteomics$Chip, "large_JAK2", "JAK2_large")
aric_aa_chip_proteomics$Chip <- str_replace(aric_aa_chip_proteomics$Chip, "^((?!_large).)*JAK2$", "JAK2_all")
aric_aa_chip_proteomics$Chip <- str_replace(aric_aa_chip_proteomics$Chip, "large_CHIP", "vf_01")
aric_aa_chip_proteomics$Chip <- str_replace(aric_aa_chip_proteomics$Chip, "^((?!_large).)*CHIP$", "chip")


aric_nopeer_ea_chip_proteomics<- read_excel_allsheets("/medpop/esp2/zyu/chip_protemoics/output/aric/covariate_adjusted_results/aric_chip_proteom_cov_results_ea_2023-01-26.xlsx")
aric_nopeer_ea_chip_proteomics=bind_rows(aric_nopeer_ea_chip_proteomics, .id = "CHIP")
aric_nopeer_ea_chip_proteomics=merge(aric_nopeer_ea_chip_proteomics, aric_map, by.x="seqid_in_sample", by.y="SeqId")
aric_nopeer_ea_chip_proteomics=aric_nopeer_ea_chip_proteomics%>%mutate(aric_ea_fdr=NA)%>%select(aric_ea_estimate=beta, aric_ea_se=se, aric_ea_p=p, aric_ea_fdr, aric_ea_n=N, Protein=SomaId, Chip=CHIP)
aric_nopeer_ea_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_chip_proteomics$Chip, "large_ASXL1", "ASXL1_large")
aric_nopeer_ea_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_chip_proteomics$Chip, "^((?!_large).)*ASXL1$", "ASXL1_all")
aric_nopeer_ea_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_chip_proteomics$Chip, "large_DNMT3A", "DNMT3A_large")
aric_nopeer_ea_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_chip_proteomics$Chip, "^((?!_large).)*DNMT3A$", "DNMT3A_all")
aric_nopeer_ea_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_chip_proteomics$Chip, "large_TET2", "TET2_large")
aric_nopeer_ea_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_chip_proteomics$Chip, "^((?!_large).)*TET2$", "TET2_all")
aric_nopeer_ea_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_chip_proteomics$Chip, "large_JAK2", "JAK2_large")
aric_nopeer_ea_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_chip_proteomics$Chip, "^((?!_large).)*JAK2$", "JAK2_all")
aric_nopeer_ea_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_chip_proteomics$Chip, "large_CHIP", "vf_01")
aric_nopeer_ea_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_chip_proteomics$Chip, "^((?!_large).)*CHIP$", "chip")


aric_nopeer_aa_chip_proteomics<- read_excel_allsheets("/medpop/esp2/zyu/chip_protemoics/output/aric/covariate_adjusted_results/aric_chip_proteom_cov_results_aa_2023-01-26.xlsx")
aric_nopeer_aa_chip_proteomics=bind_rows(aric_nopeer_aa_chip_proteomics, .id = "CHIP")
aric_nopeer_aa_chip_proteomics=merge(aric_nopeer_aa_chip_proteomics, aric_map, by.x="seqid_in_sample", by.y="SeqId")
aric_nopeer_aa_chip_proteomics=aric_nopeer_aa_chip_proteomics%>%mutate(aric_aa_fdr=NA)%>%select(aric_aa_estimate=beta, aric_aa_se=se, aric_aa_p=p, aric_aa_fdr, aric_aa_n=N, Protein=SomaId, Chip=CHIP)
aric_nopeer_aa_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_chip_proteomics$Chip, "large_ASXL1", "ASXL1_large")
aric_nopeer_aa_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_chip_proteomics$Chip, "^((?!_large).)*ASXL1$", "ASXL1_all")
aric_nopeer_aa_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_chip_proteomics$Chip, "large_DNMT3A", "DNMT3A_large")
aric_nopeer_aa_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_chip_proteomics$Chip, "^((?!_large).)*DNMT3A$", "DNMT3A_all")
aric_nopeer_aa_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_chip_proteomics$Chip, "large_TET2", "TET2_large")
aric_nopeer_aa_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_chip_proteomics$Chip, "^((?!_large).)*TET2$", "TET2_all")
aric_nopeer_aa_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_chip_proteomics$Chip, "large_JAK2", "JAK2_large")
aric_nopeer_aa_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_chip_proteomics$Chip, "^((?!_large).)*JAK2$", "JAK2_all")
aric_nopeer_aa_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_chip_proteomics$Chip, "large_CHIP", "vf_01")
aric_nopeer_aa_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_chip_proteomics$Chip, "^((?!_large).)*CHIP$", "chip")


aric_ea_female_chip_proteomics <- read_excel_allsheets("/medpop/esp2/zyu/chip_protemoics/output/aric/covariate_and_PEERS_adjusted_results/aric_chip_proteom_cov_PEER_results_ea_female_2023-01-26.xlsx")
aric_ea_female_chip_proteomics=bind_rows(aric_ea_female_chip_proteomics, .id = "CHIP")
aric_ea_female_chip_proteomics=merge(aric_ea_female_chip_proteomics, aric_map, by.x="seqid_in_sample", by.y="SeqId")
aric_ea_female_chip_proteomics=aric_ea_female_chip_proteomics%>%mutate(aric_ea_fdr=NA)%>%select(aric_ea_estimate=beta, aric_ea_se=se, aric_ea_p=p, aric_ea_fdr, aric_ea_n=N, Protein=SomaId, Chip=CHIP)
aric_ea_female_chip_proteomics$Chip <- str_replace(aric_ea_female_chip_proteomics$Chip, "large_ASXL1", "ASXL1_large")
aric_ea_female_chip_proteomics$Chip <- str_replace(aric_ea_female_chip_proteomics$Chip, "^((?!_large).)*ASXL1$", "ASXL1_all")
aric_ea_female_chip_proteomics$Chip <- str_replace(aric_ea_female_chip_proteomics$Chip, "large_DNMT3A", "DNMT3A_large")
aric_ea_female_chip_proteomics$Chip <- str_replace(aric_ea_female_chip_proteomics$Chip, "^((?!_large).)*DNMT3A$", "DNMT3A_all")
aric_ea_female_chip_proteomics$Chip <- str_replace(aric_ea_female_chip_proteomics$Chip, "large_TET2", "TET2_large")
aric_ea_female_chip_proteomics$Chip <- str_replace(aric_ea_female_chip_proteomics$Chip, "^((?!_large).)*TET2$", "TET2_all")
aric_ea_female_chip_proteomics$Chip <- str_replace(aric_ea_female_chip_proteomics$Chip, "large_JAK2", "JAK2_large")
aric_ea_female_chip_proteomics$Chip <- str_replace(aric_ea_female_chip_proteomics$Chip, "^((?!_large).)*JAK2$", "JAK2_all")
aric_ea_female_chip_proteomics$Chip <- str_replace(aric_ea_female_chip_proteomics$Chip, "large_CHIP", "vf_01")
aric_ea_female_chip_proteomics$Chip <- str_replace(aric_ea_female_chip_proteomics$Chip, "^((?!_large).)*CHIP$", "chip")


aric_aa_female_chip_proteomics <- read_excel_allsheets("/medpop/esp2/zyu/chip_protemoics/output/aric/covariate_and_PEERS_adjusted_results/aric_chip_proteom_cov_PEER_results_aa_female_2023-01-26.xlsx")
aric_aa_female_chip_proteomics=bind_rows(aric_aa_female_chip_proteomics, .id = "CHIP")
aric_aa_female_chip_proteomics=merge(aric_aa_female_chip_proteomics, aric_map, by.x="seqid_in_sample", by.y="SeqId")
aric_aa_female_chip_proteomics=aric_aa_female_chip_proteomics%>%mutate(aric_aa_fdr=NA)%>%select(aric_aa_estimate=beta, aric_aa_se=se, aric_aa_p=p, aric_aa_fdr, aric_aa_n=N, Protein=SomaId, Chip=CHIP)
aric_aa_female_chip_proteomics$Chip <- str_replace(aric_aa_female_chip_proteomics$Chip, "large_ASXL1", "ASXL1_large")
aric_aa_female_chip_proteomics$Chip <- str_replace(aric_aa_female_chip_proteomics$Chip, "^((?!_large).)*ASXL1$", "ASXL1_all")
aric_aa_female_chip_proteomics$Chip <- str_replace(aric_aa_female_chip_proteomics$Chip, "large_DNMT3A", "DNMT3A_large")
aric_aa_female_chip_proteomics$Chip <- str_replace(aric_aa_female_chip_proteomics$Chip, "^((?!_large).)*DNMT3A$", "DNMT3A_all")
aric_aa_female_chip_proteomics$Chip <- str_replace(aric_aa_female_chip_proteomics$Chip, "large_TET2", "TET2_large")
aric_aa_female_chip_proteomics$Chip <- str_replace(aric_aa_female_chip_proteomics$Chip, "^((?!_large).)*TET2$", "TET2_all")
aric_aa_female_chip_proteomics$Chip <- str_replace(aric_aa_female_chip_proteomics$Chip, "large_JAK2", "JAK2_large")
aric_aa_female_chip_proteomics$Chip <- str_replace(aric_aa_female_chip_proteomics$Chip, "^((?!_large).)*JAK2$", "JAK2_all")
aric_aa_female_chip_proteomics$Chip <- str_replace(aric_aa_female_chip_proteomics$Chip, "large_CHIP", "vf_01")
aric_aa_female_chip_proteomics$Chip <- str_replace(aric_aa_female_chip_proteomics$Chip, "^((?!_large).)*CHIP$", "chip")

aric_nopeer_ea_female_chip_proteomics<- read_excel_allsheets("/medpop/esp2/zyu/chip_protemoics/output/aric/covariate_adjusted_results/aric_chip_proteom_cov_results_ea_female_2023-01-26.xlsx")
aric_nopeer_ea_female_chip_proteomics=bind_rows(aric_nopeer_ea_female_chip_proteomics, .id = "CHIP")
aric_nopeer_ea_female_chip_proteomics=merge(aric_nopeer_ea_female_chip_proteomics, aric_map, by.x="seqid_in_sample", by.y="SeqId")
aric_nopeer_ea_female_chip_proteomics=aric_nopeer_ea_female_chip_proteomics%>%mutate(aric_ea_fdr=NA)%>%select(aric_ea_estimate=beta, aric_ea_se=se, aric_ea_p=p, aric_ea_fdr, aric_ea_n=N, Protein=SomaId, Chip=CHIP)
aric_nopeer_ea_female_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_female_chip_proteomics$Chip, "large_ASXL1", "ASXL1_large")
aric_nopeer_ea_female_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_female_chip_proteomics$Chip, "^((?!_large).)*ASXL1$", "ASXL1_all")
aric_nopeer_ea_female_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_female_chip_proteomics$Chip, "large_DNMT3A", "DNMT3A_large")
aric_nopeer_ea_female_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_female_chip_proteomics$Chip, "^((?!_large).)*DNMT3A$", "DNMT3A_all")
aric_nopeer_ea_female_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_female_chip_proteomics$Chip, "large_TET2", "TET2_large")
aric_nopeer_ea_female_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_female_chip_proteomics$Chip, "^((?!_large).)*TET2$", "TET2_all")
aric_nopeer_ea_female_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_female_chip_proteomics$Chip, "large_JAK2", "JAK2_large")
aric_nopeer_ea_female_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_female_chip_proteomics$Chip, "^((?!_large).)*JAK2$", "JAK2_all")
aric_nopeer_ea_female_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_female_chip_proteomics$Chip, "large_CHIP", "vf_01")
aric_nopeer_ea_female_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_female_chip_proteomics$Chip, "^((?!_large).)*CHIP$", "chip")

aric_nopeer_aa_female_chip_proteomics<- read_excel_allsheets("/medpop/esp2/zyu/chip_protemoics/output/aric/covariate_adjusted_results/aric_chip_proteom_cov_results_aa_female_2023-01-26.xlsx")
aric_nopeer_aa_female_chip_proteomics=bind_rows(aric_nopeer_aa_female_chip_proteomics, .id = "CHIP")
aric_nopeer_aa_female_chip_proteomics=merge(aric_nopeer_aa_female_chip_proteomics, aric_map, by.x="seqid_in_sample", by.y="SeqId")
aric_nopeer_aa_female_chip_proteomics=aric_nopeer_aa_female_chip_proteomics%>%mutate(aric_aa_fdr=NA)%>%select(aric_aa_estimate=beta, aric_aa_se=se, aric_aa_p=p, aric_aa_fdr, aric_aa_n=N, Protein=SomaId, Chip=CHIP)
aric_nopeer_aa_female_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_female_chip_proteomics$Chip, "large_ASXL1", "ASXL1_large")
aric_nopeer_aa_female_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_female_chip_proteomics$Chip, "^((?!_large).)*ASXL1$", "ASXL1_all")
aric_nopeer_aa_female_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_female_chip_proteomics$Chip, "large_DNMT3A", "DNMT3A_large")
aric_nopeer_aa_female_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_female_chip_proteomics$Chip, "^((?!_large).)*DNMT3A$", "DNMT3A_all")
aric_nopeer_aa_female_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_female_chip_proteomics$Chip, "large_TET2", "TET2_large")
aric_nopeer_aa_female_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_female_chip_proteomics$Chip, "^((?!_large).)*TET2$", "TET2_all")
aric_nopeer_aa_female_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_female_chip_proteomics$Chip, "large_JAK2", "JAK2_large")
aric_nopeer_aa_female_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_female_chip_proteomics$Chip, "^((?!_large).)*JAK2$", "JAK2_all")
aric_nopeer_aa_female_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_female_chip_proteomics$Chip, "large_CHIP", "vf_01")
aric_nopeer_aa_female_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_female_chip_proteomics$Chip, "^((?!_large).)*CHIP$", "chip")

aric_ea_male_chip_proteomics <- read_excel_allsheets("/medpop/esp2/zyu/chip_protemoics/output/aric/covariate_and_PEERS_adjusted_results/aric_chip_proteom_cov_PEER_results_ea_male_2023-01-26.xlsx")
aric_ea_male_chip_proteomics=bind_rows(aric_ea_male_chip_proteomics, .id = "CHIP")
aric_ea_male_chip_proteomics=merge(aric_ea_male_chip_proteomics, aric_map,  by.x="seqid_in_sample", by.y="SeqId")
aric_ea_male_chip_proteomics=aric_ea_male_chip_proteomics%>%mutate(aric_ea_fdr=NA)%>%select(aric_ea_estimate=beta, aric_ea_se=se, aric_ea_p=p, aric_ea_fdr, aric_ea_n=N, Protein=SomaId, Chip=CHIP)
aric_ea_male_chip_proteomics$Chip <- str_replace(aric_ea_male_chip_proteomics$Chip, "large_ASXL1", "ASXL1_large")
aric_ea_male_chip_proteomics$Chip <- str_replace(aric_ea_male_chip_proteomics$Chip, "^((?!_large).)*ASXL1$", "ASXL1_all")
aric_ea_male_chip_proteomics$Chip <- str_replace(aric_ea_male_chip_proteomics$Chip, "large_DNMT3A", "DNMT3A_large")
aric_ea_male_chip_proteomics$Chip <- str_replace(aric_ea_male_chip_proteomics$Chip, "^((?!_large).)*DNMT3A$", "DNMT3A_all")
aric_ea_male_chip_proteomics$Chip <- str_replace(aric_ea_male_chip_proteomics$Chip, "large_TET2", "TET2_large")
aric_ea_male_chip_proteomics$Chip <- str_replace(aric_ea_male_chip_proteomics$Chip, "^((?!_large).)*TET2$", "TET2_all")
aric_ea_male_chip_proteomics$Chip <- str_replace(aric_ea_male_chip_proteomics$Chip, "large_JAK2", "JAK2_large")
aric_ea_male_chip_proteomics$Chip <- str_replace(aric_ea_male_chip_proteomics$Chip, "^((?!_large).)*JAK2$", "JAK2_all")
aric_ea_male_chip_proteomics$Chip <- str_replace(aric_ea_male_chip_proteomics$Chip, "large_CHIP", "vf_01")
aric_ea_male_chip_proteomics$Chip <- str_replace(aric_ea_male_chip_proteomics$Chip, "^((?!_large).)*CHIP$", "chip")

aric_aa_male_chip_proteomics <- read_excel_allsheets("/medpop/esp2/zyu/chip_protemoics/output/aric/covariate_and_PEERS_adjusted_results/aric_chip_proteom_cov_PEER_results_aa_male_2023-01-26.xlsx")
aric_aa_male_chip_proteomics=bind_rows(aric_aa_male_chip_proteomics, .id = "CHIP")
aric_aa_male_chip_proteomics=merge(aric_aa_male_chip_proteomics, aric_map,  by.x="seqid_in_sample", by.y="SeqId")
aric_aa_male_chip_proteomics=aric_aa_male_chip_proteomics%>%mutate(aric_aa_fdr=NA)%>%select(aric_aa_estimate=beta, aric_aa_se=se, aric_aa_p=p, aric_aa_fdr, aric_aa_n=N, Protein=SomaId, Chip=CHIP)
aric_aa_male_chip_proteomics$Chip <- str_replace(aric_aa_male_chip_proteomics$Chip, "large_ASXL1", "ASXL1_large")
aric_aa_male_chip_proteomics$Chip <- str_replace(aric_aa_male_chip_proteomics$Chip, "^((?!_large).)*ASXL1$", "ASXL1_all")
aric_aa_male_chip_proteomics$Chip <- str_replace(aric_aa_male_chip_proteomics$Chip, "large_DNMT3A", "DNMT3A_large")
aric_aa_male_chip_proteomics$Chip <- str_replace(aric_aa_male_chip_proteomics$Chip, "^((?!_large).)*DNMT3A$", "DNMT3A_all")
aric_aa_male_chip_proteomics$Chip <- str_replace(aric_aa_male_chip_proteomics$Chip, "large_TET2", "TET2_large")
aric_aa_male_chip_proteomics$Chip <- str_replace(aric_aa_male_chip_proteomics$Chip, "^((?!_large).)*TET2$", "TET2_all")
aric_aa_male_chip_proteomics$Chip <- str_replace(aric_aa_male_chip_proteomics$Chip, "large_JAK2", "JAK2_large")
aric_aa_male_chip_proteomics$Chip <- str_replace(aric_aa_male_chip_proteomics$Chip, "^((?!_large).)*JAK2$", "JAK2_all")
aric_aa_male_chip_proteomics$Chip <- str_replace(aric_aa_male_chip_proteomics$Chip, "large_CHIP", "vf_01")
aric_aa_male_chip_proteomics$Chip <- str_replace(aric_aa_male_chip_proteomics$Chip, "^((?!_large).)*CHIP$", "chip")

aric_nopeer_ea_male_chip_proteomics<- read_excel_allsheets("/medpop/esp2/zyu/chip_protemoics/output/aric/covariate_adjusted_results/aric_chip_proteom_cov_results_ea_male_2023-01-26.xlsx")
aric_nopeer_ea_male_chip_proteomics=bind_rows(aric_nopeer_ea_male_chip_proteomics, .id = "CHIP")
aric_nopeer_ea_male_chip_proteomics=merge(aric_nopeer_ea_male_chip_proteomics, aric_map, by.x="seqid_in_sample", by.y="SeqId")
aric_nopeer_ea_male_chip_proteomics=aric_nopeer_ea_male_chip_proteomics%>%mutate(aric_ea_fdr=NA)%>%select(aric_ea_estimate=beta, aric_ea_se=se, aric_ea_p=p, aric_ea_fdr, aric_ea_n=N, Protein=SomaId, Chip=CHIP)
aric_nopeer_ea_male_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_male_chip_proteomics$Chip, "large_ASXL1", "ASXL1_large")
aric_nopeer_ea_male_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_male_chip_proteomics$Chip, "^((?!_large).)*ASXL1$", "ASXL1_all")
aric_nopeer_ea_male_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_male_chip_proteomics$Chip, "large_DNMT3A", "DNMT3A_large")
aric_nopeer_ea_male_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_male_chip_proteomics$Chip, "^((?!_large).)*DNMT3A$", "DNMT3A_all")
aric_nopeer_ea_male_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_male_chip_proteomics$Chip, "large_TET2", "TET2_large")
aric_nopeer_ea_male_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_male_chip_proteomics$Chip, "^((?!_large).)*TET2$", "TET2_all")
aric_nopeer_ea_male_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_male_chip_proteomics$Chip, "large_JAK2", "JAK2_large")
aric_nopeer_ea_male_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_male_chip_proteomics$Chip, "^((?!_large).)*JAK2$", "JAK2_all")
aric_nopeer_ea_male_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_male_chip_proteomics$Chip, "large_CHIP", "vf_01")
aric_nopeer_ea_male_chip_proteomics$Chip <- str_replace(aric_nopeer_ea_male_chip_proteomics$Chip, "^((?!_large).)*CHIP$", "chip")

aric_nopeer_aa_male_chip_proteomics<- read_excel_allsheets("/medpop/esp2/zyu/chip_protemoics/output/aric/covariate_adjusted_results/aric_chip_proteom_cov_results_aa_male_2023-01-26.xlsx")
aric_nopeer_aa_male_chip_proteomics=bind_rows(aric_nopeer_aa_male_chip_proteomics, .id = "CHIP")
aric_nopeer_aa_male_chip_proteomics=merge(aric_nopeer_aa_male_chip_proteomics, aric_map, by.x="seqid_in_sample", by.y="SeqId")
aric_nopeer_aa_male_chip_proteomics=aric_nopeer_aa_male_chip_proteomics%>%mutate(aric_aa_fdr=NA)%>%select(aric_aa_estimate=beta, aric_aa_se=se, aric_aa_p=p, aric_aa_fdr, aric_aa_n=N, Protein=SomaId, Chip=CHIP)
aric_nopeer_aa_male_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_male_chip_proteomics$Chip, "large_ASXL1", "ASXL1_large")
aric_nopeer_aa_male_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_male_chip_proteomics$Chip, "^((?!_large).)*ASXL1$", "ASXL1_all")
aric_nopeer_aa_male_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_male_chip_proteomics$Chip, "large_DNMT3A", "DNMT3A_large")
aric_nopeer_aa_male_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_male_chip_proteomics$Chip, "^((?!_large).)*DNMT3A$", "DNMT3A_all")
aric_nopeer_aa_male_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_male_chip_proteomics$Chip, "large_TET2", "TET2_large")
aric_nopeer_aa_male_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_male_chip_proteomics$Chip, "^((?!_large).)*TET2$", "TET2_all")
aric_nopeer_aa_male_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_male_chip_proteomics$Chip, "large_JAK2", "JAK2_large")
aric_nopeer_aa_male_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_male_chip_proteomics$Chip, "^((?!_large).)*JAK2$", "JAK2_all")
aric_nopeer_aa_male_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_male_chip_proteomics$Chip, "large_CHIP", "vf_01")
aric_nopeer_aa_male_chip_proteomics$Chip <- str_replace(aric_nopeer_aa_male_chip_proteomics$Chip, "^((?!_large).)*CHIP$", "chip")


aric_interaction_ea_chip_proteomics <- read_excel_allsheets("/medpop/esp2/zyu/chip_protemoics/output/aric/interaction/aric_chip_proteom_model_1A_ea_2023-05-19.xlsx")
aric_interaction_ea_chip_proteomics=bind_rows(aric_interaction_ea_chip_proteomics, .id = "CHIP")
aric_interaction_ea_chip_proteomics=merge(aric_interaction_ea_chip_proteomics, aric_map, by.x="seqid_in_sample", by.y="SeqId")
aric_interaction_ea_chip_proteomics=aric_interaction_ea_chip_proteomics%>%mutate(aric_ea_fdr=NA)%>%select(aric_ea_estimate=beta, aric_ea_se=se, aric_ea_p=p, aric_ea_fdr, aric_ea_n=N, Protein=SomaId, Chip=CHIP)
aric_interaction_ea_chip_proteomics$Chip <- str_replace(aric_interaction_ea_chip_proteomics$Chip, "large_ASXL1", "ASXL1_large")
aric_interaction_ea_chip_proteomics$Chip <- str_replace(aric_interaction_ea_chip_proteomics$Chip, "^((?!_large).)*ASXL1$", "ASXL1_all")
aric_interaction_ea_chip_proteomics$Chip <- str_replace(aric_interaction_ea_chip_proteomics$Chip, "large_DNMT3A", "DNMT3A_large")
aric_interaction_ea_chip_proteomics$Chip <- str_replace(aric_interaction_ea_chip_proteomics$Chip, "^((?!_large).)*DNMT3A$", "DNMT3A_all")
aric_interaction_ea_chip_proteomics$Chip <- str_replace(aric_interaction_ea_chip_proteomics$Chip, "large_TET2", "TET2_large")
aric_interaction_ea_chip_proteomics$Chip <- str_replace(aric_interaction_ea_chip_proteomics$Chip, "^((?!_large).)*TET2$", "TET2_all")
aric_interaction_ea_chip_proteomics$Chip <- str_replace(aric_interaction_ea_chip_proteomics$Chip, "large_JAK2", "JAK2_large")
aric_interaction_ea_chip_proteomics$Chip <- str_replace(aric_interaction_ea_chip_proteomics$Chip, "^((?!_large).)*JAK2$", "JAK2_all")
aric_interaction_ea_chip_proteomics$Chip <- str_replace(aric_interaction_ea_chip_proteomics$Chip, "large_CHIP", "vf_01")
aric_interaction_ea_chip_proteomics$Chip <- str_replace(aric_interaction_ea_chip_proteomics$Chip, "^((?!_large).)*CHIP$", "chip")


aric_interaction_aa_chip_proteomics <- read_excel_allsheets("/medpop/esp2/zyu/chip_protemoics/output/aric/interaction/aric_chip_proteom_model_1A_aa_2023-05-19.xlsx")
aric_interaction_aa_chip_proteomics=bind_rows(aric_interaction_aa_chip_proteomics, .id = "CHIP")
aric_interaction_aa_chip_proteomics=merge(aric_interaction_aa_chip_proteomics, aric_map, by.x="seqid_in_sample", by.y="SeqId")
aric_interaction_aa_chip_proteomics=aric_interaction_aa_chip_proteomics%>%mutate(aric_aa_fdr=NA)%>%select(aric_aa_estimate=beta, aric_aa_se=se, aric_aa_p=p, aric_aa_fdr, aric_aa_n=N, Protein=SomaId, Chip=CHIP)
aric_interaction_aa_chip_proteomics$Chip <- str_replace(aric_interaction_aa_chip_proteomics$Chip, "large_ASXL1", "ASXL1_large")
aric_interaction_aa_chip_proteomics$Chip <- str_replace(aric_interaction_aa_chip_proteomics$Chip, "^((?!_large).)*ASXL1$", "ASXL1_all")
aric_interaction_aa_chip_proteomics$Chip <- str_replace(aric_interaction_aa_chip_proteomics$Chip, "large_DNMT3A", "DNMT3A_large")
aric_interaction_aa_chip_proteomics$Chip <- str_replace(aric_interaction_aa_chip_proteomics$Chip, "^((?!_large).)*DNMT3A$", "DNMT3A_all")
aric_interaction_aa_chip_proteomics$Chip <- str_replace(aric_interaction_aa_chip_proteomics$Chip, "large_TET2", "TET2_large")
aric_interaction_aa_chip_proteomics$Chip <- str_replace(aric_interaction_aa_chip_proteomics$Chip, "^((?!_large).)*TET2$", "TET2_all")
aric_interaction_aa_chip_proteomics$Chip <- str_replace(aric_interaction_aa_chip_proteomics$Chip, "large_JAK2", "JAK2_large")
aric_interaction_aa_chip_proteomics$Chip <- str_replace(aric_interaction_aa_chip_proteomics$Chip, "^((?!_large).)*JAK2$", "JAK2_all")
aric_interaction_aa_chip_proteomics$Chip <- str_replace(aric_interaction_aa_chip_proteomics$Chip, "large_CHIP", "vf_01")
aric_interaction_aa_chip_proteomics$Chip <- str_replace(aric_interaction_aa_chip_proteomics$Chip, "^((?!_large).)*CHIP$", "chip")

aric_interaction_race_chip_proteomics <- read_excel_allsheets("/medpop/esp2/zyu/chip_protemoics/output/aric/interaction/aric_chip_proteom_model_1B_2023-05-19.xlsx")
aric_interaction_race_chip_proteomics=bind_rows(aric_interaction_race_chip_proteomics, .id = "CHIP")
aric_interaction_race_chip_proteomics=merge(aric_interaction_race_chip_proteomics, aric_map, by.x="seqid_in_sample", by.y="SeqId")
aric_interaction_race_chip_proteomics=aric_interaction_race_chip_proteomics%>%mutate(aric_race_fdr=NA)%>%select(aric_race_estimate=beta, aric_race_se=se, aric_race_p=p, aric_race_fdr, aric_race_n=N, Protein=SomaId, Chip=CHIP)
aric_interaction_race_chip_proteomics$Chip <- str_replace(aric_interaction_race_chip_proteomics$Chip, "large_ASXL1", "ASXL1_large")
aric_interaction_race_chip_proteomics$Chip <- str_replace(aric_interaction_race_chip_proteomics$Chip, "^((?!_large).)*ASXL1$", "ASXL1_all")
aric_interaction_race_chip_proteomics$Chip <- str_replace(aric_interaction_race_chip_proteomics$Chip, "large_DNMT3A", "DNMT3A_large")
aric_interaction_race_chip_proteomics$Chip <- str_replace(aric_interaction_race_chip_proteomics$Chip, "^((?!_large).)*DNMT3A$", "DNMT3A_all")
aric_interaction_race_chip_proteomics$Chip <- str_replace(aric_interaction_race_chip_proteomics$Chip, "large_TET2", "TET2_large")
aric_interaction_race_chip_proteomics$Chip <- str_replace(aric_interaction_race_chip_proteomics$Chip, "^((?!_large).)*TET2$", "TET2_all")
aric_interaction_race_chip_proteomics$Chip <- str_replace(aric_interaction_race_chip_proteomics$Chip, "large_JAK2", "JAK2_large")
aric_interaction_race_chip_proteomics$Chip <- str_replace(aric_interaction_race_chip_proteomics$Chip, "^((?!_large).)*JAK2$", "JAK2_all")
aric_interaction_race_chip_proteomics$Chip <- str_replace(aric_interaction_race_chip_proteomics$Chip, "large_CHIP", "vf_01")
aric_interaction_race_chip_proteomics$Chip <- str_replace(aric_interaction_race_chip_proteomics$Chip, "^((?!_large).)*CHIP$", "chip")


aric_ea_cad_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/aric/CAD_results/aric_ea_CHD_proteom_cov_PEER_results_2023-01-26.csv")
aric_ea_cad_proteomics=merge(aric_ea_cad_proteomics, aric_map,  by.x="seqid_in_sample", by.y="SeqId")
aric_ea_cad_proteomics=aric_ea_cad_proteomics%>%mutate(aric_ea_fdr=NA)%>%select(aric_ea_estimate=beta,aric_ea_se=se,aric_ea_p=p,aric_ea_fdr,aric_ea_n=N,Protein=SomaId)

aric_aa_cad_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/aric/CAD_results/aric_aa_CHD_proteom_cov_PEER_results_2023-01-26.csv")
aric_aa_cad_proteomics=merge(aric_aa_cad_proteomics, aric_map,  by.x="seqid_in_sample", by.y="SeqId")
aric_aa_cad_proteomics=aric_aa_cad_proteomics%>%mutate(aric_aa_fdr=NA)%>%select(aric_aa_estimate=beta,aric_aa_se=se,aric_aa_p=p,aric_aa_fdr,aric_aa_n=N,Protein=SomaId)

aric_ea_nopeer_cad_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/aric/CAD_results/aric_ea_CHD_proteom_cov_results_2023-01-26.csv")
aric_ea_nopeer_cad_proteomics=merge(aric_ea_nopeer_cad_proteomics, aric_map,  by.x="seqid_in_sample", by.y="SeqId")
aric_ea_nopeer_cad_proteomics=aric_ea_nopeer_cad_proteomics%>%mutate(aric_ea_fdr=NA)%>%select(aric_ea_estimate=beta,aric_ea_se=se,aric_ea_p=p,aric_ea_fdr,aric_ea_n=N,Protein=SomaId)

aric_aa_nopeer_cad_proteomics=fread("/medpop/esp2/zyu/chip_protemoics/output/aric/CAD_results/aric_aa_CHD_proteom_cov_results_2023-01-26.csv")
aric_aa_nopeer_cad_proteomics=merge(aric_aa_nopeer_cad_proteomics, aric_map,  by.x="seqid_in_sample", by.y="SeqId")
aric_aa_nopeer_cad_proteomics=aric_aa_nopeer_cad_proteomics%>%mutate(aric_aa_fdr=NA)%>%select(aric_aa_estimate=beta,aric_aa_se=se,aric_aa_p=p,aric_aa_fdr,aric_aa_n=N,Protein=SomaId)

########################################META step######################################
########################################
######CAD and proteomics################
###########################################

#####no peer as primary analyses###########
meta_df_prot <- data.frame()
df1 = jhs_nopeer_cad_proteomics
df2 = mesa_nopeer_cad_proteomics
df3 = chs_nopeer_cad_proteomics
df4 = aric_ea_nopeer_cad_proteomics
df5 = aric_aa_nopeer_cad_proteomics

df <- merge(df1,df2,by=c("Protein"))
df <- merge(df,df3,by=c("Protein"))
df <- merge(df,df4,by=c("Protein"))
df <- merge(df,df5,by=c("Protein"))
	
for (j in 1:nrow(df)) {
    	df$meta_beta_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$summary
    	df$meta_se_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$se
    	df$meta_p_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$test[2]
    	df$meta_het_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$het[3]

    	df$meta_beta_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$summary
    	df$meta_se_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$se
    	df$meta_p_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$test[2]
    	df$meta_het_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$het[3]
}
  	#str(df)
meta_df_prot = df

meta_df_prot$meta_p_fix=ifelse(meta_df_prot$meta_p_fix==0, 1*10^-15, meta_df_prot$meta_p_fix)
meta_df_prot=meta_df_prot%>%mutate(meta_fdr_random = p.adjust (meta_p_random, method='fdr'))
meta_df_prot=meta_df_prot%>%mutate(meta_fdr_fix = p.adjust (meta_p_fix, method='fdr'))
meta_df_prot=merge(meta_df_prot,map,by.x="Protein",by.y="SomaId",all.x=T)

write.table(meta_df_prot,"/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Finalfinal0318_meta_nopeer_aricchsjhsmesa_cad_proteomics_notchangemap.txt", sep="\t", quote=FALSE, row.names=FALSE)

#############################################
####CHIP and proteomics#########
#############################################
#####univariate################
variants_list=c("vf_01", "DNMT3A_all", "TET2_all","ASXL1_all", "DNMT3A_large", "TET2_large","ASXL1_large","chip")

meta_df_prot <- data.frame()
for (i in 1:length(variants_list)) {
  	name = variants_list[i]
  	print(name)
  	df1 = jhs_univariate_chip_proteomics[which(jhs_univariate_chip_proteomics$Chip == paste0(name)),]
  	df2 = mesa_univariate_chip_proteomics[which(mesa_univariate_chip_proteomics$Chip == paste0(name)),]
  	df3 = chs_univariate_chip_proteomics[which(chs_univariate_chip_proteomics$Chip == paste0(name)),]
  	df4 = aric_univariate_ea_chip_proteomics[which(aric_univariate_ea_chip_proteomics$Chip == paste0(name)),]
  	df5 = aric_univariate_aa_chip_proteomics[which(aric_univariate_aa_chip_proteomics$Chip == paste0(name)),]
  
  	df <- merge(df1,df2,by=c("Protein","Chip"))
  	df <- merge(df,df3,by=c("Protein","Chip"))
  	df <- merge(df,df4,by=c("Protein","Chip"))
  	df <- merge(df,df5,by=c("Protein","Chip"))

  	for (j in 1:nrow(df)) {
    		df$meta_beta_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$summary
    		df$meta_se_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$se
    		df$meta_p_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$test[2]
    		df$meta_het_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$het[3]

    		df$meta_beta_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$summary
    		df$meta_se_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$se
    		df$meta_p_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$test[2]
    		df$meta_het_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$het[3]
  	}
  	#str(df)
  	meta_df_prot <- rbind(meta_df_prot,df)
}

variants_list=c("JAK2_all", "JAK2_large")

jak2_df_prot <- data.frame()
for (i in 1:length(variants_list)) {
  	name = variants_list[i]
  	print(name)
  	df1 = chs_univariate_chip_proteomics[which(chs_univariate_chip_proteomics$Chip == paste0(name)),]
  	df2 = aric_univariate_ea_chip_proteomics[which(aric_univariate_ea_chip_proteomics$Chip == paste0(name)),]
  
  	df <- merge(df1,df2,by=c("Protein","Chip"))

  	for (j in 1:nrow(df)) {
    		df$meta_beta_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$summary
    		df$meta_se_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$se
    		df$meta_p_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$test[2]
    		df$meta_het_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$het[3]

    		df$meta_beta_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$summary
    		df$meta_se_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$se
    		df$meta_p_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$test[2]
    		df$meta_het_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$het[3]
  	}
  	#str(df)
  	jak2_df_prot <- rbind(jak2_df_prot,df)
}

meta_df_prot=bind_rows(meta_df_prot, jak2_df_prot)
meta_df_prot$meta_p_fix=ifelse(meta_df_prot$meta_p_fix==0, 1*10^-15, meta_df_prot$meta_p_fix)
meta_df_prot_bytype=meta_df_prot
meta_df_prot_bytype$Chip=ifelse(meta_df_prot_bytype$Chip=="chip", "chip_all", ifelse(meta_df_prot_bytype$Chip=="vf_01", "chip_large", meta_df_prot_bytype$Chip))
meta_df_prot_bytype$Chip_cat=paste0(substr(meta_df_prot_bytype$Chip,1,3), meta_df_prot_bytype$TargetFullName)
meta_df_prot_bytype=meta_df_prot_bytype %>% group_by(Protein, Chip_cat) %>% slice(which.min(meta_p_fix))
meta_df_prot_bytype$Chip_use=sub("_.*", "", meta_df_prot_bytype$Chip) 
meta_df_prot_bytype=merge(meta_df_prot_bytype,map,by.x="Protein",by.y="SomaId",all.x=T)
meta_df_prot_bytype_nojak2=meta_df_prot_bytype%>% filter(Chip_use != "JAK2") %>% mutate(meta_fdr_fix_fourgroup = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_bytype_nojak2=meta_df_prot_bytype_nojak2%>%select(Protein, Chip, meta_fdr_fix_fourgroup)
meta_df_prot_bytype=meta_df_prot_bytype%>% ungroup()%>% mutate(meta_fdr_fix_group_all = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_bytype=merge(meta_df_prot_bytype, meta_df_prot_bytype_nojak2, by=c("Protein","Chip"), all.x=T)

write.table(meta_df_prot_bytype,"/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Finalfinal0318_univariate_fivegroups_meta_chip_proteomics_log2_nochangemap.txt", sep="\t", quote=FALSE, row.names=FALSE)

#meta_df_prot<-meta_df_prot[order(meta_df_prot$meta_fdr_fix),]
meta_df_prot_nojak2=meta_df_prot%>% filter(Chip %in% c("vf_01", "DNMT3A_all", "TET2_all","ASXL1_all", "DNMT3A_large", "TET2_large","ASXL1_large","chip")) %>% mutate(meta_fdr_fix_eightgroup = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_nojak2=meta_df_prot_nojak2%>%select(Protein, Chip, meta_fdr_fix_eightgroup)
meta_df_prot=meta_df_prot%>% ungroup()%>% mutate(meta_fdr_fix_group_all = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot=merge(meta_df_prot, meta_df_prot_nojak2, by=c("Protein","Chip"), all.x=T)
meta_df_prot=merge(meta_df_prot,map,by.x="Protein",by.y="SomaId",all.x=T)

write.table(meta_df_prot,"/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Finalfinal0318_univariate_meta_chip_proteomics_log2_nochangemap.txt", sep="\t", quote=FALSE, row.names=FALSE)


####with peers##############
variants_list=c("vf_01", "DNMT3A_all", "TET2_all","ASXL1_all", "DNMT3A_large", "TET2_large","ASXL1_large","chip")

meta_df_prot <- data.frame()
for (i in 1:length(variants_list)) {
  	name = variants_list[i]
  	print(name)
  	df1 = jhs_chip_proteomics[which(jhs_chip_proteomics$Chip == paste0(name)),]
  	df2 = mesa_chip_proteomics[which(mesa_chip_proteomics$Chip == paste0(name)),]
  	df3 = chs_chip_proteomics[which(chs_chip_proteomics$Chip == paste0(name)),]
  	df4 = aric_ea_chip_proteomics[which(aric_ea_chip_proteomics$Chip == paste0(name)),]
  	df5 = aric_aa_chip_proteomics[which(aric_aa_chip_proteomics$Chip == paste0(name)),]
  
  	df <- merge(df1,df2,by=c("Protein","Chip"))
  	df <- merge(df,df3,by=c("Protein","Chip"))
  	df <- merge(df,df4,by=c("Protein","Chip"))
  	df <- merge(df,df5,by=c("Protein","Chip"))

  	for (j in 1:nrow(df)) {
    		df$meta_beta_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$summary
    		df$meta_se_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$se
    		df$meta_p_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$test[2]
    		df$meta_het_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$het[3]

    		df$meta_beta_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$summary
    		df$meta_se_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$se
    		df$meta_p_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$test[2]
    		df$meta_het_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$het[3]
  	}
  	#str(df)
  	meta_df_prot <- rbind(meta_df_prot,df)
}

variants_list=c("JAK2_all", "JAK2_large")

jak2_df_prot <- data.frame()
for (i in 1:length(variants_list)) {
  	name = variants_list[i]
  	print(name)
  	df1 = chs_chip_proteomics[which(chs_chip_proteomics$Chip == paste0(name)),]
  	df2 = aric_ea_chip_proteomics[which(aric_ea_chip_proteomics$Chip == paste0(name)),]
  
  	df <- merge(df1,df2,by=c("Protein","Chip"))

  	for (j in 1:nrow(df)) {
    		df$meta_beta_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$summary
    		df$meta_se_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$se
    		df$meta_p_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$test[2]
    		df$meta_het_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$het[3]

    		df$meta_beta_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$summary
    		df$meta_se_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$se
    		df$meta_p_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$test[2]
    		df$meta_het_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$het[3]
	}
  	#str(df)
  	jak2_df_prot <- rbind(jak2_df_prot,df)
}

meta_df_prot=bind_rows(meta_df_prot, jak2_df_prot)
meta_df_prot$meta_p_fix=ifelse(meta_df_prot$meta_p_fix==0, 1*10^-15, meta_df_prot$meta_p_fix)
meta_df_prot_bytype=meta_df_prot
meta_df_prot_bytype$Chip=ifelse(meta_df_prot_bytype$Chip=="chip", "chip_all", ifelse(meta_df_prot_bytype$Chip=="vf_01", "chip_large", meta_df_prot_bytype$Chip))
meta_df_prot_bytype$Chip_cat=paste0(substr(meta_df_prot_bytype$Chip,1,3), meta_df_prot_bytype$TargetFullName)
meta_df_prot_bytype=meta_df_prot_bytype %>% group_by(Protein, Chip_cat) %>% slice(which.min(meta_p_fix))
meta_df_prot_bytype$Chip_use=sub("_.*", "", meta_df_prot_bytype$Chip) 
meta_df_prot_bytype=merge(meta_df_prot_bytype,map,by.x="Protein",by.y="SomaId",all.x=T)
meta_df_prot_bytype_nojak2=meta_df_prot_bytype%>% filter(Chip_use != "JAK2") %>% mutate(meta_fdr_fix_fourgroup = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_bytype_nojak2=meta_df_prot_bytype_nojak2%>%select(Protein, Chip, meta_fdr_fix_fourgroup)
meta_df_prot_bytype=meta_df_prot_bytype%>% ungroup()%>% mutate(meta_fdr_fix_group_all = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_bytype=merge(meta_df_prot_bytype, meta_df_prot_bytype_nojak2, by=c("Protein","Chip"), all.x=T)

write.table(meta_df_prot_bytype,"/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Finalfinal0318_fivegroups_meta_chip_proteomics_log2peer50_nochangemap.txt", sep="\t", quote=FALSE, row.names=FALSE)

meta_df_prot_nojak2=meta_df_prot%>% filter(Chip %in% c("vf_01", "DNMT3A_all", "TET2_all","ASXL1_all", "DNMT3A_large", "TET2_large","ASXL1_large","chip")) %>% mutate(meta_fdr_fix_eightgroup = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_nojak2=meta_df_prot_nojak2%>%select(Protein, Chip, meta_fdr_fix_eightgroup)
meta_df_prot=meta_df_prot%>% ungroup()%>% mutate(meta_fdr_fix_group_all = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot=merge(meta_df_prot, meta_df_prot_nojak2, by=c("Protein","Chip"), all.x=T)
meta_df_prot=merge(meta_df_prot,map,by.x="Protein",by.y="SomaId",all.x=T)

write.table(meta_df_prot,"/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Finalfinal0318_meta_chip_proteomics_log2peer50_nochangemap.txt", sep="\t", quote=FALSE, row.names=FALSE)

##########no peers as sensitivy analysis
variants_list=c("vf_01", "DNMT3A_all", "TET2_all","ASXL1_all", "DNMT3A_large", "TET2_large","ASXL1_large","chip")

meta_df_prot <- data.frame()
for (i in 1:length(variants_list)) {
	name = variants_list[i]
  	print(name)
  	df1 = jhs_nopeer_chip_proteomics[which(jhs_nopeer_chip_proteomics$Chip == paste0(name)),]
  	df2 = mesa_nopeer_chip_proteomics[which(mesa_nopeer_chip_proteomics$Chip == paste0(name)),]
  	df3 = chs_nopeer_chip_proteomics[which(chs_nopeer_chip_proteomics$Chip == paste0(name)),]
  	df4 = aric_nopeer_ea_chip_proteomics[which(aric_nopeer_ea_chip_proteomics$Chip == paste0(name)),]
  	df5 = aric_nopeer_aa_chip_proteomics[which(aric_nopeer_aa_chip_proteomics$Chip == paste0(name)),]
  
  	df <- merge(df1,df2,by=c("Protein","Chip"))
  	df <- merge(df,df3,by=c("Protein","Chip"))
  	df <- merge(df,df4,by=c("Protein","Chip"))
  	df <- merge(df,df5,by=c("Protein","Chip"))

  	for (j in 1:nrow(df)) {
    		df$meta_beta_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$summary
    		df$meta_se_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$se
    		df$meta_p_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$test[2]
    		df$meta_het_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$het[3]

    		df$meta_beta_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$summary
    		df$meta_se_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$se
    		df$meta_p_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$test[2]
    		df$meta_het_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$het[3]
  	}
  	#str(df)
  	meta_df_prot <- rbind(meta_df_prot,df)
}

variants_list=c("JAK2_all", "JAK2_large")

jak2_df_prot <- data.frame()
for (i in 1:length(variants_list)) {
  	name = variants_list[i]
  	print(name)
  	df1 = chs_nopeer_chip_proteomics[which(chs_nopeer_chip_proteomics$Chip == paste0(name)),]
  	df2 = aric_nopeer_ea_chip_proteomics[which(aric_nopeer_ea_chip_proteomics$Chip == paste0(name)),]

  	df <- merge(df1,df2,by=c("Protein","Chip"))

  	for (j in 1:nrow(df)) {
    		df$meta_beta_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$summary
    		df$meta_se_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$se
    		df$meta_p_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$test[2]
    		df$meta_het_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$het[3]

    		df$meta_beta_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$summary
    		df$meta_se_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$se
    		df$meta_p_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$test[2]
    		df$meta_het_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$het[3]
	}
  	#str(df)
  	jak2_df_prot <- rbind(jak2_df_prot,df)
}

meta_df_prot=bind_rows(meta_df_prot, jak2_df_prot)
meta_df_prot$meta_p_fix=ifelse(meta_df_prot$meta_p_fix==0, 1*10^-15, meta_df_prot$meta_p_fix)
meta_df_prot_bytype=meta_df_prot
meta_df_prot_bytype$Chip=ifelse(meta_df_prot_bytype$Chip=="chip", "chip_all", ifelse(meta_df_prot_bytype$Chip=="vf_01", "chip_large", meta_df_prot_bytype$Chip))
meta_df_prot_bytype$Chip_cat=paste0(substr(meta_df_prot_bytype$Chip,1,3), meta_df_prot_bytype$TargetFullName)
meta_df_prot_bytype=meta_df_prot_bytype %>% group_by(Protein, Chip_cat) %>% slice(which.min(meta_p_fix))
meta_df_prot_bytype$Chip_use=sub("_.*", "", meta_df_prot_bytype$Chip) 
meta_df_prot_bytype=merge(meta_df_prot_bytype,map,by.x="Protein",by.y="SomaId",all.x=T)
meta_df_prot_bytype_nojak2=meta_df_prot_bytype%>% filter(Chip_use != "JAK2") %>% mutate(meta_fdr_fix_fourgroup = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_bytype_nojak2=meta_df_prot_bytype_nojak2%>%select(Protein, Chip, meta_fdr_fix_fourgroup)
meta_df_prot_bytype=meta_df_prot_bytype%>% ungroup()%>% mutate(meta_fdr_fix_group_all = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_bytype=merge(meta_df_prot_bytype, meta_df_prot_bytype_nojak2, by=c("Protein","Chip"), all.x=T)

write.table(meta_df_prot_bytype,"/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final0318_fivegroups_nopeer_meta_chip_proteomics_log2peer50_nochangemap.txt", sep="\t", quote=FALSE, row.names=FALSE)

meta_df_prot_nojak2=meta_df_prot%>% filter(Chip %in% c("vf_01", "DNMT3A_all", "TET2_all","ASXL1_all", "DNMT3A_large", "TET2_large","ASXL1_large","chip")) %>% mutate(meta_fdr_fix_eightgroup = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_nojak2=meta_df_prot_nojak2%>%select(Protein, Chip, meta_fdr_fix_eightgroup)
meta_df_prot=meta_df_prot%>% ungroup()%>% mutate(meta_fdr_fix_group_all = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot=merge(meta_df_prot, meta_df_prot_nojak2, by=c("Protein","Chip"), all.x=T)
meta_df_prot=merge(meta_df_prot,map,by.x="Protein",by.y="SomaId",all.x=T)

write.table(meta_df_prot,"/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final0318_nopeer_meta_chip_proteomics_log2peer50_nochangemap.txt", sep="\t", quote=FALSE, row.names=FALSE)

###################Female only
variants_list=c("vf_01", "DNMT3A_all", "TET2_all","ASXL1_all", "DNMT3A_large", "TET2_large","ASXL1_large","chip")

meta_df_prot <- data.frame()
for (i in 1:length(variants_list)) {
  	name = variants_list[i]
  	print(name)
  	df1 = jhs_female_chip_proteomics[which(jhs_female_chip_proteomics$Chip == paste0(name)),]
  	df2 = mesa_female_chip_proteomics[which(mesa_female_chip_proteomics$Chip == paste0(name)),]
  	df3 = chs_female_chip_proteomics[which(chs_female_chip_proteomics$Chip == paste0(name)),]
	df4 = aric_ea_female_chip_proteomics[which(aric_ea_female_chip_proteomics$Chip == paste0(name)),]
	df5 = aric_aa_female_chip_proteomics[which(aric_aa_female_chip_proteomics$Chip == paste0(name)),]
  
  	df <- merge(df1,df2,by=c("Protein","Chip"))
  	df <- merge(df,df3,by=c("Protein","Chip"))
  	df <- merge(df,df4,by=c("Protein","Chip"))
  	df <- merge(df,df5,by=c("Protein","Chip"))

  	for (j in 1:nrow(df)) {
    		df$meta_beta_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$summary
    		df$meta_se_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$se
    		df$meta_p_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$test[2]
    		df$meta_het_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$het[3]

    		df$meta_beta_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$summary
    		df$meta_se_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$se
    		df$meta_p_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$test[2]
    		df$meta_het_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$het[3]
  	}
  	#str(df)
  	meta_df_prot <- rbind(meta_df_prot,df)
}

variants_list=c("JAK2_all", "JAK2_large")

jak2_df_prot <- data.frame()
for (i in 1:length(variants_list)) {
  	name = variants_list[i]
  	print(name)
  	df1 = chs_female_chip_proteomics[which(chs_female_chip_proteomics$Chip == paste0(name)),]
	df2 = aric_ea_female_chip_proteomics[which(aric_ea_female_chip_proteomics$Chip == paste0(name)),]

  	df <- merge(df1,df2,by=c("Protein","Chip"))

  	for (j in 1:nrow(df)) {
    		df$meta_beta_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$summary
    		df$meta_se_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$se
    		df$meta_p_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$test[2]
    		df$meta_het_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$het[3]

    		df$meta_beta_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$summary
    		df$meta_se_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$se
    		df$meta_p_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$test[2]
    		df$meta_het_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$het[3]
	}
  	#str(df)
  	jak2_df_prot <- rbind(jak2_df_prot,df)
}

meta_df_prot=bind_rows(meta_df_prot, jak2_df_prot)
meta_df_prot$meta_p_fix=ifelse(meta_df_prot$meta_p_fix==0, 1*10^-15, meta_df_prot$meta_p_fix)
meta_df_prot_bytype=meta_df_prot
meta_df_prot_bytype$Chip=ifelse(meta_df_prot_bytype$Chip=="chip", "chip_all", ifelse(meta_df_prot_bytype$Chip=="vf_01", "chip_large", meta_df_prot_bytype$Chip))
meta_df_prot_bytype$Chip_cat=paste0(substr(meta_df_prot_bytype$Chip,1,3), meta_df_prot_bytype$TargetFullName)
meta_df_prot_bytype=meta_df_prot_bytype %>% group_by(Protein, Chip_cat) %>% slice(which.min(meta_p_fix))
meta_df_prot_bytype$Chip_use=sub("_.*", "", meta_df_prot_bytype$Chip) 
meta_df_prot_bytype=merge(meta_df_prot_bytype,map,by.x="Protein",by.y="SomaId",all.x=T)
meta_df_prot_bytype_nojak2=meta_df_prot_bytype%>% filter(Chip_use != "JAK2") %>% mutate(meta_fdr_fix_fourgroup = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_bytype_nojak2=meta_df_prot_bytype_nojak2%>%select(Protein, Chip, meta_fdr_fix_fourgroup)
meta_df_prot_bytype=meta_df_prot_bytype%>% ungroup()%>% mutate(meta_fdr_fix_group_all = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_bytype=merge(meta_df_prot_bytype, meta_df_prot_bytype_nojak2, by=c("Protein","Chip"), all.x=T)

write.table(meta_df_prot_bytype,"/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final0318_femaleonly_fivegroups_meta_chip_proteomics_log2peer50_nochangemap.txt", sep="\t", quote=FALSE, row.names=FALSE)

meta_df_prot_nojak2=meta_df_prot%>% filter(Chip %in% c("vf_01", "DNMT3A_all", "TET2_all","ASXL1_all", "DNMT3A_large", "TET2_large","ASXL1_large","chip")) %>% mutate(meta_fdr_fix_eightgroup = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_nojak2=meta_df_prot_nojak2%>%select(Protein, Chip, meta_fdr_fix_eightgroup)
meta_df_prot=meta_df_prot%>% ungroup()%>% mutate(meta_fdr_fix_group_all = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot=merge(meta_df_prot, meta_df_prot_nojak2, by=c("Protein","Chip"), all.x=T)
meta_df_prot=merge(meta_df_prot,map,by.x="Protein",by.y="SomaId",all.x=T)

write.table(meta_df_prot,"/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final0318_femaleonly_meta_chip_proteomics_log2peer50_nochangemap.txt", sep="\t", quote=FALSE, row.names=FALSE)


###################Male only
variants_list=c("vf_01", "DNMT3A_all", "TET2_all","ASXL1_all", "DNMT3A_large", "TET2_large","ASXL1_large","chip")

meta_df_prot <- data.frame()
for (i in 1:length(variants_list)) {
  	name = variants_list[i]
  	print(name)
  	df1 = jhs_male_chip_proteomics[which(jhs_male_chip_proteomics$Chip == paste0(name)),]
  	df2 = mesa_male_chip_proteomics[which(mesa_male_chip_proteomics$Chip == paste0(name)),]
  	df3 = chs_male_chip_proteomics[which(chs_male_chip_proteomics$Chip == paste0(name)),]
	df4 = aric_ea_male_chip_proteomics[which(aric_ea_male_chip_proteomics$Chip == paste0(name)),]
	df5 = aric_aa_male_chip_proteomics[which(aric_aa_male_chip_proteomics$Chip == paste0(name)),]
  
  	df <- merge(df1,df2,by=c("Protein","Chip"))
  	df <- merge(df,df3,by=c("Protein","Chip"))
  	df <- merge(df,df4,by=c("Protein","Chip"))
  	df <- merge(df,df5,by=c("Protein","Chip"))

  	for (j in 1:nrow(df)) {
    		df$meta_beta_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$summary
    		df$meta_se_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$se
    		df$meta_p_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$test[2]
    		df$meta_het_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$het[3]

    		df$meta_beta_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$summary
    		df$meta_se_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$se
    		df$meta_p_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$test[2]
    		df$meta_het_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$het[3]
  	}
  	#str(df)
  	meta_df_prot <- rbind(meta_df_prot,df)
}

variants_list=c("JAK2_all", "JAK2_large")

jak2_df_prot <- data.frame()
for (i in 1:length(variants_list)) {
  	name = variants_list[i]
  	print(name)
  	df1 = chs_male_chip_proteomics[which(chs_male_chip_proteomics$Chip == paste0(name)),]
	df2 = aric_ea_male_chip_proteomics[which(aric_ea_male_chip_proteomics$Chip == paste0(name)),]

  	df <- merge(df1,df2,by=c("Protein","Chip"))

  	for (j in 1:nrow(df)) {
    		df$meta_beta_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$summary
    		df$meta_se_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$se
    		df$meta_p_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$test[2]
    		df$meta_het_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$het[3]

    		df$meta_beta_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$summary
    		df$meta_se_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$se
    		df$meta_p_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$test[2]
    		df$meta_het_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$het[3]
  	}
  	#str(df)
  	jak2_df_prot <- rbind(jak2_df_prot,df)
}

meta_df_prot=bind_rows(meta_df_prot, jak2_df_prot)
meta_df_prot$meta_p_fix=ifelse(meta_df_prot$meta_p_fix==0, 1*10^-15, meta_df_prot$meta_p_fix)
meta_df_prot_bytype=meta_df_prot
meta_df_prot_bytype$Chip=ifelse(meta_df_prot_bytype$Chip=="chip", "chip_all", ifelse(meta_df_prot_bytype$Chip=="vf_01", "chip_large", meta_df_prot_bytype$Chip))
meta_df_prot_bytype$Chip_cat=paste0(substr(meta_df_prot_bytype$Chip,1,3), meta_df_prot_bytype$TargetFullName)
meta_df_prot_bytype=meta_df_prot_bytype %>% group_by(Protein, Chip_cat) %>% slice(which.min(meta_p_fix))
meta_df_prot_bytype$Chip_use=sub("_.*", "", meta_df_prot_bytype$Chip) 
meta_df_prot_bytype=merge(meta_df_prot_bytype,map,by.x="Protein",by.y="SomaId",all.x=T)
meta_df_prot_bytype_nojak2=meta_df_prot_bytype%>% filter(Chip_use != "JAK2") %>% mutate(meta_fdr_fix_fourgroup = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_bytype_nojak2=meta_df_prot_bytype_nojak2%>%select(Protein, Chip, meta_fdr_fix_fourgroup)
meta_df_prot_bytype=meta_df_prot_bytype%>% ungroup()%>% mutate(meta_fdr_fix_group_all = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_bytype=merge(meta_df_prot_bytype, meta_df_prot_bytype_nojak2, by=c("Protein","Chip"), all.x=T)

write.table(meta_df_prot_bytype,"/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final0318_maleonly_fivegroups_meta_chip_proteomics_log2peer50_nochangemap.txt", sep="\t", quote=FALSE, row.names=FALSE)

meta_df_prot_nojak2=meta_df_prot%>% filter(Chip %in% c("vf_01", "DNMT3A_all", "TET2_all","ASXL1_all", "DNMT3A_large", "TET2_large","ASXL1_large","chip")) %>% mutate(meta_fdr_fix_eightgroup = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_nojak2=meta_df_prot_nojak2%>%select(Protein, Chip, meta_fdr_fix_eightgroup)
meta_df_prot=meta_df_prot%>% ungroup()%>% mutate(meta_fdr_fix_group_all = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot=merge(meta_df_prot, meta_df_prot_nojak2, by=c("Protein","Chip"), all.x=T)
meta_df_prot=merge(meta_df_prot,map,by.x="Protein",by.y="SomaId",all.x=T)

write.table(meta_df_prot,"/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final0318_maleonly_meta_chip_proteomics_log2peer50_nochangemap.txt", sep="\t", quote=FALSE, row.names=FALSE)




###################black only
variants_list=c("vf_01", "DNMT3A_all", "TET2_all","ASXL1_all", "DNMT3A_large", "TET2_large","ASXL1_large","chip")

meta_df_prot <- data.frame()
for (i in 1:length(variants_list)) {
  	name = variants_list[i]
  	print(name)
  	df1 = jhs_chip_proteomics[which(jhs_chip_proteomics$Chip == paste0(name)),]
  	df2 = mesa_black_chip_proteomics[which(mesa_black_chip_proteomics$Chip == paste0(name)),]
  	df3 = chs_black_chip_proteomics[which(chs_black_chip_proteomics$Chip == paste0(name)),]
  	df4 = aric_aa_chip_proteomics[which(aric_aa_chip_proteomics$Chip == paste0(name)),]
  
  	df <- merge(df1,df2,by=c("Protein","Chip"))
  	df <- merge(df,df3,by=c("Protein","Chip"))
  	df <- merge(df,df4,by=c("Protein","Chip"))

  	for (j in 1:nrow(df)) {
    		df$meta_beta_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$summary
    		df$meta_se_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$se
    		df$meta_p_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$test[2]
    		df$meta_het_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$het[3]

    		df$meta_beta_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$summary
    		df$meta_se_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$se
    		df$meta_p_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$test[2]
    		df$meta_het_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$het[3]
  	}
  	#str(df)
  	meta_df_prot <- rbind(meta_df_prot,df)
}

variants_list=c("JAK2_all", "JAK2_large")

jak2_df_prot <- data.frame()
for (i in 1:length(variants_list)) {
  	name = variants_list[i]
  	print(name)
  	df = chs_black_chip_proteomics[which(chs_black_chip_proteomics$Chip == paste0(name)),]

  	for (j in 1:nrow(df)) {
    		df$meta_beta_random[j] <- meta.summaries(d=c(df$chs_estimate[j]),se=c(df$chs_se[j]), method = "random", logscale=FALSE)$summary
    		df$meta_se_random[j] <- meta.summaries(d=c(df$chs_estimate[j]),se=c(df$chs_se[j]), method = "random", logscale=FALSE)$se
    		df$meta_p_random[j] <- meta.summaries(d=c(df$chs_estimate[j]),se=c(df$chs_se[j]), method = "random", logscale=FALSE)$test[2]
    		df$meta_het_random[j] <- meta.summaries(d=c(df$chs_estimate[j]),se=c(df$chs_se[j]), method = "random", logscale=FALSE)$het[3]

    		df$meta_beta_fix[j] <- meta.summaries(d=c(df$chs_estimate[j]),se=c(df$chs_se[j]), method = "fixed", logscale=FALSE)$summary
    		df$meta_se_fix[j] <- meta.summaries(d=c(df$chs_estimate[j]),se=c(df$chs_se[j]), method = "fixed", logscale=FALSE)$se
    		df$meta_p_fix[j] <- meta.summaries(d=c(df$chs_estimate[j]),se=c(df$chs_se[j]), method = "fixed", logscale=FALSE)$test[2]
    		df$meta_het_fix[j] <- meta.summaries(d=c(df$chs_estimate[j]),se=c(df$chs_se[j]), method = "fixed", logscale=FALSE)$het[3]
  	}
  	#str(df)
  	jak2_df_prot <- rbind(jak2_df_prot,df)
}

meta_df_prot=bind_rows(meta_df_prot, jak2_df_prot)
meta_df_prot$meta_p_fix=ifelse(meta_df_prot$meta_p_fix==0, 1*10^-15, meta_df_prot$meta_p_fix)
meta_df_prot_bytype=meta_df_prot
meta_df_prot_bytype$Chip=ifelse(meta_df_prot_bytype$Chip=="chip", "chip_all", ifelse(meta_df_prot_bytype$Chip=="vf_01", "chip_large", meta_df_prot_bytype$Chip))
meta_df_prot_bytype$Chip_cat=paste0(substr(meta_df_prot_bytype$Chip,1,3), meta_df_prot_bytype$TargetFullName)
meta_df_prot_bytype=meta_df_prot_bytype %>% group_by(Protein, Chip_cat) %>% slice(which.min(meta_p_fix))
meta_df_prot_bytype$Chip_use=sub("_.*", "", meta_df_prot_bytype$Chip) 
meta_df_prot_bytype=merge(meta_df_prot_bytype,map,by.x="Protein",by.y="SomaId",all.x=T)
meta_df_prot_bytype_nojak2=meta_df_prot_bytype%>% filter(Chip_use != "JAK2") %>% mutate(meta_fdr_fix_fourgroup = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_bytype_nojak2=meta_df_prot_bytype_nojak2%>%select(Protein, Chip, meta_fdr_fix_fourgroup)
meta_df_prot_bytype=meta_df_prot_bytype%>% ungroup()%>% mutate(meta_fdr_fix_group_all = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_bytype=merge(meta_df_prot_bytype, meta_df_prot_bytype_nojak2, by=c("Protein","Chip"), all.x=T)

write.table(meta_df_prot_bytype,"/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final0318_blackonly_fivegroups_meta_chip_proteomics_log2peer50_nochangemap.txt", sep="\t", quote=FALSE, row.names=FALSE)

meta_df_prot_nojak2=meta_df_prot%>% filter(Chip %in% c("vf_01", "DNMT3A_all", "TET2_all","ASXL1_all", "DNMT3A_large", "TET2_large","ASXL1_large","chip")) %>% mutate(meta_fdr_fix_eightgroup = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_nojak2=meta_df_prot_nojak2%>%select(Protein, Chip, meta_fdr_fix_eightgroup)
meta_df_prot=meta_df_prot%>% ungroup()%>% mutate(meta_fdr_fix_group_all = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot=merge(meta_df_prot, meta_df_prot_nojak2, by=c("Protein","Chip"), all.x=T)
meta_df_prot=merge(meta_df_prot,map,by.x="Protein",by.y="SomaId",all.x=T)

write.table(meta_df_prot,"/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final0318_blackonly_meta_chip_proteomics_log2peer50_nochangemap.txt", sep="\t", quote=FALSE, row.names=FALSE)


###################white only
variants_list=c("vf_01", "DNMT3A_all", "TET2_all","ASXL1_all", "DNMT3A_large", "TET2_large","ASXL1_large","chip")

meta_df_prot <- data.frame()
for (i in 1:length(variants_list)) {
  	name = variants_list[i]
  	print(name)
  	df1 = mesa_white_chip_proteomics[which(mesa_white_chip_proteomics$Chip == paste0(name)),]
  	df2 = chs_white_chip_proteomics[which(chs_white_chip_proteomics$Chip == paste0(name)),]
  	df3 = aric_ea_chip_proteomics[which(aric_ea_chip_proteomics$Chip == paste0(name)),]
  
  	df <- merge(df1,df2,by=c("Protein","Chip"))
  	df <- merge(df,df3,by=c("Protein","Chip"))

  	for (j in 1:nrow(df)) {
    		df$meta_beta_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$summary
    		df$meta_se_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$se
    		df$meta_p_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$test[2]
    		df$meta_het_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$het[3]

    		df$meta_beta_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$summary
    		df$meta_se_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$se
    		df$meta_p_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$test[2]
    		df$meta_het_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$het[3]
  	}
  	#str(df)
  	meta_df_prot <- rbind(meta_df_prot,df)
}

variants_list=c("JAK2_all", "JAK2_large")

jak2_df_prot <- data.frame()
for (i in 1:length(variants_list)) {
  	name = variants_list[i]
  	print(name)
  	df1 = chs_white_chip_proteomics[which(chs_white_chip_proteomics$Chip == paste0(name)),]
	df2 = aric_ea_chip_proteomics[which(aric_ea_chip_proteomics$Chip == paste0(name)),]

  	df <- merge(df1,df2,by=c("Protein","Chip"))

  	for (j in 1:nrow(df)) {
    		df$meta_beta_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$summary
    		df$meta_se_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$se
    		df$meta_p_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$test[2]
    		df$meta_het_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$het[3]

    		df$meta_beta_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$summary
    		df$meta_se_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$se
    		df$meta_p_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$test[2]
    		df$meta_het_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$het[3]
  	}
  	#str(df)
  	jak2_df_prot <- rbind(jak2_df_prot,df)
}

meta_df_prot=bind_rows(meta_df_prot, jak2_df_prot)
meta_df_prot$meta_p_fix=ifelse(meta_df_prot$meta_p_fix==0, 1*10^-15, meta_df_prot$meta_p_fix)
meta_df_prot_bytype=meta_df_prot
meta_df_prot_bytype$Chip=ifelse(meta_df_prot_bytype$Chip=="chip", "chip_all", ifelse(meta_df_prot_bytype$Chip=="vf_01", "chip_large", meta_df_prot_bytype$Chip))
meta_df_prot_bytype$Chip_cat=paste0(substr(meta_df_prot_bytype$Chip,1,3), meta_df_prot_bytype$TargetFullName)
meta_df_prot_bytype=meta_df_prot_bytype %>% group_by(Protein, Chip_cat) %>% slice(which.min(meta_p_fix))
meta_df_prot_bytype$Chip_use=sub("_.*", "", meta_df_prot_bytype$Chip) 
meta_df_prot_bytype=merge(meta_df_prot_bytype,map,by.x="Protein",by.y="SomaId",all.x=T)
meta_df_prot_bytype_nojak2=meta_df_prot_bytype%>% filter(Chip_use != "JAK2") %>% mutate(meta_fdr_fix_fourgroup = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_bytype_nojak2=meta_df_prot_bytype_nojak2%>%select(Protein, Chip, meta_fdr_fix_fourgroup)
meta_df_prot_bytype=meta_df_prot_bytype%>% ungroup()%>% mutate(meta_fdr_fix_group_all = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_bytype=merge(meta_df_prot_bytype, meta_df_prot_bytype_nojak2, by=c("Protein","Chip"), all.x=T)

write.table(meta_df_prot_bytype,"/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final0318_whiteonly_fivegroups_meta_chip_proteomics_log2peer50_nochangemap.txt", sep="\t", quote=FALSE, row.names=FALSE)

meta_df_prot_nojak2=meta_df_prot%>% filter(Chip %in% c("vf_01", "DNMT3A_all", "TET2_all","ASXL1_all", "DNMT3A_large", "TET2_large","ASXL1_large","chip")) %>% mutate(meta_fdr_fix_eightgroup = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_nojak2=meta_df_prot_nojak2%>%select(Protein, Chip, meta_fdr_fix_eightgroup)
meta_df_prot=meta_df_prot%>% ungroup()%>% mutate(meta_fdr_fix_group_all = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot=merge(meta_df_prot, meta_df_prot_nojak2, by=c("Protein","Chip"), all.x=T)
meta_df_prot=merge(meta_df_prot,map,by.x="Protein",by.y="SomaId",all.x=T)

write.table(meta_df_prot,"/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final0318_whiteonly_meta_chip_proteomics_log2peer50_nochangemap.txt", sep="\t", quote=FALSE, row.names=FALSE)

######################################################################################################
###########################after this I ddin't change mapping information###################################
######################################################################################################

#####interaction################
variants_list=c("vf_01", "DNMT3A_all", "TET2_all","ASXL1_all", "DNMT3A_large", "TET2_large","ASXL1_large","chip")

meta_df_prot <- data.frame()
for (i in 1:length(variants_list)) {
  	name = variants_list[i]
  	print(name)
  	df1 = jhs_interaction_chip_proteomics[which(jhs_interaction_chip_proteomics$Chip == paste0(name)),]
  	df2 = mesa_interaction_chip_proteomics[which(mesa_interaction_chip_proteomics$Chip == paste0(name)),]
  	df3 = chs_interaction_chip_proteomics[which(chs_interaction_chip_proteomics$Chip == paste0(name)),]
  	df4 = aric_interaction_ea_chip_proteomics[which(aric_interaction_ea_chip_proteomics$Chip == paste0(name)),]
  	df5 = aric_interaction_aa_chip_proteomics[which(aric_interaction_aa_chip_proteomics$Chip == paste0(name)),]
  
  	df <- merge(df1,df2,by=c("Protein","Chip"))
  	df <- merge(df,df3,by=c("Protein","Chip"))
  	df <- merge(df,df4,by=c("Protein","Chip"))
  	df <- merge(df,df5,by=c("Protein","Chip"))

  	for (j in 1:nrow(df)) {
    		df$meta_beta_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$summary
    		df$meta_se_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$se
    		df$meta_p_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$test[2]
    		df$meta_het_random[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "random", logscale=FALSE)$het[3]

    		df$meta_beta_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$summary
    		df$meta_se_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$se
    		df$meta_p_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$test[2]
    		df$meta_het_fix[j] <- meta.summaries(d=c(df$jhs_estimate[j],df$mesa_estimate[j],df$chs_estimate[j],df$aric_ea_estimate[j],df$aric_aa_estimate[j]),se=c(df$jhs_se[j],df$mesa_se[j],df$chs_se[j], df$aric_ea_se[j],df$aric_aa_se[j]), method = "fixed", logscale=FALSE)$het[3]
  	}
  	#str(df)
  	meta_df_prot <- rbind(meta_df_prot,df)
}

variants_list=c("JAK2_all", "JAK2_large")

jak2_df_prot <- data.frame()
for (i in 1:length(variants_list)) {
  	name = variants_list[i]
  	print(name)
  	df1 = chs_interaction_chip_proteomics[which(chs_interaction_chip_proteomics$Chip == paste0(name)),]
  	df2 = aric_interaction_ea_chip_proteomics[which(aric_interaction_ea_chip_proteomics$Chip == paste0(name)),]
  
  	df <- merge(df1,df2,by=c("Protein","Chip"))

  	for (j in 1:nrow(df)) {
    		df$meta_beta_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$summary
    		df$meta_se_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$se
    		df$meta_p_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$test[2]
    		df$meta_het_random[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "random", logscale=FALSE)$het[3]

    		df$meta_beta_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$summary
    		df$meta_se_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$se
    		df$meta_p_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$test[2]
    		df$meta_het_fix[j] <- meta.summaries(d=c(df$chs_estimate[j],df$aric_ea_estimate[j]),se=c(df$chs_se[j], df$aric_ea_se[j]), method = "fixed", logscale=FALSE)$het[3]
  	}
  	#str(df)
  	jak2_df_prot <- rbind(jak2_df_prot,df)
}

meta_df_prot=bind_rows(meta_df_prot, jak2_df_prot)
#meta_df_prot$meta_p_fix=ifelse(meta_df_prot$meta_p_fix==0, 1*10^-15, meta_df_prot$meta_p_fix)
meta_df_prot_bytype=meta_df_prot
meta_df_prot_bytype$Chip=ifelse(meta_df_prot_bytype$Chip=="chip", "chip_all", ifelse(meta_df_prot_bytype$Chip=="vf_01", "chip_large", meta_df_prot_bytype$Chip))
meta_df_prot_bytype$Chip_cat=paste0(substr(meta_df_prot_bytype$Chip,1,3), meta_df_prot_bytype$TargetFullName)
meta_df_prot_bytype=meta_df_prot_bytype %>% group_by(Protein, Chip_cat) %>% slice(which.min(meta_p_fix))
meta_df_prot_bytype$Chip_use=sub("_.*", "", meta_df_prot_bytype$Chip) 
meta_df_prot_bytype=merge(meta_df_prot_bytype,map,by.x="Protein",by.y="SomaId",all.x=T)
meta_df_prot_bytype_nojak2=meta_df_prot_bytype%>% filter(Chip_use != "JAK2") %>% mutate(meta_fdr_fix_fourgroup = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_bytype_nojak2=meta_df_prot_bytype_nojak2%>%select(Protein, Chip, meta_fdr_fix_fourgroup)
meta_df_prot_bytype=meta_df_prot_bytype%>% ungroup()%>% mutate(meta_fdr_fix_group_all = p.adjust(meta_p_fix, method = "fdr"))
meta_df_prot_bytype=merge(meta_df_prot_bytype, meta_df_prot_bytype_nojak2, by=c("Protein","Chip"), all.x=T)
meta_df_prot_bytype=merge(meta_df_prot_bytype,map,by.x="Protein",by.y="SomaId",all.x=T)

write.table(meta_df_prot_bytype,"/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Finalfinal0624_interaction_fivegroups_meta_chip_proteomics_log2.txt", sep="\t", quote=FALSE, row.names=FALSE)

##match with significant ones in male or female
male=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final0318_maleonly_fivegroups_meta_chip_proteomics_log2peer50.txt")
female=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final0318_femaleonly_fivegroups_meta_chip_proteomics_log2peer50.txt")
male=male%>%filter(meta_fdr_fix_fourgroup<0.05)
male$indicator=paste(male$Protein, male$Chip_use,sep="_")
female=female%>%filter(meta_fdr_fix_fourgroup<0.05)
female$indicator=paste(female$Protein, female$Chip_use,sep="_")

male_diff <- setdiff(male$indicator, female$indicator)
female_diff <- setdiff(female$indicator, male$indicator)

# Combine the results
all_diff <- c(male_diff, female_diff)

meta_df_prot_bytype$indicator=paste(meta_df_prot_bytype$Protein, meta_df_prot_bytype$Chip_use,sep="_")
meta_df_prot_bytype_sub=meta_df_prot_bytype[which(meta_df_prot_bytype$indicator %in% all_diff ),]
meta_df_prot_bytype_sub=meta_df_prot_bytype_sub%>% ungroup()%>% mutate(meta_fdr_fix_fourgroup = p.adjust(meta_p_fix, method = "fdr"), meta_fdr_fix_group_all=NA) 

write.table(meta_df_prot_bytype_sub,"/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Finalfinal0624_sub_interaction_fivegroups_meta_chip_proteomics_log2.txt", sep="\t", quote=FALSE, row.names=FALSE)

###race
##match with significant ones in white or black
white=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final0318_whiteonly_fivegroups_meta_chip_proteomics_log2peer50.txt")
black=fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final0318_blackonly_fivegroups_meta_chip_proteomics_log2peer50.txt")
white=white%>%filter(meta_fdr_fix_fourgroup<0.05)
white$indicator=paste(white$Protein, white$Chip_use,sep="_")
black=black%>%filter(meta_fdr_fix_fourgroup<0.05)
black$indicator=paste(black$Protein, black$Chip_use,sep="_")

white_diff <- setdiff(white$indicator, black$indicator)
black_diff <- setdiff(black$indicator, white$indicator)

# Combine the results
all_diff <- c(white_diff, black_diff)

meta_df_prot_bytype=aric_interaction_race_chip_proteomics
meta_df_prot_bytype$Chip=ifelse(meta_df_prot_bytype$Chip=="chip", "chip_all", ifelse(meta_df_prot_bytype$Chip=="vf_01", "chip_large", meta_df_prot_bytype$Chip))
meta_df_prot_bytype$Chip_cat=paste0(substr(meta_df_prot_bytype$Chip,1,3), meta_df_prot_bytype$TargetFullName)
meta_df_prot_bytype=meta_df_prot_bytype %>% group_by(Protein, Chip_cat) %>% slice(which.min(aric_race_p))
meta_df_prot_bytype$Chip_use=sub("_.*", "", meta_df_prot_bytype$Chip) 
meta_df_prot_bytype=merge(meta_df_prot_bytype,map,by.x="Protein",by.y="SomaId",all.x=T)
meta_df_prot_bytype=meta_df_prot_bytype%>% filter(Chip_use != "JAK2") 

meta_df_prot_bytype$indicator=paste(meta_df_prot_bytype$Protein, meta_df_prot_bytype$Chip_use,sep="_")
meta_df_prot_bytype_sub=meta_df_prot_bytype[which(meta_df_prot_bytype$indicator %in% all_diff ),]
meta_df_prot_bytype_sub=meta_df_prot_bytype_sub%>% ungroup()%>% mutate(fdr = p.adjust(aric_race_p, method = "fdr")) 

write.table(meta_df_prot_bytype_sub,"/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Finalfinal0624_sub_raceinteraction_meta_chip_proteomics_log2.txt", sep="\t", quote=FALSE, row.names=FALSE)
