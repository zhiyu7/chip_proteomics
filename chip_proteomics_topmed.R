# Clean the environment
gc()
rm(list = ls())

# Load required libraries
library(dplyr)
library(data.table)
library(rmeta)
library(tableone)
library(survival)

# Load data files
TOPMed_chip_old <- fread("/medpop/esp2/zyu/chip_protemoics/data/TOPMed_chip/TOPMed_CHIP_Lipid_CAD.tsv")
TOPMed_chip_old <- TOPMed_chip_old %>%
  filter(study %in% c("JHS", "MESA", "CHS", "FHS", "WHI")) %>%
  select(Sample, study, Age, sex, BMI, Race, T2D, Current_smoker_baseline, Ever_smoker_baseline, contains("CAD"), contains(".norm"))

TOPMed_chip_new <- fread("/medpop/esp2/zyu/chip_protemoics/data/TOPMed_chip/012020_chip/TOPMed_100k_CHIP_variants_1_20_20_fordbgap.csv")
TOPMed_chip_new$hasCHIP <- TRUE

# Merge old and new datasets
TOPMed_chip <- merge(TOPMed_chip_old, TOPMed_chip_new, by = "Sample", all.x = TRUE)

# Extract gene-specific subsets
TET2 <- TOPMed_chip[TOPMed_chip$Gene == "TET2",]
DNMT3A <- TOPMed_chip[TOPMed_chip$Gene == "DNMT3A",]
ASXL1 <- TOPMed_chip[TOPMed_chip$Gene == "ASXL1",]
JAK2 <- TOPMed_chip[TOPMed_chip$Gene == "JAK2",]
TP53 <- TOPMed_chip[TOPMed_chip$Gene == "TP53",]

# Filter by variant allele frequency (VAF > 0.1)
vf_01 <- TOPMed_chip[TOPMed_chip$VAF > 0.1,]
chip <- TOPMed_chip[TOPMed_chip$hasCHIP == TRUE,]

# Load PCA data
pc <- fread("/medpop/esp2/projects/topmed/freeze10/freeze.10.pca.txt.gz")
pc <- pc %>%
  select(Sample = IID, PC1 = PC1_AVG, PC2 = PC2_AVG, PC3 = PC3_AVG, PC4 = PC4_AVG, PC5 = PC5_AVG, PC6 = PC6_AVG, PC7 = PC7_AVG, PC8 = PC8_AVG, PC9 = PC9_AVG, PC10 = PC10_AVG)

# Load annotation data
anno <- fread("/medpop/esp2/zyu/chip_protemoics/data/JHS_protein_ID_key.csv")

##########################################################
######################## JHS Data #######################
##########################################################

# Filter JHS data from TOPMed_chip
jhs <- TOPMed_chip[TOPMed_chip$study == "JHS",] %>%
  select(-c(CHROM:F2R1)) %>%
  distinct()

# Create binary variables for genes and large variants
jhs$DNMT3A_all <- ifelse(jhs$Sample %in% DNMT3A$Sample, 1, 0)
jhs$TET2_all <- ifelse(jhs$Sample %in% TET2$Sample, 1, 0)
jhs$ASXL1_all <- ifelse(jhs$Sample %in% ASXL1$Sample, 1, 0)
jhs$JAK2_all <- ifelse(jhs$Sample %in% JAK2$Sample, 1, 0)
jhs$DNMT3A_large <- ifelse(jhs$Sample %in% DNMT3A$Sample & jhs$Sample %in% vf_01$Sample, 1, 0)
jhs$TET2_large <- ifelse(jhs$Sample %in% TET2$Sample & jhs$Sample %in% vf_01$Sample, 1, 0)
jhs$ASXL1_large <- ifelse(jhs$Sample %in% ASXL1$Sample & jhs$Sample %in% vf_01$Sample, 1, 0)
jhs$JAK2_large <- ifelse(jhs$Sample %in% JAK2$Sample & jhs$Sample %in% vf_01$Sample, 1, 0)
jhs$vf_01 <- ifelse(jhs$Sample %in% vf_01$Sample, 1, 0)
jhs$chip <- ifelse(jhs$Sample %in% chip$Sample, 1, 0)

# Select required columns
jhs <- jhs %>%
  select(DNMT3A_all, TET2_all, ASXL1_all, JAK2_all, DNMT3A_large, TET2_large, ASXL1_large, JAK2_large, vf_01, chip, Sample, Age, Race, sex, BMI, T2D, Current_smoker_baseline, Ever_smoker_baseline, t_CAD_all, i_CAD_all, i_CAD_all_date_censor, TOTAL_ADJ.norm, HDL.norm)

# Merge with PC and ID data
jhs <- merge(jhs, pc, by = "Sample")
id <- fread("/medpop/esp2/zyu/chip_protemoics/data/JHS/JHS_WGS_ID_map_03082017_NWDID.csv")
jhs <- merge(jhs, id, by.x = "Sample", by.y = "NWDID")

# Add batch information
batch <- fread("/medpop/esp2/zyu/chip_protemoics/data/JHS/JHS_proteomics_raw.csv") %>%
  select(subjid, prot_batch)
jhs <- merge(jhs, batch, by.x = "SUBJID", by.y = "subjid")

# Load peer-adjusted proteomics data
peer_prot <- fread("/medpop/esp2/zyu/chip_protemoics/data/JHS/PEER/peernum/peerandresid/prot.JHS_log2_nopeer.txt")
prot_peercovar <- fread("/medpop/esp2/zyu/chip_protemoics/data/JHS/PEER/peernum/peerandresid/peers.JHS_log2_50.txt")

# Merge with proteomics data
jhs_prot <- merge(jhs, peer_prot, by.x = "SUBJID", by.y = "SampleId")
jhs_prot <- merge(jhs_prot, prot_peercovar, by.x = "SUBJID", by.y = "SampleId")

# Summary of age in JHS proteomics data
summary(jhs_prot$Age)
jhs_sample <- jhs_prot$Sample

# Define factor variables and variable list for creating TableOne
factorVars <- c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "DNMT3A_large", 
                "TET2_large", "ASXL1_large", "JAK2_large", "vf_01", "chip", "sex", 
                "Race", "T2D", "Current_smoker_baseline", "Ever_smoker_baseline", 
                "t_CAD_all", "i_CAD_all")

vars <- c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "DNMT3A_large", 
          "TET2_large", "ASXL1_large", "JAK2_large", "vf_01", "chip", "Age", 
          "sex", "Race", "BMI", "T2D", "Current_smoker_baseline", 
          "Ever_smoker_baseline", "t_CAD_all", "i_CAD_all")

# Create TableOne and write to file
tableOne <- CreateTableOne(vars = vars, data = jhs_prot, factorVars = factorVars)
write.table(print(tableOne), "/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/jhs_tableone.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE)

# Create TableOne with stratification by sex and race for CHIP prevalence
factorVars_chip <- c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "chip")
vars_chip <- c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "chip")

jhschipprev <- CreateTableOne(vars = vars_chip, data = jhs_prot, factorVars = factorVars_chip, strata = c("sex", "Race"))

# List proteins and variants for analysis
proteins <- colnames(jhs_prot)[36:1340]
variants_list <- colnames(jhs_prot)[3:12]

########## CAD Proteomics Analysis ##########
jhs_prot$BMI <- ifelse(is.na(jhs_prot$BMI), median(jhs_prot$BMI, na.rm = TRUE), jhs_prot$BMI)
jhs_prot$Ever_smoker_baseline <- ifelse(is.na(jhs_prot$Ever_smoker_baseline), median(jhs_prot$Ever_smoker_baseline, na.rm = TRUE), jhs_prot$Ever_smoker_baseline)
jhs_prot$HDL.norm <- ifelse(is.na(jhs_prot$HDL.norm), median(jhs_prot$HDL.norm, na.rm = TRUE), jhs_prot$HDL.norm)

# Define model adjustments
ds <- c("jhs_prot")
outcome <- c("t_CAD_all")
Model <- c("+Age+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10", 
           "+Age+sex+prot_batch+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+T2D+BMI+Ever_smoker_baseline")

model <- c("Agesexracepc_adjust", "fully_adjust")

# Loop through datasets and models to perform logistic regression for CAD outcome
for (j in 1:length(ds)) {
    d <- get(ds[j])
    k <- 1
    for (m in 1:length(Model)) {
        h <- data.frame(seq_id = character(), N = double(), beta = double(), se = double(), 
                        z = double(), p = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)

        for (i in 1:length(proteins)) {
            h[i, "seq_id"] <- proteins[i]
            fmla <- as.formula(paste(outcome[k], "~", proteins[i], Model[m], sep = ""))
            f <- glm(fmla, data = d, na.action = na.exclude, family = "binomial")
            h[i, "N"] <- summary(f)$df.null + 1
            h[i, "beta"] <- summary(f)$coefficient[2, 1]
            h[i, "se"] <- summary(f)$coefficient[2, 2]
            h[i, "z"] <- summary(f)$coefficient[2, 3]
            h[i, "p"] <- summary(f)$coefficient[2, 4]
            h[i, "lci"] <- summary(f)$coefficient[2, 1] - 1.96 * summary(f)$coefficient[2, 2]
            h[i, "uci"] <- summary(f)$coefficient[2, 1] + 1.96 * summary(f)$coefficient[2, 2]
        }
        h$p_FDR <- p.adjust(h$p, method = "fdr")
        h <- merge(anno, h, by.x = "SomaId", by.y = "seq_id")
        h <- h[order(h$p),]

        # Write results to file
        write.table(h, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/JHS_proteomics_CAD_", 
                                    outcome[k], "_", model[m], ".txt", sep = ""), sep = "\t", 
                    quote = FALSE, row.names = FALSE)
    }
}

###### CHIP Proteomics Analysis ######
# Define models for CHIP analysis with and without adjustment
Model <- c("", 
           "+Age+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+T2D+BMI+Ever_smoker_baseline+prot_batch+peer1+peer2+peer3+peer4+peer5+peer6+peer7+peer8+peer9+peer10+peer11+
           peer12+peer13+peer14+peer15+peer16+peer17+peer18+peer19+peer20+peer21+peer22+peer23+peer24+peer25+peer26+peer27+peer28+peer29+peer30+peer31+peer32+peer33+peer34+peer35+
           peer36+peer37+peer38+peer39+peer40+peer41+peer42+peer43+peer44+peer45+peer46+peer47+peer48+peer49+peer50", 
           "+Age+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+T2D+BMI+Ever_smoker_baseline+prot_batch")

model <- c("noadjust", "fully_adjust", "nopeer")
ds <- c("jhs_prot")
dss <- c("")

# Loop through datasets and models for CHIP analysis
for (j in 1:length(dss)) {
    d <- get(ds[j])
    for (m in 1:length(Model)) {
        h <- data.frame(seq_id = character(), chip = character(), N = double(), beta = double(), se = double(),
                        z = double(), p = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)
        
        for (k in 1:length(variants_list)) {
            for (i in 1:length(proteins)) {
                h[i, "seq_id"] <- proteins[i]
                h[i, "chip"] <- variants_list[k]
                fmla <- as.formula(paste(proteins[i], "~", variants_list[k], Model[m], sep = ""))
                f <- glm(fmla, data = d, na.action = na.exclude, family = "gaussian")
                h[i, "N"] <- summary(f)$df.null + 1
                h[i, "beta"] <- summary(f)$coefficient[2, 1]
                h[i, "se"] <- summary(f)$coefficient[2, 2]
                h[i, "z"] <- summary(f)$coefficient[2, 3]
                h[i, "p"] <- summary(f)$coefficient[2, 4]
                h[i, "lci"] <- summary(f)$coefficient[2, 1] - 1.96 * summary(f)$coefficient[2, 2]
                h[i, "uci"] <- summary(f)$coefficient[2, 1] + 1.96 * summary(f)$coefficient[2, 2]
            }
            if (k == 1) {
                h_merge <- h
            } else {
                h_merge <- rbind(h_merge, h)
            }
        }
        
        h_merge <- merge(anno, h_merge, by.x = "SomaId", by.y = "seq_id")
        h_merge <- h_merge[order(h_merge$p),]
        
        # Save CHIP proteomics analysis results
        write.table(h_merge, 
                    file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_JHS_proteomics_CHIP_", model[m], ".txt", sep = ""), 
                    sep = "\t", quote = FALSE, row.names = FALSE)

        # Perform FDR adjustment for specific CHIP variants
        h_merge_fdr <- h_merge[h_merge$chip %in% c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "chip"),]
        h_merge_fdr$p_FDR <- p.adjust(h_merge_fdr$p, method = "fdr")

        # Save FDR adjusted results
        write.table(h_merge_fdr, 
                    file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_FDR_JHS_proteomics_CHIP_", model[m], ".txt", sep = ""), 
                    sep = "\t", quote = FALSE, row.names = FALSE)
    }
}

## Interaction of Sex
Model <- c("*sex+Age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+T2D+BMI+Ever_smoker_baseline+prot_batch+peer1+peer2+peer3+peer4+peer5+peer6+peer7+peer8+peer9+peer10+peer11+
           peer12+peer13+peer14+peer15+peer16+peer17+peer18+peer19+peer20+peer21+peer22+peer23+peer24+peer25+peer26+peer27+peer28+peer29+peer30+peer31+peer32+peer33+peer34+
           peer35+peer36+peer37+peer38+peer39+peer40+peer41+peer42+peer43+peer44+peer45+peer46+peer47+peer48+peer49+peer50")
model <- c("fully_adjust")
ds <- c("jhs_prot")
dss <- c("")

for (j in 1:length(dss)) {
    d <- get(ds[j])
    for (m in 1:length(Model)) {
        h <- data.frame(seq_id = character(), chip = character(), N = double(), beta = double(), se = double(), 
                        z = double(), p = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)

        for (k in 1:length(variants_list)) {
            for (i in 1:length(proteins)) {
                h[i, "seq_id"] <- proteins[i]
                h[i, "chip"] <- variants_list[k]
                fmla <- as.formula(paste(proteins[i], "~", variants_list[k], Model[m], sep = ""))
                f <- glm(fmla, data = d, na.action = na.exclude, family = "gaussian")

                # Check if the interaction term is NA, if so, fill in NA values
                if (nrow(summary(f)$coefficient) < 69 || is.na(summary(f)$coefficient[69, 1])) {
                    h[i, c("N", "beta", "se", "z", "p", "lci", "uci")] <- NA
                } else {
                    h[i, "N"] <- summary(f)$df.null + 1
                    h[i, "beta"] <- summary(f)$coefficient[69, 1]
                    h[i, "se"] <- summary(f)$coefficient[69, 2]
                    h[i, "z"] <- summary(f)$coefficient[69, 3]
                    h[i, "p"] <- summary(f)$coefficient[69, 4]
                    h[i, "lci"] <- summary(f)$coefficient[69, 1] - 1.96 * summary(f)$coefficient[69, 2]
                    h[i, "uci"] <- summary(f)$coefficient[69, 1] + 1.96 * summary(f)$coefficient[69, 2]
                }
            }
            if (k == 1) {
                h_merge <- h
            } else {
                h_merge <- rbind(h_merge, h)
            }
        }

        h_merge <- merge(anno, h_merge, by.x = "SomaId", by.y = "seq_id")
        h_merge <- h_merge[order(h_merge$p),]

        # Save results for interaction of sex
        write.table(h_merge, 
                    file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_JHS_proteomics_interaction_CHIP_", model[m], ".txt", sep = ""), 
                    sep = "\t", quote = FALSE, row.names = FALSE)
    }
}

####### Sex-specific Analysis #######
Model <- c("+Age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+T2D+BMI+Ever_smoker_baseline+prot_batch+peer1+peer2+peer3+peer4+
            peer5+peer6+peer7+peer8+peer9+peer10+peer11+peer12+peer13+peer14+peer15+peer16+peer17+peer18+peer19+peer20+peer21+
            peer22+peer23+peer24+peer25+peer26+peer27+peer28+peer29+peer30+peer31+peer32+peer33+peer34+peer35+peer36+peer37+
            peer38+peer39+peer40+peer41+peer42+peer43+peer44+peer45+peer46+peer47+peer48+peer49+peer50")

model <- c("fully_adjust")

# Subset male and female data
jhs_prot_male <- jhs_prot[jhs_prot$sex == "M", ]
jhs_prot_female <- jhs_prot[jhs_prot$sex == "F", ]

ds <- c("jhs_prot_male", "jhs_prot_female")
dss <- c("male", "female")

# Loop through male and female datasets for sex-specific analysis
for (j in 1:length(dss)) {
    d <- get(ds[j])
    for (m in 1:length(Model)) {
        h <- data.frame(seq_id = character(), chip = character(), N = double(), beta = double(), se = double(), 
                        z = double(), p = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)

        for (k in 1:length(variants_list)) {
            for (i in 1:length(proteins)) {
                h[i, "seq_id"] <- proteins[i]
                h[i, "chip"] <- variants_list[k]
                fmla <- as.formula(paste(proteins[i], "~", variants_list[k], Model[m], sep = ""))
                f <- glm(fmla, data = d, na.action = na.exclude, family = "gaussian")
                h[i, "N"] <- summary(f)$df.null + 1
                h[i, "beta"] <- summary(f)$coefficient[2, 1]
                h[i, "se"] <- summary(f)$coefficient[2, 2]
                h[i, "z"] <- summary(f)$coefficient[2, 3]
                h[i, "p"] <- summary(f)$coefficient[2, 4]
                h[i, "lci"] <- summary(f)$coefficient[2, 1] - 1.96 * summary(f)$coefficient[2, 2]
                h[i, "uci"] <- summary(f)$coefficient[2, 1] + 1.96 * summary(f)$coefficient[2, 2]
            }

            if (k == 1) {
                h_merge <- h
            } else {
                h_merge <- rbind(h_merge, h)
            }
        }

        h_merge <- merge(anno, h_merge, by.x = "SomaId", by.y = "seq_id")
        h_merge <- h_merge[order(h_merge$p),]

        # Save sex-specific analysis results
        write.table(h_merge, 
                    file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_JHS_proteomics_CHIP_", dss[j], "_", model[m], ".txt", sep = ""), 
                    sep = "\t", quote = FALSE, row.names = FALSE)

        # Perform FDR adjustment for sex-specific CHIP variants
        h_merge_fdr <- h_merge[h_merge$chip %in% c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "chip"), ]
        h_merge_fdr$p_FDR <- p.adjust(h_merge_fdr$p, method = "fdr")

        # Save FDR adjusted sex-specific results
        write.table(h_merge_fdr, 
                    file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_FDR_JHS_proteomics_CHIP_", dss[j], "_", model[m], ".txt", sep = ""), 
                    sep = "\t", quote = FALSE, row.names = FALSE)
    }
}

####### Race-specific Analysis #######
Model <- c("+Age+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+T2D+BMI+Ever_smoker_baseline+peer1+peer2+peer3+peer4+
            peer5+peer6+peer7+peer8+peer9+peer10+peer11+peer12+peer13+peer14+peer15+peer16+peer17+peer18+peer19+peer20+
            peer21+peer22+peer23+peer24+peer25+peer26+peer27+peer28+peer29+peer30+peer31+peer32+peer33+peer34+peer35+
            peer36+peer37+peer38+peer39+peer40+peer41+peer42+peer43+peer44+peer45+peer46+peer47+peer48+peer49+peer50")

model <- c("fully_adjust")

# Subset data by race
jhs_prot_white <- jhs_prot[jhs_prot$Race == "White", ]
jhs_prot_black <- jhs_prot[jhs_prot$Race == "Black", ]

ds <- c("jhs_prot_white", "jhs_prot_black")
dss <- c("white", "black")

# Loop through race-specific datasets for race-specific analysis
for (j in 1:length(dss)) {
    d <- get(ds[j])
    for (m in 1:length(Model)) {
        h <- data.frame(seq_id = character(), chip = character(), N = double(), beta = double(), se = double(), 
                        z = double(), p = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)

        for (k in 1:length(variants_list)) {
            for (i in 1:length(proteins)) {
                h[i, "seq_id"] <- proteins[i]
                h[i, "chip"] <- variants_list[k]
                fmla <- as.formula(paste(proteins[i], "~", variants_list[k], Model[m], sep = ""))
                f <- glm(fmla, data = d, na.action = na.exclude, family = "gaussian")
                h[i, "N"] <- summary(f)$df.null + 1
                h[i, "beta"] <- summary(f)$coefficient[2, 1]
                h[i, "se"] <- summary(f)$coefficient[2, 2]
                h[i, "z"] <- summary(f)$coefficient[2, 3]
                h[i, "p"] <- summary(f)$coefficient[2, 4]
                h[i, "lci"] <- summary(f)$coefficient[2, 1] - 1.96 * summary(f)$coefficient[2, 2]
                h[i, "uci"] <- summary(f)$coefficient[2, 1] + 1.96 * summary(f)$coefficient[2, 2]
            }

            if (k == 1) {
                h_merge <- h
            } else {
                h_merge <- rbind(h_merge, h)
            }
        }

        h_merge <- merge(anno, h_merge, by.x = "SomaId", by.y = "seq_id")
        h_merge <- h_merge[order(h_merge$p),]

        # Save race-specific analysis results
        write.table(h_merge, 
                    file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_JHS_proteomics_CHIP_", dss[j], "_", model[m], ".txt", sep = ""), 
                    sep = "\t", quote = FALSE, row.names = FALSE)

        # Perform FDR adjustment for race-specific CHIP variants
        h_merge_fdr <- h_merge[h_merge$chip %in% c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "chip"), ]
        h_merge_fdr$p_FDR <- p.adjust(h_merge_fdr$p, method = "fdr")

        # Save FDR adjusted race-specific results
        write.table(h_merge_fdr, 
                    file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_FDR_JHS_proteomics_CHIP_", dss[j], "_", model[m], ".txt", sep = ""), 
                    sep = "\t", quote = FALSE, row.names = FALSE)
    }
}

################## PEER Testing ##################
ds <- c("jhs_prot")
dss <- c("")
outcome <- c("t_CAD_all", variants_list)
Model <- c("peer1+peer2+peer3+peer4+peer5+peer6+peer7+peer8+peer9+peer10+peer11+peer12+peer13+peer14+peer15+peer16+peer17+
           peer18+peer19+peer20+peer21+peer22+peer23+peer24+peer25+peer26+peer27+peer28+peer29+peer30+peer31+peer32+peer33+
           peer34+peer35+peer36+peer37+peer38+peer39+peer40+peer41+peer42+peer43+peer44+peer45+peer46+peer47+peer48+peer49+
           peer50")

model <- c("peer")

for (j in 1:length(dss)) {
    d <- get(ds[j])
    for (m in 1:length(Model)) {
        h <- data.frame(outcome = character(), N = double(), n_sig = double(), min_sig = double(), stringsAsFactors = FALSE)

        for (i in 1:length(outcome)) {
            h[i, "outcome"] <- outcome[i]
            fmla <- as.formula(paste(outcome[i], "~", Model[m], sep = ""))
            f <- glm(fmla, data = d, na.action = na.exclude, family = "binomial")
            h[i, "N"] <- summary(f)$df.null + 1
            h[i, "n_sig"] <- sum(summary(f)$coefficient[2:51, 4] < 0.05)
            h[i, "min_sig"] <- min(summary(f)$coefficient[2:51, 4])
        }
    }
}

##########################################################
######################### MESA ###########################
##########################################################

# Filter MESA cohort and remove specific columns
mesa <- TOPMed_chip[TOPMed_chip$study == "MESA", ] %>% select(-c(CHROM:F2R1))

# Remove duplicates
mesa <- mesa[!duplicated(mesa)]

# Create mutation indicator variables for each gene
mesa$DNMT3A_all <- ifelse(mesa$Sample %in% DNMT3A$Sample, 1, 0)
mesa$TET2_all <- ifelse(mesa$Sample %in% TET2$Sample, 1, 0)
mesa$ASXL1_all <- ifelse(mesa$Sample %in% ASXL1$Sample, 1, 0)
mesa$JAK2_all <- ifelse(mesa$Sample %in% JAK2$Sample, 1, 0)

# Create large variant indicator variables
mesa$DNMT3A_large <- ifelse(mesa$Sample %in% DNMT3A$Sample & mesa$Sample %in% vf_01$Sample, 1, 0)
mesa$TET2_large <- ifelse(mesa$Sample %in% TET2$Sample & mesa$Sample %in% vf_01$Sample, 1, 0)
mesa$ASXL1_large <- ifelse(mesa$Sample %in% ASXL1$Sample & mesa$Sample %in% vf_01$Sample, 1, 0)
mesa$JAK2_large <- ifelse(mesa$Sample %in% JAK2$Sample & mesa$Sample %in% vf_01$Sample, 1, 0)

# Create additional variables
mesa$vf_01 <- ifelse(mesa$Sample %in% vf_01$Sample, 1, 0)
mesa$chip <- ifelse(mesa$Sample %in% chip$Sample, 1, 0)

# Select specific columns
mesa <- mesa[, c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "DNMT3A_large", "TET2_large", "ASXL1_large", 
                 "JAK2_large", "vf_01", "chip", "Sample", "sex", "Age", "Race", "BMI", "T2D", "Current_smoker_baseline", 
                 "Ever_smoker_baseline", "t_CAD_all", "i_CAD_all", "i_CAD_all_date_censor", "TOTAL_ADJ.norm", "HDL.norm")]

# Merge with principal component data
mesa <- merge(mesa, pc, by = "Sample")

# Load proteomics data
peer_prot <- fread("/medpop/esp2/zyu/chip_protemoics/data/MESA/PEER/peernum/peerandresid/prot.MESA_log2_nopeer.txt")
prot_peercovar <- fread("/medpop/esp2/zyu/chip_protemoics/data/MESA/PEER/peernum/peerandresid/peers.MESA_log2_50.txt")

# Merge proteomics data with cohort data
mesa_prot <- merge(mesa, peer_prot, by.x = "Sample", by.y = "SampleId")
mesa_prot <- merge(mesa_prot, prot_peercovar, by.x = "Sample", by.y = "SampleId")

# Summary of age
summary(mesa_prot$Age)

# Extract sample list
mesa_sample <- mesa_prot$Sample

# Define factor and numeric variables for TableOne
factorVars <- c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "DNMT3A_large", "TET2_large", 
                "ASXL1_large", "JAK2_large", "vf_01", "chip", "sex", "Race", "T2D", 
                "Current_smoker_baseline", "Ever_smoker_baseline", "t_CAD_all", "i_CAD_all")

vars <- c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "DNMT3A_large", "TET2_large", "ASXL1_large", 
          "JAK2_large", "vf_01", "chip", "Age", "sex", "Race", "BMI", "T2D", 
          "Current_smoker_baseline", "Ever_smoker_baseline", "t_CAD_all", "i_CAD_all")

# Create TableOne summary
tableOne <- CreateTableOne(vars = vars, data = mesa_prot, factorVars = factorVars)

# Output TableOne to file
write.table(print(tableOne), "/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/mesa_tableone.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE)

# Create race-specific tables for Black and White populations
tableOne_black <- CreateTableOne(vars = vars, data = mesa_prot[which(mesa_prot$Race == "Black"), ], factorVars = factorVars)
write.table(print(tableOne_black), "/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/mesa_black_tableone.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE)

tableOne_white <- CreateTableOne(vars = vars, data = mesa_prot[which(mesa_prot$Race == "White"), ], factorVars = factorVars)
write.table(print(tableOne_white), "/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/mesa_white_tableone.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE)

# Stratify by sex and race
mesachipprev_sex <- CreateTableOne(vars = vars, data = mesa_prot, factorVars = factorVars, strata = c("sex"))
mesachipprev_race <- CreateTableOne(vars = vars, data = mesa_prot, factorVars = factorVars, strata = c("Race"))

# Adjust missing values in mesa_prot
mesa_prot$Age <- ifelse(is.na(mesa_prot$Age), median(mesa_prot$Age, na.rm = TRUE), mesa_prot$Age)
mesa_prot$BMI <- ifelse(is.na(mesa_prot$BMI), median(mesa_prot$BMI, na.rm = TRUE), mesa_prot$BMI)
mesa_prot$T2D <- ifelse(is.na(mesa_prot$T2D), median(mesa_prot$T2D, na.rm = TRUE), mesa_prot$T2D)
mesa_prot$Ever_smoker_baseline <- ifelse(is.na(mesa_prot$Ever_smoker_baseline), median(mesa_prot$Ever_smoker_baseline, na.rm = TRUE), mesa_prot$Ever_smoker_baseline)
mesa_prot$HDL.norm <- ifelse(is.na(mesa_prot$HDL.norm), median(mesa_prot$HDL.norm, na.rm = TRUE), mesa_prot$HDL.norm)
mesa_prot$White <- ifelse(mesa_prot$Race == "White", 1, 0)

# Variables and models for CAD proteomics analysis
ds <- c("mesa_prot")
outcome <- c("t_CAD_all")
Model <- c("+Age+sex+White+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10", 
           "+Age+sex+White+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+T2D+BMI+Ever_smoker_baseline")

for (j in 1:length(ds)) {
  d <- get(ds[j])
  for (m in 1:length(Model)) {
    h <- data.frame(seq_id = character(), N = double(), beta = double(), se = double(), z = double(), 
                    p = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)
    for (i in 1:length(proteins)) {
      h[i, "seq_id"] <- proteins[i]
      fmla <- as.formula(paste(outcome[1], "~", proteins[i], Model[m], sep = ""))
      f <- glm(fmla, data = d, na.action = na.exclude, family = "binomial")
      h[i, "N"] <- summary(f)$df.null + 1
      h[i, "beta"] <- summary(f)$coefficient[2, 1]
      h[i, "se"] <- summary(f)$coefficient[2, 2]
      h[i, "z"] <- summary(f)$coefficient[2, 3]
      h[i, "p"] <- summary(f)$coefficient[2, 4]
      h[i, "lci"] <- summary(f)$coefficient[2, 1] - 1.96 * summary(f)$coefficient[2, 2]
      h[i, "uci"] <- summary(f)$coefficient[2, 1] + 1.96 * summary(f)$coefficient[2, 2]
    }
    h$p_FDR <- p.adjust(h$p, method = "fdr")
    h <- merge(anno, h, by.x = "SomaId", by.y = "seq_id")
    h <- h[order(h$p), ]
    write.table(h, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/MESA_proteomics_CAD_", outcome[1], "_", model[m], ".txt", sep = ""), 
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

##########################################################
######################### MESA ###########################
##########################################################

# Filter MESA cohort and remove specific columns
mesa <- TOPMed_chip[TOPMed_chip$study == "MESA", ] %>% select(-c(CHROM:F2R1))

# Remove duplicates
mesa <- mesa[!duplicated(mesa)]

# Create mutation indicator variables for each gene
mesa$DNMT3A_all <- ifelse(mesa$Sample %in% DNMT3A$Sample, 1, 0)
mesa$TET2_all <- ifelse(mesa$Sample %in% TET2$Sample, 1, 0)
mesa$ASXL1_all <- ifelse(mesa$Sample %in% ASXL1$Sample, 1, 0)
mesa$JAK2_all <- ifelse(mesa$Sample %in% JAK2$Sample, 1, 0)

# Create large variant indicator variables
mesa$DNMT3A_large <- ifelse(mesa$Sample %in% DNMT3A$Sample & mesa$Sample %in% vf_01$Sample, 1, 0)
mesa$TET2_large <- ifelse(mesa$Sample %in% TET2$Sample & mesa$Sample %in% vf_01$Sample, 1, 0)
mesa$ASXL1_large <- ifelse(mesa$Sample %in% ASXL1$Sample & mesa$Sample %in% vf_01$Sample, 1, 0)
mesa$JAK2_large <- ifelse(mesa$Sample %in% JAK2$Sample & mesa$Sample %in% vf_01$Sample, 1, 0)

# Create additional variables
mesa$vf_01 <- ifelse(mesa$Sample %in% vf_01$Sample, 1, 0)
mesa$chip <- ifelse(mesa$Sample %in% chip$Sample, 1, 0)

# Select specific columns
mesa <- mesa[, c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "DNMT3A_large", "TET2_large", "ASXL1_large", 
                 "JAK2_large", "vf_01", "chip", "Sample", "sex", "Age", "Race", "BMI", "T2D", "Current_smoker_baseline", 
                 "Ever_smoker_baseline", "t_CAD_all", "i_CAD_all", "i_CAD_all_date_censor", "TOTAL_ADJ.norm", "HDL.norm")]

# Merge with principal component data
mesa <- merge(mesa, pc, by = "Sample")

# Load proteomics data
peer_prot <- fread("/medpop/esp2/zyu/chip_protemoics/data/MESA/PEER/peernum/peerandresid/prot.MESA_log2_nopeer.txt")
prot_peercovar <- fread("/medpop/esp2/zyu/chip_protemoics/data/MESA/PEER/peernum/peerandresid/peers.MESA_log2_50.txt")

# Merge proteomics data with cohort data
mesa_prot <- merge(mesa, peer_prot, by.x = "Sample", by.y = "SampleId")
mesa_prot <- merge(mesa_prot, prot_peercovar, by.x = "Sample", by.y = "SampleId")

# Summary of age
summary(mesa_prot$Age)

# Extract sample list
mesa_sample <- mesa_prot$Sample

# Define factor and numeric variables for TableOne
factorVars <- c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "DNMT3A_large", "TET2_large", 
                "ASXL1_large", "JAK2_large", "vf_01", "chip", "sex", "Race", "T2D", 
                "Current_smoker_baseline", "Ever_smoker_baseline", "t_CAD_all", "i_CAD_all")

vars <- c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "DNMT3A_large", "TET2_large", "ASXL1_large", 
          "JAK2_large", "vf_01", "chip", "Age", "sex", "Race", "BMI", "T2D", 
          "Current_smoker_baseline", "Ever_smoker_baseline", "t_CAD_all", "i_CAD_all")

# Create TableOne summary
tableOne <- CreateTableOne(vars = vars, data = mesa_prot, factorVars = factorVars)

# Output TableOne to file
write.table(print(tableOne), "/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/mesa_tableone.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE)

# Create race-specific tables for Black and White populations
tableOne_black <- CreateTableOne(vars = vars, data = mesa_prot[which(mesa_prot$Race == "Black"), ], factorVars = factorVars)
write.table(print(tableOne_black), "/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/mesa_black_tableone.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE)

tableOne_white <- CreateTableOne(vars = vars, data = mesa_prot[which(mesa_prot$Race == "White"), ], factorVars = factorVars)
write.table(print(tableOne_white), "/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/mesa_white_tableone.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE)

# Stratify by sex and race
mesachipprev_sex <- CreateTableOne(vars = vars, data = mesa_prot, factorVars = factorVars, strata = c("sex"))
mesachipprev_race <- CreateTableOne(vars = vars, data = mesa_prot, factorVars = factorVars, strata = c("Race"))

# Adjust missing values in mesa_prot
mesa_prot$Age <- ifelse(is.na(mesa_prot$Age), median(mesa_prot$Age, na.rm = TRUE), mesa_prot$Age)
mesa_prot$BMI <- ifelse(is.na(mesa_prot$BMI), median(mesa_prot$BMI, na.rm = TRUE), mesa_prot$BMI)
mesa_prot$T2D <- ifelse(is.na(mesa_prot$T2D), median(mesa_prot$T2D, na.rm = TRUE), mesa_prot$T2D)
mesa_prot$Ever_smoker_baseline <- ifelse(is.na(mesa_prot$Ever_smoker_baseline), median(mesa_prot$Ever_smoker_baseline, na.rm = TRUE), mesa_prot$Ever_smoker_baseline)
mesa_prot$HDL.norm <- ifelse(is.na(mesa_prot$HDL.norm), median(mesa_prot$HDL.norm, na.rm = TRUE), mesa_prot$HDL.norm)
mesa_prot$White <- ifelse(mesa_prot$Race == "White", 1, 0)

# Variables and models for CAD proteomics analysis
ds <- c("mesa_prot")
outcome <- c("t_CAD_all")
Model <- c("+Age+sex+White+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10", 
           "+Age+sex+White+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+T2D+BMI+Ever_smoker_baseline")

for (j in 1:length(ds)) {
  d <- get(ds[j])
  for (m in 1:length(Model)) {
    h <- data.frame(seq_id = character(), N = double(), beta = double(), se = double(), z = double(), 
                    p = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)
    for (i in 1:length(proteins)) {
      h[i, "seq_id"] <- proteins[i]
      fmla <- as.formula(paste(outcome[1], "~", proteins[i], Model[m], sep = ""))
      f <- glm(fmla, data = d, na.action = na.exclude, family = "binomial")
      h[i, "N"] <- summary(f)$df.null + 1
      h[i, "beta"] <- summary(f)$coefficient[2, 1]
      h[i, "se"] <- summary(f)$coefficient[2, 2]
      h[i, "z"] <- summary(f)$coefficient[2, 3]
      h[i, "p"] <- summary(f)$coefficient[2, 4]
      h[i, "lci"] <- summary(f)$coefficient[2, 1] - 1.96 * summary(f)$coefficient[2, 2]
      h[i, "uci"] <- summary(f)$coefficient[2, 1] + 1.96 * summary(f)$coefficient[2, 2]
    }
    h$p_FDR <- p.adjust(h$p, method = "fdr")
    h <- merge(anno, h, by.x = "SomaId", by.y = "seq_id")
    h <- h[order(h$p), ]
    write.table(h, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/MESA_proteomics_CAD_", outcome[1], "_", model[m], ".txt", sep = ""), 
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

###### CHIP proteomics ######

# Models for CHIP proteomics
Model <- c("", 
           "+Age+sex+White+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+T2D+BMI+Ever_smoker_baseline+peer1+peer2+peer3+peer4+peer5+peer6+peer7+peer8+peer9+peer10+peer11+peer12+peer13+peer14+peer15+peer16+peer17+peer18+peer19+peer20+peer21+peer22+peer23+peer24+peer25+peer26+peer27+peer28+peer29+peer30+peer31+peer32+peer33+peer34+peer35+peer36+peer37+peer38+peer39+peer40+peer41+peer42+peer43+peer44+peer45+peer46+peer47+peer48+peer49+peer50", 
           "+Age+sex+White+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+T2D+BMI+Ever_smoker_baseline")
model <- c("noadjust", "fully_adjust", "nopeer")

# Data set for CHIP proteomics
ds <- c("mesa_prot")
dss <- c("")

for (j in 1:length(dss)) {
  d <- get(ds[j])
  for (m in 1:length(Model)) {
    h <- data.frame(seq_id = character(), chip = character(), N = double(), beta = double(), se = double(), 
                    z = double(), p = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)
    for (k in 1:length(variants_list)) {
      for (i in 1:length(proteins)) {
        h[i, "seq_id"] <- proteins[i]
        h[i, "chip"] <- variants_list[k]
        fmla <- as.formula(paste(proteins[i], "~", variants_list[k], Model[m], sep = ""))
        f <- glm(fmla, data = d, na.action = na.exclude, family = "gaussian")
        h[i, "N"] <- summary(f)$df.null + 1
        h[i, "beta"] <- summary(f)$coefficient[2, 1]
        h[i, "se"] <- summary(f)$coefficient[2, 2]
        h[i, "z"] <- summary(f)$coefficient[2, 3]
        h[i, "p"] <- summary(f)$coefficient[2, 4]
        h[i, "lci"] <- summary(f)$coefficient[2, 1] - 1.96 * summary(f)$coefficient[2, 2]
        h[i, "uci"] <- summary(f)$coefficient[2, 1] + 1.96 * summary(f)$coefficient[2, 2]
      }
      if (k == 1) {
        h_merge <- h
      } else {
        h_merge <- rbind(h_merge, h)
      }
    }
    h_merge <- merge(anno, h_merge, by.x = "SomaId", by.y = "seq_id")
    h_merge <- h_merge[order(h_merge$p), ]
    write.table(h_merge, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_MESA_proteomics_CHIP_", model[m], ".txt", sep = ""), 
                sep = "\t", quote = FALSE, row.names = FALSE)

    # Filter results for FDR and specific variants
    h_merge_fdr <- h_merge[which(h_merge$chip %in% c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "chip")), ]
    h_merge_fdr$p_FDR <- p.adjust(h_merge_fdr$p, method = "fdr")
    write.table(h_merge_fdr, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_FDR_MESA_proteomics_CHIP_", model[m], ".txt", sep = ""), 
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

## Interaction with sex

# Model for interaction analysis
Model <- c("*sex+Age+White+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+T2D+BMI+Ever_smoker_baseline+peer1+peer2+peer3+peer4+
           peer5+peer6+peer7+peer8+peer9+peer10+peer11+peer12+peer13+peer14+peer15+peer16+peer17+peer18+peer19+peer20+peer21+
           peer22+peer23+peer24+peer25+peer26+peer27+peer28+peer29+peer30+peer31+peer32+peer33+peer34+peer35+peer36+peer37+
           peer38+peer39+peer40+peer41+peer42+peer43+peer44+peer45+peer46+peer47+peer48+peer49+peer50")

model <- c("fully_adjust")

for (j in 1:length(dss)) {
  d <- get(ds[j])
  for (m in 1:length(Model)) {
    h <- data.frame(seq_id = character(), chip = character(), N = double(), beta = double(), se = double(), 
                    z = double(), p = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)
    for (k in 1:length(variants_list)) {
      for (i in 1:length(proteins)) {
        h[i, "seq_id"] <- proteins[i]
        h[i, "chip"] <- variants_list[k]
        fmla <- as.formula(paste(proteins[i], "~", variants_list[k], Model[m], sep = ""))
        f <- glm(fmla, data = d, na.action = na.exclude, family = "gaussian")

        # Handle NA interaction terms
        if (nrow(summary(f)$coefficient) < 69 || is.na(summary(f)$coefficient[69, 1])) {
          h[i, c("N", "beta", "se", "z", "p", "lci", "uci")] <- NA
        } else {
          h[i, "N"] <- summary(f)$df.null + 1
          h[i, "beta"] <- summary(f)$coefficient[69, 1]
          h[i, "se"] <- summary(f)$coefficient[69, 2]
          h[i, "z"] <- summary(f)$coefficient[69, 3]
          h[i, "p"] <- summary(f)$coefficient[69, 4]
          h[i, "lci"] <- summary(f)$coefficient[69, 1] - 1.96 * summary(f)$coefficient[69, 2]
          h[i, "uci"] <- summary(f)$coefficient[69, 1] + 1.96 * summary(f)$coefficient[69, 2]
        }
      }
      if (k == 1) {
        h_merge <- h
      } else {
        h_merge <- rbind(h_merge, h)
      }
    }
    h_merge <- merge(anno, h_merge, by.x = "SomaId", by.y = "seq_id")
    h_merge <- h_merge[order(h_merge$p), ]
    write.table(h_merge, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_MESA_proteomics_interaction_CHIP_", model[m], ".txt", sep = ""), 
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

####### Sex-specific analysis #######

# Define sex-specific datasets
mesa_prot_male <- mesa_prot[which(mesa_prot$sex == "M"), ]
mesa_prot_female <- mesa_prot[which(mesa_prot$sex == "F"), ]
ds <- c("mesa_prot_male", "mesa_prot_female")
dss <- c("male", "female")

# Loop over sex-specific datasets
for (j in 1:length(dss)) {
  d <- get(ds[j])
  for (m in 1:length(Model)) {
    h <- data.frame(seq_id = character(), chip = character(), N = double(), beta = double(), se = double(), 
                    z = double(), p = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)
    for (k in 1:length(variants_list)) {
      for (i in 1:length(proteins)) {
        h[i, "seq_id"] <- proteins[i]
        h[i, "chip"] <- variants_list[k]
        fmla <- as.formula(paste(proteins[i], "~", variants_list[k], Model[m], sep = ""))
        f <- glm(fmla, data = d, na.action = na.exclude, family = "gaussian")
        h[i, "N"] <- summary(f)$df.null + 1
        h[i, "beta"] <- summary(f)$coefficient[2, 1]
        h[i, "se"] <- summary(f)$coefficient[2, 2]
        h[i, "z"] <- summary(f)$coefficient[2, 3]
        h[i, "p"] <- summary(f)$coefficient[2, 4]
        h[i, "lci"] <- summary(f)$coefficient[2, 1] - 1.96 * summary(f)$coefficient[2, 2]
        h[i, "uci"] <- summary(f)$coefficient[2, 1] + 1.96 * summary(f)$coefficient[2, 2]
      }
      if (k == 1) {
        h_merge <- h
      } else {
        h_merge <- rbind(h_merge, h)
      }
    }
    h_merge <- merge(anno, h_merge, by.x = "SomaId", by.y = "seq_id")
    h_merge <- h_merge[order(h_merge$p), ]
    write.table(h_merge, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_MESA_proteomics_CHIP_", dss[j], "_", model[m], ".txt", sep = ""), 
                sep = "\t", quote = FALSE, row.names = FALSE)

    # Filter results for FDR and specific variants
    h_merge_fdr <- h_merge[which(h_merge$chip %in% c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "chip")), ]
    h_merge_fdr$p_FDR <- p.adjust(h_merge_fdr$p, method = "fdr")
    write.table(h_merge_fdr, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_FDR_MESA_proteomics_CHIP_", dss[j], "_", model[m], ".txt", sep = ""), 
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

######################## Race-specific analysis ########################

# Define race-specific datasets
mesa_prot_white <- mesa_prot[which(mesa_prot$Race == "White"), ]
mesa_prot_black <- mesa_prot[which(mesa_prot$Race == "Black"), ]
ds <- c("mesa_prot_white", "mesa_prot_black")
dss <- c("white", "black")

# Loop over race-specific datasets
for (j in 1:length(dss)) {
  d <- get(ds[j])
  for (m in 1:length(Model)) {
    h <- data.frame(seq_id = character(), chip = character(), N = double(), beta = double(), se = double(),
                    z = double(), p = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)
    for (k in 1:length(variants_list)) {
      for (i in 1:length(proteins)) {
        h[i, "seq_id"] <- proteins[i]
        h[i, "chip"] <- variants_list[k]
        fmla <- as.formula(paste(proteins[i], "~", variants_list[k], Model[m], sep = ""))
        f <- glm(fmla, data = d, na.action = na.exclude, family = "gaussian")
        h[i, "N"] <- summary(f)$df.null + 1
        h[i, "beta"] <- summary(f)$coefficient[2, 1]
        h[i, "se"] <- summary(f)$coefficient[2, 2]
        h[i, "z"] <- summary(f)$coefficient[2, 3]
        h[i, "p"] <- summary(f)$coefficient[2, 4]
        h[i, "lci"] <- summary(f)$coefficient[2, 1] - 1.96 * summary(f)$coefficient[2, 2]
        h[i, "uci"] <- summary(f)$coefficient[2, 1] + 1.96 * summary(f)$coefficient[2, 2]
      }
      if (k == 1) {
        h_merge <- h
      } else {
        h_merge <- rbind(h_merge, h)
      }
    }
    h_merge <- merge(anno, h_merge, by.x = "SomaId", by.y = "seq_id")
    h_merge <- h_merge[order(h_merge$p), ]
    write.table(h_merge, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_MESA_proteomics_CHIP_", dss[j], "_", model[m], ".txt", sep = ""),
                sep = "\t", quote = FALSE, row.names = FALSE)

    # Filter results for FDR and specific variants
    h_merge_fdr <- h_merge[which(h_merge$chip %in% c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "chip")), ]
    h_merge_fdr$p_FDR <- p.adjust(h_merge_fdr$p, method = "fdr")
    write.table(h_merge_fdr, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_FDR_MESA_proteomics_CHIP_", dss[j], "_", model[m], ".txt", sep = ""),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}


######################## PEER testing ########################

# Define dataset and model for PEER testing
ds <- c("mesa_prot")
dss <- c("")
outcome <- c("t_CAD_all", variants_list)
Model <- c("peer1+peer2+peer3+peer4+
            peer5+peer6+peer7+peer8+peer9+peer10+peer11+peer12+peer13+peer14+peer15+peer16+peer17+peer18+peer19+peer20+
            peer21+peer22+peer23+peer24+peer25+peer26+peer27+peer28+peer29+peer30+peer31+peer32+peer33+peer34+peer35+
            peer36+peer37+peer38+peer39+peer40+peer41+peer42+peer43+peer44+peer45+peer46+peer47+peer48+peer49+peer50")

model <- c("peer")

# Loop over datasets for PEER testing
for (j in 1:length(dss)) {
  d <- get(ds[j])
  for (m in 1:length(Model)) {
    h <- data.frame(outcome = character(), N = double(), n_sig = double(), min_sig = double(), stringsAsFactors = FALSE)
    for (i in 1:length(outcome)) {
      h[i, "outcome"] <- outcome[i]
      fmla <- as.formula(paste(outcome[i], "~", Model[m], sep = ""))
      f <- glm(fmla, data = d, na.action = na.exclude, family = "binomial")
      h[i, "N"] <- summary(f)$df.null + 1
      h[i, "n_sig"] <- sum(summary(f)$coefficient[2:51, 4] < 0.05)
      h[i, "min_sig"] <- min(summary(f)$coefficient[2:51, 4])
    }
  }
}

##########################################################
###################### CHS Proteomics ####################
##########################################################

chs <- TOPMed_chip[TOPMed_chip$study == "CHS", ] %>% select(-c(CHROM: F2R1))
chs <- chs[!duplicated(chs)]
chs$DNMT3A_all <- ifelse(chs$Sample %in% DNMT3A$Sample, 1, 0)
chs$TET2_all <- ifelse(chs$Sample %in% TET2$Sample, 1, 0)
chs$ASXL1_all <- ifelse(chs$Sample %in% ASXL1$Sample, 1, 0)
chs$JAK2_all <- ifelse(chs$Sample %in% JAK2$Sample, 1, 0)
chs$DNMT3A_large <- ifelse(chs$Sample %in% DNMT3A$Sample & chs$Sample %in% vf_01$Sample, 1, 0)
chs$TET2_large <- ifelse(chs$Sample %in% TET2$Sample & chs$Sample %in% vf_01$Sample, 1, 0)
chs$ASXL1_large <- ifelse(chs$Sample %in% ASXL1$Sample & chs$Sample %in% vf_01$Sample, 1, 0)
chs$JAK2_large <- ifelse(chs$Sample %in% JAK2$Sample & chs$Sample %in% vf_01$Sample, 1, 0)
chs$vf_01 <- ifelse(chs$Sample %in% vf_01$Sample, 1, 0)
chs$chip <- ifelse(chs$Sample %in% chip$Sample, 1, 0)

chs <- chs[, c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "DNMT3A_large", 
               "TET2_large", "ASXL1_large", "JAK2_large", "vf_01", "chip", "sex", 
               "Sample", "Age", "Race", "BMI", "T2D", "Current_smoker_baseline", 
               "Ever_smoker_baseline", "t_CAD_all", "i_CAD_all", "i_CAD_all_date_censor", 
               "TOTAL_ADJ.norm", "HDL.norm")]

chs <- merge(chs, pc, by = "Sample")
chs_linker <- fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/chs_linker.txt")
peer_prot <- fread("/medpop/esp2/zyu/chip_protemoics/data/CHS/PEER/peernum/peerandresid/prot.CHS_log2_nopeer.txt")
peer_prot <- merge(chs_linker, peer_prot, by.x = "subject_id", by.y = "SampleId")
prot_peercovar <- fread("/medpop/esp2/zyu/chip_protemoics/data/CHS/PEER/peernum/peerandresid/peers.CHS_log2_50.txt")

chs_prot <- merge(chs, peer_prot, by = "Sample")
chs_prot <- merge(chs_prot, prot_peercovar, by.x = "subject_id", by.y = "SampleId")
summary(chs_prot$Age)

chs_sample <- chs_prot$Sample

factorVars <- c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "DNMT3A_large", 
                "TET2_large", "ASXL1_large", "JAK2_large", "vf_01", "chip", "sex", 
                "Race", "T2D", "Current_smoker_baseline", "Ever_smoker_baseline", 
                "t_CAD_all", "i_CAD_all")
vars <- c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "DNMT3A_large", 
          "TET2_large", "ASXL1_large", "JAK2_large", "vf_01", "chip", "Age", 
          "sex", "Race", "BMI", "T2D", "Current_smoker_baseline", "Ever_smoker_baseline", 
          "t_CAD_all", "i_CAD_all")

# Create table for full cohort
tableOne <- CreateTableOne(vars = vars, data = chs_prot, factorVars = factorVars)
write.table(print(tableOne), "/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/chs_tableone.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE)

# Create table for Black cohort
tableOne <- CreateTableOne(vars = vars, data = chs_prot[which(chs_prot$Race == "Black"), ], 
                           factorVars = factorVars)
write.table(print(tableOne), "/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/chs_black_tableone.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE)

# Create table for White cohort
tableOne <- CreateTableOne(vars = vars, data = chs_prot[which(chs_prot$Race == "White"), ], 
                           factorVars = factorVars)
write.table(print(tableOne), "/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/chs_white_tableone.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE)


###################### Prevalence Analysis ######################

# Stratify by sex for CHIP prevalence analysis
factorVars <- c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "chip")
vars <- c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "chip")
chschipprev_sex <- CreateTableOne(vars = vars, data = chs_prot, factorVars = factorVars, strata = c("sex"))
write.table(print(chschipprev_sex), "/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/chschipprev_sex.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE)

# Stratify by race for CHIP prevalence analysis
chschipprev_race <- CreateTableOne(vars = vars, data = chs_prot, factorVars = factorVars, strata = c("Race"))
write.table(print(chschipprev_race), "/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/chschipprev_race.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE)

write.table(print(mesachipprev_sex), "/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/mesachipprev_sex.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE)
write.table(print(mesachipprev_race), "/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/mesachipprev_race.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE)


###################### Adjusting Column Names ######################

# Rename specific protein
names(chs_prot)[which(names(chs_prot) == "SL003198_2728_62")] <- "SL000497_2728_62"

# Remove suffix from column names
colnames(chs_prot)[35:ncol(chs_prot)] <- sub("_.*", "", colnames(chs_prot)[35:ncol(chs_prot)])

proteins <- colnames(chs_prot)[35:1331]
variants_list <- colnames(chs_prot)[3:12]


###################### CAD Proteomics Analysis ######################

chs_prot$BMI <- ifelse(is.na(chs_prot$BMI), median(chs_prot$BMI, na.rm = TRUE), chs_prot$BMI)
chs_prot$Ever_smoker_baseline <- ifelse(is.na(chs_prot$Ever_smoker_baseline), median(chs_prot$Ever_smoker_baseline, na.rm = TRUE), chs_prot$Ever_smoker_baseline)
chs_prot$HDL.norm <- ifelse(is.na(chs_prot$HDL.norm), median(chs_prot$HDL.norm, na.rm = TRUE), chs_prot$HDL.norm)
chs_prot$White <- ifelse(chs_prot$Race == "White", 1, 0)

# Define datasets and models for CAD proteomics analysis
ds <- c("chs_prot")
dss <- c("")
outcome <- c("t_CAD_all")
Model <- c("+Age+sex+White+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10", 
           "+Age+sex+White+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+T2D+BMI+Ever_smoker_baseline")
model <- c("agesexracepc_adjust", "fully_adjust")

# Loop over datasets and models
for (j in 1:length(dss)) {
  d <- get(ds[j])
  k <- 1
  for (m in 1:length(Model)) {
    h <- data.frame(seq_id = character(), N = double(), beta = double(), se = double(),
                    z = double(), p = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)
    for (i in 1:length(proteins)) {
      h[i, "seq_id"] <- proteins[i]
      fmla <- as.formula(paste(outcome[k], "~", proteins[i], Model[m], sep = ""))
      f <- glm(fmla, data = d, na.action = na.exclude, family = "binomial")
      h[i, "N"] <- summary(f)$df.null + 1
      h[i, "beta"] <- summary(f)$coefficient[2, 1]
      h[i, "se"] <- summary(f)$coefficient[2, 2]
      h[i, "z"] <- summary(f)$coefficient[2, 3]
      h[i, "p"] <- summary(f)$coefficient[2, 4]
      h[i, "lci"] <- summary(f)$coefficient[2, 1] - 1.96 * summary(f)$coefficient[2, 2]
      h[i, "uci"] <- summary(f)$coefficient[2, 1] + 1.96 * summary(f)$coefficient[2, 2]
    }
    h$p_FDR <- p.adjust(h$p, method = "fdr")
    h <- merge(anno, h, by.x = "SomaId", by.y = "seq_id")
    h <- h[order(h$p), ]
    write.table(h, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/CHS_proteomics_CAD", outcome[k], "_", model[m], ".txt", sep = ""),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

###################### CHIP Proteomics Analysis ######################

Model <- c("", "+Age+sex+White+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+T2D+BMI+Ever_smoker_baseline+peer1+peer2+peer3+peer4+
peer5+peer6+peer7+peer8+peer9+peer10+peer11+peer12+peer13+peer14+peer15+peer16+peer17+peer18+peer19+peer20+peer21+peer22+
peer23+peer24+peer25+peer26+peer27+peer28+peer29+peer30+peer31+peer32+peer33+peer34+peer35+peer36+peer37+peer38+peer39+peer40+
peer41+peer42+peer43+peer44+peer45+peer46+peer47+peer48+peer49+peer50", 
"+Age+sex+White+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+T2D+BMI+Ever_smoker_baseline")
model <- c("noadjust", "fully_adjust", "nopeer")
ds <- c("chs_prot")
dss <- c("")

for (j in 1:length(dss)) {
  d <- get(ds[j])
  for (m in 1:length(Model)) {    
    h <- data.frame(seq_id = character(), chip = character(), N = double(), beta = double(), se = double(), 
                    z = double(), p = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)
    for (k in 1:length(variants_list)) {
      for (i in 1:length(proteins)) {
        h[i, "seq_id"] <- proteins[i]
        h[i, "chip"] <- variants_list[k]
        fmla <- as.formula(paste(proteins[i], "~", variants_list[k], Model[m], sep = ""))
        f <- glm(fmla, data = d, na.action = na.exclude, family = "gaussian")
        h[i, "N"] <- summary(f)$df.null + 1
        h[i, "beta"] <- summary(f)$coefficient[2, 1]
        h[i, "se"] <- summary(f)$coefficient[2, 2]
        h[i, "z"] <- summary(f)$coefficient[2, 3]
        h[i, "p"] <- summary(f)$coefficient[2, 4]
        h[i, "lci"] <- summary(f)$coefficient[2, 1] - 1.96 * summary(f)$coefficient[2, 2]
        h[i, "uci"] <- summary(f)$coefficient[2, 1] + 1.96 * summary(f)$coefficient[2, 2]
      }
      if (k == 1) {
        h_merge <- h
      } else {
        h_merge <- rbind(h_merge, h)
      }
    }
    h_merge <- merge(anno, h_merge, by.x = "SomaId", by.y = "seq_id")
    h_merge <- h_merge[order(h_merge$p), ]
    write.table(h_merge, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_CHS_proteomics_CHIP_", model[m], ".txt", sep = ""), 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Filter significant results and apply FDR correction
    h_merge_fdr <- h_merge[which(h_merge$chip %in% c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "chip")), ]
    h_merge_fdr$p_FDR <- p.adjust(h_merge_fdr$p, method = "fdr")
    write.table(h_merge_fdr, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_FDR_CHS_proteomics_CHIP_", model[m], ".txt", sep = ""), 
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

###################### Interaction of Sex ######################

Model <- c("*sex+Age+White+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+T2D+BMI+Ever_smoker_baseline+peer1+peer2+peer3+peer4+
peer5+peer6+peer7+peer8+peer9+peer10+peer11+peer12+peer13+peer14+peer15+peer16+peer17+peer18+peer19+peer20+peer21+peer22+
peer23+peer24+peer25+peer26+peer27+peer28+peer29+peer30+peer31+peer32+peer33+peer34+peer35+peer36+peer37+peer38+peer39+peer40+
peer41+peer42+peer43+peer44+peer45+peer46+peer47+peer48+peer49+peer50")
model <- c("fully_adjust")
ds <- c("chs_prot")
dss <- c("")

for (j in 1:length(dss)) {
  d <- get(ds[j])
  for (m in 1:length(Model)) {    
    h <- data.frame(seq_id = character(), chip = character(), N = double(), beta = double(), se = double(), 
                    z = double(), p = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)
    for (k in 1:length(variants_list)) {
      for (i in 1:length(proteins)) {
        h[i, "seq_id"] <- proteins[i]
        h[i, "chip"] <- variants_list[k]
        fmla <- as.formula(paste(proteins[i], "~", variants_list[k], Model[m], sep = ""))
        f <- glm(fmla, data = d, na.action = na.exclude, family = "gaussian")

        # Check if the interaction term is NA, if it is then fill in NA values
        if (nrow(summary(f)$coefficient) < 69 || is.na(summary(f)$coefficient[69, 1])) {
          h[i, c("N", "beta", "se", "z", "p", "lci", "uci")] <- NA
        } else {
          h[i, "N"] <- summary(f)$df.null + 1
          h[i, "beta"] <- summary(f)$coefficient[69, 1]
          h[i, "se"] <- summary(f)$coefficient[69, 2]
          h[i, "z"] <- summary(f)$coefficient[69, 3]
          h[i, "p"] <- summary(f)$coefficient[69, 4]
          h[i, "lci"] <- summary(f)$coefficient[69, 1] - 1.96 * summary(f)$coefficient[69, 2]
          h[i, "uci"] <- summary(f)$coefficient[69, 1] + 1.96 * summary(f)$coefficient[69, 2]
        }
      }
      if (k == 1) {
        h_merge <- h
      } else {
        h_merge <- rbind(h_merge, h)
      }
    }
    h_merge <- merge(anno, h_merge, by.x = "SomaId", by.y = "seq_id")
    h_merge <- h_merge[order(h_merge$p), ]
    write.table(h_merge, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_CHS_proteomics_interaction_CHIP_", model[m], ".txt", sep = ""), 
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

####### Sex Specific Analysis #######

Model <- c("+Age+White+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+T2D+BMI+Ever_smoker_baseline+peer1+peer2+peer3+peer4+
peer5+peer6+peer7+peer8+peer9+peer10+peer11+peer12+peer13+peer14+peer15+peer16+peer17+peer18+peer19+peer20+peer21+peer22+
peer23+peer24+peer25+peer26+peer27+peer28+peer29+peer30+peer31+peer32+peer33+peer34+peer35+peer36+peer37+peer38+peer39+peer40+
peer41+peer42+peer43+peer44+peer45+peer46+peer47+peer48+peer49+peer50")
model <- c("fully_adjust")

chs_prot_male <- chs_prot[which(chs_prot$sex == "M"), ]
chs_prot_female <- chs_prot[which(chs_prot$sex == "F"), ]
ds <- c("chs_prot_male", "chs_prot_female")
dss <- c("male", "female")

for (j in 1:length(dss)) {
  d <- get(ds[j])
  for (m in 1:length(Model)) {    
    h <- data.frame(seq_id = character(), chip = character(), N = double(), beta = double(), se = double(), 
                    z = double(), p = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)
    for (k in 1:length(variants_list)) {
      for (i in 1:length(proteins)) {
        h[i, "seq_id"] <- proteins[i]
        h[i, "chip"] <- variants_list[k]
        fmla <- as.formula(paste(proteins[i], "~", variants_list[k], Model[m], sep = ""))
        f <- glm(fmla, data = d, na.action = na.exclude, family = "gaussian")
        h[i, "N"] <- summary(f)$df.null + 1
        h[i, "beta"] <- summary(f)$coefficient[2, 1]
        h[i, "se"] <- summary(f)$coefficient[2, 2]
        h[i, "z"] <- summary(f)$coefficient[2, 3]
        h[i, "p"] <- summary(f)$coefficient[2, 4]
        h[i, "lci"] <- summary(f)$coefficient[2, 1] - 1.96 * summary(f)$coefficient[2, 2]
        h[i, "uci"] <- summary(f)$coefficient[2, 1] + 1.96 * summary(f)$coefficient[2, 2]
      }
      if (k == 1) {
        h_merge <- h
      } else {
        h_merge <- rbind(h_merge, h)
      }
    }
    h_merge <- merge(anno, h_merge, by.x = "SomaId", by.y = "seq_id")
    h_merge <- h_merge[order(h_merge$p), ]
    write.table(h_merge, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_CHS_proteomics_CHIP_", dss[j], "_", model[m], ".txt", sep = ""), 
                sep = "\t", quote = FALSE, row.names = FALSE)

    # Apply FDR correction for selected variants
    h_merge_fdr <- h_merge[which(h_merge$chip %in% c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "chip")), ]
    h_merge_fdr$p_FDR <- p.adjust(h_merge_fdr$p, method = "fdr")
    write.table(h_merge_fdr, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_FDR_CHS_proteomics_CHIP_", dss[j], "_", model[m], ".txt", sep = ""), 
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

######################## Race Specific Analysis ########################

Model <- c("+Age+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+T2D+BMI+Ever_smoker_baseline+peer1+peer2+peer3+peer4+
peer5+peer6+peer7+peer8+peer9+peer10+peer11+peer12+peer13+peer14+peer15+peer16+peer17+peer18+peer19+peer20+peer21+peer22+
peer23+peer24+peer25+peer26+peer27+peer28+peer29+peer30+peer31+peer32+peer33+peer34+peer35+peer36+peer37+peer38+peer39+peer40+
peer41+peer42+peer43+peer44+peer45+peer46+peer47+peer48+peer49+peer50")
model <- c("fully_adjust")

chs_prot_white <- chs_prot[which(chs_prot$Race == "White"), ]
chs_prot_black <- chs_prot[which(chs_prot$Race == "Black"), ]
ds <- c("chs_prot_white", "chs_prot_black")
dss <- c("white", "black")

for (j in 1:length(dss)) {
  d <- get(ds[j])
  for (m in 1:length(Model)) {    
    h <- data.frame(seq_id = character(), chip = character(), N = double(), beta = double(), se = double(), 
                    z = double(), p = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)
    for (k in 1:length(variants_list)) {
      for (i in 1:length(proteins)) {
        h[i, "seq_id"] <- proteins[i]
        h[i, "chip"] <- variants_list[k]
        fmla <- as.formula(paste(proteins[i], "~", variants_list[k], Model[m], sep = ""))
        f <- glm(fmla, data = d, na.action = na.exclude, family = "gaussian")
        h[i, "N"] <- summary(f)$df.null + 1
        h[i, "beta"] <- summary(f)$coefficient[2, 1]
        h[i, "se"] <- summary(f)$coefficient[2, 2]
        h[i, "z"] <- summary(f)$coefficient[2, 3]
        h[i, "p"] <- summary(f)$coefficient[2, 4]
        h[i, "lci"] <- summary(f)$coefficient[2, 1] - 1.96 * summary(f)$coefficient[2, 2]
        h[i, "uci"] <- summary(f)$coefficient[2, 1] + 1.96 * summary(f)$coefficient[2, 2]
      }
      if (k == 1) {
        h_merge <- h
      } else {
        h_merge <- rbind(h_merge, h)
      }
    }
    h_merge <- merge(anno, h_merge, by.x = "SomaId", by.y = "seq_id")
    h_merge <- h_merge[order(h_merge$p), ]
    write.table(h_merge, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_CHS_proteomics_CHIP_", dss[j], "_", model[m], ".txt", sep = ""), 
                sep = "\t", quote = FALSE, row.names = FALSE)

    # Apply FDR correction for selected variants
    h_merge_fdr <- h_merge[which(h_merge$chip %in% c("DNMT3A_all", "TET2_all", "ASXL1_all", "JAK2_all", "chip")), ]
    h_merge_fdr$p_FDR <- p.adjust(h_merge_fdr$p, method = "fdr")
    write.table(h_merge_fdr, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_FDR_CHS_proteomics_CHIP_", dss[j], "_", model[m], ".txt", sep = ""), 
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

######################## PEER Testing ########################

ds <- c("chs_prot")
dss <- c("")
outcome <- c("t_CAD_all", variants_list)

Model <- c("peer1+peer2+peer3+peer4+peer5+peer6+peer7+peer8+peer9+peer10+peer11+peer12+peer13+peer14+peer15+peer16+
peer17+peer18+peer19+peer20+peer21+peer22+peer23+peer24+peer25+peer26+peer27+peer28+peer29+peer30+peer31+peer32+
peer33+peer34+peer35+peer36+peer37+peer38+peer39+peer40+peer41+peer42+peer43+peer44+peer45+peer46+peer47+peer48+
peer49+peer50")
model <- c("peer")

for (j in 1:length(dss)) {
  d <- get(ds[j])
  for (m in 1:length(Model)) {
    h <- data.frame(outcome = character(), N = double(), n_sig = double(), min_sig = double(), stringsAsFactors = FALSE)
    for (i in 1:length(outcome)) {
      h[i, "outcome"] <- outcome[i]
      fmla <- as.formula(paste(outcome[i], "~", Model[m], sep = ""))
      f <- glm(fmla, data = d, na.action = na.exclude, family = "binomial")
      h[i, "N"] <- summary(f)$df.null + 1
      h[i, "n_sig"] <- sum(summary(f)$coefficient[2:51, 4] < 0.05)
      h[i, "min_sig"] <- min(summary(f)$coefficient[2:51, 4])
    }
  }
}

############## Output IDs ##############
sample <- c(jhs_sample, mesa_sample, chs_sample)

# Save sample IDs to output
write.table(sample, 
            file = "/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/topmed_id.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



