# Clear environment and garbage collection
gc()
rm(list = ls())

# Load required libraries
library(dplyr)
library(tidyverse)
library(data.table)
library(ggrepel)
library(gridExtra)
library(survival)
library(tableone) # Library for table creation

########################################### Read in Data ###################################################
###### Read in Proteins ############
proteins_wide <- fread("/medpop/esp2/projects/UK_Biobank/baskets/ukb672544/olink_data3k_pivot.tsv.gz")
proteins_wide <- proteins_wide %>% filter(ins_index == 0)

linker <- fread("/medpop/esp2/projects/UK_Biobank/baskets/ukb672544/coding143.tsv")
linker$meaning <- sapply(strsplit(linker$meaning, ";"), `[`, 1)

# Create a named vector for the new column names
new_colnames <- setNames(linker$meaning, linker$coding)

# Get numeric column names and convert to numeric values
num_colnames <- as.numeric(names(proteins_wide)[3:ncol(proteins_wide)])

# Match the column names with the linker codes
matched_colnames <- new_colnames[match(num_colnames, as.numeric(names(new_colnames)))]

# Replace NAs with original column names
matched_colnames[is.na(matched_colnames)] <- as.character(num_colnames)[is.na(matched_colnames)]

# Replace the column names in the main data frame
names(proteins_wide)[3:ncol(proteins_wide)] <- matched_colnames

####### Read in CHIP Data ############
load("/medpop/esp2/mesbah/projects/Meta_GWAS/n650k/chip_call/ukb450k_mgbb53k.CHIP_vars.for_phewas.rda")
newchip_overall <- ukb450k_ch

# Replace NA in CHIP positive cases with 0
newchip_overall[which(newchip_overall$CHIP == 1),][is.na(newchip_overall[which(newchip_overall$CHIP == 1),])] <- 0

# Load CHIP VAF data and rename columns
newchip_vaf <- fread("/medpop/esp2/zyu/chip_protemoics/data/200K_CHIP/UKB200k_CHIP_Variants.11Aug2021.csv")
newchip_250k_vaf <- fread("/medpop/esp2/mesbah/datasets/CHIP/UKBB/450k/chipVars_CV_AB.ukb250k_all.csv.gz")
newchip_250k_vaf <- newchip_250k_vaf %>% rename(eid_7089 = Broad_ID)

# Summarize CHIP VAF
newchip_vaf_CHIP <- newchip_vaf %>% group_by(eid_7089) %>% summarise(CHIP_VAF = sum(AF))
newchip_250k_vaf_CHIP <- newchip_250k_vaf %>% group_by(eid_7089) %>% summarise(CHIP_VAF = sum(AF))

# Combine CHIP data
newchip_vaf_CHIP <- bind_rows(newchip_vaf_CHIP, newchip_250k_vaf_CHIP)

# Merge and mutate to calculate CHIP VAF
newchip_overall <- merge(newchip_overall, newchip_vaf_CHIP, by = "eid_7089", all.x = TRUE)
newchip_overall <- newchip_overall %>%
  mutate(CHIP_VAF = ifelse(CHIP == 0, 0, ifelse(CHIP == 1, CHIP_VAF, NA)))

# Process JAK2 Gene
newchip_vaf_JAK2 <- newchip_vaf %>% filter(Gene == "JAK2") %>%
  group_by(eid_7089) %>% summarise(JAK2_VAF = sum(AF))
newchip_250k_vaf_JAK2 <- newchip_250k_vaf %>% filter(Gene.refGene == "JAK2") %>%
  group_by(eid_7089) %>% summarise(JAK2_VAF = sum(AF))
newchip_vaf_JAK2 <- bind_rows(newchip_vaf_JAK2, newchip_250k_vaf_JAK2)

# Merge and update JAK2 VAF
newchip_overall <- merge(newchip_overall, newchip_vaf_JAK2, by = "eid_7089", all.x = TRUE)
newchip_overall <- newchip_overall %>%
  mutate(JAK2_VAF = ifelse(JAK2 == 0, 0, ifelse(JAK2 == 1, JAK2_VAF, NA)))

# Process DNMT3A Gene
newchip_vaf_DNMT3A <- newchip_vaf %>% filter(Gene == "DNMT3A") %>%
  group_by(eid_7089) %>% summarise(DNMT3A_VAF = sum(AF))
newchip_250k_vaf_DNMT3A <- newchip_250k_vaf %>% filter(Gene.refGene == "DNMT3A") %>%
  group_by(eid_7089) %>% summarise(DNMT3A_VAF = sum(AF))
newchip_vaf_DNMT3A <- bind_rows(newchip_vaf_DNMT3A, newchip_250k_vaf_DNMT3A)

# Merge and update DNMT3A VAF
newchip_overall <- merge(newchip_overall, newchip_vaf_DNMT3A, by = "eid_7089", all.x = TRUE)
newchip_overall <- newchip_overall %>%
  mutate(DNMT3A_VAF = ifelse(DNMT3A == 0, 0, ifelse(DNMT3A == 1, DNMT3A_VAF, NA)))

# Process TET2 Gene
newchip_vaf_TET2 <- newchip_vaf %>% filter(Gene == "TET2") %>%
  group_by(eid_7089) %>% summarise(TET2_VAF = sum(AF))
newchip_250k_vaf_TET2 <- newchip_250k_vaf %>% filter(Gene.refGene == "TET2") %>%
  group_by(eid_7089) %>% summarise(TET2_VAF = sum(AF))
newchip_vaf_TET2 <- bind_rows(newchip_vaf_TET2, newchip_250k_vaf_TET2)

# Merge and update TET2 VAF
newchip_overall <- merge(newchip_overall, newchip_vaf_TET2, by = "eid_7089", all.x = TRUE)
newchip_overall <- newchip_overall %>%
  mutate(TET2_VAF = ifelse(TET2 == 0, 0, ifelse(TET2 == 1, TET2_VAF, NA)))

# Process ASXL1 Gene
newchip_vaf_ASXL1 <- newchip_vaf %>% filter(Gene == "ASXL1") %>%
  group_by(eid_7089) %>% summarise(ASXL1_VAF = sum(AF))
newchip_250k_vaf_ASXL1 <- newchip_250k_vaf %>% filter(Gene.refGene == "ASXL1") %>%
  group_by(eid_7089) %>% summarise(ASXL1_VAF = sum(AF))
newchip_vaf_ASXL1 <- bind_rows(newchip_vaf_ASXL1, newchip_250k_vaf_ASXL1)

# Merge and update ASXL1 VAF
newchip_overall <- merge(newchip_overall, newchip_vaf_ASXL1, by = "eid_7089", all.x = TRUE)
newchip_overall <- newchip_overall %>%
  mutate(ASXL1_VAF = ifelse(ASXL1 == 0, 0, ifelse(ASXL1 == 1, ASXL1_VAF, NA)))

# Process TP53 Gene
newchip_vaf_TP53 <- newchip_vaf %>% filter(Gene == "TP53") %>%
  group_by(eid_7089) %>% summarise(TP53_VAF = sum(AF))
newchip_250k_vaf_TP53 <- newchip_250k_vaf %>% filter(Gene.refGene == "TP53") %>%
  group_by(eid_7089) %>% summarise(TP53_VAF = sum(AF))
newchip_vaf_TP53 <- bind_rows(newchip_vaf_TP53, newchip_250k_vaf_TP53)

# Merge and update TP53 VAF
newchip_overall <- merge(newchip_overall, newchip_vaf_TP53, by = "eid_7089", all.x = TRUE)
newchip_overall <- newchip_overall %>%
  mutate(TP53_VAF = ifelse(TP53 == 0, 0, ifelse(TP53 == 1, TP53_VAF, NA)))

# Process PPM1D Gene
newchip_vaf_PPM1D <- newchip_vaf %>% filter(Gene == "PPM1D") %>%
  group_by(eid_7089) %>% summarise(PPM1D_VAF = sum(AF))
newchip_250k_vaf_PPM1D <- newchip_250k_vaf %>% filter(Gene.refGene == "PPM1D") %>%
  group_by(eid_7089) %>% summarise(PPM1D_VAF = sum(AF))
newchip_vaf_PPM1D <- bind_rows(newchip_vaf_PPM1D, newchip_250k_vaf_PPM1D)

# Merge and update PPM1D VAF
newchip_overall <- merge(newchip_overall, newchip_vaf_PPM1D, by = "eid_7089", all.x = TRUE)
newchip_overall <- newchip_overall %>%
  mutate(PPM1D_VAF = ifelse(PPM1D == 0, 0, ifelse(PPM1D == 1, PPM1D_VAF, NA)))

# Process SF3B1 Gene
newchip_vaf_SF3B1 <- newchip_vaf %>% filter(Gene == "SF3B1") %>%
  group_by(eid_7089) %>% summarise(SF3B1_VAF = sum(AF))
newchip_250k_vaf_SF3B1 <- newchip_250k_vaf %>% filter(Gene.refGene == "SF3B1") %>%
  group_by(eid_7089) %>% summarise(SF3B1_VAF = sum(AF))
newchip_vaf_SF3B1 <- bind_rows(newchip_vaf_SF3B1, newchip_250k_vaf_SF3B1)

# Merge and update SF3B1 VAF
newchip_overall <- merge(newchip_overall, newchip_vaf_SF3B1, by = "eid_7089", all.x = TRUE)
newchip_overall <- newchip_overall %>%
  mutate(SF3B1_VAF = ifelse(SF3B1 == 0, 0, ifelse(SF3B1 == 1, SF3B1_VAF, NA)))

# Process SRSF2 Gene
newchip_vaf_SRSF2 <- newchip_vaf %>% filter(Gene == "SRSF2") %>%
  group_by(eid_7089) %>% summarise(SRSF2_VAF = sum(AF))
newchip_250k_vaf_SRSF2 <- newchip_250k_vaf %>% filter(Gene.refGene == "SRSF2") %>%
  group_by(eid_7089) %>% summarise(SRSF2_VAF = sum(AF))
newchip_vaf_SRSF2 <- bind_rows(newchip_vaf_SRSF2, newchip_250k_vaf_SRSF2)

# Merge and update SRSF2 VAF
newchip_overall <- merge(newchip_overall, newchip_vaf_SRSF2, by = "eid_7089", all.x = TRUE)
newchip_overall <- newchip_overall %>%
  mutate(SRSF2_VAF = ifelse(SRSF2 == 0, 0, ifelse(SRSF2 == 1, SRSF2_VAF, NA)))

# Update column names to append "_CHIPvar" to relevant columns
current_colnames <- colnames(newchip_overall)
current_colnames[-(1:4)] <- paste(current_colnames[-(1:4)], "_CHIPvar", sep = "")
colnames(newchip_overall) <- current_colnames

######### Covariates ############
pheno <- fread("/medpop/esp2/mzekavat/UKBB/ukbb_PhenoFile.ALL_500k.UpdatedIncdPhenos_202020.REAL.txt")
pheno <- data.frame(pheno) %>%
  select(id, Sex, genotyping_array, PC1:PC10, BMI, in_white_British_ancestry_subset, BMI, ever_smoked, Prev_Diabetes_Type_2, Creatinine.x)

ethnic <- fread("/medpop/esp2/zyu/chip_protemoics/data/ethnic.txt.gz")

print("Covariate loading complete")

######### Combine Phenotype Data with Ethnicity ##############
pheno <- merge(pheno, ethnic, all.x = TRUE, by = "id")

remove <- fread("/medpop/esp2/projects/UK_Biobank/withdrawn_samples/w7089_2023-04-25.csv", header = FALSE)
print("Sample removal list loaded")

######### Combine CHIP Proteomics Data with Phenotype ############
chip_proteomics <- merge(newchip_overall, proteins_wide, by.x = "eid_7089", by.y = "eid")
chip_proteomics <- chip_proteomics %>%
  filter(!(eid_7089 %in% remove$V1))

chip_proteomics_pheno <- merge(chip_proteomics, pheno, by.x = "eid_7089", by.y = "id")

########## Generate Additional Covariates (eGFR calculation) #############
calculate_eGFRcr <- function(gender, age, scr_umol_L) {
  # Convert Scr from umol/L to mg/dL
  scr <- scr_umol_L / 88.4
  
  # Define gender-based parameters for the eGFR formula
  kappa <- ifelse(gender == "Female", 0.7, 0.9)
  alpha <- ifelse(gender == "Female", -0.329, -0.411)
  
  # Calculate eGFR using the CKD-EPI formula
  eGFRcr <- 141 * (min(scr / kappa, 1) ^ alpha) * (max(scr / kappa, 1) ^ -1.209) * (0.993 ^ age)
  
  # Adjust eGFR for females
  if (gender == "Female") {
    eGFRcr <- eGFRcr * 1.018
  }
  
  return(eGFRcr)
}

# Apply eGFR calculation to the dataset
chip_proteomics_pheno$egfr_cr <- mapply(calculate_eGFRcr, chip_proteomics_pheno$Sex, chip_proteomics_pheno$AGE_assessment, chip_proteomics_pheno$Creatinine.x)

print("eGFR calculations complete")

########################################### Table 1 ###################################################
factorVars <- c("CHIP_CHIPvar", "DNMT3A_CHIPvar", "TET2_CHIPvar", "ASXL1_CHIPvar", "JAK2_CHIPvar", "in_white_British_ancestry_subset", "Sex", "ever_smoked", 
"Prev_Diabetes_Type_2")
vars <- c("AGE_assessment", "CHIP_CHIPvar", "DNMT3A_CHIPvar", "TET2_CHIPvar", "ASXL1_CHIPvar", "JAK2_CHIPvar", "in_white_British_ancestry_subset", "Sex", "ever_smoked", 
"Prev_Diabetes_Type_2", "BMI", "whiteblack")

tableOne <- CreateTableOne(vars = vars, data = chip_proteomics_pheno, factorVars = factorVars)
write.table(print(tableOne),"/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/ukb_table1_updatedwith3kprotein.txt", sep="\t", quote=FALSE, row.names=TRUE)

factorVars <- c("CHIP_CHIPvar", "DNMT3A_CHIPvar", "TET2_CHIPvar", "ASXL1_CHIPvar", "JAK2_CHIPvar")
vars <- c("CHIP_CHIPvar", "DNMT3A_CHIPvar", "TET2_CHIPvar", "ASXL1_CHIPvar", "JAK2_CHIPvar")

tableOne <- CreateTableOne(vars = vars, data = chip_proteomics_pheno, strata="Sex", factorVars = factorVars)
write.table(print(tableOne),"/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/ukb_chip_sex_updatedwith3kprotein.txt", sep="\t", quote=FALSE, row.names=TRUE)

tableOne <- CreateTableOne(vars = vars, data = chip_proteomics_pheno, strata="whiteblack", factorVars = factorVars)
write.table(print(tableOne),"/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/ukb_chip_race_updatedwith3kprotein.txt", sep="\t", quote=FALSE, row.names=TRUE)

########################################### Running Analysis ###################################################
########## Filter Out Missing Data for Proteomics #############
# Get the column names
col_names <- names(chip_proteomics_pheno)

# Replace "-" with "_" in the column names
new_col_names <- gsub("-", "_", col_names)

# Set the new column names
names(chip_proteomics_pheno) <- new_col_names

# Define protein columns
proteins <- colnames(chip_proteomics_pheno)[39:2961]

# Select relevant columns for analysis
proteomics <- chip_proteomics_pheno[, c("eid_7089", proteins)]

# Set row names to eid_7089 and remove the identifier column
rownames(proteomics) <- proteomics[, 1]
proteomics <- proteomics[, -1]

# Remove participants with more than 20% missing data and proteins with more than 10% missing data
proteomics <- proteomics[rowMeans(is.na(proteomics)) < 0.2, ]
proteomics <- proteomics[, colMeans(is.na(proteomics)) < 0.1]

print("Missing data filtering complete")

########## Combine CHIP Proteomics and Phenotype Data after Missingness Filter ############
# Identify protein columns that were kept after filtering
protein_columns_indices <- 39:2961
protein_names_in_proteomics <- setdiff(colnames(proteomics), "eid_7089")

# Keep columns relevant to analysis, excluding dropped proteins
all_columns_to_keep <- c(colnames(chip_proteomics_pheno)[-protein_columns_indices], protein_names_in_proteomics)

# Filter the master dataset for participants and select appropriate columns
chip_proteomics_pheno_remove1020missing <- chip_proteomics_pheno %>%
  filter(eid_7089 %in% row.names(proteomics)) %>%
  select(all_of(all_columns_to_keep))

print("CHIP proteomics and phenotype data cleaned and ready")

########## Running Analysis Models for Different Variants and Proteins ############
variants_list <- c("CHIP_CHIPvar", "DNMT3A_CHIPvar", "TET2_CHIPvar", "ASXL1_CHIPvar", "JAK2_CHIPvar", 
                   "expandedCHIP_CHIPvar", "expandedDNMT3A_CHIPvar", "expandedTET2_CHIPvar", 
                   "expandedASXL1_CHIPvar", "expandedJAK2_CHIPvar")

# Define adjustment models
Model <- c("+AGE_assessment+Sex+in_white_British_ancestry_subset+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10", 
           "+AGE_assessment+Sex+in_white_British_ancestry_subset+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+BMI+ever_smoked+Prev_Diabetes_Type_2", 
           "+AGE_assessment+Sex+in_white_British_ancestry_subset+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+BMI+ever_smoked+Prev_Diabetes_Type_2+egfr_cr")

# Define the corresponding model names
model <- c("minimal_adjust", "fully_adjust", "additionally_egfr")

# Start running the models
print("Start running models")

# Iterate through each dataset (starting with "chip_proteomics_pheno")
ds <- c("chip_proteomics_pheno")

for (j in seq_along(ds)) {
  d <- get(ds[j])
  for (m in seq_along(Model)) {    
    h <- data.frame(seq_id = character(), chip = character(), N = double(), beta = double(), se = double(),
                    z = double(), p = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)
    for (k in seq_along(variants_list)) {
      for (i in seq_along(proteins_to_keep)) {
        h[i, "seq_id"] <- proteins_to_keep[i]
        h[i, "chip"] <- variants_list[k]
        fmla <- as.formula(paste(proteins_to_keep[i], "~", variants_list[k], Model[m], sep = ""))
        f <- glm(fmla, data = d, na.action = na.exclude, family = "gaussian")
        h[i, "N"] <- summary(f)$df.null + 1
        h[i, "beta"] <- summary(f)$coefficient[2, 1]
        h[i, "se"] <- summary(f)$coefficient[2, 2]
        h[i, "z"] <- summary(f)$coefficient[2, 3]
        h[i, "p"] <- summary(f)$coefficient[2, 4]
        h[i, "lci"] <- h[i, "beta"] - 1.96 * h[i, "se"]
        h[i, "uci"] <- h[i, "beta"] + 1.96 * h[i, "se"]
      }
      
      if (k == 1) {
        h_merge <- h
      } else {
        h_merge <- rbind(h_merge, h)
      }
    }
    
    # Sort by p-values and write results to file
    h_merge <- h_merge[order(h_merge$p), ]
    write.table(h_merge, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_updatedwith3kprotein_UKB_proteomics_CHIP_", model[m], ".txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

########## Process Fully Adjusted Models and Apply FDR ############
h_merge <- fread("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_updatedwith3kprotein_UKB_proteomics_CHIP_fully_adjust.txt")

# Adjust the CHIP categories to match with Soma analysis categories
h_merge <- h_merge %>%
  mutate(chip = case_when(
    chip == "CHIP_CHIPvar" ~ "chip_all",
    chip == "DNMT3A_CHIPvar" ~ "DNMT3A_all",
    chip == "TET2_CHIPvar" ~ "TET2_all",
    chip == "ASXL1_CHIPvar" ~ "ASXL1_all",
    chip == "JAK2_CHIPvar" ~ "JAK2_all",
    chip == "expandedCHIP_CHIPvar" ~ "chip_large",
    chip == "expandedDNMT3A_CHIPvar" ~ "DNMT3A_large",
    chip == "expandedTET2_CHIPvar" ~ "TET2_large",
    chip == "expandedASXL1_CHIPvar" ~ "ASXL1_large",
    chip == "expandedJAK2_CHIPvar" ~ "JAK2_large",
    TRUE ~ chip # Ensures values not matched above remain unchanged
  ))

# Concatenate chip and seq_id to create unique identifiers for grouping
h_merge$Chip_cat <- paste0(substr(h_merge$chip, 1, 3), h_merge$seq_id)

# Select the minimum p-value for each Chip_cat group
h_merge <- h_merge %>%
  group_by(Chip_cat) %>%
  slice(which.min(p))

# Separate JAK2 variants for additional processing
h_merge$Chip_use <- sub("_.*", "", h_merge$chip)
h_merge_nojak2 <- h_merge %>%
  ungroup() %>%
  filter(Chip_use != "JAK2") %>%
  mutate(fdr_fourgroup = p.adjust(p, method = "fdr"))

# Keep the fdr for four groups of CHIP and merge back
h_merge_nojak2 <- h_merge_nojak2 %>%
  select(seq_id, chip, fdr_fourgroup)

# Apply FDR adjustment across all CHIP variants
h_merge <- h_merge %>%
  ungroup() %>%
  mutate(fdr_all = p.adjust(p, method = "fdr"))

# Merge the FDR results back
h_merge <- merge(h_merge, h_merge_nojak2, by = c("seq_id", "chip"), all.x = TRUE)

# Write the processed results to a file
write.table(h_merge, file = "/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_calculatedforfinaluse_updatedwith3kprotein_UKB_proteomics_CHIP_fully_adjust.txt", sep = "\t", quote = FALSE, row.names = FALSE)

print("FDR adjustments complete and results saved")

########## Peer Adjustment Model Analysis ############
# Load PEER factors data
peer_150 <- fread("/medpop/esp2/zyu/chip_protemoics/data/UKB/PEER/peernum/peerandresid/peers.UKB_150.txt")

# Merge PEER factors with the cleaned CHIP-proteomics data
chip_proteomics_pheno_remove1020missing <- merge(chip_proteomics_pheno_remove1020missing, peer_150, by.x = "eid_7089", by.y = "SampleId")

# Generate a list of PEER factor covariates
peer_factors <- paste0("peer", 1:150)

# Combine base covariates with PEER factors
base_covariates <- c("AGE_assessment", "Sex", "in_white_British_ancestry_subset", "PC1", "PC2", "PC3", "PC4", "PC5", 
                     "PC6", "PC7", "PC8", "PC9", "PC10", "BMI", "ever_smoked", "Prev_Diabetes_Type_2")
full_model_covariates <- c(base_covariates, peer_factors)

# Convert the full covariate list to formula syntax
model_formula <- paste("+", paste(full_model_covariates, collapse = "+"))

# Define models and variants for analysis
Model <- c(model_formula)
model <- c("peer_adjust")

# Run the models
print("Start running Peer-adjusted models")

ds <- c("chip_proteomics_pheno_remove1020missing")
for (j in 1:length(ds)) {
  d <- get(ds[j])
  for (m in 1:length(Model)) {    
    h <- data.frame(seq_id = character(), chip = character(), N = double(), beta = double(), se = double(),
                    z = double(), p = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)
    for (k in 1:length(variants_list)) {
      for (i in 1:length(proteins_to_keep)) {
        h[i, "seq_id"] <- proteins_to_keep[i]
        h[i, "chip"] <- variants_list[k]
        fmla <- as.formula(paste(proteins_to_keep[i], "~", variants_list[k], Model[m], sep = ""))
        f <- glm(fmla, data = d, na.action = na.exclude, family = "gaussian")
        h[i, "N"] <- summary(f)$df.null + 1
        h[i, "beta"] <- summary(f)$coefficient[2, 1]
        h[i, "se"] <- summary(f)$coefficient[2, 2]
        h[i, "z"] <- summary(f)$coefficient[2, 3]
        h[i, "p"] <- summary(f)$coefficient[2, 4]
        h[i, "lci"] <- h[i, "beta"] - 1.96 * h[i, "se"]
        h[i, "uci"] <- h[i, "beta"] + 1.96 * h[i, "se"]
      }
      
      if (k == 1) {
        h_merge <- h
      } else {
        h_merge <- rbind(h_merge, h)
      }
    }
    
    # Sort and write results to file
    h_merge <- h_merge[order(h_merge$p), ]
    write.table(h_merge, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_updatedwith3kprotein_withPEER150_UKB_proteomics_CHIP_", model[m], ".txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

print("Peer adjustment analysis complete")

########## CHIP Variant Allele Frequency (VAF) Analysis ############
variants_list <- c("CHIP_VAF_CHIPvar", "DNMT3A_VAF_CHIPvar", "TET2_VAF_CHIPvar", "ASXL1_VAF_CHIPvar", "JAK2_VAF_CHIPvar")

# Run the analysis for CHIP VAF variants
for (j in 1:length(ds)) {
  d <- get(ds[j])
  for (m in 1:length(Model)) {    
    h <- data.frame(seq_id = character(), chip = character(), N = double(), beta = double(), se = double(),
                    z = double(), p = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)
    for (k in 1:length(variants_list)) {
      for (i in 1:length(proteins_to_keep)) {
        h[i, "seq_id"] <- proteins_to_keep[i]
        h[i, "chip"] <- variants_list[k]
        fmla <- as.formula(paste(proteins_to_keep[i], "~", variants_list[k], Model[m], sep = ""))
        f <- glm(fmla, data = d, na.action = na.exclude, family = "gaussian")
        h[i, "N"] <- summary(f)$df.null + 1
        h[i, "beta"] <- summary(f)$coefficient[2, 1]
        h[i, "se"] <- summary(f)$coefficient[2, 2]
        h[i, "z"] <- summary(f)$coefficient[2, 3]
        h[i, "p"] <- summary(f)$coefficient[2, 4]
        h[i, "lci"] <- h[i, "beta"] - 1.96 * h[i, "se"]
        h[i, "uci"] <- h[i, "beta"] + 1.96 * h[i, "se"]
      }
      
      if (k == 1) {
        h_merge <- h
      } else {
        h_merge <- rbind(h_merge, h)
      }
    }
    
    # Sort and write results to file
    h_merge <- h_merge[order(h_merge$p), ]
    write.table(h_merge, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_VAF_updatedwith3kprotein_withPEER150_UKB_proteomics_CHIP_", model[m], ".txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

print("VAF analysis complete")

########## Caucasian-Only Analysis ############
# Subset data to white participants
chip_proteomics_pheno_remove1020missing_white <- chip_proteomics_pheno_remove1020missing[chip_proteomics_pheno_remove1020missing$whiteblack == "white", ]

# Specify dataset
ds <- c("chip_proteomics_pheno_remove1020missing_white")

# Run the model for Caucasian-only participants
for (j in 1:length(ds)) {
  d <- get(ds[j])
  for (m in 1:length(Model)) {    
    h <- data.frame(seq_id = character(), chip = character(), N = double(), beta = double(), se = double(),
                    z = double(), p = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)
    for (k in 1:length(variants_list)) {
      for (i in 1:length(proteins_to_keep)) {
        h[i, "seq_id"] <- proteins_to_keep[i]
        h[i, "chip"] <- variants_list[k]
        fmla <- as.formula(paste(proteins_to_keep[i], "~", variants_list[k], Model[m], sep = ""))
        f <- glm(fmla, data = d, na.action = na.exclude, family = "gaussian")
        h[i, "N"] <- summary(f)$df.null + 1
        h[i, "beta"] <- summary(f)$coefficient[2, 1]
        h[i, "se"] <- summary(f)$coefficient[2, 2]
        h[i, "z"] <- summary(f)$coefficient[2, 3]
        h[i, "p"] <- summary(f)$coefficient[2, 4]
        h[i, "lci"] <- h[i, "beta"] - 1.96 * h[i, "se"]
        h[i, "uci"] <- h[i, "beta"] + 1.96 * h[i, "se"]
      }
      
      if (k == 1) {
        h_merge <- h
      } else {
        h_merge <- rbind(h_merge, h)
      }
    }
    
    # Apply FDR corrections and save the results
    h_merge <- h_merge[order(h_merge$p), ]
    h_merge <- h_merge %>%
      mutate(chip = case_when(
        chip == "CHIP_VAF_CHIPvar" ~ "chip_all",
        chip == "DNMT3A_VAF_CHIPvar" ~ "DNMT3A_all",
        chip == "TET2_VAF_CHIPvar" ~ "TET2_all",
        chip == "ASXL1_VAF_CHIPvar" ~ "ASXL1_all",
        chip == "JAK2_VAF_CHIPvar" ~ "JAK2_all",
        TRUE ~ chip
      ))
    h_merge$Chip_cat <- paste0(substr(h_merge$chip, 1, 3), h_merge$seq_id)
    h_merge <- h_merge %>%
      group_by(Chip_cat) %>%
      slice(which.min(p))
    h_merge$Chip_use <- sub("_.*", "", h_merge$chip)
    h_merge_nojak2 <- h_merge %>%
      ungroup() %>%
      filter(Chip_use != "JAK2") %>%
      mutate(fdr_fourgroup = p.adjust(p, method = "fdr"))
    h_merge_nojak2 <- h_merge_nojak2 %>%
      select(seq_id, chip, fdr_fourgroup)
    h_merge <- h_merge %>%
      ungroup() %>%
      mutate(fdr_all = p.adjust(p, method = "fdr"))
    h_merge <- merge(h_merge, h_merge_nojak2, by = c("seq_id", "chip"), all.x = TRUE)

    # Save Caucasian-only results
    write.table(h_merge, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_VAF_updatedwith3kprotein_withPEER150_UKB_whiteonly_proteomics_CHIP_", model[m], ".txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

print("Caucasian-only analysis complete")

########## Male and Female-Specific Analysis ############
# Define datasets for male and female participants
ds <- c("Female", "Male")

# Loop over male and female datasets
for (j in seq_along(ds)) {
  if (ds[j] == "Female") {
    d <- chip_proteomics_pheno_remove1020missing[chip_proteomics_pheno_remove1020missing$Sex == "Female", ]
  } else if (ds[j] == "Male") {
    d <- chip_proteomics_pheno_remove1020missing[chip_proteomics_pheno_remove1020missing$Sex == "Male", ]
  }
  
  for (m in 1:length(Model)) {    
    h <- data.frame(seq_id = character(), chip = character(), N = double(), beta = double(), se = double(),
                    z = double(), p = double(), lci = double(), uci = double(), stringsAsFactors = FALSE)
    for (k in 1:length(variants_list)) {
      for (i in 1:length(proteins_to_keep)) {
        h[i, "seq_id"] <- proteins_to_keep[i]
        h[i, "chip"] <- variants_list[k]
        fmla <- as.formula(paste(proteins_to_keep[i], "~", variants_list[k], Model[m], sep = ""))
        f <- glm(fmla, data = d, na.action = na.exclude, family = "gaussian")
        h[i, "N"] <- summary(f)$df.null + 1
        h[i, "beta"] <- summary(f)$coefficient[2, 1]
        h[i, "se"] <- summary(f)$coefficient[2, 2]
        h[i, "z"] <- summary(f)$coefficient[2, 3]
        h[i, "p"] <- summary(f)$coefficient[2, 4]
        h[i, "lci"] <- h[i, "beta"] - 1.96 * h[i, "se"]
        h[i, "uci"] <- h[i, "beta"] + 1.96 * h[i, "se"]
      }
      
      if (k == 1) {
        h_merge <- h
      } else {
        h_merge <- rbind(h_merge, h)
      }
    }
    
    # Apply FDR corrections and save the results
    h_merge <- h_merge[order(h_merge$p), ]
    h_merge <- h_merge %>%
      mutate(chip = case_when(
        chip == "CHIP_CHIPvar" ~ "chip_all",
        chip == "DNMT3A_CHIPvar" ~ "DNMT3A_all",
        chip == "TET2_CHIPvar" ~ "TET2_all",
        chip == "ASXL1_CHIPvar" ~ "ASXL1_all",
        chip == "JAK2_CHIPvar" ~ "JAK2_all",
        chip == "expandedCHIP_CHIPvar" ~ "chip_large",
        chip == "expandedDNMT3A_CHIPvar" ~ "DNMT3A_large",
        chip == "expandedTET2_CHIPvar" ~ "TET2_large",
        chip == "expandedASXL1_CHIPvar" ~ "ASXL1_large",
        chip == "expandedJAK2_CHIPvar" ~ "JAK2_large",
        TRUE ~ chip
      ))
    h_merge$Chip_cat <- paste0(substr(h_merge$chip, 1, 3), h_merge$seq_id)
    h_merge <- h_merge %>%
      group_by(Chip_cat) %>%
      slice(which.min(p))
    h_merge$Chip_use <- sub("_.*", "", h_merge$chip)
    h_merge_nojak2 <- h_merge %>%
      ungroup() %>%
      filter(Chip_use != "JAK2") %>%
      mutate(fdr_fourgroup = p.adjust(p, method = "fdr"))
    h_merge_nojak2 <- h_merge_nojak2 %>%
      select(seq_id, chip, fdr_fourgroup)
    h_merge <- h_merge %>%
      ungroup() %>%
      mutate(fdr_all = p.adjust(p, method = "fdr"))
    h_merge <- merge(h_merge, h_merge_nojak2, by = c("seq_id", "chip"), all.x = TRUE)

    # Save male/female-specific results
    write.table(h_merge, file = paste("/medpop/esp2/zyu/chip_protemoics/output/chipproteomicscad/Final_updatedwith3kprotein_withPEER150_UKB_proteomics_CHIP_", model[m], "_", ds[j], ".txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

print("Male and Female-specific analysis complete")


