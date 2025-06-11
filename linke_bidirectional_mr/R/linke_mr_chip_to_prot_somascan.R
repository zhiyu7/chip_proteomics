######## prep ########
library(R.utils)
library(TwoSampleMR)
library(data.table)
library(ieugwasr)
library(stringr)
library(mr.raps)
library(tidyr)
library(dplyr)
library(data.table)
library(optparse)



gc()
rm(list = ls())


# inputting file for shell script
option_list = list(
    make_option(c("-p", "--protein"), type="character", default=NULL,
                help="the protein", metavar="character")
    
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# # the option for running for gene using bash script in g_ shell script
protein <- opt$protein




#### 02 reading in clumped summary stats #####
# mutating chr and bp position for data formatting, then remove duplicates based on the exact marker (duplicated snps)
chip <- data.frame(fread("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/clumped_chip_ss/CHIP_2022ss_clumped.tsv"))
tet2 <- data.frame(fread("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/clumped_chip_ss/TET2_2022ss_clumped.tsv"))
dnmt3a <- data.frame(fread("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/clumped_chip_ss/DNMT3A_2022ss_clumped.tsv"))





######### 03 - data prep ##########
df_sum <- data.frame(file=NA, exp=NA, outc=NA, nsnp=NA, method=NA, b=NA, se=NA, pval=NA)[-1,]
df_instr <- data.frame(file =NA,pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                       other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                       palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                       outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                       pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                       action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA, marker_ld=NA)[-1,]



# all the files located for pQTL
all_files <- list.files("/broad/hptmp/lli", pattern = "\\.txt\\.gz$", full.names = TRUE)


######## looping #######
protein_tag = protein
# protein_tag <- strsplit(protein, "_")[[1]][3]
# fetch the file from all the file in the pathway
file_match <- all_files[grepl(paste0("_",protein_tag,"_"), all_files)]
# incase there are multiple files for the same protein
for (file in file_match){
    # reading in the file
    chrom_u <- data.frame(fread(file))
    
    # selecting column and renaming them
    chrom_u <- chrom_u %>%
        select(Chrom, Pos, otherAllele, effectAllele, rsids, Pval, Beta, SE, ImpMAF) %>%
        rename(CHROM = Chrom,
               GENPOS = Pos,
               ALLELE0 = otherAllele, ALLELE1 = effectAllele,
               SNP = rsids, P = Pval, BETA = Beta, SE = SE, A1FREQ = ImpMAF) %>%
        mutate(CHROM = as.numeric(gsub("chr", "", CHROM)),
               marker = paste(CHROM, GENPOS, ALLELE0, ALLELE1, sep = ":"),
               marker1 = paste(CHROM, GENPOS, ALLELE1, ALLELE0, sep = ":")
               )
    
    if (is.null(chrom_u) || nrow(chrom_u) == 0) {
        print(paste0("Skipping ", protein, "(1)"))
    } 
    else {
        for (j in c("chip", "dnmt3a", "tet2")) {
            # get the summary statistics
            outcome <- get(j)
            # Only selecting the chromosome of interest to speed stuff up downstream from here
            outcome_overlap <- outcome[outcome$SNP %in% chrom_u$SNP,]
            nrow(outcome_overlap)
            
            if (is.null(outcome_overlap) || nrow(outcome_overlap) == 0) {
                print(paste0("Skipping ", protein, "(2)"))
            } else {
                # if it's not empty
                # setting up recording data frame
                outcome_overlap$phen <- paste(j)
                # formatting exposure data
                outcome_overlap <- format_data(outcome_overlap, type="exposure", phenotype_col="phen", 
                                               snp_col="SNP", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                               effect_allele_col="alt", other_allele_col="ref", pval_col="pval", 
                                               chr_col="chr", pos_col="bp")
                
                # reducing dimension of outcome ss(pQTL) and adding rsID as SNP column
                chrom_overlap <- chrom_u %>%
                    filter(SNP %in% outcome_overlap$SNP) %>%
                    as.data.frame()
                chrom_overlap$phen <- protein
                chrom_overlap <- format_data(chrom_overlap, type="outcome", phenotype_col="phen", 
                                             snp_col="SNP", beta_col="BETA", se_col="SE", eaf_col="A1FREQ",
                                             effect_allele_col="ALLELE1", other_allele_col="ALLELE0", 
                                             pval_col="P", chr_col="CHROM", pos_col="GENPOS")
                # harmonize data
                dat <- harmonise_data(exposure_dat=outcome_overlap, outcome_dat=chrom_overlap)                                             # This is where the matching happens
                
                dat <- dat[order(dat$pval.exposure),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
                dat <- dat[!duplicated(dat$SNP),]
                dat <- dat[dat$mr_keep,]
                
                if (nrow(dat[dat$mr_keep,])==0) {
                    print(paste0("Skipping ",  protein, "(3)"))
                    results <- NULL
                } else {
                    # if there's only one SNP
                    if (nrow(dat)==1) {                                                                                                        # This is where the magic happens: if you have 1 variant, you use the Wald ratio as your method
                        results_mr <- mr(dat, method_list=c("mr_wald_ratio"))
                        results <- data.frame(file=file, exp=paste(j), outc=protein_tag, nsnp=results_mr$nsnp, method="Wald_ratio", b=results_mr$b, 
                                              se=results_mr$se, pval=results_mr$pval)
                    } else if (nrow(dat)==2) {                                                                                                # If you have 2 variants, you can use the classic IVW method but not the MR-Egger method
                        # if there are 2
                        output_mr_ivw <- mr(dat, method_list=c("mr_ivw"))
                        results_1 <- data.frame(file=file, exp=paste(j), outc=protein_tag, nsnp=output_mr_ivw$nsnp, 
                                                method="IVW", b=output_mr_ivw$b, 
                                                se=output_mr_ivw$se, pval=output_mr_ivw$pval)
                        output_mr_raps <- mr.raps(data.frame(beta.exposure=dat$beta.exposure, beta.outcome=dat$beta.outcome, 
                                                             se.exposure=dat$se.exposure, se.outcome=dat$se.outcome))
                        results_2 <- data.frame(file=file, exp=paste(j), outc=protein_tag, nsnp=output_mr_ivw$nsnp, 
                                                method="MR-RAPS", b=output_mr_raps$beta.hat, 
                                                se=output_mr_raps$beta.se, pval=2 * pnorm(-abs(output_mr_raps$beta.hat / output_mr_raps$beta.se)))
                        results <- rbind(results_1, results_2)
                        
                    } else {                                                                                                                  # If you have more than 2 variants, you can do anything (including IVW and MR-Egger)
                        output_mr_ivw <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression"))
                        results_1 <- data.frame(file=file, exp=paste(j), outc=protein_tag, nsnp=output_mr_ivw$nsnp, 
                                                method=c("IVW", "Egger"), b=output_mr_ivw$b, 
                                                se=output_mr_ivw$se, pval=output_mr_ivw$pval)
                        output_mr_egger <- mr_egger_regression(b_exp=dat$beta.exposure, b_out=dat$beta.outcome, se_exp=dat$se.exposure, se_out=dat$se.outcome)
                        results_2 <- data.frame(file=file, exp=paste(j), outc=protein_tag, nsnp=output_mr_ivw$nsnp[1], 
                                                method=c("Egger_int"), b=output_mr_egger$b_i, 
                                                se=output_mr_egger$se_i, pval=output_mr_egger$pval_i)
                        output_mr_raps <- mr.raps(data.frame(beta.exposure=dat$beta.exposure, beta.outcome=dat$beta.outcome, 
                                                             se.exposure=dat$se.exposure, se.outcome=dat$se.outcome))
                        results_3 <- data.frame(file=file, exp=paste(j), outc=protein_tag, nsnp=output_mr_ivw$nsnp[1], 
                                                method="MR-RAPS", b=output_mr_raps$beta.hat, 
                                                se=output_mr_raps$beta.se, pval=2 * pnorm(-abs(output_mr_raps$beta.hat / output_mr_raps$beta.se)))
                        results <- rbind(results_1, results_2, results_3)
                        
                    }
                }
                if (is.null(results) || nrow(results) == 0) {
                    print(paste0("Skipping ",  protein, "(4)"))
                } else {
                    df_sum <- rbind(df_sum, results)
                    df_instr <- rbind(df_instr, dat)
                }
            }
        }
    }
    # writing the output once per protein
    write.csv(df_sum, paste0("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/somascan_chip_to_22prot_raw/", protein,"_3chip_raw.csv"), row.names = F, quote = F) 
    write.csv(df_instr, paste0("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/somascan_chip_to_22prot_raw/",protein,"_3chip_instr_raw.csv"), row.names = F, quote = F) 
    print(paste(protein, "done"))
}

