######## 01 - prep ########
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
    make_option(c("-s", "--sumstats"), type="character", default=NULL, 
                help="the summary statistics of CHIP/TET2/DNMT3A", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# # the option for running for gene using bash script in g_ shell script
sumstats <- opt$sumstats


## testing section ##
# sumstats <- "CHIP"

##### 02 - files ######
# a mini loop to check smaller protein upon request
protein_list <- fread("/medpop/esp2/lli/Zhi-mr_protein/Data/somascan_protein_ss/new_protlist", header = F)$V1
all_files <- list.files("/broad/hptmp/lli", pattern = "\\.txt\\.gz$", full.names = TRUE)


# reading in summary stats for CHIPs
ss_df <- data.frame(fread(paste0("/medpop/esp2/mesbah/projects/Meta_GWAS/n650k/eur_only_summary/lifted_hg37.eur_metaGWAS.", sumstats, ".GWAMA.hg37_dbSNP.eaf001_min2Studies.tsv.gz"))) %>%
    rename(SNP = MarkerID,
           beta = BETA,
           sebeta = SE,
           pval = P,
           ref = REF,
           alt = ATL,
           af_alt = EAF) %>%
    select(SNP, beta, sebeta, pval, ref,alt,af_alt)
print(paste0("Successfully reading in summary statistics of ", sumstats))



###### 03 - loops ######
# preparing empty dataframe to store the results
df_sum <- data.frame(file=NA, exp=NA, outc=NA, nsnp=NA, method=NA, b=NA, se=NA, pval=NA)[-1,]
df_instr <- data.frame(file =NA,pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                       other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                       palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                       outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                       pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                       action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA, marker_ld=NA)[-1,]


# the protein tag
######## looping #######
for (protein in protein_list){
    
    protein_tag <- protein
    # protein_tag <- strsplit(protein, "_")[[1]][3]
    
    # fetch the file from all the file in the pathway
    file_match <- all_files[grepl(paste0("_",protein_tag,"_"), all_files)]
    
    # incase there are multiple files for the same protein
    for (file in file_match){
        chrom_u <- data.frame(fread(file))
        
        # selecting column and renaming them
        chrom_u <- chrom_u %>%
            filter(Pval < 5e-6) %>%
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
            print(paste0("Skipping ", protein, ", pQTL file is empty after clumping"))
        } else {
            
            # it already has rsID so can be clumped directly
            clump <- ld_clump(dplyr::tibble(rsid=chrom_u$SNP, pval=chrom_u$P),
                              plink_bin = "/medpop/esp2/lli/Mike-METAP1/TwoSampleMR/plink",
                              clump_kb=10000, clump_r2 = 0.1,
                              bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
            # filtering to obtain the clumped data
            chrom_u <- chrom_u[chrom_u$SNP %in% clump$rsid, ]
            
            # Only selecting the chromosome of interest to speed stuff up downstream from here
            outcome_overlap <- ss_df[ss_df$SNP %in% chrom_u$SNP,]
            nrow(outcome_overlap)
            
            if (is.null(outcome_overlap) || nrow(outcome_overlap) == 0) {
                print(paste0("Skipping ", protein, ", no overlap between exposure and outcome data"))
            } else {
                # if it's not empty
                # setting up recording data frame
                outcome_overlap$phen <- paste(sumstats)
                # formatting exposure data
                outcome_overlap <- format_data(outcome_overlap, type="outcome", phenotype_col="phen", 
                                               snp_col="SNP", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                               effect_allele_col="alt", other_allele_col="ref", pval_col="pval", 
                                               chr_col="chr", pos_col="bp")
                
                # reducing dimension of outcome ss(pQTL) and adding rsID as SNP column
                # same here
                chrom_overlap <- chrom_u %>%
                    filter(SNP %in% outcome_overlap$SNP) %>%
                    as.data.frame()
                # adding the phenotype
                chrom_overlap$phen <- protein
                chrom_overlap <- format_data(chrom_overlap, type="exposure", phenotype_col="phen", 
                                             snp_col="SNP", beta_col="BETA", se_col="SE", eaf_col="A1FREQ",
                                             effect_allele_col="ALLELE1", other_allele_col="ALLELE0", 
                                             pval_col="P", chr_col="CHROM", pos_col="GENPOS")
                # harmonize data, handles strandflipping
                dat <- harmonise_data(exposure_dat=chrom_overlap, outcome_dat=outcome_overlap)                                  
                
                dat <- dat[order(dat$pval.exposure),]
                dat <- dat[!duplicated(dat$SNP),]
                dat <- dat[dat$mr_keep,]
                
                # if theres no data left after harmonizing
                if (nrow(dat[dat$mr_keep,])==0) {
                    print(paste0("Skipping ",  protein, ", no data after harmonizing"))
                    results <- NULL
                } else {
                    # if there's only one SNP
                    if (nrow(dat)==1) {                                                                                                        # This is where the magic happens: if you have 1 variant, you use the Wald ratio as your method
                        results_mr <- mr(dat, method_list=c("mr_wald_ratio"))
                        results <- data.frame(file=file, exp=protein, outc=protein_tag, nsnp=results_mr$nsnp, method="Wald_ratio", b=results_mr$b, 
                                              se=results_mr$se, pval=results_mr$pval)
                    } else if (nrow(dat)==2) {                                                                                                # If you have 2 variants, you can use the classic IVW method but not the MR-Egger method
                        # if there are 2
                        output_mr_ivw <- mr(dat, method_list=c("mr_ivw"))
                        results_1 <- data.frame(file=file, exp=protein, outc=protein_tag, nsnp=output_mr_ivw$nsnp, 
                                                method="IVW", b=output_mr_ivw$b, 
                                                se=output_mr_ivw$se, pval=output_mr_ivw$pval)
                        # output_mr_raps <- mr.raps(data.frame(beta.exposure=dat$beta.exposure, beta.outcome=dat$beta.outcome, 
                        #                                      se.exposure=dat$se.exposure, se.outcome=dat$se.outcome))
                        # results_2 <- data.frame(file=file, exp=protein, outc=protein_tag, nsnp=output_mr_ivw$nsnp, 
                        #                         method="MR-RAPS", b=output_mr_raps$beta.hat, 
                        #                         se=output_mr_raps$beta.se, pval=2 * pnorm(-abs(output_mr_raps$beta.hat / output_mr_raps$beta.se)))
                        results <- rbind(results_1)
                        
                    } else {
                        # If you have more than 2 variants, you can do anything (including IVW and MR-Egger)
                        output_mr_ivw <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression"))
                        results_1 <- data.frame(file=file, exp=protein_tag, outc=sumstats, nsnp=output_mr_ivw$nsnp, 
                                                method=c("IVW", "Egger"), b=output_mr_ivw$b, 
                                                se=output_mr_ivw$se, pval=output_mr_ivw$pval)
                        output_mr_egger <- mr_egger_regression(b_exp=dat$beta.exposure, b_out=dat$beta.outcome, se_exp=dat$se.exposure, se_out=dat$se.outcome)
                        results_2 <- data.frame(file=file, exp=protein_tag, outc=sumstats, nsnp=output_mr_ivw$nsnp[1], 
                                                method=c("Egger_int"), b=output_mr_egger$b_i, 
                                                se=output_mr_egger$se_i, pval=output_mr_egger$pval_i)
                        output_mr_raps <- mr.raps(data.frame(beta.exposure=dat$beta.exposure, beta.outcome=dat$beta.outcome,
                                                             se.exposure=dat$se.exposure, se.outcome=dat$se.outcome))
                        results_3 <- data.frame(file=file, exp=protein_tag, outc=sumstats, nsnp=output_mr_ivw$nsnp[1],
                                                method="MR-RAPS", b=output_mr_raps$beta.hat,
                                                se=output_mr_raps$beta.se, pval=2 * pnorm(-abs(output_mr_raps$beta.hat / output_mr_raps$beta.se)))
                        results <- rbind(results_1, results_2, results_3)
                        
                    }
                }
                if (is.null(results) || nrow(results) == 0) {
                    print(paste0("Skipping ",  protein, ", results are empty"))
                } else {
                    df_sum <- rbind(df_sum, results)
                    df_instr <- rbind(df_instr, dat)
                }
            }
        }
        # writing the output once per protein
        write.csv(df_sum, paste0("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/somascan_22prot_to_chip_raw/mr_prot_to_", sumstats,"_somascan.csv"), row.names = F, quote = F) 
        write.csv(df_instr, paste0("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/somascan_22prot_to_chip_raw/mr_prot_to_",sumstats,"_somascan_instruments.csv"), row.names = F, quote = F)
        print(paste(protein, "done"))
    }
}