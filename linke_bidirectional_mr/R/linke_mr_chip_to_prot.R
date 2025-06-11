######### 01 - prep ########
library(reticulate)
use_condaenv("synapser_env", required = TRUE) ### My own conda environment configured to use all the packages
library(synapser) 
library(R.utils)
library(TwoSampleMR)
library(tidyr)
library(dplyr)
library(data.table)
library(stringr)
library(ieugwasr)
library(MendelianRandomization)
library(mr.raps)
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


## testing area
# protein <- "FLT3LG"


# login information for synapse, needs email and PAT
synLogin(email=, 
         authToken=)



######## 02 - reading in the summary stats ##########
# mutating chr and bp position for data formatting, then remove duplicates based on the exact marker (duplicated snps)
# when CHIP function as instruments/exposure data, use preclumped version with rsID
chip <- data.frame(fread("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/clumped_chip_ss/CHIP_2022ss_clumped.tsv"))
tet2 <- data.frame(fread("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/clumped_chip_ss/TET2_2022ss_clumped.tsv"))
dnmt3a <- data.frame(fread("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/clumped_chip_ss/DNMT3A_2022ss_clumped.tsv"))




############ 03 - looping ###########
# preparing files
sumstats_info = fread("/medpop/esp2/aschuerm/ukb_proteomics_cvd/input_files/olink_protein_map_3k_v1.tsv") %>%
    filter(Assay %in% protein)
## !!!!!!! this is only because I checked the 2 proteins with multiple rows were the ones having the same docname !!!!!!!! ###
sumstats_info <- sumstats_info[!duplicated(sumstats_info$Docname),]

# obtaining the synapser ID for all the proteins

# empty dataframes to store the data
# MR reuslts
df_sum <- data.frame(exp=NA, outc=NA, nsnp=NA, method=NA, b=NA, se=NA, pval=NA)[-1,]
# instrument used
df_instr <- data.frame(pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                       other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                       palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                       outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                       pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                       action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA)[-1,]


#### for each of the protein in the list
# Downloading the summary statistics for the protein of interest
path = "/medpop/esp2/aschuerm/misc/others/zhi/"
# create unique unpack directory for files, 
# as some proteins shares the same SS and can conflict when running parallele jobs
unique_unpack_dir <- paste0(path, "sumstats_", protein, "_", Sys.getpid())
dir.create(unique_unpack_dir, showWarnings = FALSE, recursive = TRUE)
# download the file into unique directory, and unpack in that origin
syn_code <- synGet(entity = sumstats_info$Code, downloadLocation = unique_unpack_dir)
untar(syn_code$path, exdir = unique_unpack_dir)

# appending 22 chromosomes as CHIP signals comes from all chromosomes
bm_whole <- data.frame(matrix(ncol = 12, nrow = 0))
for (nr in 1:22){
    chrom_part <- fread(paste0(syn_code$cacheDir, "/", gsub(".tar", "", sumstats_info$Docname), "/", 
                               "discovery_chr", nr, "_", sumstats_info$UKBPPP_ProteinID, 
                               ":", sumstats_info$Panel, ".gz"))
    bm_whole <- rbind(bm_whole, chrom_part)
    print(paste("chromosome", nr, "bound"))
}

# processing data
bm_whole <- bm_whole %>%
    mutate(P = 10^-LOG10P, 
           marker = paste(CHROM, GENPOS, ALLELE0, ALLELE1, sep=":"),
           marker1 = paste(CHROM, GENPOS, ALLELE1, ALLELE0, sep=":")
    )


# if the protein file isn't empty
if (is.null(bm_whole) || nrow(bm_whole) == 0) {
    # print an error message
    print(paste0("Skipping ", sumstats_info$Assay, "protein summary statistics file is empty"))
    } else {
        
        # looping over the 3 summary statistics file
        for (j in c("chip", "dnmt3a", "tet2")) {
            # Only selecting the chromosome of interest to speed stuff up downstream from here
            outcome <- get(j)
            outcome_overlap <- outcome[outcome$marker %in% bm_whole$marker |
                                       outcome$marker %in% bm_whole$marker1,]
            nrow(outcome_overlap)

        if (is.null(outcome_overlap) || nrow(outcome_overlap) == 0) {
            print(paste0("Skipping ", protein, "no overlap between CHIP and protein SS"))
        } else {
            # setting up recording data frame
            outcome_overlap$phen <- paste(j)
            rsid <- outcome_overlap[,c("marker", "SNP")] # dictionary for SNP rsID
            # formatting exposure data
            outcome_overlap <- format_data(outcome_overlap, type="exposure", phenotype_col="phen", 
                                           snp_col="SNP", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                           effect_allele_col="alt", other_allele_col="ref", pval_col="pval", 
                                           chr_col="chr", pos_col="bp")
            
            # reducing dimention of outcome ss(pQTL) and adding rsID as SNP column
            chrom_overlap <- bm_whole %>%
                filter(marker %in% rsid$marker | marker1 %in% rsid$marker) %>%
                left_join(rsid, by = "marker") %>%                
                mutate(SNP_alt = SNP) %>%                                   
                select(-SNP) %>%
                left_join(rsid, by = c("marker1" = "marker")) %>%               
                mutate(
                    SNP = coalesce(SNP, SNP_alt)                              
                ) %>%
                select(-SNP_alt) %>%
                as.data.frame()
            # format outcome data
            chrom_overlap$phen <- protein
            chrom_overlap <- format_data(chrom_overlap, type="outcome", phenotype_col="phen", 
                                         snp_col="SNP", beta_col="BETA", se_col="SE", eaf_col="A1FREQ",
                                         effect_allele_col="ALLELE1", other_allele_col="ALLELE0", 
                                         pval_col="LOG10P", chr_col="CHROM", samplesize_col="N", pos_col="GENPOS", log_pval=T)
            
            # harmonizing data
            dat <- harmonise_data(exposure_dat=outcome_overlap, outcome_dat=chrom_overlap)                                             
            dat <- dat[order(dat$pval.exposure),]                                                                                   
            dat <- dat[!duplicated(dat$SNP),]
            dat <- dat[dat$mr_keep,]
 
            # if there's no harmonized data,skipping the result of the loop and make the results as null
            if (nrow(dat[dat$mr_keep,])==0) {
                print(paste0("Skipping ", protein, ", no data left after harmonising"))
                results <- NULL
            } else {
            
                # if there is only one SNP/row
                if (nrow(dat)==1) { 
                    # one SNP can only run wald ratio
                    results_mr <- mr(dat, method_list=c("mr_wald_ratio"))
                    results <- data.frame(exp=protein, outc=paste("outcome"), 
                                          nsnp=results_mr$nsnp, method=results_mr$method, b=results_mr$b, 
                                          se=results_mr$se, pval=results_mr$pval)
                    
                } else if (nrow(dat)==2) {
                    # If you have 2 variants, you can use the classic IVW method but not the MR-Egger method
                    output_mr_ivw <- mr(dat, method_list=c("mr_ivw"))
                    results_1 <- data.frame(exp=paste(j), outc=protein, nsnp=output_mr_ivw$nsnp, 
                                            method="IVW", b=output_mr_ivw$b, 
                                            se=output_mr_ivw$se, pval=output_mr_ivw$pval)
                    output_mr_raps <- mr.raps(data.frame(beta.exposure=dat$beta.exposure, beta.outcome=dat$beta.outcome, 
                                                         se.exposure=dat$se.exposure, se.outcome=dat$se.outcome))
                    results_2 <- data.frame(exp=paste(j), outc=protein, nsnp=output_mr_ivw$nsnp, 
                                            method="MR-RAPS", b=output_mr_raps$beta.hat, 
                                            se=output_mr_raps$beta.se, pval=2 * pnorm(-abs(output_mr_raps$beta.hat / output_mr_raps$beta.se)))
                    # binding the results of MR-RAPS and MR IVW
                    results <- rbind(results_1, results_2)
                    
                } else {               
                # If you have more than 2 variants, you can do anything (including IVW and MR-Egger) 
                    output_mr_ivw <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression"))
                    results_1 <- data.frame(exp=paste(j), outc=protein, nsnp=output_mr_ivw$nsnp, 
                                            method=c("IVW", "Egger"), b=output_mr_ivw$b, 
                                            se=output_mr_ivw$se, pval=output_mr_ivw$pval)
                    # egger
                    output_mr_egger <- mr_egger_regression(b_exp=dat$beta.exposure, b_out=dat$beta.outcome, se_exp=dat$se.exposure, se_out=dat$se.outcome)
                    results_2 <- data.frame(exp=paste(j), outc=protein, nsnp=output_mr_ivw$nsnp[1], 
                                            method=c("Egger_int"), b=output_mr_egger$b_i, 
                                            se=output_mr_egger$se_i, pval=output_mr_egger$pval_i)
                    # mr raps
                    output_mr_raps <- mr.raps(data.frame(beta.exposure=dat$beta.exposure, beta.outcome=dat$beta.outcome, 
                                                         se.exposure=dat$se.exposure, se.outcome=dat$se.outcome))
                    results_3 <- data.frame(exp=paste(j), outc=protein, nsnp=output_mr_ivw$nsnp[1], 
                                            method="MR-RAPS", b=output_mr_raps$beta.hat, 
                                            se=output_mr_raps$beta.se, pval=2 * pnorm(-abs(output_mr_raps$beta.hat / output_mr_raps$beta.se)))
                    # appending results
                    results <- rbind(results_1, results_2, results_3)
                }
            }
            
            if (is.null(results) || nrow(results) == 0) {
                print(paste0("Skipping ", protein, ", no results"))
            } else {
                # appending results to bigger dataframe
                df_sum <- rbind(df_sum, results)
                df_instr <- rbind(df_instr, dat)
            }
        }
    }
}


write.csv(df_sum, paste0("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/chip_to_prot_raw/", protein,"_3chip_raw.csv"), row.names = F, quote = F)
write.csv(df_instr, paste0("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/chip_to_prot_raw/", protein,"_3chip_instr_raw.csv"), row.names = F, quote = F)
print(paste(protein, "done"))
unlink(unique_unpack_dir, recursive = TRUE)


