######### 00 - prep ########
library(reticulate)
use_condaenv("synapser_env", required = TRUE)
library(synapser) 
library(R.utils)
library(optparse)
library(TwoSampleMR)
library(stringr)
library(ieugwasr)
library(genetics.binaRies)
library(mr.raps)
library(tidyr)
library(dplyr)
library(data.table)

gc()
rm(list = ls())

# inputting file for shell script
option_list = list(
    make_option(c("-s", "--sumstats"), type="character", default=NULL, 
                help="the summary statistics of CHIP/TET2/DNMT3A", metavar="character"),
    make_option(c("-p", "--protein"), type="character", default=NULL,
                help="the protein", metavar="character")
    
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# # the option for running for gene using bash script in g_ shell script
sumstats <- opt$sumstats
protein <- opt$protein


## testing section ##
# sumstats <- "CHIP"
# protein <- "OTUD7B"

# loggin information for synapse, needs email and PAT
synLogin(email=, 
         authToken=)



###### 02 - reading in files #######
# the full CHIP/TET2/DNMT3A summary statistics are required
# pQTL as instruments
ss_df <- data.frame(fread(paste0("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/clumped_chip_ss/", sumstats, "_2022ss_short.tsv.gz")))
print(paste0("Successfully reading in summary statistics of ", sumstats))
ss_df <- ss_df[!duplicated(ss_df$marker),]


# preparing files
sumstats_info = fread("/medpop/esp2/aschuerm/ukb_proteomics_cvd/input_files/olink_protein_map_3k_v1.tsv") %>%
    filter(Assay %in% protein)
## !!!!!!! this is only because I checked the 2 proteins with multiple rows were the ones having the same docname !!!!!!!! ###
sumstats_info <- sumstats_info[!duplicated(sumstats_info$Docname),]
print(paste0("Successfully reading in info needed for downloading of ", protein))



###### 03 - preparing to loop #####
# two empty dataframe for storage
df_sum <- data.frame(exp=NA, outc=NA, nsnp=NA, method=NA, b=NA, se=NA, pval=NA)[-1,]
df_instr <- data.frame(pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                       other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                       palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                       outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                       pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                       action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA, marker_ld=NA)[-1,]



### writing the loop
# temporary cache 
# Downloading the summary statistics for the protein of interest
path = "/medpop/esp2/aschuerm/misc/others/zhi/"
# create unique unpack directory for files, 
# as some proteins shares the same SS and can conflict when running parallele jobs
unique_unpack_dir <- paste0(path, "sumstats_", protein, "_", Sys.getpid())
dir.create(unique_unpack_dir, showWarnings = FALSE, recursive = TRUE)
# download the file into unique directory, and unpack in that origin
syn_code <- synGet(entity = sumstats_info$Code, downloadLocation = unique_unpack_dir)
untar(syn_code$path, exdir = unique_unpack_dir)

# after downloading 
chr <- sumstats_info$chr
chrom_u <- data.frame(fread(paste0(syn_code$cacheDir, "/", gsub(".tar", "", sumstats_info$Docname), "/", 
                                   "discovery_chr", chr, "_", sumstats_info$UKBPPP_ProteinID, 
                                   ":", sumstats_info$Panel, ".gz")))

# selecting regions # Selecting the cis-region only (here defined as 1Mb before or after the protein-encoding region)
# window remains debatable, and filter level is also subject to change
chrom_u <- chrom_u %>%
    filter(GENPOS > sumstats_info$gene_start - 5e5 & GENPOS < sumstats_info$gene_end + 5e5) %>%
    mutate(P = 10^-LOG10P,
           marker = paste(CHROM, GENPOS, ALLELE0, ALLELE1, sep=":"),
           marker1 = paste(CHROM, GENPOS, ALLELE1, ALLELE0, sep=":")) %>%
    filter(P < 5e-6)



# if there's no significant resutls from the pQTL data
if (is.null(chrom_u) || nrow(chrom_u) == 0) {
    print(paste0("Skipping ", sumstats_info$Assay, ", because pQTL data is empty after filtering for significance"))
} else {
    # at this point there's still no rsID for the downloaded protein data
    # reducing dimention of outcome ss(pQTL) and adding rsID as SNP column
    rsid <- ss_df[,c("marker", "SNP")] # dictionary for SNP rsID
    outcome_overlap <- ss_df[ss_df$marker %in% chrom_u$marker | ss_df$marker %in% chrom_u$marker1,]
    nrow(outcome_overlap)
    
    # if there's no overlap between exposure and the outcome data, skip MR analysis completely
    if (is.null(outcome_overlap) || nrow(outcome_overlap) == 0) {
        print(paste0("Skipping ", sumstats_info$Assay, ", because no overlap between exposure pQTL and outcome CHIP SS data"))
    } else {
        # mutating the data for SNP rsID
        chrom_overlap <- chrom_u %>%
            left_join(rsid, by = "marker") %>%                
            mutate(SNP_alt = SNP) %>%                                   
            select(-SNP) %>%
            left_join(rsid, by = c("marker1" = "marker")) %>%               
            mutate(
                SNP = coalesce(SNP, SNP_alt)                              
            ) %>%
            select(-SNP_alt) %>%
            as.data.frame()
        
        # as the data has rsID, clump the data
        clump <- ld_clump(dplyr::tibble(rsid=chrom_overlap$SNP, pval=chrom_overlap$P),
                          plink_bin = "/medpop/esp2/lli/Mike-METAP1/TwoSampleMR/plink",
                          clump_kb=10000, clump_r2 = 0.1,
                          bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
        
        # filter original file based on the clump
        chrom_clumped <- chrom_overlap[chrom_overlap$SNP %in% clump$rsid,] %>%
            mutate(phen = protein)
        # formatting the overlapped data
        chrom_clumped <- format_data(chrom_clumped, type="exposure", phenotype_col="phen", snp_col="SNP", 
                                     beta_col="BETA", se_col="SE", eaf_col="A1FREQ",
                                     effect_allele_col="ALLELE1", other_allele_col="ALLELE0", 
                                     pval_col="P", chr_col="CHROM", samplesize_col="N", pos_col="GENPOS")
        
        # outcome data
        outcome_overlap <- ss_df[ss_df$SNP %in% chrom_clumped$SNP,]
        nrow(outcome_overlap)
        outcome_overlap$phen <- paste(sumstats)
        # formatting the data as outcome (chip ss)
        outcome_overlap <- as.data.frame(outcome_overlap)
        outcome_overlap <- format_data(outcome_overlap, type="outcome", phenotype_col="phen", 
                                       snp_col="SNP", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                       effect_allele_col="alt", other_allele_col="ref", pval_col="pval", 
                                       chr_col="chr", pos_col="bp")
        # this keeps the first instance of duplicated SNPs/multiallelic sites
        
        # harmonize data (this handles the strand orientation and remove palindromic ones)
        dat <- harmonise_data(exposure_dat=chrom_clumped, outcome_dat=outcome_overlap)
        dat <- dat[order(dat$pval.exposure),]                                                                                   
        dat <- dat[!duplicated(dat$SNP),]
        dat <- dat[dat$mr_keep,]
        
        # if after clumping the summary statistics is empty, skip the MR analysis
        if (nrow(dat)==0) {
            print(paste0("Skipping ", sumstats_info$Assay, "harmonized data has 0 rows"))
            results <- NULL
        } else {
            # if non zero, there are 3 different situations: 1 snp, 2 snps, 3+ snps
            if (nrow(dat)==1) {                                                                                                        
                results_mr <- mr(dat, method_list=c("mr_wald_ratio"))
                results <- data.frame(exp=protein, outc=sumstats, 
                                      nsnp=results_mr$nsnp, method=results_mr$method, b=results_mr$b,
                                      se=results_mr$se, pval=results_mr$pval)
                
            } else if (nrow(dat)==2) {
                # if you have 2 SNPs, IVW, EGGER and RAPS
                output_mr_ivw <- mr(dat, method_list=c("mr_ivw"))
                results_1 <- data.frame(exp=protein, outc=sumstats, nsnp=output_mr_ivw$nsnp, 
                                        method="IVW", b=output_mr_ivw$b, 
                                        se=output_mr_ivw$se, pval=output_mr_ivw$pval)
                output_mr_raps <- mr.raps(data.frame(beta.exposure=dat$beta.exposure, beta.outcome=dat$beta.outcome, 
                                                     se.exposure=dat$se.exposure, se.outcome=dat$se.outcome))
                results_2 <- data.frame(exp=protein, outc=sumstats, nsnp=output_mr_ivw$nsnp, 
                                        method="MR-RAPS", b=output_mr_raps$beta.hat, 
                                        se=output_mr_raps$beta.se, pval=2 * pnorm(-abs(output_mr_raps$beta.hat / output_mr_raps$beta.se)))
                # binding the results of MR-RAPS and MR IVW
                results <- rbind(results_1, results_2)
                
            } else {
                # if you have 3+SNPs, IVW, EGGER and RAPS
                output_mr_ivw <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression"))
                results_1 <- data.frame(exp=protein, outc=sumstats, nsnp=output_mr_ivw$nsnp, 
                                        method=c("Inverse variance weighted", "Egger"), b=output_mr_ivw$b, 
                                        se=output_mr_ivw$se, pval=output_mr_ivw$pval)
                output_mr_egger <- mr_egger_regression(b_exp=dat$beta.exposure, b_out=dat$beta.outcome, se_exp=dat$se.exposure, se_out=dat$se.outcome)
                results_2 <- data.frame(exp=protein, outc=sumstats, nsnp=output_mr_ivw$nsnp[1], 
                                        method=c("Egger (intercept)"), b=output_mr_egger$b_i, 
                                        se=output_mr_egger$se_i, pval=output_mr_egger$pval_i)
                output_mr_raps <- mr.raps(data.frame(beta.exposure=dat$beta.exposure, beta.outcome=dat$beta.outcome, 
                                                     se.exposure=dat$se.exposure, se.outcome=dat$se.outcome))
                results_3 <- data.frame(exp=protein, outc=sumstats, nsnp=output_mr_ivw$nsnp[1], 
                                        method="MR-RAPS", b=output_mr_raps$beta.hat, 
                                        se=output_mr_raps$beta.se, pval=2 * pnorm(-abs(output_mr_raps$beta.hat / output_mr_raps$beta.se)))
                results <- rbind(results_1, results_2, results_3)
            }
        }
        if (is.null(results) || nrow(results) == 0) {
            print(paste0("Skipping ", protein, "because results are empty"))
        } else {
            df_sum <- rbind(df_sum, results)
            df_instr <- rbind(df_instr, dat)
        }
    }
}


write.csv(df_sum, paste0("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/prot_to_chip_raw/", protein, "_", sumstats,"_nocorr.csv"), row.names = F, quote = F)
write.csv(df_instr, paste0("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/prot_to_chip_raw/", protein, "_", sumstats,"_nocorr_instruments.csv"), row.names = F, quote = F) 
print(paste0(protein," & ", sumstats, " done"))
unlink(unique_unpack_dir, recursive = TRUE)






















