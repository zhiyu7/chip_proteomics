
  library(synapser) 
  library(R.utils)
  library(TwoSampleMR)
  library(tidyr)
  library(dplyr)
  library(data.table)
  library(ieugwasr)
  library(genetics.binaRies)
  library(MendelianRandomization)
  library(mr.raps)
  
  synLogin('XXX','XXX') 
  path = "/medpop/esp2/aschuerm/misc/others/zhi/"
  
  protein_list = fread(paste0(path, "updated_olink.txt"), header = FALSE)
  protein_list = protein_list$V1
  
  chip_sumstats <- fread("/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary/lifted_hg37.eur_metaGWAS.CHIP.GWAMA.hg37_dbSNP.eaf001_min2Studies.tsv.gz")
    chip_sumstats$POS_hg38 <- as.numeric(gsub("chr\\d+:(\\d+):.*", "\\1", chip_sumstats$varID_hg38))
    chip = chip_sumstats[,c("CHR", "POS_hg38", "OtherAllele", "EffectAllele", "MarkerID", "P", "BETA", "SE", "EAF")]
    colnames(chip) = c("#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt")
    chip$marker <- paste(chip$"#chrom", chip$pos, chip$ref, chip$alt, sep=":")
    chip <- chip[chip$pval < 5e-6,]
    clump <- ld_clump(dplyr::tibble(rsid=chip$rsids, pval=chip$pval),                          
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = 0.1,
                      bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
    chip <- chip[chip$rsids %in% clump$rsid,]
  tet2_sumstats <- fread("/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary/lifted_hg37.eur_metaGWAS.TET2.GWAMA.hg37_dbSNP.eaf001_min2Studies.tsv.gz")
    tet2_sumstats$POS_hg38 <- as.numeric(gsub("chr\\d+:(\\d+):.*", "\\1", tet2_sumstats$varID_hg38))
    tet2 = tet2_sumstats[,c("CHR", "POS_hg38", "OtherAllele", "EffectAllele", "MarkerID", "P", "BETA", "SE", "EAF")]
    colnames(tet2) = c("#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt")
    tet2$marker <- paste(tet2$"#chrom", tet2$pos, tet2$ref, tet2$alt, sep=":")
    tet2 <- tet2[tet2$pval < 5e-6,]
    clump <- ld_clump(dplyr::tibble(rsid=tet2$rsids, pval=tet2$pval),              
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = 0.1,
                      bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
    tet2 <- tet2[tet2$rsids %in% clump$rsid,]
  dnmt3a_sumstats <- fread("/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary/lifted_hg37.eur_metaGWAS.DNMT3A.GWAMA.hg37_dbSNP.eaf001_min2Studies.tsv.gz")
    dnmt3a_sumstats$POS_hg38 <- as.numeric(gsub("chr\\d+:(\\d+):.*", "\\1", dnmt3a_sumstats$varID_hg38))
    dnmt3a = dnmt3a_sumstats[,c("CHR", "POS_hg38", "OtherAllele", "EffectAllele", "MarkerID", "P", "BETA", "SE", "EAF")]
    colnames(dnmt3a) = c("#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt")
    dnmt3a$marker <- paste(dnmt3a$"#chrom", dnmt3a$pos, dnmt3a$ref, dnmt3a$alt, sep=":")
    dnmt3a <- dnmt3a[dnmt3a$pval < 5e-6,]
    clump <- ld_clump(dplyr::tibble(rsid=dnmt3a$rsids, pval=dnmt3a$pval),           
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = 0.001,
                      bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
    dnmt3a <- dnmt3a[dnmt3a$rsids %in% clump$rsid,]
 
  sumstats_info = fread("/medpop/esp2/aschuerm/ukb_proteomics_cvd/input_files/olink_protein_map_3k_v1.tsv")    
  sumstats_info = sumstats_info[which(sumstats_info$Assay %in% protein_list),]
  
  df_sum <- data.frame(exp=NA, outc=NA, nsnp=NA, method=NA, b=NA, se=NA, pval=NA)[-1,]
  df_instr <- data.frame(pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                     other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                     palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                     outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                     pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                     action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA)[-1,]

for (i in sumstats_info$Code){

  syn_code <- synGet(entity=i, downloadLocation = paste0(path, "sumstats/")) 
    untar(paste(syn_code$path), list=F, exdir=paste(syn_code$cacheDir))
    
    
   bm_whole <- data.frame(matrix(ncol = 12, nrow = 0))
    for (nr in 1:22){
      chrom_part <- fread(paste0(syn_code$cacheDir, "/", gsub(".tar", "", sumstats_info[sumstats_info$Code==i,]$Docname[1]), "/", 
                    "discovery_chr", nr, "_", sumstats_info[sumstats_info$Code==i,]$UKBPPP_ProteinID[1], 
                    ":", sumstats_info[sumstats_info$Code==i,]$Panel[1], ".gz"))
      bm_whole <- rbind(bm_whole, chrom_part)
      print(paste("chromosome", nr, "bound"))
    }
    bm_whole$P <- 10^-bm_whole$LOG10P    
    bm_whole$marker <- paste(bm_whole$CHROM, bm_whole$GENPOS, bm_whole$ALLELE0, bm_whole$ALLELE1, sep=":")

     
  if (is.null(bm_whole) || nrow(bm_whole) == 0) {
      print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1], "(1)"))
  } else {
  
  for (j in c("chip", "dnmt3a", "tet2")) {
  
   outcome <- get(j)
   outcome_overlap <- outcome[outcome$marker %in% bm_whole$marker,]                     

    if (is.null(outcome_overlap) || nrow(outcome_overlap) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1], "(2)"))
    } else {
      
     outcome_overlap$phen <- paste(j)
     rsid <- outcome_overlap[,c("marker", "rsids")]
     rsid$marker <- tolower(rsid$marker)
     outcome_overlap <- format_data(outcome_overlap, type="exposure", phenotype_col="phen", snp_col="marker", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     
     chrom <- bm_whole
     chrom$marker <- tolower(chrom$marker)
     chrom_overlap <- chrom[chrom$marker %in% outcome_overlap$SNP,]                                                
     chrom_overlap$phen <- sumstats_info[sumstats_info$Code==i,]$Assay[1]
     chrom_overlap <- format_data(chrom_overlap, type="outcome", phenotype_col="phen", snp_col="marker", beta_col="BETA", se_col="SE", eaf_col="A1FREQ",
                                effect_allele_col="ALLELE1", other_allele_col="ALLELE0", pval_col="LOG10P", chr_col="CHROM", samplesize_col="N", pos_col="GENPOS", log_pval=T)
     dat <- harmonise_data(exposure_dat=outcome_overlap, outcome_dat=chrom_overlap)                                        
     
     dat <- merge(dat, rsid, all.x=T, all.y=F, by.x="SNP", by.y="marker")
     dat <- dat[order(dat$pval.exposure),]                                                                                 
     names(dat)[names(dat) %in% c('SNP', 'rsids')] <- c('marker', 'SNP')
     dat <- dat[!duplicated(dat$SNP),]
     dat <- dat[dat$mr_keep,]
     rm(chrom_overlap, outcome_overlap)
     
 

   if (nrow(dat[dat$mr_keep,])==0) {
      print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1], "(3)"))
     results <- NULL
    } else {
  
   if (nrow(dat)==1) {                                                                                           
     results_mr <- mr(dat, method_list=c("mr_wald_ratio"))
     results <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste("outcome"), nsnp=results_mr$nsnp, method=results_mr$method, b=results_mr$b, 
                         se=results_mr$se, pval=results_mr$pval)
    } else if (nrow(dat)==2) {                                                                                            

     output_mr_ivw <- mr(dat, method_list=c("mr_ivw"))
     results_1 <- data.frame(exp=paste(j), outc=sumstats_info[sumstats_info$Code==i,]$Assay[1], nsnp=output_mr_ivw$nsnp, 
                           method="Inverse variance weighted", b=output_mr_ivw$b, 
                           se=output_mr_ivw$se, pval=output_mr_ivw$pval)
     output_mr_raps <- mr.raps(b_exp=dat$beta.exposure, b_out=dat$beta.outcome, se_exp=dat$se.exposure, se_out=dat$se.outcome)
     results_2 <- data.frame(exp=paste(j), outc=sumstats_info[sumstats_info$Code==i,]$Assay[1], nsnp=output_mr_ivw$nsnp, 
                           method="MR-RAPS", b=output_mr_raps$beta.hat, 
                           se=output_mr_raps$beta.se, pval=output_mr_raps$beta.p.value)
     results <- rbind(results_1, results_2)
      
    } else {                                                                                                                  # If you have more than 2 variants, you can do anything (including IVW and MR-Egger)

     output_mr_ivw <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression"))
     results_1 <- data.frame(exp=paste(j), outc=sumstats_info[sumstats_info$Code==i,]$Assay[1], nsnp=output_mr_ivw$nsnp, 
                           method=c("Inverse variance weighted", "Egger"), b=output_mr_ivw$b, 
                           se=output_mr_ivw$se, pval=output_mr_ivw$pval)
     output_mr_egger <- mr_egger_regression(b_exp=dat$beta.exposure, b_out=dat$beta.outcome, se_exp=dat$se.exposure, se_out=dat$se.outcome)
     results_2 <- data.frame(exp=paste(j), outc=sumstats_info[sumstats_info$Code==i,]$Assay[1], nsnp=output_mr_ivw$nsnp[1], 
                           method=c("Egger (intercept)"), b=output_mr_egger$b_i, 
                           se=output_mr_egger$se_i, pval=output_mr_egger$pval_i)
     output_mr_raps <- mr.raps(b_exp=dat$beta.exposure, b_out=dat$beta.outcome, se_exp=dat$se.exposure, se_out=dat$se.outcome)
     results_3 <- data.frame(exp=paste(j), outc=sumstats_info[sumstats_info$Code==i,]$Assay[1], nsnp=output_mr_ivw$nsnp[1], 
                           method="MR-RAPS", b=output_mr_raps$beta.hat, 
                           se=output_mr_raps$beta.se, pval=output_mr_raps$beta.p.value)
     results <- rbind(results_1, results_2, results_3)
      
      
    }
    }
   
  if (is.null(results) || nrow(results) == 0) {
      print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1], "(4)"))
  } else {
    df_sum <- rbind(df_sum, results)
    df_instr <- rbind(df_instr, dat)
    rm(dat, results, ld, dat2, output_mr_ivw_corr, output_mr_egger_corr)
  
  }
  }
  }
  }
  
  write.csv(df_sum, "/medpop/esp2/aschuerm/misc/others/zhi/mr_chip_to_prot.csv") 
  write.csv(df_instr, "/medpop/esp2/aschuerm/misc/others/zhi/mr_chip_to_prot_instruments.csv") 
  print(paste(sumstats_info[sumstats_info$Code==i,]$Assay[1], "done"))
  unlink(paste0(gsub("\\\\[^\\\\]*$", "", syn_code$cacheDir)), recursive=T)
}
  
