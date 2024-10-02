
  library(synapser) 
  library(R.utils)
  library(TwoSampleMR)
  library(tidyr)
  library(dplyr)
  library(data.table)
  library(ieugwasr)
  library(genetics.binaRies)
  library(MendelianRandomization)
  
  synLogin(email='XXX', authToken="XXX")
  path = "/medpop/esp2/aschuerm/misc/others/zhi/"
  
  protein_list = fread(paste0(path, "updated_olink.txt"), header = FALSE)
  protein_list = protein_list$V1
  
  chip_sumstats <- fread("/medpop/esp2/mesbah/projects/Meta_GWAS/n650k/eur_only_summary/lifted_hg37.eur_metaGWAS.CHIP.GWAMA.hg37_dbSNP.eaf001_min2Studies.tsv.gz")
    chip_sumstats$POS_hg38 <- as.numeric(gsub("chr\\d+:(\\d+):.*", "\\1", chip_sumstats$varID_hg38))
    chip = chip_sumstats[,c("CHR", "POS_hg38", "OtherAllele", "EffectAllele", "MarkerID", "P", "BETA", "SE", "EAF")]
    colnames(chip) = c("#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt")
    chip$marker <- paste(chip$"#chrom", chip$pos, chip$ref, chip$alt, sep=":")
  tet2_sumstats <- fread("/medpop/esp2/mesbah/projects/Meta_GWAS/n650k/eur_only_summary/lifted_hg37.eur_metaGWAS.TET2.GWAMA.hg37_dbSNP.eaf001_min2Studies.tsv.gz")
    tet2_sumstats$POS_hg38 <- as.numeric(gsub("chr\\d+:(\\d+):.*", "\\1", tet2_sumstats$varID_hg38))
    tet2 = tet2_sumstats[,c("CHR", "POS_hg38", "OtherAllele", "EffectAllele", "MarkerID", "P", "BETA", "SE", "EAF")]
    colnames(tet2) = c("#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt")
    tet2$marker <- paste(tet2$"#chrom", tet2$pos, tet2$ref, tet2$alt, sep=":")
  dnmt3a_sumstats <- fread("/medpop/esp2/mesbah/projects/Meta_GWAS/n650k/eur_only_summary/lifted_hg37.eur_metaGWAS.DNMT3A.GWAMA.hg37_dbSNP.eaf001_min2Studies.tsv.gz")
    dnmt3a_sumstats$POS_hg38 <- as.numeric(gsub("chr\\d+:(\\d+):.*", "\\1", dnmt3a_sumstats$varID_hg38))
    dnmt3a = dnmt3a_sumstats[,c("CHR", "POS_hg38", "OtherAllele", "EffectAllele", "MarkerID", "P", "BETA", "SE", "EAF")]
    colnames(dnmt3a) = c("#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt")
    dnmt3a$marker <- paste(dnmt3a$"#chrom", dnmt3a$pos, dnmt3a$ref, dnmt3a$alt, sep=":")
  
  sumstats_info = fread("/medpop/esp2/aschuerm/ukb_proteomics_cvd/input_files/olink_protein_map_3k_v1.tsv")    
  sumstats_info = sumstats_info[which(sumstats_info$Assay %in% protein_list),]
  
  df_sum <- data.frame(exp=NA, outc=NA, nsnp=NA, method=NA, b=NA, se=NA, pval=NA)[-1,]
  df_instr <- data.frame(pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                     other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                     palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                     outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                     pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                     action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA, marker_ld=NA)[-1,]

for (i in sumstats_info$Code){

  syn_code <- synGet(entity=i, downloadLocation = paste0(path, "sumstats/")) 
    untar(paste(syn_code$path), list=F, exdir=paste(syn_code$cacheDir))
    chrom_u <- fread(paste0(syn_code$cacheDir, "/", gsub(".tar", "", sumstats_info[sumstats_info$Code==i,]$Docname[1]), "/", 
                  "discovery_chr", sumstats_info[sumstats_info$Code==i,]$chr[1], "_", sumstats_info[sumstats_info$Code==i,]$UKBPPP_ProteinID[1], 
                  ":", sumstats_info[sumstats_info$Code==i,]$Panel[1], ".gz"))
  
  chrom_u <- chrom_u[chrom_u$GENPOS > (sumstats_info[sumstats_info$Code==i,]$gene_start[1] - 200000) &                    
                  chrom_u$GENPOS < (sumstats_info[sumstats_info$Code==i,]$gene_end[1] + 200000), ]
    chrom_u$P <- 10^-chrom_u$LOG10P    
    chrom_u <- chrom_u[chrom_u$P<5e-6,]                                                                                        
    chrom_u2 <- chrom_u
    chrom_u2$BETA <- -chrom_u2$BETA
    names(chrom_u2)[names(chrom_u2) %in% c('ALLELE0', 'ALLELE1')] <- c('ALLELE1', 'ALLELE0')
    chrom_u$marker <- paste(chrom_u$CHROM, chrom_u$GENPOS, chrom_u$ALLELE0, chrom_u$ALLELE1, sep=":")
    chrom_u2$marker <- paste(chrom_u2$CHROM, chrom_u2$GENPOS, chrom_u2$ALLELE0, chrom_u2$ALLELE1, sep=":")
    chrom_u <- rbind(chrom_u, chrom_u2)
    rm(chrom_u2)
    
  if (is.null(chrom_u) || nrow(chrom_u) == 0) {
      print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1], "(1)"))
  } else {
  
  for (j in c("chip", "dnmt3a", "tet2")) {
  
   outcome <- get(j)
   outcome_overlap <- outcome[outcome$marker %in% chrom_u$marker,]                    

    if (is.null(outcome_overlap) || nrow(outcome_overlap) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1], "(2)"))
    } else {
      
     outcome_overlap$phen <- paste(j)
     rsid <- outcome_overlap[,c("marker", "rsids")]
     rsid$marker <- tolower(rsid$marker)
     outcome_overlap <- as.data.frame(outcome_overlap)
     outcome_overlap <- format_data(outcome_overlap, type="outcome", phenotype_col="phen", snp_col="marker", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     
     chrom <- chrom_u
     chrom$marker <- tolower(chrom$marker)
     chrom_overlap <- chrom[chrom$marker %in% outcome_overlap$SNP,]                                
     chrom_overlap$phen <- sumstats_info[sumstats_info$Code==i,]$Assay[1]
     chrom_overlap <- as.data.frame(chrom_overlap)
     chrom_overlap <- format_data(chrom_overlap, type="exposure", phenotype_col="phen", snp_col="marker", beta_col="BETA", se_col="SE", eaf_col="A1FREQ",
                                effect_allele_col="ALLELE1", other_allele_col="ALLELE0", pval_col="LOG10P", chr_col="CHROM", samplesize_col="N", pos_col="GENPOS", log_pval=T)
     dat <- harmonise_data(exposure_dat=chrom_overlap, outcome_dat=outcome_overlap)                                    
     
     dat <- merge(dat, rsid, all.x=T, all.y=F, by.x="SNP", by.y="marker")
     dat <- dat[order(dat$pval.exposure),]                                                            
     names(dat)[names(dat) %in% c('SNP', 'rsids')] <- c('marker', 'SNP')
     dat <- dat[!duplicated(dat$SNP),]
     dat <- dat[dat$mr_keep,]
     clump <- ld_clump(dplyr::tibble(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$id.exposure),   
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = 0.1,
                      bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
     dat <- dat[dat$SNP %in% clump$rsid,]
     rm(chrom_overlap, outcome_overlap, clump)
     
   
   if (nrow(dat[dat$mr_keep,])==0) {
      print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1], "(3)"))
     results <- NULL
    } else {
  
   if (nrow(dat)==1) {                                                                                                    
     results_mr <- mr(dat, method_list=c("mr_wald_ratio"))
     results <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste(j), nsnp=results_mr$nsnp, method=results_mr$method, b=results_mr$b, 
                         se=results_mr$se, pval=results_mr$pval)
    } else if (nrow(dat)==2) {                                                                                             
      ld <- ld_matrix(dat$SNP, bfile="/medpop/esp2/aschuerm/tools/g1000_eur", plink_bin=genetics.binaRies::get_plink_binary())
       rownames(ld) <- gsub("\\_.*", "", rownames(ld))
       dat_dupl <- dat[!(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_") %in% colnames(ld)),]
        colnames(dat_dupl)[which(colnames(dat_dupl) %in% c("effect_allele.exposure", "other_allele.exposure", "effect_allele.outcome", "other_allele.outcome"))] <- c("other_allele.exposure", "effect_allele.exposure", "other_allele.outcome", "effect_allele.outcome")
        dat_dupl$beta.exposure <- dat_dupl$beta.exposure*-1
        dat_dupl$beta.outcome <- dat_dupl$beta.outcome*-1
        dat_dupl$eaf.exposure <- 1-dat_dupl$eaf.exposure
        dat_dupl$eaf.outcome <- 1-dat_dupl$eaf.outcome
        dat <- dat[(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_") %in% colnames(ld)),]
        dat <- rbind(dat, dat_dupl)
        dat <- dat[ order(match(dat$SNP, rownames(ld))), ]
      dat2 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                               by = dat$beta.outcome, byse = dat$se.outcome,
                                               correlation = ld)
     output_mr_ivw_corr <- MendelianRandomization::mr_ivw(dat2, correl = TRUE)
     results <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste(j), nsnp=output_mr_ivw_corr@SNPs, 
                           method="Inverse variance weighted (correlation inc)", b=output_mr_ivw_corr@Estimate, 
                           se=output_mr_ivw_corr@StdError, pval=output_mr_ivw_corr@Pvalue)
    } else {                                                                                                                  
      ld <- ld_matrix(dat$SNP, bfile="/medpop/esp2/aschuerm/tools/g1000_eur", plink_bin=genetics.binaRies::get_plink_binary())
       rownames(ld) <- gsub("\\_.*", "", rownames(ld))
       dat_dupl <- dat[!(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_") %in% colnames(ld)),]
        colnames(dat_dupl)[which(colnames(dat_dupl) %in% c("effect_allele.exposure", "other_allele.exposure", "effect_allele.outcome", "other_allele.outcome"))] <- c("other_allele.exposure", "effect_allele.exposure", "other_allele.outcome", "effect_allele.outcome")
        dat_dupl$beta.exposure <- dat_dupl$beta.exposure*-1
        dat_dupl$beta.outcome <- dat_dupl$beta.outcome*-1
        dat_dupl$eaf.exposure <- 1-dat_dupl$eaf.exposure
        dat_dupl$eaf.outcome <- 1-dat_dupl$eaf.outcome
        dat <- dat[(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_") %in% colnames(ld)),]
        dat <- rbind(dat, dat_dupl)
        dat <- dat[ order(match(dat$SNP, rownames(ld))), ]
      dat2 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                               by = dat$beta.outcome, byse = dat$se.outcome,
                                               correlation = ld)
     output_mr_ivw_corr <- MendelianRandomization::mr_ivw(dat2, correl = TRUE)
     output_mr_egger_corr <- MendelianRandomization::mr_egger(dat2, correl = TRUE)
     results1 <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste(j), nsnp=output_mr_ivw_corr@SNPs, 
                           method="Inverse variance weighted (correlation inc)", b=output_mr_ivw_corr@Estimate, 
                           se=output_mr_ivw_corr@StdError, pval=output_mr_ivw_corr@Pvalue)
     results2 <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste(j), nsnp=output_mr_egger_corr@SNPs, 
                           method="Egger (correlation inc)", b=output_mr_egger_corr@Estimate, 
                           se=output_mr_egger_corr@StdError.Est, pval=output_mr_egger_corr@Pvalue.Est)
     results3 <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste(j), nsnp=output_mr_egger_corr@SNPs, 
                           method="Egger intercept (correlation inc)", b=output_mr_egger_corr@Intercept, 
                           se=output_mr_egger_corr@StdError.Int, pval=output_mr_egger_corr@Pvalue.Int)
     results <- rbind(results1, results2, results3)
     rm(results1, results2, results3)
     
    }
    }
   
  if (is.null(results) || nrow(results) == 0) {
      print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1], "(4)"))
  } else {
    df_sum <- rbind(df_sum, results)
    df_instr <- rbind(df_instr, dat)
    rm(dat, results, ld, dat2, output_mr_ivw_corr)
  
  }
  }
  }
  }
  
  write.csv(df_sum, "/medpop/esp2/aschuerm/misc/others/zhi/mr_prot_to_chip_rev.csv") 
  write.csv(df_instr, "/medpop/esp2/aschuerm/misc/others/zhi/mr_prot_to_chip_instruments_rev.csv") 
  print(paste(sumstats_info[sumstats_info$Code==i,]$Assay[1], "done"))
  unlink(paste0(gsub("\\\\[^\\\\]*$", "", syn_code$cacheDir)), recursive=T)
}
  
