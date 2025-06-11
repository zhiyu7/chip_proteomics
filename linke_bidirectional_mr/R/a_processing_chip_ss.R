######## prep ########
library(R.utils)
library(data.table)
library(stringr)
library(tidyr)
library(dplyr)
library(ieugwasr)

gc()
rm(list = ls())


path = "/medpop/esp2/aschuerm/misc/others/zhi/"

setwd("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr")

##### dictionary for rsID in build 38
# dictionary <- data.frame(fread("/medpop/esp2/btruong/Tools/hg38_common_1_chrpos.txt")) %>%
#     dplyr::select(V3, V5) %>%
#     dplyr::rename(rsids = V3,
#                   chr_bp = V5)

###### old GWAS summary statistics
# /medpop/esp2/mesbah/projects/Meta_GWAS/n650k/eur_only_summary/lifted_hg37.eur_metaGWAS.DNMT3A.GWAMA.hg37_dbSNP.eaf001_min2Studies.tsv.gz

######### CHIP ##########
chip <- data.frame(fread("/medpop/esp2/mesbah/projects/Meta_GWAS/n650k/eur_only_summary/lifted_hg37.eur_metaGWAS.CHIP.GWAMA.hg37_dbSNP.eaf001_min2Studies.tsv.gz")) %>%
    rename(SNP = MarkerID,
           beta = BETA,
           sebeta = SE,
           pval = P,
           ref = REF,
           alt = ATL,
           af_alt = EAF,
           marker = varID_hg38) %>%
    select(SNP, beta, sebeta, pval, ref, alt, af_alt, marker) %>%
    # extract the chr string from the marker, make a chr_bp, then create 2 separate column for chr and bp
    mutate(chr_bp = str_extract(marker, "(?<=chr)\\d+:\\d+"),
           chr = as.integer(str_remove(str_split_fixed(marker, ":", 4)[, 1], "chr")),
           bp  = as.integer(str_split_fixed(marker, ":", 4)[, 2]),
           marker = str_remove(marker, "^chr")
           ) 

# now filter
chip_filtered <- chip %>%
    filter(pval < 5e-6)

clump <- ld_clump(dplyr::tibble(rsid=chip_filtered$SNP, pval=chip_filtered$pval), 
                  plink_bin = "/medpop/esp2/lli/Mike-METAP1/TwoSampleMR/plink",
                  clump_kb=10000, clump_r2 = 0.1,
                  bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
chip_clumped <- chip[chip$SNP %in% clump$rsid,]

write.table(chip, gzfile("clumped_chip_ss/CHIP_2022ss_short.tsv.gz"), row.names = F, quote = F, sep = "\t")
write.table(chip_clumped, "clumped_chip_ss/CHIP_2022ss_clumped.tsv", row.names = F, quote = F, sep = "\t")
rm(chip, chip_clumped)


# chip <- data.frame(fread("/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.hasCH.EUR.ukbb450k_AoU250k_TOPMed72k_MGBB53k_BioVU54k.out.gz")) %>%
#     dplyr::select(rs_number, reference_allele, other_allele, p.value, beta, se, eaf) %>%
#     dplyr::rename(marker = rs_number,
#                   ref = reference_allele,
#                   alt = other_allele,
#                   pval = p.value,
#                   af_alt = eaf,
#                   beta = beta,
#                   sebeta=se) %>%
#     mutate(chr_bp = str_extract(marker, "(?<=chr)\\d+:\\d+")) %>%
#     left_join(dictionary, by = "chr_bp") %>% # appending and making the rsids column
#     dplyr::select(-chr_bp) %>%
#     # create the marker required for MR
    # mutate(chr = as.integer(str_remove(str_split_fixed(marker, ":", 4)[, 1], "chr")),
    #        bp  = as.integer(str_split_fixed(marker, ":", 4)[, 2]),
    #        marker = str_remove(marker, "^chr"))



# chip_filtered <- chip %>%
#     filter(pval < 5e-8)
# clump <- ld_clump(dplyr::tibble(rsid=chip_filtered$rsids, pval=chip_filtered$pval),
#                   plink_bin = "/medpop/esp2/lli/Mike-METAP1/TwoSampleMR/plink",
#                   clump_kb=10000, clump_r2 = 0.1,
#                   bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
# # final chip file
# chip_clumped <- chip[chip$rsids %in% clump$rsid,]
# 
# # write.table(chip, gzfile("clumped_chip_ss/CHIP_ss_short.tsv.gz"), row.names = F, quote = F, sep = "\t")
# write.table(chip_clumped, "clumped_chip_ss/CHIP_ss_clumped.tsv", row.names = F, quote = F, sep = "\t")
# 
# rm(chip, chip_filtered, chip_clumped)



############ TET2 ###########
tet2 <- data.frame(fread("/medpop/esp2/mesbah/projects/Meta_GWAS/n650k/eur_only_summary/lifted_hg37.eur_metaGWAS.TET2.GWAMA.hg37_dbSNP.eaf001_min2Studies.tsv.gz")) %>%
    rename(SNP = MarkerID,
           beta = BETA,
           sebeta = SE,
           pval = P,
           ref = REF,
           alt = ATL,
           af_alt = EAF,
           marker = varID_hg38) %>%
    select(SNP, beta, sebeta, pval, ref, alt, af_alt, marker) %>%
    # extract the chr string from the marker, make a chr_bp, then create 2 separate column for chr and bp
    mutate(chr_bp = str_extract(marker, "(?<=chr)\\d+:\\d+"),
           chr = as.integer(str_remove(str_split_fixed(marker, ":", 4)[, 1], "chr")),
           bp  = as.integer(str_split_fixed(marker, ":", 4)[, 2]),
           marker = str_remove(marker, "^chr")
    ) 


# now filter
tet2_filtered <- tet2 %>%
    filter(pval < 5e-6)

clump <- ld_clump(dplyr::tibble(rsid=tet2_filtered$SNP, pval=tet2_filtered$pval),                                 # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
                  plink_bin = "/medpop/esp2/lli/Mike-METAP1/TwoSampleMR/plink",
                  clump_kb=10000, clump_r2 = 0.1,
                  bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
tet2_clumped <- tet2[tet2$SNP %in% clump$rsid,]

write.table(tet2, gzfile("clumped_chip_ss/TET2_2022ss_short.tsv.gz"), row.names = F, quote = F, sep = "\t")
write.table(tet2_clumped, "clumped_chip_ss/TET2_2022ss_clumped.tsv", row.names = F, quote = F, sep = "\t")
rm(tet2, tet2_clumped)
    
# tet2 <- data.frame(fread("/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.hasTET2.EUR.ukbb450k_AoU250k_TOPMed72k_MGBB53k_BioVU54k.out.gz")) %>%
#     dplyr::select(rs_number, reference_allele, other_allele, p.value, beta, se, eaf) %>%
#     dplyr::rename(marker = rs_number,
#                   ref = reference_allele,
#                   alt = other_allele,
#                   pval = p.value,
#                   af_alt = eaf,
#                   beta = beta,
#                   sebeta=se) %>%
#     mutate(chr_bp = str_extract(marker, "(?<=chr)\\d+:\\d+")) %>%
#     left_join(dictionary, by = "chr_bp") %>%
#     dplyr::select(-chr_bp) %>%
#     # create the marker required for MR
#     mutate(chr = as.integer(str_remove(str_split_fixed(marker, ":", 4)[, 1], "chr")),
#            bp  = as.integer(str_split_fixed(marker, ":", 4)[, 2]),
#            marker = str_remove(marker, "^chr"))


# tet2_filtered <- tet2 %>%
#     filter(pval < 5e-8)
# clump <- ld_clump(dplyr::tibble(rsid=tet2_filtered$rsids, pval=tet2_filtered$pval),
#                   plink_bin = "/medpop/esp2/lli/Mike-METAP1/TwoSampleMR/plink",
#                   clump_kb=10000, clump_r2 = 0.1,
#                   bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
# tet2_clumped <- tet2[tet2$rsids %in% clump$rsid,]

# write.table(tet2, gzfile("clumped_chip_ss/TET2_ss_short.tsv.gz"), row.names = F, quote = F, sep = "\t")
# write.table(tet2_clumped, "clumped_chip_ss/TET2_ss_clumped.tsv", row.names = F, quote = F, sep = "\t")





########### DNMT3A #########
dnmt3a <- data.frame(fread("/medpop/esp2/mesbah/projects/Meta_GWAS/n650k/eur_only_summary/lifted_hg37.eur_metaGWAS.DNMT3A.GWAMA.hg37_dbSNP.eaf001_min2Studies.tsv.gz")) %>%
    rename(SNP = MarkerID,
           beta = BETA,
           sebeta = SE,
           pval = P,
           ref = REF,
           alt = ATL,
           af_alt = EAF,
           marker = varID_hg38) %>%
    select(SNP, beta, sebeta, pval, ref, alt, af_alt, marker) %>%
    # extract the chr string from the marker, make a chr_bp, then create 2 separate column for chr and bp
    mutate(chr_bp = str_extract(marker, "(?<=chr)\\d+:\\d+"),
           chr = as.integer(str_remove(str_split_fixed(marker, ":", 4)[, 1], "chr")),
           bp  = as.integer(str_split_fixed(marker, ":", 4)[, 2]),
           marker = str_remove(marker, "^chr")
    ) 

# now filter
dnmt3a_filtered <- dnmt3a %>%
    filter(pval < 5e-6)

clump <- ld_clump(dplyr::tibble(rsid=dnmt3a_filtered$SNP, pval=dnmt3a_filtered$pval),                                 # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
                  plink_bin = "/medpop/esp2/lli/Mike-METAP1/TwoSampleMR/plink",
                  clump_kb=10000, clump_r2 = 0.1,
                  bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
dnmt3a_clumped <- dnmt3a[dnmt3a$SNP %in% clump$rsid,]


write.table(dnmt3a, gzfile("clumped_chip_ss/DNMT3A_2022ss_short.tsv.gz"), row.names = F, quote = F, sep = "\t")
write.table(dnmt3a_clumped, "clumped_chip_ss/DNMT3A_2022ss_clumped.tsv", row.names = F, quote = F, sep = "\t")


# dnmt3a <- data.frame(fread("/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.hasDNMT3A.EUR.ukbb450k_AoU250k_TOPMed72k_MGBB53k_BioVU54k.out.gz")) %>%
#     dplyr::select(rs_number, reference_allele, other_allele, p.value, beta, se, eaf) %>%
#     dplyr::rename(marker = rs_number,
#                   ref = reference_allele,
#                   alt = other_allele,
#                   pval = p.value,
#                   af_alt = eaf,
#                   beta = beta,
#                   sebeta=se) %>%
#     mutate(chr_bp = str_extract(marker, "(?<=chr)\\d+:\\d+")) %>%
#     left_join(dictionary, by = "chr_bp") %>%
#     dplyr::select(-chr_bp) %>%
#     # create the marker required for MR
#     mutate(chr = as.integer(str_remove(str_split_fixed(marker, ":", 4)[, 1], "chr")),
#            bp  = as.integer(str_split_fixed(marker, ":", 4)[, 2]),
#            marker = str_remove(marker, "^chr"))
# 
# 
# dnmt3a_filtered <- dnmt3a %>%
#     filter(pval < 5e-8)
# clump <- ld_clump(dplyr::tibble(rsid=dnmt3a_filtered$rsids, pval=dnmt3a_filtered$pval),                                 # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
#                   plink_bin = "/medpop/esp2/lli/Mike-METAP1/TwoSampleMR/plink",
#                   clump_kb=10000, clump_r2 = 0.001,
#                   bfile = "/medpop/esp2/aschuerm/tools/g1000_eur")
# dnmt3a_clumped <- dnmt3a[dnmt3a$rsids %in% clump$rsid,]
# 
# # write.table(dnmt3a, gzfile("clumped_chip_ss/DNMT3A_ss_short.tsv.gz"), row.names = F, quote = F, sep = "\t")
# write.table(dnmt3a_clumped, "clumped_chip_ss/DNMT3A_ss_clumped.tsv", row.names = F, quote = F, sep = "\t")








