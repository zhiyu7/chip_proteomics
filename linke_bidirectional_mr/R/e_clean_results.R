##### 01 - prep #####
library(dplyr)
library(data.table)


gc()
rm(list = ls())


##### 02 - reading in files #####
all_files <- list.files("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/prot_to_chip_raw/")
protlist <- fread("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/protlist_clean", header = F)$V1
SS <- c("CHIP", "TET2", "DNMT3A")


##### CHIP TO PROT #####
results <- data.frame()
log <- data.frame()
# gather results
for (protein in protlist){
    file <- paste0("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/chip_to_prot_raw/",protein, "_3chip_raw.csv")
    if (file.exists(file)){
        df <- data.frame(fread(file, header = T))
        results <- rbind(results, df)
        log <- rbind(log, data.frame(protein = protein, status = "Y"))
    } else {
        print(paste0(protein, " file doesn't exist"))
        log <- rbind(log, data.frame(protein = protein, status = "N"))
        next
    }
}

n_distinct(results$exp, results$outc)

# output results
write.csv(results, "/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/mr_chip_to_prot.csv", quote = F, row.names = F)
write.csv(log, "/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/mr_chip_to_prot_statuslog.csv", quote = F, row.names = F)




####### PROT TO CHIP #######
results <- data.frame()
log <- data.frame()

# # reading in MR results, yes corr
# for (ss in SS){
#     for (protein in protlist){
#         file <- paste0("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/prot_to_chip_raw/",protein, "_", ss,".csv")
#         if (file.exists(file)){
#             df <- data.frame(fread(file, header = T))
#             results <- rbind(results, df)
#             log <- rbind(log, data.frame(protein = protein, sumstats = ss, status = "Y"))
#         } else {
#             print(paste0(ss, " & ", protein, " file doesn't exist"))
#             log <- rbind(log, data.frame(protein = protein, sumstats = ss, status = "N"))
#             next
#         }
#     }
# }
# 
# write.csv(results, "/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/mr_prot_to_chip.csv", quote = F, row.names = F)
# write.csv(log, "/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/mr_prot_to_chip_statuslog.csv", quote = F, row.names = F)





##### no corr
# reading in MR results
results <- data.frame()
log <- data.frame()

for (ss in SS){
    for (protein in protlist){
        file <- paste0("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/prot_to_chip_raw/",protein, "_", ss,"_nocorr.csv")
        if (file.exists(file)){
            df <- data.frame(fread(file, header = T))
            results <- rbind(results, df)
            log <- rbind(log, data.frame(protein = protein, sumstats = ss, status = "Y"))
        } else {
            print(paste0(ss, " & ", protein, " file doesn't exist"))
            log <- rbind(log, data.frame(protein = protein, sumstats = ss, status = "N"))
            next
        }
    }
}
n_distinct(results$exp, results$outc)

write.csv(results, "/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/mr_prot_to_chip_nocorr.csv", quote = F, row.names = F)
write.csv(log, "/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/mr_prot_to_chip_nocorr_statuslog.csv", quote = F, row.names = F)


####### CHIP TO PROT, SOMASCAN #######
# ***waiting for the full protein list but here are the selected results between CHIP SSs and certain proteins
# 22 proteins of interest
protlist <- unlist(fread("/medpop/esp2/lli/Zhi-mr_protein/Data/somascan_protein_ss/new_protlist", header = F))
log <- data.frame()
results <- data.frame()

for (protein in protlist){
    file <- paste0("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/somascan_chip_to_22prot_raw/", protein, "_3chip_raw.csv")
    if (file.exists(file)){
        df <- data.frame(fread(file, header = T))
        results <- rbind(results, df)
        log <- rbind(log, data.frame(protein = protein, status = "Y"))
    } else {
        print(paste0(protein, " file doesn't exist"))
        log <- rbind(log, data.frame(protein = protein, status = "N"))
        next
    }
}

n_distinct(results$exp, results$outc)

# output results
write.csv(results, "/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/mr_chip_to_22prot_somascan.csv", quote = F, row.names = F)
write.csv(log, "/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/mr_chip_to_22prot_somascan_statuslog.csv", quote = F, row.names = F)




####### PROT TO CHIP, SOMASCAN ########
results <- data.frame()
for (ss in SS){
    file <- paste0("/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/somascan_22prot_to_chip_raw/", "mr_prot_to_", ss, "_somascan.csv")
    if (file.exists(file)){
        df <- data.frame(fread(file, header = T))
        results <- rbind(results, df)
    } else {
        print(paste0(ss, " & ", protein, " file doesn't exist"))
        next
    }
}

n_distinct(results$exp, results$outc)
write.csv(results, "/medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/results/mr_22prot_to_chip_somascan.csv", quote = F, row.names = F)



