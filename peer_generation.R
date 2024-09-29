#######################################
# PEER Analysis
gc()
rm(list=ls()) # Clear all objects in memory
#######################################

# Load required libraries
library(data.table)
library(dplyr)
library(peer)
library(doParallel)

# Set up parallel processing
num_cores <- detectCores() - 1 # Reserve 1 core for other tasks
registerDoParallel(num_cores)

# Define parameters
n_peer <- 150 # Number of PEER factors, adjust based on server capabilities
cohort <- "UKB" # Name of the cohort

# Define file paths
data_path <- '/medpop/esp2/zyu/chip_protemoics/data/UKB/UKB_proteomics_imputed.csv'
result_path <- '/medpop/esp2/zyu/chip_protemoics/data/UKB/PEER/'

# Create necessary directories if they don't exist
dir.create(file.path(result_path, 'peernum', 'peerandresid'), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(result_path, 'peernum', 'other'), showWarnings = FALSE, recursive = TRUE)

# Load proteomics data
olink <- fread(data_path)

# Prepare proteomics data for PEER (excluding the first column if it's Sample IDs)
prot <- as.matrix(olink[, -1, with = FALSE]) # Ensure 'with = FALSE' for matrix conversion

# Configure and run PEER model
model <- PEER()

# Set phenotype data (proteomics data)
PEER_setPhenoMean(model, prot)

# Set the number of factors to be inferred (Nk)
PEER_setNk(model, n_peer)

# Enable mean factor to adjust for mean expression level
PEER_setAdd_mean(model, TRUE)

# Run the PEER model (factor inference and update)
PEER_update(model)

# Extract and save the inferred factors (PEER factors)
factors <- PEER_getX(model)
colnames(factors) <- c('MeanFactor', paste('peer', 1:n_peer, sep=''))
peer_factors <- data.table(SampleId = olink[[1]], factors[, -1, drop = FALSE]) # Drop 'MeanFactor'
fwrite(peer_factors, file.path(result_path, 'peernum', 'peerandresid', paste('peers_', cohort, '_', n_peer, '.txt')))

# Extract and save residuals (adjusted expression levels)
residuals <- PEER_getResiduals(model)
colnames(residuals) <- colnames(prot)
residuals_data <- data.table(SampleId = olink[[1]], residuals)
fwrite(residuals_data, file.path(result_path, 'peernum', 'peerandresid', paste('resid_', cohort, '_', n_peer, '.txt')))

# Extract and save precision values (weights of each factor)
precision <- PEER_getAlpha(model)
write.table(precision, file.path(result_path, 'peernum', 'other', paste('precision_', n_peer, '.txt')), col.names = FALSE, row.names = FALSE)

# Extract and save the weights (W matrix)
weights <- PEER_getW(model)
write.table(weights, file.path(result_path, 'peernum', 'other', paste('weights_', n_peer, '.txt')), col.names = FALSE, row.names = FALSE)

# Cleanup and close parallel processing
stopImplicitCluster()

# Final garbage collection to free memory
gc()
