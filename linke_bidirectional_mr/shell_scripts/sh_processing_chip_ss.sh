#!/bin/bash -l
#$ -wd /medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr
#$ -N shorter_ss
#$ -l h_vmem=64G
#$ -o /medpop/esp2/lli/Zhi-mr_protein/logs/shorter_ss.log
#$ -e /medpop/esp2/lli/Zhi-mr_protein/logs/shorter_ss.log
#$ -t 1

source /broad/software/scripts/useuse
use Anaconda3
source activate synapser_env

Rscript /medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/a_processing_chip_ss.R
# Rscript /medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/a_processing_chip_ss_noukb.R