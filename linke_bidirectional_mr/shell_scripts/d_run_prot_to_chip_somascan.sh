#!/bin/bash -l
#$ -wd /medpop/esp2/lli/Zhi-mr_protein/
#$ -N new_soma_prot_to_chip
#$ -l h_vmem=30G
#$ -o /medpop/esp2/lli/Zhi-mr_protein/logs/new_soma_prot_to_chip
#$ -e /medpop/esp2/lli/Zhi-mr_protein/logs/new_soma_prot_to_chip
#$ -t 1


source /broad/software/scripts/useuse
use Anaconda3
source activate synapser_env

# for SS in CHIP TET2 DNMT3A; do
for SS in DNMT3A; do
    Rscript /medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/linke_mr_prot_to_chip_somascan.R -s ${SS}
done;