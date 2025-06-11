#!/bin/bash -l
#$ -wd /medpop/esp2/lli/Zhi-mr_protein/
#$ -N soma_chip_to_prot
#$ -l h_vmem=30G
#$ -o /medpop/esp2/lli/Zhi-mr_protein/logs/new_soma_chip_to_prot.log
#$ -e /medpop/esp2/lli/Zhi-mr_protein/logs/new_soma_chip_to_prot.log
#$ -t 1-22

i=$(expr ${SGE_TASK_ID} - 1)
# create looping parameters
# protlist=($(cat /medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/protlist_clean))
protlist=($(cat /medpop/esp2/lli/Zhi-mr_protein/Data/somascan_protein_ss/new_protlist))
protein=${protlist[$i]}

source /broad/software/scripts/useuse
use Anaconda3
source activate synapser_env

Rscript /medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/linke_mr_chip_to_prot_somascan.R -p ${protein}
