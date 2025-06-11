#!/bin/bash -l
#$ -wd /medpop/esp2/lli/Zhi-mr_protein/
#$ -N prot_to_chip
#$ -l h_vmem=40G
#$ -o /medpop/esp2/lli/Zhi-mr_protein/logs/prot_to_chip.log
#$ -e /medpop/esp2/lli/Zhi-mr_protein/logs/prot_to_chip.log
#$ -t 1-242

i=$(expr ${SGE_TASK_ID} - 1)
# create looping parameters
protlist=($(cat /medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/protlist_clean))
protein=${protlist[$i]}

source /broad/software/scripts/useuse
use Anaconda3
source activate synapser_env

# run R
for SS in CHIP TET2 DNMT3A; do
echo "Running MR for protein: $protein and summary stat: $SS"
  Rscript /medpop/esp2/zyu/chip_protemoics/code/CHIP_proteomics_CAD/mr/linke_mr_prot_to_chip.R -p ${protein} -s ${SS}
done
