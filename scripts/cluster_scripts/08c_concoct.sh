#!/bin/bash -l
#
#SBATCH -n 48 #number cores
#SBATCH --mem 248G
#SBATCH -t 167:00:00 #time in hours:min:sec
#SBATCH --partition=production
#SBATCH -D /share/eisenlab/casett/crabs_eukrep/megahit_coassembly_mapping/
#SBATCH -e anvio_concoct_euk_2020.txt
#SBATCH -J crabs_concoct



module load anvio/6.2
module load concoct/1.1.0

source activate concoct-1.1.0
conda activate --stack anvio-6.2
#source activate anvio-6.1

anvi-cluster-contigs -p Crab_merged_kaiju_1000bp_EUK/PROFILE.db -c Crab_Euk.db -C CONCOCT --driver concoct -T 48  --just-do-it





