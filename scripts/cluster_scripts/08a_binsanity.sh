#!/bin/bash -l
#
#SBATCH -n 48 #number cores
#SBATCH --mem 980G
#SBATCH -t 167:00:00 #time in hours:min:sec
#SBATCH --partition=production
#SBATCH -D /share/eisenlab/casett/crabs_eukrep/megahit_coassembly_mapping/
#SBATCH -e anvio_binsanity_euk_2020.txt
#SBATCH -J crabs_binsanity



module load anvio/6.2
module load binsanity/0.3.3

source activate anvio-6.2

anvi-cluster-contigs -p Crab_merged_kaiju_1000bp_EUK/PROFILE.db -c Crab_Euk.db -C BINSANITY --driver binsanity -T 48  --just-do-it








