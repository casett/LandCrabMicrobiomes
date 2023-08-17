#!/bin/bash -l
#
#SBATCH -n 48 #number cores
#SBATCH --mem 490G
#SBATCH -t 100:00:00 #time in hours:min:sec
#SBATCH --partition=production
#SBATCH -D /share/eisenlab/casett/crabs_eukrep/megahit_coassembly_mapping/
#SBATCH -e anvio_merge_err.txt
#SBATCH -J crabs_merge


module load anvio/6.2
module load prodigal



source activate anvio-6.2







anvi-merge *.er.bam-ANVIO_PROFILE/PROFILE.db -o Crab_merged_kaiju_1000bp_EUK/ -c Crab_Euk.db --enforce-hierarchical-clustering

