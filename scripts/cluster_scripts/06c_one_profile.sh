#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH --mem 180G
#SBATCH -t 48:00:00 #time in hours:min:sec
#SBATCH --partition=production
#SBATCH -D /share/eisenlab/casett/crabs_eukrep/
#SBATCH -e anvio_mapping_err_failedsamp.txt
#SBATCH -J crabs_profile_failedsample

#source /share/eisenlab/gjospin/.profile

module load anvio/6.2
module load prodigal

source activate anvio-6.2



anvi-init-bam /share/eisenlab/casett/crabs_eukrep/megahit_coassembly_mapping/140_S26.er-RAW.bam -o /share/eisenlab/casett/crabs_eukrep/megahit_coassembly_mapping/140_S26.er.bam

anvi-profile -i megahit_coassembly_mapping/140_S26.er.bam -c Crab_Euk.db  --num-threads 24 --min-contig-length 1000
