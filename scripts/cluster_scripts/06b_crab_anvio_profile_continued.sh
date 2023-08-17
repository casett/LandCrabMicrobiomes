#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH --mem 180G
#SBATCH -t 144:00:00 #time in hours:min:sec
#SBATCH --partition=production
#SBATCH -D /share/eisenlab/casett/crabs_eukrep/
#SBATCH -e anvio_mapping_err_2ndhalf.txt
#SBATCH -J crabs_profile_2ndhalf

#source /share/eisenlab/gjospin/.profile

module load anvio/6.2
module load prodigal

source activate anvio-6.2


#make text file of all sample ids
for sample in $(cat File_names_2ndHalf.txt); 
do anvi-profile -i megahit_coassembly_mapping/$sample'.er.bam' -c Crab_Euk.db  --num-threads 24 --min-contig-length 1000; 
done


