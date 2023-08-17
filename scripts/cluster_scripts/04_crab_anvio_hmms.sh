#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH --mem 180G
#SBATCH -t 144:00:00 #time in hours:min:sec
#SBATCH --partition=production
#SBATCH -D /share/eisenlab/casett/crabs_eukrep/
#SBATCH -e anvio_mapping_err_euk.txt
#SBATCH -J crabs_anvio_hmms_2020_euk

#source /share/eisenlab/gjospin/.profile

module load anvio/6.2
module load prodigal

source activate anvio-6.2


anvi-gen-contigs-database -f Crab_euk_contigs_fixed.fa -o Crab_Euk.db

anvi-run-hmms -c Crab_Euk.db --num-threads 24

anvi-run-hmms -T 24 -c Crab_Euk.db -H /share/eisenlab/gjospin/software/anvio/anvio/data/hmm/Bacteria_71
anvi-run-hmms -T 24 -c Crab_Euk.db -H /share/eisenlab/gjospin/software/anvio/anvio/data/hmm/Archaea_76
anvi-run-hmms -T 24 -c Crab_Euk.db -H /share/eisenlab/gjospin/software/anvio/anvio/data/hmm/Ribosomal_RNAs
anvi-run-hmms -T 24 -c Crab_Euk.db -H /share/eisenlab/gjospin/software/anvio/anvio/data/hmm/Protista_83


#can't do COGs b/c of permissions issues
#anvi-setup-ncbi-cogs --num-threads 24 --reset

#anvi-run-ncbi-cogs -c SGChyt.db --num-threads 24

anvi-get-sequences-for-gene-calls -c Crab_Euk.db -o Crab_Euk_gene_calls.fa

#make text file of all sample ids
#for sample in $(cat File_names_euk.txt); 
#do anvi-profile -i $sample'.er.bam' -c SGChyt_Euk.db --num-threads 48 --min-contig-length 1000; 
#done




