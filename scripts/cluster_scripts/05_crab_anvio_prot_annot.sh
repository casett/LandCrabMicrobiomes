#!/bin/bash -l
#
#SBATCH --ntasks 16 #number cores
#SBATCH -J crab_anvio_annot
#SBATCH --mem=350G #memory
#SBATCH -p intel,batch
#SBATCH -D /rhome/cassande/bigdata/eisenlab/crabs_eukrep/
#SBATCH -o logs/5_annot.log
#SBATCH -e logs/5_annot.log


module unload miniconda2
module unload anaconda3
module load miniconda3
module load eggnog-mapper/1.0.3 
#module load interproscan/5.52-86.0

#source activate anvio-7


#anvi-get-sequences-for-gene-calls -c Crab_Euk.db --get-aa-sequences -o Crab_amino_acid_sequences.fa

#anvi-setup-ncbi-cogs --num-threads 24

#anvi-run-ncbi-cogs -c Crab_Euk.db --num-threads 24


#emapper.py -i Crab_amino_acid_sequences_fix.fa --output Crab_prot_fix -m diamond --cpu 24

source activate anvio-7
#anvi-script-run-eggnog-mapper -c Crab_Euk.db --annotation Crab_prot_fix.emapper.annotations --use-version 1.0.3

module load interproscan/5.52-86.0
interproscan.sh -i Crab_amino_acid_sequences_fix.fa -f tsv -o Crab_interpro-output.tsv


anvi-import-functions -c Crab_Euk.db -i Crab_interpro-output.tsv -p interproscan


