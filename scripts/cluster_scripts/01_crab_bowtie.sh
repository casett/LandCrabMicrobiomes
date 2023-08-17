#!/bin/bash -l
#
#SBATCH --ntasks 8 #number cores
#SBATCH -D /share/eisenlab/casett/crabs_eukrep/
#SBATCH -e bowtie_stderr_eukrep.txt
#SBATCH -J crab_bowtie_coassembly_eukrep
#SBATCH --mem=48000 #memory
#SBATCH -t 12:00:00 #time in hours:min:sec
#SBATCH --partition=production


module load bowtie2
#module load anvio/6.2
#source activate anvio-6.2

#anvi-script-reformat-fasta Crab_euk_contigs.fa -o Crab_euk_contigs_fixed.fa -l 0 --simplify-names

#cd megahit_coassembly_mapping

bowtie2-build Crab_euk_contigs_fixed.fa megahit_coassembly_mapping/Crab_eukrep_contigs

