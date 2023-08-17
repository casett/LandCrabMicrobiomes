#!/bin/bash -l
#
#SBATCH -n 16 #number cores
#SBATCH --mem 98000
#SBATCH -t 144:00:00 #time in hours:min:sec
#SBATCH --partition=production
#SBATCH -D /share/eisenlab/casett/crabs_eukrep/fastqs/
#SBATCH -e bowtie_mapping_err.txt
#SBATCH -J crab_bam_coassembly


module load bowtie2
module load samtools
module load anvio/6.2

source activate anvio-6.2

for file in $(cat File_names.txt); 
do bowtie2 --threads 16 -x /share/eisenlab/casett/crabs_eukrep/megahit_coassembly_mapping/Crab_eukrep_contigs -1 $file'_L002_R1_001.fastq.gz.reads.filtered_1.fastq' -2 $file'_L002_R2_001.fastq.gz.reads.filtered_2.fastq' -S /share/eisenlab/casett/crabs_eukrep/megahit_coassembly_mapping/$file'.er.sam'; 
done

for file in $(cat File_names.txt); 
do samtools view -F 4 -bS /share/eisenlab/casett/crabs_eukrep/megahit_coassembly_mapping/$file'.er.sam' > /share/eisenlab/casett/crabs_eukrep/megahit_coassembly_mapping/$file'.er-RAW.bam'; 
done

for file in $(cat File_names.txt); 
do anvi-init-bam /share/eisenlab/casett/crabs_eukrep/megahit_coassembly_mapping/$file'.er-RAW.bam' -o /share/eisenlab/casett/crabs_eukrep/megahit_coassembly_mapping/$file'.er.bam'; 
done

