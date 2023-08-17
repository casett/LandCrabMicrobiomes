#!/bin/bash -l
#
#SBATCH --ntasks 8 #number cores
#SBATCH -J crab_motu
#SBATCH --mem=50G #memory
#SBATCH -p intel,batch
#SBATCH -D /rhome/cassande/bigdata/eisenlab/crabs_eukrep/fastqs/
#SBATCH -o logs/14_motu.log
#SBATCH -e logs/14_motu.log


module unload miniconda2
module load miniconda3


source activate mOTUs



for file in $(cat File_names.txt);
do 
	motus profile -f $file'_L002_R1_001.fastq.gz.reads.filtered_1.fastq' -r $file'_L002_R2_001.fastq.gz.reads.filtered_2.fastq' -g 1 -l 30 > $file.motu.txt;
done

motus merge -i 120_S25.motu.txt,140_S26.motu.txt,147_S27.motu.txt,177_S28.motu.txt,178_S29.motu.txt,179_S31.motu.txt,209_S30.motu.txt,214_S32.motu.txt,82_S21.motu.txt,88_S22.motu.txt,99_S23.motu.txt,BLB_S18.motu.txt,BLB_WG4_S19.motu.txt,BLR_S17.motu.txt,BLR_WG4_S20.motu.txt,BLS_S16.motu.txt,DC06_S1.motu.txt,DC10_S2.motu.txt,DC11_S3.motu.txt,DC17_S4.motu.txt,DM13_S9.motu.txt,DM14_S10.motu.txt,DM15_S11.motu.txt,DM16_S12.motu.txt,GG01_S14.motu.txt,GH01_S13.motu.txt,GN13_S5.motu.txt,GN15_S6.motu.txt,GN18_S7.motu.txt,GN19_S8.motu.txt,NEG_S33.motu.txt,OC01_S15.motu.txt > all_sample_profiles.txt
