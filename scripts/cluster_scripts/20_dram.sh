#!/bin/bash -l
#
#SBATCH -n 48 #number cores
#SBATCH -e logs/20_dramv.log
#SBATCH -o logs/20_dramv.log
#SBATCH --mem 400G #memory per node in Gb
#SBATCH -p batch
#SBATCH -J crab_mag_dram

CPU=48
DIR=bac_mags/fin_bins
OUT=dram


mkdir $OUT


conda activate DRAM

#DRAM-setup.py import_config --config_loc  my_old_config.txt
#DRAM-setup.py prepare_databases --skip_uniref --output_dir /rhome/cassande/bigdata/software/db-dramv

# step 1 annotate
DRAM.py annotate -i 'bac_mags/fin_bins/*.fa' -o annotation --threads $CPU

#step 2 summarize anntotations
DRAM.py distill -i annotation/annotations.tsv -o genome_summaries --trna_path annotation/trnas.tsv --rrna_path annotation/rrnas.tsv