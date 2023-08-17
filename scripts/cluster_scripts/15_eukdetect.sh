#!/bin/bash -l
#
#SBATCH --ntasks 16 #number cores
#SBATCH -J crab_eukdetect
#SBATCH --mem=100G #memory
#SBATCH -p batch
#SBATCH -D /rhome/cassande/bigdata/eisenlab/crabs_eukrep/
#SBATCH -o logs/15_eukdetect.log
#SBATCH -e logs/15_eukdetect.log



module load samtools
module load miniconda3

source activate eukdetect-r3

CONFIG=eukdetect_config.yml
CPU=16

eukdetect --mode runall --configfile $CONFIG --cores $CPU


