#!/bin/bash -l
#
#SBATCH -e logs/23_metabolic.log
#SBATCH -o logs/23_metabolic.log
#SBATCH --nodes 1 --ntasks 24 -p intel
#SBATCH -J crab_metabolic --mem 250G

conda activate /rhome/cassande/bigdata/.conda/envs/METABOLIC_v4.0

BINS=whokaryote_results/who_bac
CPU=24
OUT=metabolic_output_samples
META=/rhome/cassande/bigdata/eisenlab/sg_phor/metabolic_input/METABOLIC

mkdir $OUT

perl $META/METABOLIC-G.pl -t $CPU -in-gn $BINS -p meta -o $OUT