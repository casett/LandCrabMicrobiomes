#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH --mem 96G #memory per node in Gb
#SBATCH -t 96:00:00 #time in hours:min:sec
#SBATCH --partition=production
#SBATCH -D /share/eisenlab/casett/crabs_eukrep/
#SBATCH -e eukrep_err.txt
#SBATCH -J crabs_eukrep



module load eukrep/0.6.6

EukRep -i invasionSea_megahit_NOTinterleaved.contigs.fa -o Crab_euk_contigs.fa --min 1000
