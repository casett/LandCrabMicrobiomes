#!/bin/bash -l
#
#SBATCH -n 8 #number cores
#SBATCH -e logs/pigz.log
#SBATCH -o logs/pigz.log
#SBATCH --mem 20G #memory per node in Gb
#SBATCH -p intel,batch
#SBATCH -D /rhome/cassande/bigdata/eisenlab/crabs_eukrep/
#SBATCH -J crab_zip


pigz fastqs/*.fastq
