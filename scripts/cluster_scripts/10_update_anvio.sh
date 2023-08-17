#!/bin/bash -l
#
#SBATCH --ntasks 8 #number cores
#SBATCH -J crab_anvio_update
#SBATCH --mem=50G #memory
#SBATCH -p intel,batch
#SBATCH -D /rhome/cassande/bigdata/eisenlab/crabs_eukrep/
#SBATCH -o 10_update.log
#SBATCH -e 10_update.log




module unload miniconda2
module unload anaconda3
module load miniconda3

source activate anvio-7

anvi-migrate --migrate-dbs-safely *.db

anvi-migrate --migrate-dbs-safely megahit_coassembly_mapping/Crab_merged_kaiju_1000bp_EUK/*.db
