#!/bin/bash -l
#
#SBATCH --ntasks 8 #number cores
#SBATCH -J crab_eukelele
#SBATCH --mem=150G #memory
#SBATCH -p intel,batch
#SBATCH -D /rhome/cassande/bigdata/eisenlab/crabs_eukrep/
#SBATCH -o logs/13_eukelele.log
#SBATCH -e logs/13_eukelele.log



module unload miniconda2
module load miniconda3


conda activate EUKulele

BINFOLDER=megahit_coassembly_mapping/Crab_merged_kaiju_1000bp_EUK/all_bins_prot
#BINFOLDER contain proteins in .faa for each mag
EUKulele --sample_dir $BINFOLDER -m mags
