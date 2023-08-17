#!/bin/bash -l
#
#SBATCH --ntasks 8 #number cores
#SBATCH -J crab_eukCC
#SBATCH --mem=50G #memory
#SBATCH -p intel,batch
#SBATCH -D /rhome/cassande/bigdata/eisenlab/crabs_eukrep/
#SBATCH -o logs/11_eukcc.log
#SBATCH -e logs/11_eukcc.log



module unload miniconda2
module load miniconda3

conda activate eukcc2

export EUKCC2_DB=/rhome/cassande/bigdata/software/eukccdb/eukcc2_db_ver_1

BINFOLDER=/rhome/cassande/bigdata/eisenlab/crabs_eukrep/megahit_coassembly_mapping/Crab_merged_kaiju_1000bp_EUK/all_bins
OUTFOLDER=/rhome/cassande/bigdata/eisenlab/crabs_eukrep/megahit_coassembly_mapping/Crab_merged_kaiju_1000bp_EUK/eukcc_results
CPU=8

#mkdir $OUTFOLDER
#files in folder need to end in .fa
eukcc folder --out $OUTFOLDER --threads $CPU $BINFOLDER



