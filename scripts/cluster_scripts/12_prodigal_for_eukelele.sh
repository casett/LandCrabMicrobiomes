#!/bin/bash -l
#
#SBATCH --ntasks 8 #number cores
#SBATCH -J crab_prodigal
#SBATCH --mem=50G #memory
#SBATCH -p intel,batch
#SBATCH -D /rhome/cassande/bigdata/eisenlab/crabs_eukrep/
#SBATCH -o logs/12_prodigal.log
#SBATCH -e logs/12_prodigal.log


module unload miniconda2
module load miniconda3

conda activate anvio-7

NTBINS=megahit_coassembly_mapping/Crab_merged_kaiju_1000bp_EUK/all_bins

cd $NTBINS
for f in *.fa;
do prodigal -i $f -o $f'_genes' -a $f'_proteins'.faa;
done

cd ..
mkdir all_bins_prot
mv all_bins/*.faa all_bins_prot

conda deactivate

