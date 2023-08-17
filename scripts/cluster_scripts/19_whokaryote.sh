#!/bin/bash -l
#SBATCH --ntasks=16 # Number of cores
#SBATCH -J crab_whokaryote
#SBATCH --mem=50G #memory
#SBATCH -p batch
#SBATCH -D /rhome/cassande/bigdata/eisenlab/crabs_eukrep/
#SBATCH -o logs/19_whokaryote.log
#SBATCH -e logs/19_whokaryote.log


conda activate whokaryote

INDIR=assemblies_for_fun
CPU=16
MIN=1000

#mkdir $OUT

#default minlen = 5000
for ASSEMBLY in $(ls $INDIR);
do 
	BASE=$(basename $ASSEMBLY .fasta)

	whokaryote.py --contigs $INDIR/$ASSEMBLY --outdir $BASE.whokaryote --f --minsize $MIN 
done