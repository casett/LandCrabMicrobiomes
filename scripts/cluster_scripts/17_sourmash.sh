#!/bin/bash -l
#
#SBATCH --ntasks 16 #number cores
#SBATCH -J sourmash_test
#SBATCH --mem=250G #memory
#SBATCH -p intel
#SBATCH -D /rhome/cassande/bigdata/eisenlab/crabs_eukrep/fastqs/
#SBATCH -o ../logs/17_sourmash.log
#SBATCH -e ../logs/17_sourmash.log


conda activate sourmash_to_phylo


#for infile in *_L002_R1_001.fastq.gz.reads.filtered_1.fastq.gz
#do
#    bn=$(basename ${infile} _L002_R1_001.fastq.gz.reads.filtered_1.fastq.gz)
#    sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --merge ${bn} -o ${bn}.sig ${infile} ${bn}_L002_R2_001.fastq.gz.reads.filtered_2.fastq.gz
#done


#curl -JLO https://osf.io/k2u8s/download
#curl -JLO https://osf.io/v3zmg/download
#gunzip gtdb-rs207.taxonomy.csv.gz


for infile in *sig
do
    bn=$(basename $infile .sig)
    sourmash gather ${infile} gtdb-rs207.genomic.k31.zip -o ${bn}_gather_gtdbrs207_reps.csv
done


sourmash tax prepare -t gtdb-rs207.taxonomy.csv -o gtdb-rs207.taxonomy.sqldb -F sql


for infile in *_gather_gtdbrs207_reps.csv
do
    sourmash tax annotate -g ${infile} -t gtdb-rs207.taxonomy.sqldb 
done

