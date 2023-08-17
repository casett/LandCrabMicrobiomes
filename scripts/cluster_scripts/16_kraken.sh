#!/bin/bash -l
#
#SBATCH --ntasks 16 #number cores
#SBATCH -J krakendb
#SBATCH --mem=700G #memory
#SBATCH -p highmem
#SBATCH -D /rhome/cassande/bigdata/eisenlab/crabs_eukrep/
#SBATCH -o logs/16_krak.log
#SBATCH -e logs/16_krak.log
#SBATCH -t 14-00:00

#module load bracken

conda activate kraken2

CPU=16

export KRAKEN2_DB_PATH="/rhome/cassande/bigdata/eisenlab/crabs_eukrep/"
#kraken2-build --standard --db stan --threads 32
#wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20210517.tar.gz
#tar -zxvf k2_standard_20210517.tar.gz -C stan
#bracken-build -d stan -t $CPU -k 35 -l 150

#kraken2-build --download-taxonomy --db nt
#kraken2-build --download-library nt --db nt
#kraken2-build --build --threads $CPU --db nt --fast-build
#bracken-build -d nt -t $CPU -k 35 -l 250


for prefix in `ls fastqs/*.gz | cut -f1,2 -d'_' | sort -u`;
do
	echo $prefix
	read1=( ${prefix}*_L002_R1_001.fastq.gz.reads.filtered_1.fastq.gz ) #the parentheses assign the globbed filename to an array (of length 1)
	read2=( ${prefix}*_L002_R2_001.fastq.gz.reads.filtered_2.fastq.gz )

#	kraken2 --db nt --threads $CPU --report ${prefix}_kraken2_report_paired.tsv --gzip-compressed --paired ${read1} ${read2} > ${prefix}.log

#	bracken -d nt -i ${prefix}_kraken2_report_paired.tsv -o ${prefix}.S.bracken -r 250 -l 'S' -t 10
#	bracken -d nt -i ${prefix}_kraken2_report_paired.tsv -o ${prefix}.F.bracken -r 250 -l 'F' -t 10
#	bracken -d nt -i ${prefix}_kraken2_report_paired.tsv -o ${prefix}.G.bracken -r 250 -l 'G' -t 10

	bracken -d nt -i ${prefix}_kraken2_report_paired.tsv -o ${prefix}.K.bracken -r 250 -l 'K' -t 10
	bracken -d nt -i ${prefix}_kraken2_report_paired.tsv -o ${prefix}.P.bracken -r 250 -l 'P' -t 10
	bracken -d nt -i ${prefix}_kraken2_report_paired.tsv -o ${prefix}.C.bracken -r 250 -l 'C' -t 10
	bracken -d nt -i ${prefix}_kraken2_report_paired.tsv -o ${prefix}.O.bracken -r 250 -l 'O' -t 10
done


