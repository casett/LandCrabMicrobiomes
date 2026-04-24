{\rtf1\ansi\ansicpg1252\cocoartf2822
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\froman\fcharset0 Times-Roman;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;}
\margl1440\margr1440\vieww20400\viewh13900\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs24 \cf0 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 #!/bin/bash\
#\
#SBATCH --job-name=MitoZ \
#SBATCH --cpus-per-task=16\
#SBATCH --mem=100gb\
#SBATCH --time=05:30:00 \
#SBATCH -C scratch\
#SBATCH -p medium\
#SBATCH --output=mitoZ_slurm_%j.log\
\
### usage ###\
#sbatch 01_mitoz.sh sample_name \
\
#Split interleaved reads\
#source /usr/users/laetitia.wilkins/.bashrc\
#conda activate /usr/users/laetitia.wilkins/anaconda3/envs/bbmap\
\
#rsync -v -L $1 $TMP_SCRATCH\
\
pigz -d -p 16 *reads.filtered_1.fastq.gz\
pigz -d -p 16 *reads.filtered_2.fastq.gz\
\
#reformat.sh in1=$TMP_SCRATCH/$\{1\} out1=$TMP_SCRATCH/$\{1%%.*\}_R1.fastq out2=$TMP_SCRATCH/$\{1%%.*\}_R2.fastq\
\
\
mkdir $PWD/$\{1\}_mitoz_skmer\
\
#Move split read pairs into output directory\
rsync -v -L *_R*fastq $PWD/$\{1\}_mitoz_skmer\
\
module purge\
module load singularity\
\
\
singularity exec -B $PWD \\\
/scratch/projects/eei/software/MitoZ_v3.4.sif \\\
mitoz all \\\
--workdir $PWD/$\{1\}_mitoz_skmer \\\
--outprefix $\{1\} \\\
--thread_number 16 \\\
--clade Arthropoda \\\
--genetic_code 5 \\\
--fq1 $PWD/$\{1\}_mitoz_skmer/*filtered_1.fastq \\\
--fq2 $PWD/$\{1\}_mitoz_skmer/*filtered_2.fastq \\\
--fastq_read_length 151 \\\
--data_size_for_mt_assembly 0 \\\
--assembler megahit \\\
--kmers_megahit 39 59 79 99 119 141 \\\
--memory 100 \\\
--requiring_taxa Arthropoda\
\
\
rm -rf $PWD/$\{1\}_mitoz_skmer/*fastq}