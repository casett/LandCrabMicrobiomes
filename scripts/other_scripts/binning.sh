{\rtf1\ansi\ansicpg1252\cocoartf2822
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\froman\fcharset0 Times-Roman;\f1\froman\fcharset0 Times-Bold;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs24 \cf0 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 #!/bin/bash\
#\
#SBATCH --job-name binner\
#SBATCH --cpus-per-task=24\
#SBATCH --mem=512GB\
#SBATCH \strike \strikec0 -output=binner\strike0\striked0 %j.out\
#SBATCH \strike \strikec0 -error=binner\strike0\striked0 %j.err\
#SBATCH --partition=fat\
#SBATCH --qos=long\
#SBATCH --time=120:00:00 \
#sbatch .sh <samplename>\
#do sh #run in assembly folder with reads in ../*fastq.gz\
\
\
#save current working directory\
WD="$PWD"\
SAMPLE_ID=$1\
\
module load singularity\
module load anaconda3\
source $ANACONDA3_ROOT/etc/profile.d/conda.sh \
\
\
############################################## Processing contigs ############################################## \
\
conda activate /usr/users/laetitia.wilkins/anaconda3/envs/bbmap\
#remove all contigs below 1000\
\
reformat.sh in=porcelains.contigs.fa out=contigs-fixed.fa minlength=1000\
\
conda deactivate\
\
\
#simplify names\
sed -i 's/NODE_/c/g' contigs-fixed.fa\
sed -i 's/_length.*//g' contigs-fixed.fa\
\
rsync -v -L contigs-fixed.fa ../*fastq.gz $TMP_SCRATCH\
\
\
############################################## Mapping ############################################## \
\
#map reads to scaffolds and convert to bams, then sort and index the bam file\
\
module load bwa/0.7.12\
module load samtools/1.12\
\
bwa index $TMP_SCRATCH/contigs-fixed.fa\
\
for f in $TMP_SCRATCH/
\f1\b fastq.gz; do bwa mem -t $SLURM_CPUS_PER_TASK $TMP_SCRATCH/contigs-fixed.fa $TMP_SCRATCH/$\{f##
\f0\b0 /\} | samtools sort -o $TMP_SCRATCH/$\{1\}_$\{f##*/\}_sorted.bam; done\
\
for f in $TMP_SCRATCH/
\f1\b fastq.gz; do samtools index $TMP_SCRATCH/$\{SAMPLE_ID\}_$\{f##
\f0\b0 /\}_sorted.bam; done \
\
module unload bwa\
module unload samtools\
\
echo 'bams made, sorted, and indexed'\
\
\
\
\
############################################## Binning ############################################## \
\
###################### Metabat ######################\
\
#run metabat wrapper script and metabat without coverage info\
\
echo metabat with coverage\
\
singularity exec -B $TMP_SCRATCH:/data,$WD /scratch/projects/eei/software/metabat2_2.15--h4da6f23_2.sif \\\
runMetaBat.sh --minContig 1500 -t $SLURM_CPUS_PER_TASK /data/contigs-fixed.fa /data/*sorted.bam}