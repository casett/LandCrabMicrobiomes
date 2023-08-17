#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH --mem 96G
#SBATCH -t 96:00:00 #time in hours:min:sec
#SBATCH --partition=production
#SBATCH -D /share/eisenlab/casett/crabs_eukrep/megahit_coassembly_mapping/
#SBATCH -e anvio_bin_summary.txt
#SBATCH -J crab_euk_summary



module load anvio/6.2

source activate anvio-6.2

anvi-summarize  -p Crab_merged_kaiju_1000bp_EUK/PROFILE.db -c Crab_Euk.db   -o  Crab_merged_kaiju_1000bp_EUK/sample_summary_MAXBIN -C MAXBIN

anvi-summarize -p Crab_merged_kaiju_1000bp_EUK/PROFILE.db -c Crab_Euk.db   -o  Crab_merged_kaiju_1000bp_EUK/sample_summary_METABAT -C METABAT

anvi-summarize  -p Crab_merged_kaiju_1000bp_EUK/PROFILE.db -c Crab_Euk.db   -o  Crab_merged_kaiju_1000bp_EUK/sample_summary_CONCOCT -C CONCOCT


