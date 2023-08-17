#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH --mem 96G #memory per node in Gb
#SBATCH -t 96:00:00 #time in hours:min:sec
#SBATCH --partition=production
#SBATCH -D /share/eisenlab/casett/crabs_eukrep/
#SBATCH -e kaiju_euk_err.txt
#SBATCH -J crab_kaiju_euk



/share/eisenlab/gjospin/software/kaiju/src/kaiju -z 24 -t /share/eisenlab/casett/database/kaiju_nr_euk/nodes.dmp -f /share/eisenlab/casett/database/kaiju_nr_euk/kaiju_db_nr_euk.fmi -i Crab_Euk_gene_calls.fa -o Crab_gene_calls.kaiju.out.updated.euk -v

/share/eisenlab/gjospin/software/kaiju/src/kaijuReport -t /share/eisenlab/casett/database/kaiju_nr_euk/nodes.dmp -n /share/eisenlab/casett/database/kaiju_nr_euk/names.dmp -i Crab_gene_calls.kaiju.out.updated.euk -r genus -p -o Crab_gene_calls.summary.updated.euk

/share/eisenlab/gjospin/software/kaiju/src/addTaxonNames -t /share/eisenlab/casett/database/kaiju_nr_euk/nodes.dmp -n /share/eisenlab/casett/database/kaiju_nr_euk/names.dmp -i Crab_gene_calls.kaiju.out.updated.euk -o Crab_gene_calls.kaiju.names.out.updated.samelevels_v2.euk -r superkingdom,phylum,class,order,family,genus,species

module load anvio/6.2

source activate anvio-6.2

anvi-import-taxonomy-for-genes -i Crab_gene_calls.kaiju.names.out.updated.samelevels_v2.euk -c Crab_Euk.db -p kaiju --just-do-it
