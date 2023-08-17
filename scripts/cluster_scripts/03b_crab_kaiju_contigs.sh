#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH --partition=production
#SBATCH -D /share/eisenlab/casett/crabs_eukrep/
#SBATCH -e crab_kaiju_stderr.txt
#SBATCH -o crab_kaiju_stderr.txt
#SBATCH -J crab_kaiju
#SBATCH --mem=96G #memory in Gb
#SBATCH --time=12:00:00 #time in hours:min:sec

/share/eisenlab/gjospin/software/kaiju/src/kaiju -z 24 -t /share/eisenlab/casett/database/kaiju_nr_euk/nodes.dmp -f /share/eisenlab/casett/database/kaiju_nr_euk/kaiju_db_nr_euk.fmi -i Crab_euk_contigs.fa -o kaiju.out

/share/eisenlab/gjospin/software/kaiju/src/kaijuReport -t /share/eisenlab/casett/database/kaiju_nr_euk/nodes.dmp -n /share/eisenlab/casett/database/kaiju_nr_euk/names.dmp -i kaiju.out -r genus -p -o kaiju.out.summary

/share/eisenlab/gjospin/software/kaiju/src/addTaxonNames -t /share/eisenlab/casett/database/kaiju_nr_euk/nodes.dmp -n /share/eisenlab/casett/database/kaiju_nr_euk/names.dmp -i kaiju.out -o kaiju.names.out.samelevels -r phylum,class,order,family,genus,species

/share/eisenlab/gjospin/software/kaiju/src/kaiju2krona -t /share/eisenlab/casett/database/kaiju_nr_euk/nodes.dmp -n /share/eisenlab/casett/database/kaiju_nr_euk/names.dmp -u -i kaiju.out -o kaiju.out.krona
