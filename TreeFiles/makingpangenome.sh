#!/bin/bash
#SBATCH -n 10 #Request 4 tasks (cores)
#SBATCH -N 1 #Request 1 node
#SBATCH -t 2-2:00 #Request runtime of 30 minutes
#SBATCH -C centos7 #Request only Centos7 nodes
#SBATCH -p sched_mit_chisholm #Run on sched_engaging_default partition
#SBATCH --mem-per-cpu=10000 #Request 4G of memory per CPU
#SBATCH -o output_%j.txt #redirect output to output_JOBID.txt
#SBATCH -e error_%j.txt #redirect errors to error_JOBID.txt
#SBATCH --mail-type=BEGIN,END #Mail when job starts and ends
#SBATCH --mail-user=demers@mit.edu #email recipient

#anvi-gen-genomes-storage -e external-genomes-v2.txt -o Alteromonas-v2-GENOMES.db

#anvi-pan-genome -g Alteromonas-v2-GENOMES.db --project-name Alteromonas_Pangenome_v2 --output-dir Alteromonas_Pangenome_v2 --num-threads 10 --use-ncbi-blast --I-know-this-is-not-a-good-idea

#anvi-script-compute-bayesian-pan-core -p Alteromonas_Pangenome_v2/Alteromonas_Pangenome_v2-PAN.db -g Alteromonas-v2-GENOMES.db --store-in-db

#anvi-script-add-default-collection -p Alteromonas_Pangenome_v2/Alteromonas_Pangenome_v2-PAN.db -g Alteromonas-v2-GENOMES.db
anvi-summarize -g Alteromonas-v2-GENOMES.db -p Alteromonas_Pangenome_v2/Alteromonas_Pangenome_v2-PAN.db -C DEFAULT --force-overwrite -o summary-v2

