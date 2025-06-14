#!/bin/bash
#SBATCH --job-name=rrap_codon_mini
#SBATCH --output=logs/%x.%A.%a.out
#SBATCH --error=logs/%x.%A.%a.err
#SBATCH --time=14-00:00:00 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --export=ALL

nextflow run /hps/nobackup/rdf/metagenomics/service-team/users/timrozday/rrap --samplesheet /hps/nobackup/rdf/metagenomics/service-team/users/timrozday/rrap/assets/PRJEB51728_mini_mini_codon.csv -profile codon,internal
