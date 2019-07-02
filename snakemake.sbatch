#!/bin/bash

# sbatch submission script to run main snakemake process. It then submits
# individual jobs from the compute node.

#SBATCH --job-name=snakemake
#SBATCH --output=snakemake.sbatch.log
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=4G
#SBATCH --tasks-per-node=4

echo "hello world3"

#activate environment
# source activate my_PrimerDesign_env

snakemake --jobs 200 -p --ri --cluster-config cluster-config.json --cluster "sbatch --partition={cluster.partition} --job-name={cluster.name} --output=/dev/null --job-name={cluster.name} --nodes={cluster.n} --mem={cluster.mem}"
