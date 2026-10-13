#!/bin/bash
#SBATCH --job-name=checkm_beluga_red
#SBATCH --partition=scicore
#SBATCH --qos=6hours
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=06:00:00
#SBATCH --output=checkm_beluga_reduced.out
#SBATCH --error=checkm_beluga_reduced.err

echo "Running CheckM lineage_wf (reduced tree) on directory: /scicore/home/holmp/carill0001/Beluga_MAGs/all_MAGs"
date

# Load conda so that 'conda activate' works in the batch job
source ~/.bashrc

# Activate the environment where CheckM is installed
conda activate checkm-env

# Run CheckM lineage workflow on all MAGs
checkm lineage_wf \
    -x fna \
    -t 8 \
    --reduced_tree \
    /scicore/home/holmp/carill0001/Beluga_MAGs/all_MAGs \
    /scicore/home/holmp/carill0001/Beluga_MAGs/checkm_beluga

echo "CheckM finished"
date
