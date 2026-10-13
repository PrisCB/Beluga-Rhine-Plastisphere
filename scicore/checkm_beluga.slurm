#!/bin/bash
#SBATCH --job-name=checkm_beluga
#SBATCH --partition=scicore
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=checkm_beluga.out
#SBATCH --error=checkm_beluga.err

# Load conda
module load Miniconda3

# Activate environment that contains CheckM
conda activate dram-env

# Directory containing your MAGs
MAG_DIR=/scicore/home/holmp/carill0001/Beluga_MAGs/all_MAGs

echo "Running CheckM lineage_wf on directory: ${MAG_DIR}"
date

checkm lineage_wf \
    "${MAG_DIR}" \
    checkm_out \
    -x fna \
    -t 16

echo "CheckM finished"
date

