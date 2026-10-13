#!/bin/bash
#SBATCH --job-name=bt2_index_beluga
#SBATCH --partition=scicore
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/scicore/home/holmp/carill0001/Beluga_MAGs/logs/bt2_index_beluga.out
#SBATCH --error=/scicore/home/holmp/carill0001/Beluga_MAGs/logs/bt2_index_beluga.err

module load Miniconda3
eval "$(/scicore/soft/easybuild/apps/Miniconda3/24.7.1-0/bin/conda shell.bash hook)"
conda activate mapping-env

cd /scicore/home/holmp/carill0001/Beluga_MAGs

# 1) Concatenate all MAGs (already done, but keep the check)
if [ ! -f beluga_MAGs_all.fna ]; then
    echo "$(date) Concatenating MAGs into beluga_MAGs_all.fna"
    cat all_MAGs/*.fna > beluga_MAGs_all.fna
else
    echo "$(date) beluga_MAGs_all.fna already exists, skipping concatenation."
fi

# 2) Build Bowtie2 index if it does not exist
if [ ! -f beluga_MAGs_all.1.bt2 ]; then
    echo "$(date) Building Bowtie2 index"
    bowtie2-build beluga_MAGs_all.fna beluga_MAGs_all
    echo "$(date) Bowtie2 index finished"
else
    echo "$(date) Bowtie2 index already exists, nothing to do."
fi
