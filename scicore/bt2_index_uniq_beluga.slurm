#!/bin/bash
#SBATCH --job-name=bt2_index_uniq
#SBATCH --partition=scicore
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --output=/scicore/home/holmp/carill0001/Beluga_MAGs/logs/bt2_index_uniq_%j.out
#SBATCH --error=/scicore/home/holmp/carill0001/Beluga_MAGs/logs/bt2_index_uniq_%j.err

# 1. Load conda and activate mapping-env
module load Miniconda3
eval "$(/scicore/soft/easybuild/apps/Miniconda3/24.7.1-0/bin/conda shell.bash hook)"
conda activate mapping-env

BASE_DIR="/scicore/home/holmp/carill0001/Beluga_MAGs"
cd "${BASE_DIR}" || exit 1

# 2. Build combined FASTA with unique headers if it does not yet exist
if [ ! -f "beluga_MAGs_all_uniq.fna" ]; then
    echo "$(date) Building beluga_MAGs_all_uniq.fna from all_MAGs/*.fna"

    > beluga_MAGs_all_uniq.fna

    for MAG in all_MAGs/*.fna; do
        MAG_BASE=$(basename "${MAG}" .fna)
        echo "$(date)  Processing ${MAG_BASE}"

        awk -v prefix="${MAG_BASE}" '
            /^>/ {
                sub(/^>/, "", $0);
                print ">" prefix "|" $0;
                next;
            }
            { print }
        ' "${MAG}" >> beluga_MAGs_all_uniq.fna
    done

    echo "$(date) Finished building beluga_MAGs_all_uniq.fna"
else
    echo "$(date) beluga_MAGs_all_uniq.fna already exists, skipping rebuild."
fi

# 3. Build Bowtie2 index for the unique FASTA
echo "$(date) Building Bowtie2 index beluga_MAGs_all_uniq"

bowtie2-build \
    -f \
    -p "${SLURM_CPUS_PER_TASK}" \
    beluga_MAGs_all_uniq.fna \
    beluga_MAGs_all_uniq

BT2_EXIT=$?

if [ ${BT2_EXIT} -ne 0 ]; then
    echo "ERROR: bowtie2-build failed with exit code ${BT2_EXIT}" >&2
    exit ${BT2_EXIT}
fi

echo "$(date) Finished Bowtie2 index for beluga_MAGs_all_uniq"
