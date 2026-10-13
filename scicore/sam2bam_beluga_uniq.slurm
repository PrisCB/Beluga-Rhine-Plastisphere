#!/bin/bash
#SBATCH --job-name=sam2bam_beluga
#SBATCH --partition=scicore
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --array=0-14
#SBATCH --output=/scicore/home/holmp/carill0001/Beluga_MAGs/logs/sam2bam_beluga_%A_%a.out
#SBATCH --error=/scicore/home/holmp/carill0001/Beluga_MAGs/logs/sam2bam_beluga_%A_%a.err

# 1. Load conda and activate mapping env
module load Miniconda3
eval "$(/scicore/soft/easybuild/apps/Miniconda3/24.7.1-0/bin/conda shell.bash hook)"
conda activate mapping-env

# 2. Define paths
BASE_DIR="/scicore/home/holmp/carill0001/Beluga_MAGs"
MAP_DIR="${BASE_DIR}/mapping"
SAMPLE_LIST="${BASE_DIR}/sample_ids.txt"
LOG_DIR="${BASE_DIR}/logs"

mkdir -p "${LOG_DIR}"

# 3. Get this task's sample ID
LINE_NUM=$((SLURM_ARRAY_TASK_ID + 1))
SAMPLE_ID=$(sed -n "${LINE_NUM}p" "${SAMPLE_LIST}")

echo "$(date) Array task ${SLURM_ARRAY_TASK_ID} (line ${LINE_NUM}): SAMPLE_ID=${SAMPLE_ID}"

if [ -z "${SAMPLE_ID}" ]; then
    echo "No SAMPLE_ID found for this array index. Exiting." >&2
    exit 1
fi

cd "${MAP_DIR}" || exit 1

SAM_FILE="${SAMPLE_ID}.sam"
BAM_FILE="${SAMPLE_ID}.bam"

# 4. Sanity checks
if [ ! -f "${SAM_FILE}" ]; then
    echo "SAM file not found: ${SAM_FILE}" >&2
    exit 1
fi

# If a tiny BAM exists from previous attempts, remove it so we overwrite cleanly
if [ -f "${BAM_FILE}" ]; then
    echo "Removing existing BAM file: ${BAM_FILE}"
    rm -f "${BAM_FILE}"
fi

echo "$(date) Converting ${SAM_FILE} -> ${BAM_FILE}"

# 5. Sort SAM -> BAM and index
samtools sort \
    -@ "${SLURM_CPUS_PER_TASK}" \
    -o "${BAM_FILE}" \
    "${SAM_FILE}"

if [ $? -ne 0 ]; then
    echo "samtools sort failed for ${SAM_FILE}" >&2
    exit 1
fi

samtools index "${BAM_FILE}"
if [ $? -ne 0 ]; then
    echo "samtools index failed for ${BAM_FILE}" >&2
    exit 1
fi

# 6. Remove SAM to save space
rm -f "${SAM_FILE}"

echo "$(date) Finished ${SAMPLE_ID}"
