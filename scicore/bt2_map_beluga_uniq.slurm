#!/bin/bash
#SBATCH --job-name=bt2_map_beluga_uniq
#SBATCH --partition=scicore
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --array=0-14
#SBATCH --output=/scicore/home/holmp/carill0001/Beluga_MAGs/logs/bt2_map_beluga_uniq_%A_%a.out
#SBATCH --error=/scicore/home/holmp/carill0001/Beluga_MAGs/logs/bt2_map_beluga_uniq_%A_%a.err

# 1. Load conda and activate mapping env
module load Miniconda3
eval "$(/scicore/soft/easybuild/apps/Miniconda3/24.7.1-0/bin/conda shell.bash hook)"
conda activate mapping-env

# 2. Paths and settings
BASE_DIR="/scicore/home/holmp/carill0001/Beluga_MAGs"
MAP_DIR="${BASE_DIR}/mapping"
LOG_DIR="${BASE_DIR}/logs"
SAMPLE_LIST="${BASE_DIR}/sample_ids.txt"

# Where the trimmed reads live (from the tar you unpacked)
READ_BASE="/scicore/home/holmp/carill0001/release_AA"

# Bowtie2 index prefix (WITHOUT .bt2 suffixes)
REF_PREFIX="${BASE_DIR}/beluga_MAGs_all_uniq"

mkdir -p "${MAP_DIR}" "${LOG_DIR}"

# 3. Get this task's sample ID
LINE_NUM=$((SLURM_ARRAY_TASK_ID + 1))
SAMPLE_ID=$(sed -n "${LINE_NUM}p" "${SAMPLE_LIST}")

echo "$(date) Array task ${SLURM_ARRAY_TASK_ID} (line ${LINE_NUM}): SAMPLE_ID=${SAMPLE_ID}"

if [ -z "${SAMPLE_ID}" ]; then
    echo "No SAMPLE_ID found for this array index. Exiting." >&2
    exit 1
fi

# 4. Define read files
R1="${READ_BASE}/${SAMPLE_ID}/trimmed_reads/read1_paired.fq.gz"
R2="${READ_BASE}/${SAMPLE_ID}/trimmed_reads/read2_paired.fq.gz"

if [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then
    echo "Missing reads for ${SAMPLE_ID}: ${R1} or ${R2} not found." >&2
    exit 1
fi

# 5. Run Bowtie2 mapping
cd "${MAP_DIR}" || exit 1

SAM_OUT="${MAP_DIR}/${SAMPLE_ID}.sam"

echo "$(date) Mapping ${SAMPLE_ID} against ${REF_PREFIX}"
bowtie2 \
    -x "${REF_PREFIX}" \
    -1 "${R1}" \
    -2 "${R2}" \
    -S "${SAM_OUT}" \
    -p "${SLURM_CPUS_PER_TASK}" \
    --no-unal

BT2_EXIT=$?

if [ ${BT2_EXIT} -ne 0 ]; then
    echo "ERROR: bowtie2 mapping failed for ${SAMPLE_ID} with exit code ${BT2_EXIT}" >&2
    exit ${BT2_EXIT}
fi

echo "$(date) Finished mapping ${SAMPLE_ID}"
