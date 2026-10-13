#!/bin/bash
#SBATCH --job-name=bt2_map_beluga
#SBATCH --partition=scicore
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --array=0-14
#SBATCH --output=/scicore/home/holmp/carill0001/Beluga_MAGs/logs/bt2_map_%A_%a.out
#SBATCH --error=/scicore/home/holmp/carill0001/Beluga_MAGs/logs/bt2_map_%A_%a.err

echo "$(date) Starting array task ${SLURM_ARRAY_TASK_ID}"

# 1. Load conda and activate mapping environment
module load Miniconda3
eval "$( /scicore/soft/easybuild/apps/Miniconda3/24.7.1-0/bin/conda shell.bash hook )"
conda activate mapping-env

# 2. Define paths
BASE_DIR="/scicore/home/holmp/carill0001/Beluga_MAGs"
HOME_DIR="/scicore/home/holmp/carill0001"
REF_PREFIX="${BASE_DIR}/beluga_MAGs_all"   # Bowtie2 index prefix
SAMPLE_LIST="${BASE_DIR}/sample_ids.txt"
READS_BASE="${HOME_DIR}/release_AA"
MAP_DIR="${BASE_DIR}/mapping"

mkdir -p "${MAP_DIR}"

# 3. Select sample for this array task
SAMPLE_ID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${SAMPLE_LIST}")
echo "$(date) Task ${SLURM_ARRAY_TASK_ID} using SAMPLE_ID=${SAMPLE_ID}"

# 4. Define read files
R1="${READS_BASE}/${SAMPLE_ID}/trimmed_reads/read1_paired.fq.gz"
R2="${READS_BASE}/${SAMPLE_ID}/trimmed_reads/read2_paired.fq.gz"

echo "$(date) R1 = ${R1}"
echo "$(date) R2 = ${R2}"

# 5. Check that reads exist
if [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then
    echo "$(date) ERROR: One or both read files missing for ${SAMPLE_ID}"
    exit 1
fi

# 6. Define output files
SAM_OUT="${MAP_DIR}/${SAMPLE_ID}.sam"
BAM_OUT="${MAP_DIR}/${SAMPLE_ID}.bam"
BAM_SORTED="${MAP_DIR}/${SAMPLE_ID}.sorted.bam"

# Skip if sorted BAM already exists
if [ -f "${BAM_SORTED}" ]; then
    echo "$(date) ${SAMPLE_ID}: sorted BAM already exists, skipping."
    exit 0
fi

# 7. Run Bowtie2
echo "$(date) Running Bowtie2 for ${SAMPLE_ID}"
bowtie2 \
    -x "${REF_PREFIX}" \
    -1 "${R1}" \
    -2 "${R2}" \
    --threads ${SLURM_CPUS_PER_TASK} \
    -S "${SAM_OUT}"

if [ $? -ne 0 ]; then
    echo "$(date) ERROR: Bowtie2 failed for ${SAMPLE_ID}"
    exit 1
fi

# 8. Convert SAM to BAM and sort
echo "$(date) Converting SAM to BAM for ${SAMPLE_ID}"
samtools view -@ ${SLURM_CPUS_PER_TASK} -bS "${SAM_OUT}" > "${BAM_OUT}"

if [ $? -ne 0 ]; then
    echo "$(date) ERROR: samtools view failed for ${SAMPLE_ID}"
    exit 1
fi

echo "$(date) Sorting BAM for ${SAMPLE_ID}"
samtools sort -@ ${SLURM_CPUS_PER_TASK} -o "${BAM_SORTED}" "${BAM_OUT}"

if [ $? -ne 0 ]; then
    echo "$(date) ERROR: samtools sort failed for ${SAMPLE_ID}"
    exit 1
fi

# 9. Index BAM
echo "$(date) Indexing sorted BAM for ${SAMPLE_ID}"
samtools index "${BAM_SORTED}"

# 10. Clean up intermediate files
echo "$(date) Cleaning up SAM and unsorted BAM for ${SAMPLE_ID}"
rm -f "${SAM_OUT}" "${BAM_OUT}"

echo "$(date) Finished ${SAMPLE_ID}"
EOF
