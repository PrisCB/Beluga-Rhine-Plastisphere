#!/usr/bin/env python

import os
import csv
from collections import defaultdict

# Base paths
BASE_DIR = "/scicore/home/holmp/carill0001"
BELUGA_DIR = os.path.join(BASE_DIR, "Beluga_MAGs")

# Input files
LOOKUP_PATH = os.path.join(BELUGA_DIR, "lookup", "contig_to_MAG.tsv")
SAMPLE_ID_PATH = os.path.join(BELUGA_DIR, "sample_ids.txt")

# Output
OUT_DIR = os.path.join(BELUGA_DIR, "coverage")
os.makedirs(OUT_DIR, exist_ok=True)
OUT_PATH = os.path.join(OUT_DIR, "mag_coverage_byMAG_v2.tsv")

# ----------------------------------------------------------------------
# 1. Load contig -> MAG mapping
# ----------------------------------------------------------------------
contig_to_mag = {}

with open(LOOKUP_PATH, "r") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) != 2:
            continue
        contig, mag = parts
        contig_to_mag[contig] = mag

mag_ids = sorted(set(contig_to_mag.values()))
print(f"Loaded {len(contig_to_mag)} contig->MAG mappings")
print(f"Number of distinct MAGs in mapping: {len(mag_ids)}")

# ----------------------------------------------------------------------
# 2. Load sample IDs
# ----------------------------------------------------------------------
sample_ids = []
with open(SAMPLE_ID_PATH, "r") as f:
    for line in f:
        s = line.strip()
        if not s:
            continue
        sample_ids.append(s)

print(f"Found {len(sample_ids)} samples: {', '.join(sample_ids)}")

# ----------------------------------------------------------------------
# 3. Accumulate length-weighted coverage per MAG per sample
#    For each (MAG, sample) we keep:
#       cov_sum = sum(depth * contigLen)
#       len_sum = sum(contigLen)
#    and compute coverage = cov_sum / len_sum
# ----------------------------------------------------------------------
cov_sum = defaultdict(float)   # key: (mag, sample)
len_sum = defaultdict(float)   # key: (mag, sample)

for sample in sample_ids:
    depth_path = os.path.join(
        BASE_DIR,
        "release_AA",
        sample,
        "bwa_coverage",
        "CONTIGS_MEGAHIT_depth.txt"
    )
    print(f"Reading depth file for {sample}: {depth_path}")

    if not os.path.exists(depth_path):
        print(f"WARNING: depth file not found for sample {sample}, skipping.")
        continue

    with open(depth_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        # Expect columns:
        # contigName, contigLen, totalAvgDepth, bamDepth, bamVar
        for row in reader:
            if len(row) < 3:
                continue
            contig = row[0]
            try:
                contig_len = float(row[1])
                depth = float(row[2])
            except ValueError:
                continue

            mag = contig_to_mag.get(contig)
            if mag is None:
                continue  # contig not in any MAG

            key = (mag, sample)
            cov_sum[key] += depth * contig_len
            len_sum[key] += contig_len

# ----------------------------------------------------------------------
# 4. Build list of MAGs that actually have coverage
# ----------------------------------------------------------------------
mags_with_cov = sorted(
    {mag for (mag, sample) in cov_sum.keys() if len_sum[(mag, sample)] > 0.0}
)

print(f"Number of MAGs with at least one mapped contig: {len(mags_with_cov)}")
if mags_with_cov:
    print(f"DEBUG first MAG with coverage: {mags_with_cov[0]}")

# ----------------------------------------------------------------------
# 5. Write output table
# ----------------------------------------------------------------------
with open(OUT_PATH, "w", newline="") as out_f:
    writer = csv.writer(out_f, delimiter="\t")

    # Header
    header = ["MAG_ID"] + sample_ids
    writer.writerow(header)

    # Rows per MAG
    for mag in mags_with_cov:
        row = [mag]
        for sample in sample_ids:
            key = (mag, sample)
            if len_sum[key] > 0:
                cov = cov_sum[key] / len_sum[key]
            else:
                cov = 0.0
            row.append(f"{cov:.6f}")
        writer.writerow(row)

print(f"Wrote MAG coverage table: {OUT_PATH}")
