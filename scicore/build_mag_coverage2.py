#!/usr/bin/env python3
import os
import csv
from collections import defaultdict

BASE_DIR = "/scicore/home/holmp/carill0001"
PROJECT_DIR = os.path.join(BASE_DIR, "Beluga_MAGs")
RELEASE_DIR = os.path.join(BASE_DIR, "release_AA")

LOOKUP_PATH = os.path.join(PROJECT_DIR, "lookup", "contig_to_MAG.tsv")
SAMPLE_IDS_PATH = os.path.join(PROJECT_DIR, "sample_ids.txt")
COVERAGE_DIR = os.path.join(PROJECT_DIR, "coverage")
os.makedirs(COVERAGE_DIR, exist_ok=True)
OUT_PATH = os.path.join(COVERAGE_DIR, "mag_coverage.tsv")

# ---- 1. Read contig -> MAG mapping ----
contig_to_mag = {}
with open(LOOKUP_PATH, "r") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        contig, mag = parts[0], parts[1]
        contig_to_mag[contig] = mag

if not contig_to_mag:
    raise RuntimeError(f"No contig->MAG mappings found in {LOOKUP_PATH}")

# ---- 2. Read sample IDs ----
samples = []
with open(SAMPLE_IDS_PATH, "r") as f:
    for line in f:
        s = line.strip()
        if not s:
            continue
        samples.append(s)

if not samples:
    raise RuntimeError(f"No sample IDs found in {SAMPLE_IDS_PATH}")

# ---- 3. Data structures to accumulate coverage ----
mag_cov_sum = defaultdict(lambda: defaultdict(float))  # depth*len
mag_len_sum = defaultdict(lambda: defaultdict(float))  # len

global_cov_sum = defaultdict(float)
global_len_sum = defaultdict(float)

# ---- 4. Process each sample's depth file ----
for sample in samples:
    depth_file = os.path.join(
        RELEASE_DIR,
        sample,
        "bwa_coverage",
        "CONTIGS_MEGAHIT_depth.txt"
    )
    print(f"Reading depth file for {sample}: {depth_file}")

    if not os.path.exists(depth_file):
        raise FileNotFoundError(f"Depth file not found for sample {sample}: {depth_file}")

    with open(depth_file, "r") as f:
        header = f.readline().rstrip("\n").split("\t")

        try:
            contig_idx = header.index("contigName")
            len_idx = header.index("contigLen")
            depth_idx = header.index("totalAvgDepth")
        except ValueError:
            contig_idx, len_idx, depth_idx = 0, 1, 2

        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) <= depth_idx:
                continue

            contig = parts[contig_idx]
            try:
                length = float(parts[len_idx])
                depth = float(parts[depth_idx])
            except ValueError:
                continue

            # Global coverage
            global_cov_sum[sample] += depth * length
            global_len_sum[sample] += length

            # MAG-level coverage (only if contig is in a MAG)
            mag = contig_to_mag.get(contig)
            if mag is None:
                continue

            mag_cov_sum[mag][sample] += depth * length
            mag_len_sum[mag][sample] += length

# ---- 5. Write output table ----
with open(OUT_PATH, "w", newline="") as out_f:
    writer = csv.writer(out_f, delimiter="\t")

    header = ["MAG_ID"] + samples
    writer.writerow(header)

    # Per-MAG rows
    for mag in sorted(mag_cov_sum.keys()):
        row = [mag]
        for sample in samples:
            len_sum = mag_len_sum[mag].get(sample, 0.0)
            if len_sum == 0.0:
                cov = 0.0
            else:
                cov = mag_cov_sum[mag][sample] / len_sum
            row.append(f"{cov:.6f}")
        writer.writerow(row)

    # Global summary (all MAGs)
    if global_len_sum:
        row = ["all_MAGs_concat"]
        for sample in samples:
            glen = global_len_sum.get(sample, 0.0)
            if glen == 0.0:
                gcov = 0.0
            else:
                gcov = global_cov_sum[sample] / glen
            row.append(f"{gcov:.6f}")
        writer.writerow(row)

print(f"Wrote MAG coverage table: {OUT_PATH}")
