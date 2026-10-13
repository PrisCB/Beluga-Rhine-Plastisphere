#!/usr/bin/env python3
import os
from collections import defaultdict

BASE_DIR = "/scicore/home/holmp/carill0001"
PROJECT_DIR = os.path.join(BASE_DIR, "Beluga_MAGs")
RELEASE_DIR = os.path.join(BASE_DIR, "release_AA")

LOOKUP_PATH = os.path.join(PROJECT_DIR, "lookup", "contig_to_MAG.tsv")
SAMPLE_IDS_PATH = os.path.join(PROJECT_DIR, "sample_ids.txt")
COVERAGE_DIR = os.path.join(PROJECT_DIR, "coverage")
os.makedirs(COVERAGE_DIR, exist_ok=True)
OUT_PATH = os.path.join(COVERAGE_DIR, "mag_coverage_byMAG_simple.tsv")

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
print(f"Loaded {len(contig_to_mag)} contig->MAG mappings")

# For curiosity, how many unique MAG IDs?
mag_ids = sorted(set(contig_to_mag.values()))
print(f"Unique MAG IDs in lookup: {len(mag_ids)}")

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
print(f"Found {len(samples)} samples: {', '.join(samples)}")

# ---- 3. Data structures to accumulate coverage ----
# coverage * length per MAG and sample
mag_cov_sum = defaultdict(lambda: defaultdict(float))
mag_len_sum = defaultdict(lambda: defaultdict(float))

# global across all contigs per sample
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

            # global coverage for this sample
            global_cov_sum[sample] += depth * length
            global_len_sum[sample] += length

            # MAG-level only if contig in lookup
            mag = contig_to_mag.get(contig)
            if mag is None:
                continue
            mag_cov_sum[mag][sample] += depth * length
            mag_len_sum[mag][sample] += length

print(f"Number of MAGs with at least one mapped contig: {len(mag_cov_sum)}")

# For debug: print first 5 MAG IDs that actually have coverage
for i, mag in enumerate(sorted(mag_cov_sum.keys())):
    print(f"DEBUG first MAG with coverage {i+1}: {mag}")
    if i >= 4:
        break

# ---- 5. Write output table very explicitly ----
with open(OUT_PATH, "w") as out:
    # header
    out.write("MAG_ID\t" + "\t".join(samples) + "\n")

    # per-MAG rows
    for mag in sorted(mag_cov_sum.keys()):
        row = [mag]
        for sample in samples:
            len_sum = mag_len_sum[mag].get(sample, 0.0)
            if len_sum == 0.0:
                cov = 0.0
            else:
                cov = mag_cov_sum[mag][sample] / len_sum
            row.append(f"{cov:.6f}")
        out.write("\t".join(row) + "\n")

    # single global row
    row = ["all_MAGs_concat"]
    for sample in samples:
        glen = global_len_sum.get(sample, 0.0)
        if glen == 0.0:
            gcov = 0.0
        else:
            gcov = global_cov_sum[sample] / glen
        row.append(f"{gcov:.6f}")
    out.write("\t".join(row) + "\n")

print(f"Wrote MAG coverage table: {OUT_PATH}")
