# sciCORE: MAG quality assessment and coverage processing

This directory documents upstream metagenomic processing performed on the University of Basel sciCORE HPC cluster for the Beluga Rhine Plastisphere study. Downstream R analyses and manuscript figures are documented in the [main README](../README.md).

## Input samples

The MAG coverage workflow used 14 metagenomic datasets representing seven water–plastic station pairs:

```text
GP01STb_0323  GP01T2_0323
GP06STa_0328  GP06T2_0328
GP07STb_0331  GP07T1_0331
GP08STa_0331  GP08T3_0331
GP09STb_0401  GP09T1_0401
GP13STa_0404  GP13T2_0404
GP14STb_0405  GP14T2_0405
```

`ST` denotes the water fraction; `T` denotes the plastic-associated fraction.

## Genome quality

CheckM quality statistics were exported to:

```text
checkm_beluga/checkm_qa.tsv
```

Downstream MAG selection used completeness ≥70% and contamination ≤10%, together with abundance and detection filters, yielding 120 MAGs from the original set of 302.

## Read mapping and coverage

The sciCORE working directory included MAG FASTA sequences, Bowtie2 indexes and mapping outputs, together with Python and SLURM scripts for indexing, mapping and summarisation.

Relevant paths in the HPC project were:

```text
all_MAGs/                         # Individual MAG FASTA files
lookup/contig_to_MAG.tsv         # Contig-to-MAG association table
coverage/mag_coverage_byMAG_v2.tsv
scripts/mag_cov_byMAG_v2.py
checkm_beluga/checkm_qa.tsv
```

The coverage summarisation script was executed as:

```bash
cd /scicore/home/holmp/carill0001/Beluga_MAGs
python scripts/mag_cov_byMAG_v2.py
```

Its output contained 302 MAG rows and 14 sample columns. Coverage calculations and their reference-sequence compatibility should be evaluated against the Bowtie2 alignment outputs when reproducing the analysis from raw reads.

## Downstream integration

MAG quality metrics, MAG relative abundances and DRAM annotations are integrated in R for community-weighted functional analyses. Relative abundance normalisation across 302 recovered MAGs precedes selection of the 120 quality-filtered MAGs.

For the sequencing methods, software versions, statistical procedures and accession details, see the [associated manuscript](https://doi.org/10.64898/2026.09.06.749157).
