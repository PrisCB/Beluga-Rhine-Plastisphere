# Beluga Rhine Plastisphere

Analysis and figure-generation code accompanying:

**Carrillo-Barragán P. et al. (2026).** *Riverine plastic litter restructures microbial vitamin B12 metabolism and generates functional heterogeneity.* bioRxiv. https://doi.org/10.64898/2026.09.06.749157

## Study overview

The 2019 Greenpeace *Beluga II* expedition sampled 14 stations along approximately 820 km of the River Rhine, from Rotterdam to Basel (23 March–5 April 2019). The study integrates ATR-FTIR polymer characterisation, **independently sequenced 16S and 18S rRNA gene amplicons**, and shotgun metagenomics to compare plastic-associated and water-column microbial communities.

The analyses examine community composition, relationships with dissolved nutrient gradients, functional heterogeneity, and vitamin B12 (corrinoid) biosynthetic potential. Amplicon and shotgun libraries were prepared from the same DNA extracts but constitute separate sequencing datasets.

## Repository structure

```text
Beluga-Rhine-Plastisphere/
├── README.md
├── scripts/
│   ├── beluga_16S_ASV_analysis.R
│   ├── beluga_18S_ASV_analysis.R
│   ├── beluga_FunctionalAnalysis
│   ├── figure1_MapNutrientsMPs
│   ├── figure2_NMDS&envfit
│   ├── figure3_functional_heatmap_corrinoid.R
│   └── mag_coverage_filtering.R
├── data/
│   └── metadata/
│       ├── README.md
│       ├── mag_qc.csv
│       ├── mag_metadata_120.csv
│       ├── sample_metadata.csv
│       ├── sra_amplicon_accessions.tsv
│       └── Beluga_MIMARKS_MIMS_MIMAGS_03.09.2026.xlsx
└── scicore/
    ├── README.md
    ├── mag_coverage_byMAG_v2.tsv
    └── contig_to_MAG.tsv
```

The `scripts/` directory contains the downstream R analyses and figure-generation code. The `scicore/` documentation describes upstream MAG quality assessment, read mapping and coverage summarisation on the University of Basel sciCORE HPC cluster. The repository's available input data and intermediate objects determine which analyses can be executed directly.

## Sequencing and bioinformatic workflows

### 16S and 18S rRNA gene amplicons

16S (515F/926R) and 18S (515F/951R) rRNA gene amplicons were generated in technical triplicate from the DNA extracts also used for shotgun metagenomics. Libraries were sequenced on an Illumina NextSeq 2000. The workflow used FastQC, PEAR, Cascabel v3.0 and DADA2 for sequence processing and ASV inference. Taxonomic classification used SILVA v138.1 for bacterial/archaeal ASVs and PR2 v5 for eukaryotic ASVs. Contamination correction included microDecon.

Relevant scripts:

- `scripts/beluga_16S_ASV_analysis.R`
- `scripts/beluga_18S_ASV_analysis.R`

### Shotgun metagenomics and MAG reconstruction

Shotgun reads were processed using Cutadapt v3.1, fastp v0.20.1 and FastQC v0.11.9. De novo assembly used MEGAHIT v1.2.9, with assembly quality assessed using QUAST v4.6.3. MAGs were reconstructed using MetaBAT2 v2.12.1, MaxBin2 v2.2.4 and CONCOCT v1.1.0, with bin integration and dereplication using DAS Tool v1.1.0.

CheckM v1.2.4 assessed genome completeness and contamination. From 302 recovered MAGs, 120 were retained for downstream analyses following quality filtering (≥70% completeness and ≤10% contamination) and the study's abundance/detection criteria. Taxonomy was assigned using GTDB-Tk (Galaxy wrapper v2.6.1+galaxy0; GTDB release 226).

Reads were mapped to the 302 recovered MAGs with Bowtie2 v2.5.1; samtools and custom Python scripts were used to process alignments and summarise MAG coverage. Relative abundances were normalised across the full set of 302 MAGs before selecting the 120 MAGs for downstream analyses. The upstream sciCORE workflow is described in [`scicore/README.md`](scicore/README.md).

### Functional annotation

DRAM v1.4.5 was used to annotate protein-coding genes. Community-weighted functional profiles were calculated from MAG relative abundances and KEGG Orthology annotations. The analyses compared habitat-level functional composition and dispersion, including aerobic corrin ring synthesis and other corrinoid-related pathways.

Relevant script: `scripts/beluga_FunctionalAnalysis`.

## Figures and downstream analyses

| Manuscript figure | Script | Analysis |
|---|---|---|
| Figure 1 | `scripts/figure1_MapNutrientsMPs` | Sampling stations, dissolved nutrients and ATR-FTIR polymer composition |
| Figure 2 | `scripts/figure2_NMDS&envfit` | ASV community ordination and environmental-gradient associations |
| Figure 3 | `scripts/figure3_functional_heatmap_corrinoid.R` | Community-weighted functional profiles and corrinoid biosynthetic potential |

The figure scripts use derived tables and intermediate R objects from the corresponding upstream analyses. Examples include `atr_mps_composition.rds`, bootstrap outputs for environmental fitting, and functional-profile `.rds` files. Required objects and file paths are specified in the scripts.

## Data availability

Amplicon sequencing run accessions and associated BioSample records are documented in `data/metadata/sra_amplicon_accessions.tsv` under **NCBI BioProject [PRJNA1506697](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1506697)**. The 120 selected MAG assemblies were submitted under genome submission `SUB16390494`; validation and genome accession assignment are pending. Shotgun raw-read SRA accession assignment is not established by the included amplicon accession table.

## Citation

Carrillo-Barragán P. et al. (2026). *Riverine plastic litter restructures microbial vitamin B12 metabolism and generates functional heterogeneity.* bioRxiv. https://doi.org/10.64898/2026.09.06.749157


## Reproducible MAG abundance processing (archived input)

From the repository root, run `Rscript scripts/mag_coverage_filtering.R` (base R only). This imports `scicore/mag_coverage_byMAG_v2.tsv` and `data/metadata/mag_qc.csv`, normalises each sample over **all 302 MAGs**, and exports the archived 120-MAG subset and an eligibility report to `results/mag_processing/`. The 120-MAG identifiers are supplied by `data/metadata/mag_metadata_120.csv`. This script is a reconstruction of verified downstream calculations, not an original sciCORE script.

