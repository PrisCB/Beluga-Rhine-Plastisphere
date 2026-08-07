# Beluga-Rhine-Plastisphere

Code and analysis scripts for a study of plastisphere and free-living (water) microbial communities along the River Rhine, based on samples from the 2019 Greenpeace *Beluga II* expedition.

## Study overview

Samples were collected at 14 stations spanning ~820 km of the Rhine, from Rotterdam, Netherlands to Basel, Switzerland, during the March–April 2019 *Beluga II* expedition. The study combines:

- ATR-FTIR polymer characterisation of microplastics
- 16S/18S rRNA amplicon sequence variants (ASVs), extracted *in silico* from shotgun metagenomic reads via Cascabel
- 120 quality-filtered metagenome-assembled genomes (MAGs), recovered from co-assembly and binning of the same shotgun data
to characterise how plastic surfaces restructure microbial community composition and function — with a particular focus on corrinoid (vitamin B12) biosynthesis — relative to the surrounding river water.

## Repository structure

```
├── figures/
│   ├── figure1_map.R                       # Station map (Fig. 1A)
│   ├── figure1_panelC_standalone.R         # Polymer composition panel (Fig. 1C)
│   ├── figure1_composed.R                  # Composed multi-panel Figure 1
│   ├── figure2_three_panel.R               # Community composition/NMDS (Fig. 2)
│   ├── figure_functional_heatmap.R         # Functional heatmap (general)
│   └── figure3_functional_heatmap_corrinoid.R  # Corrinoid/B12 functional heatmap (Fig. 3)
├── analysis/
│   ├── beluga_16S_ASV_analysis.R           # 16S community analysis (water + plastisphere)
│   ├── beluga_18S_eukaryotic_analysis.R    # 18S eukaryotic/diatom community analysis
│   ├── beluga_MAG_functional_analysis.R    # MAG-level functional (DRAM) analysis
│   └── Beluga_reload_checklist.R           # Session/object reload checklist for the R workspace
├── coverage/
│   ├── mag_cov_byMAG_v2.py                 # Final MAG-level coverage aggregation script (used)
│   ├── mag_cov_byMAG.py                    # Earlier version
│   ├── mag_cov_byMAG_simple.py             # Simplified version
│   └── build_mag_coverage2.py              # Coverage matrix construction
└── README.md
```

## Data availability

Raw shotgun metagenomic sequence data and BioSample metadata are deposited at NCBI under **BioProject [PRJNA1506697](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1506697)**.

- 14 water (planktonic) BioSamples: SAMN62142635–SAMN62143646
- 11 plastisphere BioSamples (7 ATR-FTIR-characterised stations plus additional stations with plastisphere sequencing)
- 120 metagenome-assembled genome (MAG) BioSamples: SAMN62216404–SAMN62216534, classified using GTDB-Tk (database release 226)
- Corresponding genome assemblies: GenBank accessions *in process*

Supplementary tables referenced in the manuscript (station metadata, nutrient data, polymer composition) are provided in the manuscript's Supplementary Material.

## Methods summary

Full methodological detail is provided in the manuscript (Carrillo-Barragán et al., *ISME Communications*, in preparation). Briefly:

- **Assembly:** MEGAHIT v1.2.9
- **Binning:** MetaBAT2, MaxBin2, CONCOCT; dereplicated with DAS Tool
- **Genome quality:** CheckM v1.2.4 (≥70% completeness, ≤10% contamination)
- **Taxonomy:** GTDB-Tk (Galaxy wrapper v2.6.1+galaxy0), GTDB reference database release 226
- **Coverage:** Bowtie2 v2.5.1 → samtools → custom Python aggregation (`mag_cov_byMAG_v2.py`)
- **Functional annotation:** DRAM v1.4.5

## Citation

If you use this code or data, please cite:

> Carrillo-Barragán P, et al. Riverine plastic substrates restructure microbial vitamin B12 metabolism and generate functional heterogeneity. *ISME Communications* (in preparation).

## Contact

Priscilla Carrillo-Barragán — University of Zurich, previously University of Basel
