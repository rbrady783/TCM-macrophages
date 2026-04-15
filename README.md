# TCM-Macrophages

Analysis code for: "Distinct tumor genomic signatures underlie canine macrophage polarization"

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18499208.svg)](https://doi.org/10.5281/zenodo.18499208)

## Data Availability

- **Analysis code and raw data:** This GitHub repository, archived on Zenodo (DOI: [10.5281/zenodo.18499208](https://doi.org/10.5281/zenodo.18499208))

## Repository Structure

The repository is organized into three top-level folders, with subfolders matching the figure or table each analysis supports.

### `01_cytokine_analysis/`

Cytokine secretion analysis from macrophages treated with tumor-conditioned media.

- **`S1Fig_cell_counts/`** — Cancer cell counts at TCM harvest (S1 Fig)
- **`Table1_ICCs/`** — Intraclass correlation coefficients for donor variability (Table 1)
- **`Fig1_modz_scores/`** — Modified z-score calculation and plotting (Fig 1)
- **`Fig2_S4Table_cytokine correlations/`** — Spearman correlations between cytokines + heatmap (Fig 2, S4 Table)
- **`S2Fig_raw_donor_values/`** — Raw cytokine values by donor (S2 Fig)
- **`S3Table_cell_count_histo_effects/`** — Effects of cell count and histology on cytokine levels (S3 Table A & B)
- **`Text_permutation_testing/`** — Permutation tests for histology and cell line effects
- **`Text_mutation_phosphorylation/`** — Mutation and phosphorylation status analysis
- **`Fig6_S7Table_VEGF_exosomes/`** — VEGF / exosome fractionation analyses (Fig 6, S7 Table), with subfolders for each panel:
  - `Fig6A_S7Table_TCM_fractionation/` — Whole vs depleted vs exosome-only fractions
  - `Fig6C_MVB12A_VEGF_stimulation/` — MVB12A high vs low validation (with/without CIN, with/without HSA)
  - `Fig6D_exosomal_VEGF_ratios/` — Exosomal vs free VEGF ratio analysis
  - `TCM_VEGF_correlation/` — TCM VEGF content vs macrophage VEGF secretion
- **`Fig7_S3Fig_S8Table_CCL3/`** — CCL3-related analyses (Fig 7, S3 Fig, S8 Table), with subfolders for:
  - `Fig7A_CCL3_levels/` — CCL3 protein/mRNA across cell lines
  - `Fig7B_S8Table_CCL3_dose_response/` — Recombinant CCL3 dose-response
  - `S3Fig_CCL3_knockdown/` — CCL3 siRNA knockdown validation in DH82

### `02_RNAseq_analysis/`

RNA-seq correlation and differential expression analyses.

- **`gene_names.csv`** — Shared lookup table mapping ENSCAF gene IDs to human-readable gene names. Used to interpret outputs in the subfolders below.
- **`Fig3_S5Table_VEGF_correlations/`** — Spearman correlations between RNA-seq expression and cytokine modified z-scores (Fig 3, S5 Table), with per-cytokine analyses and the Fig 3 plot generation
- **`Fig4_S6Table_DEGs/`** — Differential expression analysis between top/bottom modified z-score quartiles (Fig 4, S6 Table), with subfolders for:
  - `quartile_groups/` — Quartile group assignment and statistical comparisons
  - `per_cytokine_DESeq2/` — DESeq2 pipeline run separately for each cytokine
  - `Fig4_volcano_plots/` — Volcano plot generation
- **`Fig5_GSEA/`** — Gene set enrichment analysis (Fig 5), with subfolders for:
  - `per_cytokine_GSEA/` — GSEA against MSigDB Hallmark, C2, C6, and C7 collections per cytokine
  - `Fig5_swimlane_plot/` — Swimlane plot generation

### `03_prism_files/`

GraphPad Prism files containing source data and analyses for figures generated in Prism rather than R. Requires GraphPad Prism to open.

## Requirements

R (≥4.0). Required packages are listed at the top of each script.

## License

MIT License

## Other Data

- **RNA-seq data:** Normalized RNA-seq expression data (TPM and RUVg-normalized counts) for the FACC canine cancer cell line panel were previously generated and are described in Farrell et al. (2024) PLoS ONE 19(5): e0303470. https://doi.org/10.1371/journal.pone.0303470. All scripts and their output files are fully reproduced in this repository, so the analyses can be examined and verified in detail.
- **Whole-exome sequencing data:** Previously published in Das et al. (2019) Mol Cancer Ther 18(8):1460-1471. https://doi.org/10.1158/1535-7163.MCT-18-1346
- **CCL3 in histiocytic sarcoma:** NCBI BioProject PRJDB11462, PRJEB36828, PRJDB17594
