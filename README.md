#TCM-Macrophages

Analysis code for: "Tumor-conditioned media from canine cancer cell lines induces cytokine secretion by macrophages"

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18499208.svg)](https://doi.org/10.5281/zenodo.18499208)

## Data Availability

- **Analysis code:** This GitHub repository, archived on Zenodo (DOI: 10.5281/zenodo.18499208)
- **Raw cytokine data:** Dryad repository (DOI: to be added)
- **RNA-seq data:** Normalized RNA-seq expression data (TPM and RUVg-normalized counts) for the FACC canine cancer cell line panel were previously generated and are described in Farrell et al. (2024) PLoS ONE 19(5): e0303470. https://doi.org/10.1371/journal.pone.0303470
- **Whole-exome sequencing data:** Previously published in Das et al. (2019) Mol Cancer Ther 18(8):1460-1471. https://doi.org/10.1158/1535-7163.MCT-18-1346
- **CCL3 in histiocytic sarcoma:** NCBI BioProject PRJDB11462, PRJEB36828, PRJDB17594

## Repository Structure

### 01_cytokine_analysis/
Cytokine secretion analysis from macrophages treated with tumor-conditioned media.

| Figure/Table | Script |
|--------------|--------|
| Fig 2 (Modified z-scores) | `modz_plots.R`, `get modified z scores from raw data.R` |
| S1 Fig (Donor variability) | `Plos_one_figs/S1 Fig for PLOS one.R` |
| S2 Fig (Raw values) | `Plos_one_figs/S2 Fig for plos one/raw_for_plos_one_plot.R` |
| S3 Fig (Histology) | `Plos_one_figs/s3 fig for plos one.R` |
| S1 Table (ICCs) | `donor variability- ICCs and stats.R`, `updated_ICCs_for_PLOS.R` |
| S2 Table (Cytokine correlations) | `correlations between cytokines with heatmap.R` |
| S3 Table (Cell count effects) | `impact of cell count on cytokine levels.R` |
| S4 Table (Histology effects) | `impact of histo on cytokine levels.R` |

**CCL3/ (Fig 7):**
- `ccl3_stim_data_inverse__log_analysis.R` - Dose-response analysis
- `describing ccl3 levels.R` - CCL3 protein/mRNA levels
- `plot_ccl3_stim.R` - Figure generation

**VEGF_exos/ (Fig 6):**
- `exo_conditions_vegf.R` - Exosome fractionation analysis
- `mvb12a_vegf_analysis_without_CIN.R` - MVB12A validation
- `analyze ratios of high vs low mvb12a.R` - Exosomal vs free VEGF

### 02_RNAseq_analysis/
RNA-seq correlation and differential expression analysis.

**RUV_TPM code/ (Fig 3, S5 Table):**
- `*_updated_with_both_RUVg_TPM.R` - Spearman correlations for each cytokine
- `Make vegf_TPM_plot/vegf_tpm_spearmans_plots.R` - Fig 3 generation

**Pairwise code/new, good pairwise results/ (Fig 4, S6 Table):**
- `*_pairwise_final.R` - DESeq2 analysis for each cytokine
- `PLos one volcano/plos one volcano plots.R` - Volcano plot generation
- `comparing_top&bottom_quartiles_stats.R` - S6 Table statistics

**GSEA/ (Fig 5):**
- `*_gsea.R` - GSEA for each cytokine
- `gsea_swimline_for_plos.R` - Swimlane plot

## Requirements

R (≥4.0). Required packages are listed at the top of each script.

## License

MIT License
```