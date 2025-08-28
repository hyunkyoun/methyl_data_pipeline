# Methylation Data Analysis Pipeline - README

## Overview

This R workflow provides a comprehensive pipeline for analyzing DNA methylation data from Illumina mouse methylation arrays. The pipeline processes raw IDAT files through quality control, batch correction, differential methylation analysis, and motif enrichment preparation.

## Pipeline Flow
```
Raw IDAT Files → QC Filtering → ComBat Normalization → PCA Analysis → 
DMP Analysis → Visualization → FASTA Extraction → AME-Ready Output
```

---

## Detailed Step-by-Step Explanation

### Setup and Library Loading
```r
.libPaths("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/packages")
setwd("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline")
```

**Purpose**: Configures the R environment with custom package library path and working directory.

**Required Libraries**:
- **sesame**: Modern methylation array analysis
- **minfi**: Traditional methylation array processing
- **IlluminaMouseMethylation packages**: Mouse-specific annotation and manifest
- **limma**: Linear modeling for differential analysis
- **sva**: Surrogate variable analysis and ComBat batch correction
- **Visualization**: ggplot2, pheatmap, RColorBrewer, gridExtra
- **Genomics**: GenomicRanges, GenomicFeatures, BSgenome.Mmusculus.UCSC.mm10
- **Utilities**: Biostrings, parallel processing packages

---

## Step 1: Read IDAT Files

```r
betas <- openSesame(idir, BPPARAM = MulticoreParam(workers = cores), collapseToPfx = TRUE)
detP <- openSesame(idir, func = pOOBAH, return.pval = TRUE, BPPARAM = MulticoreParam(workers = cores))
```

**What it does**:
- Reads raw .idat files from the specified directory
- Extracts beta values (methylation levels, 0-1 scale)
- Calculates detection p-values using pOOBAH method
- Uses parallel processing for efficiency

**Output**:
- `betas`: Matrix of methylation beta values
- `detP`: Detection p-values for quality assessment
- `betas_raw_idat.csv`: Saved raw data

**Key Points**:
- Beta values represent the proportion of methylated cytosines
- Detection p-values indicate probe reliability
- Parallel processing speeds up large datasets

---

## Step 2: Per-Run Quality Control Filtering

```r
qc_filter_run <- function(betas_run, detP_run, targets_run, run_id) {
  # Remove samples with high mean detection p-values
  # Remove probes with high detection p-values
  # Remove incomplete cases and sex chromosomes
}
```

**Quality Control Filters Applied**:

1. **Sample-level filtering**: Removes samples with mean detection p-value > 0.05
2. **Probe-level filtering**: Keeps only probes with detection p-value < 0.05 in ALL samples
3. **Complete case filtering**: Removes probes with any missing values
4. **CpG filtering**: Keeps only CpG probes (removes control probes)
5. **Sex chromosome removal**: Excludes X and Y chromosome probes

**Why per-run filtering?**:
- Different sequencing runs may have different quality characteristics
- Ensures consistent quality standards across batches
- Prevents removal of valid probes that may fail in only one run

**Output**:
- `betas_qc_filtered.csv`: Quality-controlled beta values
- Substantially reduced probe and sample counts

---

## Step 3: Sample Sheet Matching

```r
targets$Basename <- paste0(targets$Beadchip, "_", targets$Position)
targets$Samples <- paste0(targets$Sample.ID, "_", targets$Experiment)
```

**Purpose**:
- Links methylation data to sample metadata
- Creates standardized sample identifiers
- Ensures proper sample tracking through the pipeline

**Key Variables Created**:
- `Basename`: Technical identifier matching IDAT files
- `Samples`: Biological sample identifier for analysis
- `BATCH`: Batch information for correction

---

## Step 4: ComBat Batch Correction

```r
if (skip_combat) {
  betas_adjusted <- betas_combined
} else {
  mvals_combat <- ComBat(dat = mvals_combined, batch = batch, mod = mod)
  betas_adjusted <- 2^mvals_combat / (1 + 2^mvals_combat)
}
```

**What ComBat does**:
- Removes systematic differences between batches/runs
- Preserves biological variation while removing technical variation
- Uses empirical Bayes methods for batch effect removal

**Process**:
1. Convert beta values to M-values (log2 ratio scale)
2. Apply ComBat correction in M-value space
3. Convert back to beta values
4. Handle tissue-specific design matrices when applicable

**When it's skipped**: If only one batch is present

**Output**: `betas_combat_adjusted.csv`

---

## Step 5 & 6: PCA Analysis and Top Variable CpGs

```r
create_simple_pca <- function(betas_data, targets_data, ntop = 100, output_prefix) {
  # Select most variable CpGs
  # Standardize data
  # Perform PCA
  # Create visualizations
}
```

**Process**:
1. **Variable CpG Selection**: Identifies top 10,000 most variable CpGs
2. **Data Standardization**: Centers and scales the data
3. **PCA Computation**: Performs principal component analysis
4. **Visualization Creation**: Generates plots colored by:
   - Batch (technical variation)
   - Tissue type (biological variation)
   - Genotype (experimental variation)

**Output Files**:
- `COMBAT_ADJUSTED_PCA_Batch_PC1_PC2.pdf`
- `COMBAT_ADJUSTED_PCA_Tissue_PC1_PC2.pdf`
- `COMBAT_ADJUSTED_PCA_Genotype_PC1_PC2.pdf`
- `COMBAT_ADJUSTED_betas_top10000_cpgs.csv`

**Purpose**: 
- Quality assessment (batch effects, outliers)
- Dimensionality reduction
- Data exploration

---

## Step 7: Differential Methylation Analysis

### Contrast Definitions
The pipeline analyzes 16 different biological contrasts:

**Thymus Comparisons**:
- `dp10_vs_WT_thymus1`: dp10 vs wildtype in thymus tissue
- `dp16_vs_WT_thymus1`: dp16 vs wildtype in thymus tissue  
- `dp16_dp10_vs_WT_thymus1`: double knockout vs wildtype

**T-cell Comparisons**:
- Similar contrasts in T-cell populations

**RUNX1 Studies**:
- `runx1_vs_WT_thymus3`: RUNX1 overexpression studies
- `runx1_dp16_vs_WT_thymus3`: Combined RUNX1 and dp16

**Trisomy Study**:
- `Ts65Dn_vs_WT_thymocytes5`: Down syndrome model

### Statistical Analysis
```r
mod <- model.matrix(~ 0 + group, data = targets)
fit <- lmFit(betas, mod)
cont_matrix <- makeContrasts(contrasts = contrast_formula, levels = mod)
fit2 <- eBayes(fit_cont)
```

**Method**: Linear modeling with empirical Bayes moderation (limma)

**Significance Thresholds**:
- **P-value**: < 1e-4 (very stringent)
- **Effect size**: |logFC| > 0.1 (10% difference in beta values)
- **FDR**: < 1e-4 for final DMP calls

### Visualizations Generated

#### Volcano Plots
- X-axis: Effect size (logFC)
- Y-axis: -log10(p-value)
- Colors: Red (hypermethylated), Blue (hypomethylated), Gray (not significant)
- Files: `volcano_[contrast_name].pdf`

#### Z-score Heatmaps
- Rows: Significant DMPs (differentially methylated positions)
- Columns: Samples
- Colors: Z-scored beta values (standardized methylation levels)
- Clustering: Hierarchical clustering of DMPs
- Files: `heatmap_zscore_[contrast_name].pdf`

---

## Step 8: FASTA Sequence Extraction for Motif Analysis

```r
extract_dmp_sequences <- function(DMPs, contrast_name, window_size = 100) {
  # Add genomic coordinates
  # Separate GOM (gain) vs LOM (loss) of methylation
  # Create GRanges with ±100bp windows
  # Extract sequences from mouse genome
  # Write FASTA files
}
```

### Process Flow:

1. **Coordinate Mapping**: Links CpG probe IDs to genomic coordinates using mouse annotation
2. **Directional Separation**:
   - **GOM (Gain of Methylation)**: DMPs with logFC > 0 (hypermethylated)
   - **LOM (Loss of Methylation)**: DMPs with logFC < 0 (hypomethylated)
3. **Window Creation**: Extends each DMP position by ±100bp (200bp total)
4. **Overlap Resolution**: Merges overlapping regions to avoid redundancy
5. **Sequence Extraction**: Retrieves DNA sequences from mm10 mouse genome
6. **FASTA Generation**: Creates properly formatted FASTA files

### Output Structure:
```
FASTA_sequences/
├── dp10_vs_WT_thymus1_GOM.fasta          # Hypermethylated regions
├── dp10_vs_WT_thymus1_LOM.fasta          # Hypomethylated regions
├── dp10_vs_WT_thymus1_GOM_details.csv    # Detailed DMP information
├── dp10_vs_WT_thymus1_LOM_details.csv
├── [other contrasts...]
```

### Purpose for Motif Analysis:
- **GOM sequences**: May be enriched for repressor binding sites
- **LOM sequences**: May be enriched for activator binding sites
- **AME input**: Ready for Analysis of Motif Enrichment

---

## Output Files Summary

### Data Files:
- `betas_raw_idat.csv`: Raw methylation data
- `betas_qc_filtered.csv`: Quality-controlled data
- `betas_combat_adjusted.csv`: Batch-corrected data
- `COMBAT_ADJUSTED_betas_top10000_cpgs.csv`: Top variable CpGs
- `DM_[contrast].csv`: Differential methylation results (16 files)

### Visualization Files:
- `COMBAT_ADJUSTED_PCA_*.pdf`: PCA plots (3 files)
- `volcano_[contrast].pdf`: Volcano plots (16 files)
- `heatmap_zscore_[contrast].pdf`: Heatmaps (16 files)

### Motif Analysis Files:
- `FASTA_sequences/[contrast]_GOM.fasta`: Hypermethylated regions
- `FASTA_sequences/[contrast]_LOM.fasta`: Hypomethylated regions
- `FASTA_sequences/[contrast]_*_details.csv`: DMP annotations

---

## Next Steps: AME Analysis

### Recommended AME Settings:
```bash
ame --control [background_sequences] --evalue-report-threshold 1e-4 \
    [input_fasta] [motif_database]
```

**Parameters**:
- **E-value cutoff**: 1e-4 (stringent threshold)
- **Motif database**: JASPAR or HOCOMOCO mouse motifs
- **Background**: Consider using non-significant regions

### Analysis Strategy:
1. **Compare GOM vs LOM**: Different motif enrichments indicate direction-specific regulation
2. **Cross-contrast comparison**: Identify common vs unique regulatory patterns
3. **Focus on high-DMP contrasts**: More DMPs = more statistical power
4. **Validate findings**: Cross-reference with known biology of the experimental conditions

---

## Troubleshooting

### Common Issues:
1. **Memory limitations**: Reduce `ntop` parameter or use more RAM
2. **No significant DMPs**: Relax p-value or effect size thresholds
3. **Batch effects remaining**: Check PCA plots, consider additional covariates
4. **FASTA extraction failures**: Verify genome annotation compatibility

### Quality Checks:
1. **Sample clustering in PCA**: Should group by biology, not batch
2. **Detection p-value distributions**: Should be mostly < 0.05
3. **Volcano plot patterns**: Should show clear significant hits
4. **DMP counts**: Should be reasonable (hundreds to thousands, not tens of thousands)

---

## Dependencies and Requirements

### System Requirements:
- R ≥ 4.0.0
- Minimum 16GB RAM (32GB+ recommended for large datasets)
- 50GB+ free disk space

### Critical Files Needed:
- `data/samplesheet.csv`: Sample metadata
- `idat/`: Directory containing .idat files
- `bmiq/DoBMIQ.R`: BMIQ normalization script (if used)

### Mouse Genome:
- BSgenome.Mmusculus.UCSC.mm10 package
- IlluminaMouseMethylationanno.12.v1.mm10 annotation

---

This pipeline provides a complete, publication-ready analysis of mouse methylation array data with integrated motif analysis preparation. The modular design allows for easy modification and extension based on specific research needs.