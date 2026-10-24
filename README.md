# Methylation Data Analysis Pipeline - README

## Overview

This R workflow provides a comprehensive pipeline for analyzing DNA methylation data from Illumina mouse methylation arrays. The pipeline processes raw IDAT files through quality control, batch correction, differential methylation analysis, and motif enrichment preparation.

Each execution saves its outputs under a new `results/analysis_YYYY-MM-DD_HH-MM-SS/` folder. Output paths below are relative to that run folder; see [Output Organization](#output-organization) for the complete layout.

Brain samples are excluded before IDAT selection, QC, ComBat, and all downstream analyses. Experiment numbers are preserved, so Experiment 6 has no output folder. The current analysis uses 161 samples from Experiments 1, 2, 3, 4, 5, 7, and 8.

## Pipeline Flow
```
Raw IDAT Files → QC Filtering → ComBat Normalization → DMP Analysis →
PCA Analysis → Visualization → FASTA Extraction → AME-Ready Output
```

---

## Detailed Step-by-Step Explanation

### Setup and Library Loading

Set R's working directory to the `methyl_data_pipeline` project folder. Install the required packages once:

```r
source("install_packages.R")
```

The installer creates `packages/` inside the project and passes that folder explicitly to each CRAN, Bioconductor, and GitHub installation call. Packages already installed in that folder are skipped; packages installed elsewhere do not count as local installations. Newly installed dependencies also go into `packages/`, while existing dependencies may still be used from other R libraries. The built-in `parallel` package is supplied by R. The `packages/` folder is excluded from Git.

The workflow looks in the project library first whenever it runs:

```r
package_dir <- file.path(getwd(), "packages")
.libPaths(c(package_dir, .libPaths()))
source("workflow_new.R")
```

On a new computer, rerun the installer to obtain packages for that computer's R version and operating system. SeSAMe's reference data are separate from the installed packages; download them once with `sesameData::sesameDataCache()`, which uses SeSAMe's data cache.

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
idir <- "idats"
samplesheet_path <- file.path(idir, "samplesheet_169_samples.csv")
# targets_raw has already been read and filtered to exclude brain.
idat_prefixes <- searchIDATprefixes(idir)
idat_prefixes <- idat_prefixes[names(idat_prefixes) %in% targets_raw$Basename]
betas <- openSesame(idat_prefixes, BPPARAM = MulticoreParam(workers = cores), collapseToPfx = TRUE)
detP <- openSesame(idat_prefixes, func = pOOBAH, return.pval = TRUE, BPPARAM = MulticoreParam(workers = cores))
```

**What it does**:
- Reads only eligible non-brain IDAT pairs from the specified directory
- Extracts beta values (methylation levels, 0-1 scale)
- Calculates detection p-values using pOOBAH method
- Uses parallel processing for efficiency

**Output**:
- `betas`: Matrix of methylation beta values
- `detP`: Detection p-values for quality assessment
- `tables/betas/betas_raw_idat.csv`: Saved raw data

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
- `tables/betas/betas_qc_filtered.csv`: Quality-controlled beta values
- Substantially reduced probe and sample counts

---

## Step 3: Sample Sheet Matching

```r
targets <- targets_raw  # Reuse the samplesheet validated before QC.
matched_indices <- match(colnames(betas), targets$Basename)
targets <- targets[matched_indices, , drop = FALSE]
```

**Purpose**:
- Links methylation data to sample metadata
- Creates standardized sample identifiers
- Ensures proper sample tracking through the pipeline

The samplesheet is read once from `idats/samplesheet_169_samples.csv`, with UTF-8 BOM support. Matching uses its existing `Basename` text exactly, such as `205243950018_R01C01`. The chip number is not reconstructed from `Beadchip`, which may contain rounded values after spreadsheet editing.

Missing, blank, or duplicate `Basename` identifiers cause a clear error. Every eligible sample must have a complete IDAT pair before processing. Brain and unlisted IDATs are excluded before reading intensities. After QC, the same metadata is reordered to match the surviving beta columns before they are renamed to `Samples`.

The source samplesheet retains all 169 records as an inventory. Sample `12926_8` is corrected to `True.Genotype = DP10_DP16`, `Genotype = DP10_DP16`, `cross.type = Double`, and `group = double_tcells_dp10_dp16`. Its sample and IDAT identifiers stay the same. Experiment 8 now has 7 DP10, 5 DP16, 7 DP10_DP16, and 5 WT samples.

**Key Variables Used**:
- `Basename`: Full technical identifier supplied by the samplesheet and matching IDAT filenames
- `Samples`: Biological sample identifier for analysis
- `BATCH`: Batch information for correction

---

## Step 4: ComBat Batch Correction

```r
batch <- factor(matching_samples$BATCH)
group <- interaction(matching_samples$Experiment, matching_samples$Tissue.Cell.type,
                     matching_samples$True.Genotype, drop = TRUE)
skip_combat <- nlevels(batch) < 2

if (skip_combat) {
  betas_adjusted <- betas_combined
} else {
  mod <- if (nlevels(group) > 1) {
    model.matrix(~ group)
  } else {
    model.matrix(~ 1, data = matching_samples)
  }
  # The workflow checks biological metadata, confounding, and finite M-values first.
  mvals_combat <- ComBat(dat = mvals_combined, batch = batch, mod = mod)
  betas_adjusted <- plogis(mvals_combat * log(2))
}
```

**What ComBat does**:
- Adjusts estimated technology effects on the mean and variability of M-values using empirical Bayes.
- Retains the fitted biological group mean effects supplied in `mod`; this does not guarantee that all biological variation or raw group averages remain unchanged.
- Uses the two technology categories: `BATCH1` for run `1_2` (older technology), and `BATCH2` for runs `3_4`, `5`, and `6` (newer technology). It does not estimate separate effects for the three newer-technology run labels.

**Process**:
1. Require non-missing `Experiment`, `Tissue.Cell.type`, and `True.Genotype` after QC.
2. Protect each observed experiment/tissue/genotype mean in the ComBat design. Technical correction still uses the full QC matrix; experiment is protected biology, not a technical batch to remove. For one remaining group, use an intercept-only design.
3. Check that the combined group and technology design has full rank and leaves residual degrees of freedom. Stop if the effects cannot be separated; do not drop group protection to force correction.
4. Convert beta values to M-values (log2 ratio scale). Retain the existing policy of excluding CpGs with any non-finite M-value, including beta values exactly 0 or 1.
5. Apply ComBat in M-value space and reject non-finite results.
6. Convert back to beta values using a numerically stable inverse transformation.

**When it's skipped**: If only one technology is present; the input beta matrix is passed through unchanged.

Group labels must be correct. Full rank establishes that the model can be fitted, but does not prove that residual technology effects are absent or that correction is perfect. Whether additional ComBat adjustment is needed after SeSAMe should be assessed from the data.

**Output**: `tables/betas/betas_combat_adjusted.csv`; the aligned, retained sample metadata is saved as `metadata/samples_after_qc.csv`.

---

## Step 5: Differential Methylation Analysis

Comparisons are generated automatically from `Experiment`, `Tissue.Cell.type`, and `True.Genotype`. Every non-WT genotype is compared with WT animals from **the same experiment and tissue**. In Experiments 1 and 2, WT controls must have `True.Genotype = WT` and the exact `Genotype` label: `WT_DP10` for DP10, `WT_DP16` for DP16, or `WT_DP10_WT_DP16` for DP10_DP16. All Dp animals are selected by `True.Genotype` regardless of `Genotype` or `cross.type`; for example, `11157_1` remains DP10 despite its `DP10_WT_DP16` label and Double cross. The `cross.type` column does not determine comparison membership. Each arm needs at least two retained samples; otherwise the comparison is skipped with a warning, without substituting WT from another group.

After excluding brain, the current samplesheet supports 19 comparisons:

| Experiment | Tissue | Genotypes compared with WT | Same-cohort WT samples |
| --- | --- | --- | --- |
| 1 | thymus | DP10, DP16, DP10_DP16 | 6, 4, 5 respectively |
| 2 | Tcells | DP10, DP16, DP10_DP16 | 4, 4, 3 respectively |
| 3 | thymus | DP16, RUNX1-DEL, RUNX1-DEL_DP16 | 5 |
| 4 | thymocytes | DP16, RUNX1-DEL, RUNX1-DEL_DP16 | 4 |
| 5 | thymocytes | Ts65Dn | 4 |
| 7 | Tcells | DP10, DP16, DP10_DP16 | 10 |
| 8 | Tcells | DP10, DP16, DP10_DP16 | 5 |

Thus T cells have three separate DP10-vs-WT comparisons, three DP16-vs-WT comparisons, and three DP10_DP16-vs-WT comparisons. Sample `11432_4` is the only Experiment 4 thymus sample and has no matching WT controls; it is retained in metadata but skipped for DMP testing.

### Statistical Analysis

`R/dmp_analysis.R` fits each Experiment 1/2 comparison independently using its complete genotype arm and selected WT controls. The other experiments retain a separate limma model per experiment/tissue cohort. Each fit uses every retained CpG after ComBat, with equal weight per sample. Each contrast subtracts its selected WT mean from the genotype mean, and moderation uses only samples in that fit. The unique `Samples` identifier remains the analysis unit; repeated `Sample.ID` values are not deduplicated.

Experiment 1's Dp-versus-WT counts are 9 vs 6 (DP10), 10 vs 4 (DP16), and 4 vs 5 (DP10_DP16). Experiment 2's counts are 6 vs 4, 7 vs 4, and 3 vs 3. These results use `Experiment_1_sameparent/` and `Experiment_2_sameparent/`; earlier `Experiment_1/` and `Experiment_2/` results remain available. The comparison filenames stay the same. The `WT.Cross` field in DMP tables, summaries, and contrast weights records the control rule. PCA membership and preprocessing are unaffected by this WT-only selection rule.

Full result tables and significant subsets are stored once in `results_list` and reused by volcano plots, pairwise and combined RUNX heatmaps, and FASTA extraction. No separate pooled RUNX analysis is fitted.

**Shared significance thresholds**, set in `workflow_new.R`:

- **FDR**: `adj.P.Val < dmp_fdr_cutoff`, default `0.05`, adjusted by BH within each contrast.
- **Effect size**: `abs(diff) > dmp_delta_beta_cutoff`, default `0.1` (10 percentage points in beta).
- Tables, volcano colors, heatmap rows, and FASTA selection use the same rules.

`tables/summaries/DMP_summary.csv` records experiment, tissue, sample counts, and DMP counts. `DMP_contrast_weights.csv` records each comparison's two arms, sample counts, and +1/-1 coefficients. Complete tables use short names, for example `Experiment_2_sameparent/tables/dmp/DP10_vs_WT_Exp2_Tcells.csv`. Earlier run folders retain their original results. An incremental same-parent rerun also saves its summaries and exact comparison sample lists inside each new experiment folder, preserving the earlier run-wide summaries.

---

## Step 6: PCA Analysis and Top Variable CpGs

The workflow produces both pooled overviews across experiments and a separate PCA for each experiment/tissue cohort. Each PCA selects up to 10,000 CpGs with the highest positive variance within its own sample subset, centers and scales the data, and computes PC1/PC2.

The pooled overviews exclude brain, all Ts65Dn samples, and the four Experiment 5 WT controls (`11984_5`, `11985_5`, `11988_5`, `11989_5`). This leaves 153 samples: 75 Tcells, 60 thymus, and 18 thymocytes. The all-tissue overview has genotype-, batch-, and tissue-colored versions; each pooled tissue PCA has genotype- and batch-colored versions. They are saved under `All_Experiments/figures/pca/`, with selected betas, CpGs, scores, and variance under `All_Experiments/tables/pca/`. These plots reuse the saved ComBat-adjusted values and refit feature selection and PCA after the exclusions. The exclusions apply to pooled PCA only; Experiment 5 retains its separate PCA and DMP comparison.

Each cohort gets genotype- and batch-colored plots, with consistent colors across cohorts. Titles identify the experiment, tissue, sample count, selected CpG count, and variance explained. The independent PCA axes are not directly comparable between cohorts.

Experiments 1, 2, 3, 4, 5, 7, and 8 each have their own PCA outputs. In addition, thymocytes with `age == 17.5` get a separately fitted PCA with the exact title **Age17.5**, under that experiment's `figures/pca/Age17.5/{genotype,batch}/` and `tables/pca/Age17.5/` folders. Selection uses tissue and age metadata, not a point's position on an earlier graph.

Currently, Age17.5 contains the 18 Experiment 4 thymocytes: 4 WT, 5 DP16, 4 RUNX1-DEL, and 5 RUNX1-DEL_DP16. The four WT and four Ts65Dn animals from Experiment 5 are age 2 and are excluded from this additional graph. Their separate Experiment 5 PCA and DMP comparison remain available. Because all Experiment 4 thymocytes are age 17.5, its regular and Age17.5 PCAs use the same samples.

`R/pca_analysis.R` skips subsets with fewer than three samples or fewer than two variable CpGs, with a warning. The one-sample Experiment 4 thymus cohort is therefore skipped for its separate PCA; that sample remains in the pooled thymus and all-tissue overviews.

For example, Experiment 2 T-cell outputs are:

- `Experiment_2_sameparent/figures/pca/genotype/Genotypes_Exp2_Tcells.pdf`
- `Experiment_2_sameparent/figures/pca/batch/Genotypes_Exp2_Tcells.pdf`
- `Experiment_2_sameparent/tables/pca/betas/Genotypes_Exp2_Tcells.csv`: selected beta values.
- `Experiment_2_sameparent/tables/pca/cpgs/Genotypes_Exp2_Tcells.csv`: CpG IDs and variances.
- `Experiment_2_sameparent/tables/pca/scores/Genotypes_Exp2_Tcells.csv`: PC1/PC2 and sample metadata.
- `Experiment_2_sameparent/tables/pca/variance/Genotypes_Exp2_Tcells.csv`: variance explained by each PC.

The age-specific genotype figure is `Experiment_4/figures/pca/Age17.5/genotype/Genotypes_Exp4_thymocytes.pdf`.

---

## Step 7: Visualize DMP Results

Plot titles show the comparison, experiment, tissue, and sample counts. Heatmaps use the exact samples used for DMP testing: WT animals from that experiment/tissue on the left and experimental animals on the right. Black brackets above the columns show each actual `True.Genotype` and its sample count; individual sample labels and colored genotype strips are hidden.

#### Volcano Plots
- X-axis: Beta difference (experimental minus control)
- X-axis limits stay symmetric around zero, spanning at least -0.5 to +0.5 and expanding to include all observed differences and both effect-size cutoff lines.
- Y-axis: -log10(FDR)
- Colors: Red (hypermethylated), Blue (hypomethylated), Gray (not significant)
- Files: `Experiment_<N>/figures/volcano/<comparison>_Exp<N>_<tissue>.pdf`

#### Mean-centered Heatmaps
- Rows: Significant DMPs (differentially methylated positions)
- Columns: Samples
- Colors: Each beta value minus that DMP's mean across all samples in the comparison, with no division by standard deviation. Centering is across the combined experimental and WT samples, not separately within each group. Values remain in beta units; white is zero, red is above the row mean, and blue is below it.
- Direction labels: GOM CpGs (positive experimental-minus-WT DMP effect) appear at the top with a red bracket on the right; LOM CpGs (negative effect) appear at the bottom with a blue bracket. Matching colored CpG counts appear below the heatmap, including zero for an absent direction. The horizontal color scale shows mean-centered beta values.
- Clustering: Hierarchical clustering within each direction block when it contains at least two DMPs. This keeps each bracket aligned with the correct CpGs. Direction comes from the existing significant DMP calls, not the mean-centered tile colors.
- Layout code: `R/heatmap_plots.R` prepares the direction blocks and draws the brackets, counts, and color scale.
- Files: `Experiment_<N>/figures/heatmap/<comparison>_Exp<N>_<tissue>.pdf`
- Page size: 11 inches wide, with height growing from a 4-inch minimum to an 8.5-inch maximum (US Letter landscape). All DMP rows fit within that page; larger heatmaps use thinner rows.

#### Combined RUNX Heatmaps

RUNX panels also stay within one experiment and tissue and display WT, Dp16, and Dp16Runx1 in that order. Dp16Runx1 is the display label for `RUNX1-DEL_DP16`; the `RUNX1-DEL` group is omitted from these combined panels. The current data support Experiment 3 thymus (15 samples: 5 WT, 4 Dp16, 6 Dp16Runx1) and Experiment 4 thymocytes (14 samples: 4 WT, 5 Dp16, 5 Dp16Runx1). The unmatched Experiment 4 thymus sample is excluded from these panels.

Full heatmaps show the union of significant DMPs from the cohort's Dp16-versus-WT and Dp16Runx1-versus-WT comparisons, once per CpG. Rows are mean-centered across the displayed samples. Opposite significant directions across these two comparisons produce a gray Mixed block between red GOM and blue LOM, with counts for all three blocks.

The full panel is saved as `Experiment_<N>/figures/heatmap/combined/RUNX_vs_WT_Exp<N>_<tissue>.pdf`. A focused panel uses the same basename under `figures/heatmap/top1000/`: select up to 1,000 significant DP16-versus-WT CpGs, ranked by FDR with CpG ID as tie-breaker, then show the same three genotypes. Its direction labels use only the DP16-versus-WT effects. A panel with no eligible DMPs is skipped.

---

## Step 8: FASTA Sequence Extraction for Motif Analysis

```r
extract_and_write_FASTA <- function(DMPs, contrast_name, annoMouse, genome,
                                    fasta_dir, details_dir) {
  # Add genomic coordinates
  # Separate GOM (gain) vs LOM (loss) of methylation
  # Create GRanges with ±100bp windows
  # Extract sequences from mouse genome
  # Write FASTA files
}
```

The stored significant DMPs from each requested comparison are exported; no differential model is fitted during FASTA extraction.

Gains and losses are handled independently. If one direction has no annotated DMPs, extraction and sequence naming are skipped for that direction and its count is reported as zero. Its FASTA file contains zero sequences and its details CSV contains only column headers. Repeating extraction within the same run clears any previous files for that empty direction; a fresh workflow execution uses a new run folder. The other direction is exported normally. If neither direction has valid coordinates, or the comparison has no significant DMPs, extraction is skipped and no FASTA files are written for it.

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
Experiment_1/fasta/GOM/DP10_vs_WT_Exp1_thymus.fasta
Experiment_1/fasta/LOM/DP10_vs_WT_Exp1_thymus.fasta
Experiment_1/tables/fasta_details/GOM/DP10_vs_WT_Exp1_thymus.csv
Experiment_1/tables/fasta_details/LOM/DP10_vs_WT_Exp1_thymus.csv
```

`tables/summaries/FASTA_summary.csv` records annotated gain/loss CpG counts and the total number of sequences written per comparison. Sequence counts can be lower than CpG counts because overlapping windows are merged.

### Purpose for Motif Analysis:
- **GOM sequences**: May be enriched for repressor binding sites
- **LOM sequences**: May be enriched for activator binding sites
- **AME input**: Ready for Analysis of Motif Enrichment

---

## Output Organization

Set `results_root <- "results"` near the top of `workflow_new.R` to choose the destination. Each execution creates a dated run folder; repeated names get `_01`, `_02`, and so on. Existing runs are preserved. Each experiment gets its own folder within that run.

Experiments 1 and 2 append `_sameparent` to their folder names for the WT parental-cross rule. The generic experiment layout below also applies to those folders.

```
results/
└── analysis_YYYY-MM-DD_HH-MM-SS/
    ├── All_Experiments/
    │   ├── figures/pca/{genotype,batch,tissue}/
    │   └── tables/pca/{betas,cpgs,scores,variance}/
    ├── Experiment_1/                 # repeated for every experiment
    │   ├── figures/
    │   │   ├── heatmap/              # pairwise PDFs; combined/ and top1000/ for RUNX
    │   │   ├── volcano/
    │   │   └── pca/{genotype,batch}/
    │   ├── tables/
    │   │   ├── dmp/
    │   │   ├── pca/{betas,cpgs,scores,variance}/
    │   │   └── fasta_details/{GOM,LOM}/
    │   ├── fasta/{GOM,LOM}/
    │   └── metadata/samples.csv
    ├── tables/{betas,summaries}/
    ├── metadata/
    └── output_manifest.csv
```

Comparison files contain only the comparison, experiment, and tissue, for example `DP10_vs_WT_Exp2_Tcells.pdf`. Plot type, FASTA direction, and PCA table type are distinguished by their folders.

Start with `tables/summaries/DMP_summary.csv` for comparison counts and `FASTA_summary.csv` for sequence counts. Run-wide beta matrices stay under `tables/betas/`. `metadata/` stores the input samplesheet, retained samples, workflow/helper copies, settings, and R/package versions. Each experiment also has its own sample metadata.

`output_manifest.csv` lists every output by relative path and file size. A failed run can contain partial outputs without a final manifest. Skipped comparisons have no DMP tables or plots.

Validation scripts run with `Rscript --vanilla tests/<test>.R`; set `METHYL_DMP_TEST_LIBRARY` to the project's `packages/` path if needed:

- `test_dmp_analysis.R`: exact experiment/tissue controls, independent fits, short names, missing metadata, and minimum replication.
- `test_sameparent_dmp.R`: exact WT `Genotype` selection, retained single-Dp animals from double crosses, independence from `cross.type`, independent pairwise fits, unchanged other experiments, and no fallback to unmatched WT.
- `test_sample_matching.R`: brain excluded before IDAT selection, exact identifiers, corrected genotype, and QC sample ordering.
- `test_output_layout.R`: real CSV/PDF/FASTA output, pooled/experiment/tissue and Age17.5 PCA membership, pooled-only exclusions, combined RUNX direction calls, manifest, and preservation of previous runs.
- `test_pca_analysis.R`: independent feature selection, sample alignment, consistent colors, and rejected mixed-experiment inputs.
- `test_combat.R`: protected experiment/tissue/genotype groups and estimable batch correction.
- `test_heatmap_plots.R`: direction brackets, counts, mean centering, and PDF rendering.
- `test_fasta_export.R`: genomic sequence extraction, merged regions, and empty directions.

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
- `idats/samplesheet_169_samples.csv`: Sample metadata, including complete `Basename` identifiers
- `idats/`: Directory containing .idat files

### Mouse Genome:
- BSgenome.Mmusculus.UCSC.mm10 package
- IlluminaMouseMethylationanno.12.v1.mm10 annotation

---

This pipeline provides a complete, publication-ready analysis of mouse methylation array data with integrated motif analysis preparation. The modular design allows for easy modification and extension based on specific research needs.
