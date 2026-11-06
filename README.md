# Covariate-adjusted co-expression does not establish gene–immune checkpoint specificity in glioblastoma

Code for the article *Covariate-adjusted co-expression does not establish gene–immune checkpoint specificity in glioblastoma: a genome-wide negative-control analysis with EZH2 as a case study* (Siraj Ud Din, Sartaj Ud Din, Asif Ali Shah).

   Archived version: https://doi.org/10.5281/zenodo.24015878

## What the code does

Three public cohorts of primary IDH-wildtype glioblastoma (TCGA n = 160, CGGA_325 n = 74, CGGA_693 n = 109) are analysed. The standard analysis (partial Spearman correlations adjusted for proliferation and ESTIMATE scores, pooled across cohorts) is applied to *EZH2* and eight checkpoint genes. The same pipeline is then applied to matched control genes, to all 12,130 expressed genes (genome-wide null), to positive controls, to spiked-in signals and to genes from published glioma studies.

## Files

| Script | Content | Main outputs |
|---|---|---|
| `00_setup.R` | Packages, constants and helper functions (partial Spearman correlation, Fisher-z pooling). Sourced by the other scripts; not run on its own. | – |
| `01_data_and_cohorts.R` | Downloads TCGA-GBM RNA-seq and clinical data, prepares the CGGA files, computes ESTIMATE and MCP-counter scores, builds the three cohorts | `share/*.rds`, `share/scores_*.csv`, `results/cohorts.rds` |
| `02_main_analyses.R` | Expression, PRC2 and proliferation; survival; checkpoint partial correlations; microenvironment; single-cell expression; Hallmark sensitivity | Tables S1–S5, Supplementary Figs. S1–S3 and S5, Fig. 1 |
| `03_benchmark_and_robustness.R` | Measurement-error simulation; matched negative-control benchmark; survival, robustness, pooling and cohort sensitivity analyses; benchmark, simulation and forest-plot figures | Tables S6–S7, Fig. 2, Supplementary Figs. S4 and S6 |
| `04_genomewide_controls_audit.R` | Comparator-rule sensitivity, genome-wide calibration and marker genes; genome-wide null and principal-component adjustment; positive controls, spike-in power and audit of published claims; copies all figures under their article numbers | Tables S9–S11, Figs. 3–5, `figures/article/` |

## Requirements

R ≥ 4.3 (analyses were run with R 4.6.0). Packages:

- CRAN: `survival`, `survminer`, `ggplot2`, `patchwork`, `ggrepel`, `data.table`, `msigdbr`, `dplyr`, `stringr`, `BiocManager`, `remotes`
- Bioconductor: `TCGAbiolinks`, `SummarizedExperiment`, `GSVA`
- `estimate`: `install.packages("estimate", repos = "http://r-forge.r-project.org")`
- `immunedeconv`: `remotes::install_github("omnideconv/immunedeconv")`

Step 1 needs internet access and about 3 GB of free disk space. Steps 2–4 run offline.

## Input data you must download yourself

1. **CGGA** (http://www.cgga.org.cn): put these four files in `Data/geo/` inside the project folder:
   `CGGA.mRNAseq_325.RSEM-genes.20200506.txt`, `CGGA.mRNAseq_325_clinical.20200506.txt`,
   `CGGA.mRNAseq_693.RSEM-genes.20200506.txt`, `CGGA.mRNAseq_693_clinical.20200506.txt`
2. **TISCH2** (http://tisch.compbio.cn): put the heatmap exports for the dataset Glioma_GSE131928 in `share/`:
   `TISCH_EZH2_heatmap.csv`, `TISCH_LAG3_heatmap.csv`, `TISCH_PDCD1_heatmap.csv`

TCGA data are downloaded automatically from the Genomic Data Commons in step 1. A later GDC data release may differ slightly; step 1 prints the release number and warns if the cohort sizes differ from those in the article.

## How to run

Put the scripts in a folder `scripts/` inside your project folder, set the working directory to the project folder and run the four steps in order:

```r
setwd("path/to/project")                    # the folder that contains (or will contain) share/
source("scripts/01_data_and_cohorts.R")
source("scripts/02_main_analyses.R")
source("scripts/03_benchmark_and_robustness.R")
source("scripts/04_genomewide_controls_audit.R")
```

Each step writes a log to `results/log_0X_*.txt`, including the R session information. All analyses are deterministic except the simulation, which uses `set.seed(2026)`. Steps 3 and 4 take several minutes each.

## Article figures and tables

After step 4, `figures/article/` contains every figure under its article number (Fig1–Fig5, FigS1–FigS6). The result tables in `results/` correspond to the supplementary spreadsheet as follows:

| Spreadsheet | File in `results/` |
|---|---|
| S1a, S1b | `Table_S1a_*`, `Table_S1b_*` |
| S2a, S2b; S2c (Supplementary Information) | `Table_S2a_*`, `Table_S2b_*`; `Table2_Cox_models.csv` |
| S3a, S3b, S3c | `Table3a_*`, `Table3b_*`, `Table3c_*` |
| S3d, S3e | `Table_S3a_immune_by_cohort.csv`, `Table_S3b_immune_meta.csv` |
| S4, S5 | `Table_S4_*`, `Table_S5_*` |
| S6a–S6c | `Table_S6a_*`, `Table_S6b_*`, `Table_S6c_*` |
| S7a–S7i | `Table_S7a_*` … `Table_S7i_*` |
| S8 | literature survey, compiled manually (Supplementary Methods M8) |
| S9a–S9c, S10a–S10b, S11a–S11c | `Table_S9a_*` … `Table_S11c_*` (the empty strictly non-immune rows and columns of S10 are omitted in the spreadsheet) |

Intermediate figures written by individual analyses (for example `FigureS2_specificity_benchmark`, `FigureS3_simulation`) are superseded by the article figures.

## License

MIT. Please cite the article and the archived code.
