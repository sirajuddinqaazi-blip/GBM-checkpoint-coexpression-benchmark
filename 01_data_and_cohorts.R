# =====================================================================
# 01_data_and_cohorts.R
# Step 1: download and prepare data, composition scores, analysis cohorts
# Part of: Covariate-adjusted co-expression does not establish gene-immune
# checkpoint specificity in glioblastoma (EZH2 case study).
#
# HOW TO RUN: set the working directory to the PROJECT folder (the folder
# that contains share/), then source this script from wherever it is, e.g.
#   setwd("D:/Bioinformatics/EZH2_GBM_Bioinfo_2026")
#   source("scripts/01_data_and_cohorts.R")
# Outputs go to results/ and figures/ in the project folder; a log is written
# to results/log_01_data_and_cohorts.txt.
# Before running: place the four CGGA files in Data/geo/ and the TISCH2 heatmap
# CSVs in share/ (see README). Needs internet access for the TCGA download.
# =====================================================================
.this <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/"), error = function(e) NA_character_)
source(file.path(if (is.na(.this)) getwd() else dirname(.this), "00_setup.R"))
start_log("log_01_data_and_cohorts.txt")
tryCatch({

# =====================================================================
# 1a Download and prepare TCGA and CGGA data
# =====================================================================
step("1a Download and prepare TCGA and CGGA data")
local({
  set.seed(2026)
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  for (p in c("TCGAbiolinks", "SummarizedExperiment"))
    if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
  for (p in c("dplyr", "stringr"))
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

  suppressPackageStartupMessages({
    library(TCGAbiolinks); library(SummarizedExperiment)
    library(dplyr); library(stringr)
  })

  RAW_DIR  <- "Data/raw"          # GDC download cache
  PROC_DIR <- "Data/processed"
  CGGA_DIR <- "Data/geo"          # place the four CGGA files here (Part C)
  SHARE    <- "share"
  for (d in c(RAW_DIR, PROC_DIR, CGGA_DIR, SHARE)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  # Numbers obtained in the published analysis (GDC data release used for
  # the manuscript). If a later GDC release changes the data, the script
  # warns but continues; report the release number printed below.
  EXPECTED <- c(patients = 289, primary = 185, idh_wt = 162, final = 160, normals = 5)
  check_n <- function(label, observed) {
    exp <- EXPECTED[[label]]
    if (observed == exp) cat(sprintf("  ✓ %-9s n = %d (as expected)\n", label, observed))
    else warning(sprintf("%s: n = %d, expected %d. The GDC data may have changed since the published analysis.",
                         label, observed, exp), call. = FALSE)
  }

  # ---------------------------------------------------------------------
  # Part A. Download and prepare TCGA-GBM RNA-seq
  # ---------------------------------------------------------------------
  gbm_file <- file.path(PROC_DIR, "gbm_data.rds")
  if (file.exists(gbm_file)) {
    cat("Part A: found", gbm_file, "- skipping download\n")
    gbm_data <- readRDS(gbm_file)
  } else {
    cat("Part A: querying GDC ...\n")
    info <- tryCatch(getGDCInfo(), error = function(e) NULL)
    if (!is.null(info)) cat("  GDC data release:", info$data_release, "\n")

    query <- GDCquery(
      project       = "TCGA-GBM",
      data.category = "Transcriptome Profiling",
      data.type     = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      sample.type   = c("Primary Tumor", "Solid Tissue Normal"))
    print(table(getResults(query)$sample_type))

    GDCdownload(query, method = "api", files.per.chunk = 10, directory = RAW_DIR)
    gbm_data <- GDCprepare(query, directory = RAW_DIR, summarizedExperiment = TRUE)
    saveRDS(gbm_data, gbm_file)
    if (!is.null(info)) writeLines(paste("GDC data release:", info$data_release),
                                   file.path(PROC_DIR, "GDC_release.txt"))
  }
  stopifnot(all(c("tpm_unstrand", "unstranded") %in% assayNames(gbm_data)))
  cat("  Samples:", ncol(gbm_data), "| genes:", nrow(gbm_data), "\n")

  # ---------------------------------------------------------------------
  # Part B. Build the TCGA cohort (identical rules to the published analysis)
  # ---------------------------------------------------------------------
  cat("\nPart B: building the IDH-wildtype primary GBM cohort\n")
  clin <- as.data.frame(colData(gbm_data))
  needed <- c("barcode", "patient", "sample_type", "classification_of_tumor", "paper_IDH.status",
              "vital_status", "days_to_death", "days_to_last_follow_up")
  miss <- setdiff(needed, colnames(clin))
  if (length(miss)) stop("Missing clinical columns: ", paste(miss, collapse = ", "))

  # B1. One sample per patient: prefer Primary Tumor, then the newer
  #     GDC re-processing plate (plate code starting with "A96").
  clin_unique <- clin %>%
    mutate(plate_code = str_extract(barcode, "(?<=-)[A-Z0-9]{4}(?=-\\d{2}$)")) %>%
    arrange(patient,
            desc(sample_type == "Primary Tumor"),
            desc(str_starts(plate_code, "A96"))) %>%
    distinct(patient, .keep_all = TRUE)
  check_n("patients", nrow(clin_unique))

  # B2. Primary (not recurrent) tumors
  gbm_primary <- clin_unique %>%
    filter(sample_type == "Primary Tumor", classification_of_tumor == "primary")
  check_n("primary", nrow(gbm_primary))

  # B3. IDH-wildtype (curated pan-glioma annotation, Ceccarelli et al. 2016)
  gbm_idhwt <- gbm_primary %>% filter(paper_IDH.status == "WT")
  check_n("idh_wt", nrow(gbm_idhwt))

  # B4. Complete survival information
  cohort <- gbm_idhwt %>%
    filter(!is.na(vital_status), vital_status %in% c("Alive", "Dead"),
           !(vital_status == "Dead"  & is.na(days_to_death)),
           !(vital_status == "Alive" & is.na(days_to_last_follow_up)))
  check_n("final", nrow(cohort))

  # B5. EZH2 expression (TPM) and median split (kept for compatibility)
  ezh2_row <- grep("^ENSG00000106462(\\.|$)", rownames(gbm_data))
  stopifnot(length(ezh2_row) == 1)
  cohort$EZH2_expression <- as.numeric(assay(gbm_data, "tpm_unstrand")[ezh2_row, rownames(cohort)])
  cohort$EZH2_group <- ifelse(cohort$EZH2_expression >= median(cohort$EZH2_expression), "High_EZH2", "Low_EZH2")
  saveRDS(cohort, file.path(PROC_DIR, "tumor_final_EZH2_analysis_ready.rds"))
  cat(sprintf("  EZH2 median = %.2f TPM (published: 17.38)\n", median(cohort$EZH2_expression)))

  # ---------------------------------------------------------------------
  # Part C. CGGA files (manual download)
  # ---------------------------------------------------------------------
  # From http://www.cgga.org.cn  (Download -> mRNAseq_325 and mRNAseq_693)
  # download and UNZIP these four files into Data/geo/ :
  #   CGGA.mRNAseq_325.RSEM-genes.20200506.txt
  #   CGGA.mRNAseq_325_clinical.20200506.txt
  #   CGGA.mRNAseq_693.RSEM-genes.20200506.txt
  #   CGGA.mRNAseq_693_clinical.20200506.txt
  cat("\nPart C: checking CGGA files\n")
  cgga_files <- file.path(CGGA_DIR, c(
    "CGGA.mRNAseq_325.RSEM-genes.20200506.txt", "CGGA.mRNAseq_325_clinical.20200506.txt",
    "CGGA.mRNAseq_693.RSEM-genes.20200506.txt", "CGGA.mRNAseq_693_clinical.20200506.txt"))
  absent <- cgga_files[!file.exists(cgga_files)]
  if (length(absent)) stop("Download and unzip the CGGA files first (see Part C comments). Missing:\n  ",
                           paste(absent, collapse = "\n  "))
  cat("  ✓ all four CGGA files present\n")

  # ---------------------------------------------------------------------
  # Part D. Export compact input files to share/
  # ---------------------------------------------------------------------
  cat("\nPart D: exporting to", SHARE, "\n")
  normal_bc <- colnames(gbm_data)[colData(gbm_data)$sample_type == "Solid Tissue Normal"]
  check_n("normals", length(normal_bc))
  keep_bc <- c(rownames(cohort), normal_bc)
  stopifnot(all(keep_bc %in% colnames(gbm_data)))

  gi <- as.data.frame(rowData(gbm_data))
  gi <- gi[, intersect(c("gene_id", "gene_name", "gene_type"), colnames(gi))]
  saveRDS(gi, file.path(SHARE, "tcga_gene_info.rds"), compress = "xz")

  tpm <- round(assay(gbm_data, "tpm_unstrand")[, keep_bc], 3)
  tpm <- tpm[rowSums(tpm) > 0, ]
  saveRDS(tpm, file.path(SHARE, "tcga_tpm.rds"), compress = "xz")

  cnt <- assay(gbm_data, "unstranded")[, rownames(cohort)]
  cnt <- cnt[rowSums(cnt) > 0, ]
  saveRDS(cnt, file.path(SHARE, "tcga_counts.rds"), compress = "xz")

  saveRDS(cohort, file.path(SHARE, "tcga_cohort_clinical.rds"), compress = "xz")
  saveRDS(as.data.frame(colData(gbm_data))[normal_bc, ], file.path(SHARE, "tcga_normals_clinical.rds"))

  for (id in c("325", "693")) {
    expr <- read.delim(file.path(CGGA_DIR, paste0("CGGA.mRNAseq_", id, ".RSEM-genes.20200506.txt")), check.names = FALSE)
    clin_c <- read.delim(file.path(CGGA_DIR, paste0("CGGA.mRNAseq_", id, "_clinical.20200506.txt")), check.names = FALSE)
    colnames(clin_c) <- sub(" \\(.*\\)$", "", colnames(clin_c))
    gbm_ids <- clin_c$CGGA_ID[trimws(clin_c$PRS_type) == "Primary" & trimws(clin_c$Histology) == "GBM"]
    expr_sub <- expr[, c(colnames(expr)[1], intersect(gbm_ids, colnames(expr)))]
    saveRDS(expr_sub, file.path(SHARE, paste0("cgga", id, "_expr_gbm.rds")), compress = "xz")
    saveRDS(clin_c,   file.path(SHARE, paste0("cgga", id, "_clinical.rds")))
    cat(sprintf("  ✓ CGGA_%s: %d primary GBM samples (published: %s)\n",
                id, ncol(expr_sub) - 1, ifelse(id == "325", "85", "140")))
  }

  writeLines(capture.output(sessionInfo()), file.path(SHARE, "sessionInfo_data_preparation.txt"))
  cat("\n✓ Data preparation complete.\n")
})

# =====================================================================
# 1b ESTIMATE and MCP-counter scores
# =====================================================================
step("1b ESTIMATE and MCP-counter scores")
local({
  set.seed(2026)
  library(estimate)
  library(immunedeconv)
  library(SummarizedExperiment)   # not strictly needed; keeps environment consistent

  dir.create("share/tmp_estimate", recursive = TRUE, showWarnings = FALSE)

  # ---------- helpers ----------
  collapse_by_symbol <- function(mat, symbols) {
    keep <- !is.na(symbols) & symbols != ""
    mat <- mat[keep, , drop = FALSE]; symbols <- symbols[keep]
    mean_expr <- rowMeans(mat)
    idx <- tapply(seq_len(nrow(mat)), symbols, function(i) i[which.max(mean_expr[i])])
    out <- mat[unlist(idx), , drop = FALSE]
    rownames(out) <- names(idx)
    out
  }

  run_estimate <- function(mat, label) {
    in_f  <- file.path("share/tmp_estimate", paste0(label, "_input.txt"))
    cg_f  <- file.path("share/tmp_estimate", paste0(label, "_common.gct"))
    out_f <- file.path("share/tmp_estimate", paste0(label, "_estimate.gct"))
    df <- data.frame(NAME = rownames(mat), log2(mat + 1), check.names = FALSE)
    write.table(df, in_f, sep = "\t", quote = FALSE, row.names = FALSE)
    filterCommonGenes(input.f = in_f, output.f = cg_f, id = "GeneSymbol")
    estimateScore(input.ds = cg_f, output.ds = out_f, platform = "illumina")
    sc <- read.table(out_f, skip = 2, header = TRUE, sep = "\t", check.names = FALSE)
    rownames(sc) <- sc$NAME
    sc <- t(sc[, -(1:2)])
    data.frame(sample = colnames(mat),        # keep original IDs (read.table may alter them)
               ESTIMATE_stromal = as.numeric(sc[, "StromalScore"]),
               ESTIMATE_immune  = as.numeric(sc[, "ImmuneScore"]),
               ESTIMATE_score   = as.numeric(sc[, "ESTIMATEScore"]),
               check.names = FALSE)
  }

  run_mcp <- function(mat) {
    res <- immunedeconv::deconvolute(mat, "mcp_counter")
    m <- as.data.frame(t(as.matrix(res[, -1])))
    colnames(m) <- paste0("MCP_", gsub("[^A-Za-z0-9]+", "_", res$cell_type))
    m$sample <- rownames(m)
    m
  }

  process <- function(mat, label) {
    cat("\n=== ", label, ": ", nrow(mat), " genes x ", ncol(mat), " samples ===\n", sep = "")
    est <- run_estimate(mat, label)
    mcp <- run_mcp(mat)
    out <- merge(est, mcp, by = "sample")
    stopifnot(nrow(out) == ncol(mat))
    write.csv(out, paste0("share/scores_", label, ".csv"), row.names = FALSE)
    cat("✓ Saved share/scores_", label, ".csv\n", sep = "")
    print(summary(out$ESTIMATE_immune))
  }

  # ---------- TCGA ----------
  tpm    <- readRDS("share/tcga_tpm.rds")
  gi     <- readRDS("share/tcga_gene_info.rds")
  cohort <- readRDS("share/tcga_cohort_clinical.rds")
  sym    <- gi$gene_name[match(rownames(tpm), rownames(gi))]
  tcga   <- collapse_by_symbol(tpm[, rownames(cohort)], sym)
  process(tcga, "TCGA")

  # ---------- CGGA ----------
  for (id in c("325", "693")) {
    ex <- readRDS(paste0("share/cgga", id, "_expr_gbm.rds"))
    mat <- as.matrix(ex[, -1]); mode(mat) <- "numeric"
    mat <- collapse_by_symbol(mat, ex[[1]])
    process(mat, paste0("CGGA", id))
  }

  cat("\n✓ All done. Upload: share/scores_TCGA.csv, share/scores_CGGA325.csv, share/scores_CGGA693.csv\n")
})

# =====================================================================
# 1c Build the three analysis cohorts
# =====================================================================
step("1c Build the three analysis cohorts")
local({
  set.seed(2026)
  num <- function(x) suppressWarnings(as.numeric(as.character(x)))

  # ------------------------- TCGA -------------------------
  tpm <- readRDS(file.path(DATA_DIR, "tcga_tpm.rds"))
  gi  <- readRDS(file.path(DATA_DIR, "tcga_gene_info.rds"))
  cl  <- readRDS(file.path(DATA_DIR, "tcga_cohort_clinical.rds"))
  nc  <- readRDS(file.path(DATA_DIR, "tcga_normals_clinical.rds"))
  sym <- gi$gene_name[match(rownames(tpm), rownames(gi))]

  get_tcga <- function(genes, samples) {
    sapply(genes, function(g) {
      i <- which(sym == g)
      if (length(i) == 0) stop("Gene not found in TCGA: ", g)
      i <- i[which.max(rowMeans(tpm[i, , drop = FALSE]))]   # highest-expressed Ensembl ID
      as.numeric(tpm[i, samples])
    })
  }

  tcga <- data.frame(sample = rownames(cl), get_tcga(GENES_NEEDED, rownames(cl)), check.names = FALSE)
  os_days <- ifelse(cl$vital_status == "Dead", cl$days_to_death, cl$days_to_last_follow_up)
  tcga <- within(tcga, {
    OS_months <- os_days / 30.44
    event     <- as.integer(cl$vital_status == "Dead")
    age       <- num(cl$age_at_index)
    sex       <- factor(tolower(cl$paper_Gender), levels = c("female", "male"))
    MGMT      <- factor(cl$paper_MGMT.promoter.status, levels = c("Methylated", "Unmethylated"))
    KPS       <- num(cl$paper_Karnofsky.Performance.Score)
    subtype   <- factor(cl$paper_Transcriptome.Subtype, levels = c("CL", "ME", "NE", "PN"),
                        labels = c("Classical", "Mesenchymal", "Neural", "Proneural"))
    ABSOLUTE_purity <- num(cl$paper_ABSOLUTE.purity)
  })
  tcga$cohort <- "TCGA"

  normals <- data.frame(sample = rownames(nc), get_tcga("EZH2", rownames(nc)))
  colnames(normals)[2] <- "EZH2"

  # ------------------------- CGGA -------------------------
  build_cgga <- function(id) {
    ex <- readRDS(file.path(DATA_DIR, paste0("cgga", id, "_expr_gbm.rds")))
    c2 <- readRDS(file.path(DATA_DIR, paste0("cgga", id, "_clinical.rds")))
    for (cc in names(c2)) if (is.character(c2[[cc]])) c2[[cc]] <- trimws(c2[[cc]])
    keep <- c2$PRS_type == "Primary" & c2$Histology == "GBM" & c2$Grade == "WHO IV" &
            c2$IDH_mutation_status %in% "Wildtype"
    ids <- intersect(c2$CGGA_ID[keep], colnames(ex))
    miss <- setdiff(GENES_NEEDED, ex[[1]]); if (length(miss)) stop("Missing in CGGA_", id, ": ", paste(miss, collapse = ","))
    m <- ex[match(GENES_NEEDED, ex[[1]]), ids]
    d <- data.frame(sample = ids, t(apply(m, 2, as.numeric)), check.names = FALSE)
    colnames(d)[-1] <- GENES_NEEDED
    c2 <- c2[match(ids, c2$CGGA_ID), ]
    d$OS_months <- c2$OS / 30.44
    d$event     <- c2$Censor               # CGGA: 1 = dead, 0 = alive
    d$age       <- c2$Age
    d$sex       <- factor(tolower(c2$Gender), levels = c("female", "male"))
    d$MGMT      <- factor(c2$MGMTp_methylation_status, levels = c("methylated", "un-methylated"),
                          labels = c("Methylated", "Unmethylated"))
    d$radio     <- c2$Radio_status; d$chemo <- c2$Chemo_status
    d$cohort    <- paste0("CGGA_", id)
    d
  }
  cohorts <- list(TCGA = tcga, CGGA_325 = build_cgga("325"), CGGA_693 = build_cgga("693"))

  # ------------------- scores, derived variables -------------------
  score_files <- c(TCGA = "scores_TCGA.csv", CGGA_325 = "scores_CGGA325.csv", CGGA_693 = "scores_CGGA693.csv")
  for (k in names(cohorts)) {
    d  <- cohorts[[k]]
    sc <- read.csv(file.path(DATA_DIR, score_files[[k]]), check.names = FALSE)
    d  <- merge(d, sc, by = "sample", all.x = TRUE, sort = FALSE)
    if (any(is.na(d$ESTIMATE_immune))) stop("Missing ESTIMATE scores in ", k)
    d$prolif     <- rowMeans(scale(log2(d[, PROLIF_GENES] + 1)))   # 10-gene proliferation score
    d$log2EZH2   <- log2(d$EZH2 + 1)
    d$EZH2_group <- factor(ifelse(d$EZH2 >= median(d$EZH2), "High", "Low"), levels = c("Low", "High"))
    cohorts[[k]] <- d
  }

  saveRDS(list(cohorts = cohorts, tcga_normals = normals), file.path(RES_DIR, "cohorts.rds"))
  for (k in names(cohorts)) cat(sprintf("%-9s n = %3d | events = %3d | EZH2 median = %.2f\n",
    k, nrow(cohorts[[k]]), sum(cohorts[[k]]$event, na.rm = TRUE), median(cohorts[[k]]$EZH2)))
  cat("TCGA normal brain n =", nrow(normals), "\n✓ results/cohorts.rds saved\n")
})

cat("\n\u2713 01_data_and_cohorts.R complete.\n")
}, finally = end_log())
