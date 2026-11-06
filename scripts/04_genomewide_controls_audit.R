# =====================================================================
# 04_genomewide_controls_audit.R
# Step 4: generalization, genome-wide null, positive controls, spike-in power and audit
# Part of: Covariate-adjusted co-expression does not establish gene-immune
# checkpoint specificity in glioblastoma (EZH2 case study).
#
# HOW TO RUN: set the working directory to the PROJECT folder (the folder
# that contains share/), then source this script from wherever it is, e.g.
#   setwd("D:/Bioinformatics/EZH2_GBM_Bioinfo_2026")
#   source("scripts/04_genomewide_controls_audit.R")
# Outputs go to results/ and figures/ in the project folder; a log is written
# to results/log_04_genomewide_controls_audit.txt.
# =====================================================================
.this <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/"), error = function(e) NA_character_)
source(file.path(if (is.na(.this)) getwd() else dirname(.this), "00_setup.R"))
check_project(c("results/cohorts.rds"))
start_log("log_04_genomewide_controls_audit.txt")
tryCatch({

# =====================================================================
# 4a Comparator rules, genome-wide calibration, marker genes
# =====================================================================
step("4a Comparator rules, genome-wide calibration, marker genes")
local({
  set.seed(2026)
  say <- function(...) cat(sprintf(...), "\n", sep = "")
  hdr <- function(x) cat("\n\n==================== ", x, " ====================\n", sep = "")
  r3  <- function(x) round(x, 3)
  co  <- readRDS(file.path(RES_DIR, "cohorts.rds"))$cohorts
  EST <- c("ESTIMATE_immune", "ESTIMATE_stromal"); M3 <- c("prolif", EST); M4 <- c(M3, "MCP_T_cell")
  MIN_MEDIAN <- 1

  # ---------------- expression matrices (same handling as scripts 07/08) ----------------
  collapse_symbols <- function(mat, sym) {
    ok <- !is.na(sym) & sym != ""; mat <- mat[ok, , drop = FALSE]; sym <- sym[ok]
    o <- order(rowMeans(mat), decreasing = TRUE); mat <- mat[o, , drop = FALSE]; sym <- sym[o]
    keep <- !duplicated(sym); mat <- mat[keep, , drop = FALSE]; rownames(mat) <- sym[keep]; mat
  }
  tpm <- readRDS(file.path(DATA_DIR, "tcga_tpm.rds")); gi <- readRDS(file.path(DATA_DIR, "tcga_gene_info.rds"))
  mats <- list(TCGA = collapse_symbols(as.matrix(tpm[, co$TCGA$sample, drop = FALSE]), gi$gene_name[match(rownames(tpm), rownames(gi))]))
  rm(tpm); invisible(gc())
  for (id in c("325", "693")) {
    ex <- readRDS(file.path(DATA_DIR, paste0("cgga", id, "_expr_gbm.rds")))
    m <- as.matrix(ex[, -1]); storage.mode(m) <- "numeric"
    mats[[paste0("CGGA_", id)]] <- collapse_symbols(m, as.character(ex[[1]]))
  }
  cc_genes <- tryCatch({
    h <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
                  error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
    setdiff(unique(h$gene_symbol[h$gs_name %in% c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")]), c("EZH2", CHECKPOINTS))
  }, error = function(e) NULL)

  # ---------------- vectorised partial Spearman + pooling ----------------
  resid_cols <- function(Y, Z) { Z1 <- cbind(1, Z); Y - Z1 %*% solve(crossprod(Z1), crossprod(Z1, Y)) }
  cor_mat <- function(A, B) { A <- scale(A); B <- scale(B); crossprod(A, B) / (nrow(A) - 1) }
  pool <- function(r3d, n, k) {                     # r3d: genes x cohorts
    z <- atanh(pmin(pmax(r3d, -0.999999), 0.999999)); w <- n - 3 - k
    zp <- as.vector(z %*% w) / sum(w); se <- 1 / sqrt(sum(w))
    data.frame(r = tanh(zp), p = 2 * pnorm(-abs(zp / se)))
  }
  prep <- list()
  for (k in names(co)) {
    d <- co[[k]]; m <- mats[[k]][, d$sample, drop = FALSE]
    lg <- log2(m + 1); lg <- lg[apply(lg, 1, var) > 0, , drop = FALSE]
    med <- apply(m[rownames(lg), , drop = FALSE], 1, median)
    R  <- apply(lg[names(med)[med >= MIN_MEDIAN], , drop = FALSE], 1, rank)          # samples x genes
    R  <- cbind(R[, setdiff(colnames(R), "EZH2"), drop = FALSE], EZH2 = rank(d$EZH2))
    Tk <- sapply(CHECKPOINTS, function(t) rank(d[[t]]))
    S <- ng <- Zcc <- NULL
    if (!is.null(cc_genes)) { g <- intersect(cc_genes, rownames(lg)); Zcc <- scale(t(lg[g, , drop = FALSE])); S <- rowSums(Zcc); ng <- length(g) }
    prep[[k]] <- list(d = d, R = R, Tk = Tk, Zcc = Zcc, S = S, ng = ng,
                      rho = drop(cor(R, rank(d$prolif))), rho_imm = drop(cor(R, rank(d$ESTIMATE_immune))))
  }
  n_c <- sapply(prep, function(p) nrow(p$d))
  genes <- Reduce(intersect, lapply(prep, function(p) colnames(p$R)))
  rho_mat <- sapply(prep, function(p) p$rho[genes]); rownames(rho_mat) <- genes
  imm_mat <- sapply(prep, function(p) p$rho_imm[genes]); rownames(imm_mat) <- genes
  excl <- unique(c(PROLIF_GENES, CHECKPOINTS))
  say("Expressed genes common to all cohorts: %d", length(genes))

  partial_all <- function(gs, Zvars, cc = FALSE) {   # returns genes x targets x cohorts
    arr <- array(NA_real_, c(length(gs), length(CHECKPOINTS), length(co)), dimnames = list(gs, CHECKPOINTS, names(co)))
    for (k in names(co)) {
      P <- prep[[k]]; d <- P$d
      Zr <- if (length(Zvars)) sapply(Zvars, function(z) rank(d[[z]])) else matrix(0, nrow(d), 0)
      if (cc) Zr <- cbind(Zr, cc = rank(P$S / P$ng))
      arr[, , k] <- cor_mat(resid_cols(P$R[, gs, drop = FALSE], Zr), resid_cols(P$Tk, Zr))
      if (cc) for (g in intersect(gs, colnames(P$Zcc))) {
        Zg <- cbind(Zr[, setdiff(colnames(Zr), "cc"), drop = FALSE], cc = rank((P$S - P$Zcc[, g]) / (P$ng - 1)))
        arr[g, , k] <- cor_mat(resid_cols(P$R[, g, drop = FALSE], Zg), resid_cols(P$Tk, Zg))
      }
    }
    arr
  }
  # check against pspear for EZH2
  chk <- partial_all("EZH2", M3)
  for (k in names(co)) stopifnot(abs(chk["EZH2", "LAG3", k] - pspear("EZH2", "LAG3", M3, co[[k]])[["r"]]) < 1e-8)
  say("Vectorised estimates match pspear() for EZH2 (|difference| < 1e-8).")

  MODELS <- list(M0 = list(Z = character(0)), M3 = list(Z = M3), M4 = list(Z = M4))
  if (!is.null(cc_genes)) MODELS$M3h <- list(Z = M3, cc = TRUE)
  pooled <- list()
  for (mn in names(MODELS)) {
    arr <- partial_all(genes, MODELS[[mn]]$Z, isTRUE(MODELS[[mn]]$cc))
    kk <- length(MODELS[[mn]]$Z) + isTRUE(MODELS[[mn]]$cc)
    pooled[[mn]] <- lapply(CHECKPOINTS, function(t) { x <- pool(arr[, t, ], n_c, kk); rownames(x) <- genes; x })
    names(pooled[[mn]]) <- CHECKPOINTS
  }
  rho_mean <- rowMeans(rho_mat)

  # =====================================================================
  # A. Sensitivity of the EZH2 benchmark to comparator rules
  # =====================================================================
  hdr("A. EZH2 benchmark under alternative comparator rules")
  rowsA <- list()
  for (thr in c(0.70, 0.75, 0.80)) {
    cand <- genes[apply(rho_mat >= thr, 1, all) & !(genes %in% c("EZH2", excl))]
    for (tol in c(0.03, 0.05, 0.10, Inf)) {
      set <- cand[abs(rho_mean[cand] - rho_mean["EZH2"]) <= tol]
      for (mn in intersect(c("M3", "M3h", "M4"), names(MODELS))) for (t in c("LAG3", "PDCD1", "CD276")) {
        x <- pooled[[mn]][[t]]; e <- x["EZH2", "r"]; c_r <- x[set, "r"]
        rowsA[[length(rowsA) + 1]] <- data.frame(rho_min = thr, tolerance = ifelse(is.finite(tol), tol, NA), model = mn, target = t,
          N = length(set), EZH2_r = e, control_median = median(c_r), n_ge = sum(c_r >= e),
          empirical_p = (1 + sum(c_r >= e)) / (1 + length(set)),
          frac_pos_sig = mean(c_r > 0 & x[set, "p"] < 0.05))
      }
    }
  }
  A <- do.call(rbind, rowsA)
  write.csv(A, file.path(RES_DIR, "Table_S9a_benchmark_rule_sensitivity.csv"), row.names = FALSE)
  print(transform(subset(A, model == "M3"), EZH2_r = r3(EZH2_r), control_median = r3(control_median),
                  empirical_p = r3(empirical_p), frac_pos_sig = r3(frac_pos_sig)), row.names = FALSE)
  cat("(tolerance NA = all coupled genes, no matching)\n")

  # =====================================================================
  # B. Genome-wide calibration by coupling to proliferation
  # =====================================================================
  hdr("B. Genome-wide calibration (all expressed genes x 8 checkpoints)")
  bins <- c(-1, 0.2, 0.4, 0.6, 0.75, 1)
  lab  <- c("<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.75", ">=0.75")
  gene_bin <- cut(rho_mean, bins, labels = lab, include.lowest = TRUE)
  non_immune <- apply(abs(imm_mat) < 0.3, 1, all)
  rowsB <- list()
  for (mn in intersect(c("M0", "M3", "M4"), names(MODELS))) for (t in CHECKPOINTS) for (sub in c("all genes", "non-immune genes")) {
    x <- pooled[[mn]][[t]]
    keep <- !(genes %in% c("EZH2", excl)) & (if (sub == "non-immune genes") non_immune else TRUE)
    q <- p.adjust(x$p[keep], "BH"); names(q) <- genes[keep]
    for (b in lab) {
      g <- genes[keep & gene_bin == b]
      if (!length(g)) next
      rowsB[[length(rowsB) + 1]] <- data.frame(model = mn, target = t, gene_set = sub, coupling_bin = b, n_genes = length(g),
        median_r = median(x[g, "r"]), frac_pos_nominal = mean(x[g, "r"] > 0 & x[g, "p"] < 0.05),
        frac_neg_nominal = mean(x[g, "r"] < 0 & x[g, "p"] < 0.05), frac_pos_BH = mean(x[g, "r"] > 0 & q[g] < 0.05))
    }
  }
  B <- do.call(rbind, rowsB)
  B$coupling_bin <- factor(B$coupling_bin, levels = lab)
  write.csv(B, file.path(RES_DIR, "Table_S9b_genomewide_calibration.csv"), row.names = FALSE)
  cat("--- Share of genes with a positive, nominally significant pooled association (M3), by coupling bin ---\n")
  w <- reshape(subset(B, model == "M3" & gene_set == "all genes")[, c("target", "coupling_bin", "frac_pos_nominal")],
               idvar = "target", timevar = "coupling_bin", direction = "wide")
  names(w) <- sub("frac_pos_nominal.", "", names(w)); print(transform(w, check.names = FALSE), row.names = FALSE, digits = 2)
  cat("--- Same, non-immune genes only (|rho with ESTIMATE immune| < 0.3 in all cohorts) ---\n")
  w2 <- reshape(subset(B, model == "M3" & gene_set == "non-immune genes")[, c("target", "coupling_bin", "frac_pos_nominal")],
                idvar = "target", timevar = "coupling_bin", direction = "wide")
  names(w2) <- sub("frac_pos_nominal.", "", names(w2)); print(w2, row.names = FALSE, digits = 2)
  say("Genes per coupling bin: %s", paste(names(table(gene_bin)), table(gene_bin), sep = "=", collapse = " | "))

  pB <- ggplot(subset(B, model %in% c("M3", "M4")),
               aes(coupling_bin, frac_pos_nominal, colour = model, group = interaction(model, gene_set), linetype = gene_set)) +
    geom_hline(yintercept = 0.025, linetype = 3, colour = "grey50") +
    geom_line() + geom_point(size = 1.4) + facet_wrap(~target, nrow = 2) +
    scale_colour_manual(values = c(M3 = "#C0392B", M4 = "#2C3E50"),
                        labels = c(M3 = "Proliferation + ESTIMATE (M3)", M4 = "M3 + T cells (M4)")) +
    labs(x = "Mean correlation of the gene with the ten-gene proliferation score", y = "Share of genes with positive pooled P < 0.05",
         colour = NULL, linetype = NULL, title = "Share of genes with a positive adjusted checkpoint association, by coupling to proliferation",
         subtitle = "Dotted line: 0.025, the expected share under no association (one-sided half of P < 0.05)") +
    theme_pub + theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1), strip.text = element_text(face = "bold.italic"))
  save_fig(pB, "FigureS6_calibration", 11, 6.5)

  # =====================================================================
  # C. Published proliferation-coupled glioma markers
  # =====================================================================
  hdr("C. Proliferation-coupled marker genes, each against its own matched controls")
  MARKERS <- c("EZH2", "FOXM1", "CENPA", "E2F1", "MELK", "UBE2T", "KIF4A", "HELLS", "MCM8", "BUB1B", "CDK4", "AURKB", "CDC20", "PLK1")
  MARKERS <- intersect(MARKERS, genes)
  pool_set <- genes[!(genes %in% excl) | genes %in% MARKERS]
  rowsC <- list()
  for (g in MARKERS) {
    ctrl <- setdiff(pool_set[abs(rho_mean[pool_set] - rho_mean[g]) <= 0.05], c(g, MARKERS, excl))
    for (mn in intersect(c("M3", "M4"), names(MODELS))) for (t in CHECKPOINTS) {
      x <- pooled[[mn]][[t]]; e <- x[g, "r"]; c_r <- x[ctrl, "r"]
      rowsC[[length(rowsC) + 1]] <- data.frame(gene = g, score_gene = g %in% PROLIF_GENES, coupling = rho_mean[[g]], model = mn, target = t,
        r_pooled = e, p_pooled = x[g, "p"], n_controls = length(ctrl), control_median = median(c_r),
        rank = 1 + sum(c_r > e), empirical_p = (1 + sum(c_r >= e)) / (1 + length(ctrl)))
    }
  }
  C <- do.call(rbind, rowsC)
  C$standard_test_positive <- C$r_pooled > 0 & C$p_pooled < 0.05
  C$benchmark_specific <- C$empirical_p < 0.05
  write.csv(C, file.path(RES_DIR, "Table_S9c_marker_genes.csv"), row.names = FALSE)
  cs <- subset(C, model == "M3")
  cat("--- M3: genes passing the standard test (positive, P < 0.05) vs genes specific in the benchmark (empirical P < 0.05) ---\n")
  print(aggregate(cbind(standard_test_positive, benchmark_specific) ~ gene, data = cs, FUN = sum), row.names = FALSE)
  say("Over all marker-checkpoint pairs (M3): %d of %d pass the standard test; %d of those are specific in the benchmark.",
      sum(cs$standard_test_positive), nrow(cs), sum(cs$standard_test_positive & cs$benchmark_specific))
  print(transform(subset(cs, target %in% c("LAG3", "PDCD1"))[, c("gene", "coupling", "target", "r_pooled", "p_pooled", "control_median", "rank", "n_controls", "empirical_p")],
                  coupling = r3(coupling), r_pooled = r3(r_pooled), p_pooled = fmt_p(p_pooled), control_median = r3(control_median), empirical_p = r3(empirical_p)), row.names = FALSE)
})

# =====================================================================
# 4b Genome-wide null and principal-component adjustment
# =====================================================================
step("4b Genome-wide null and principal-component adjustment")
local({
  set.seed(2026)
  suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
  say <- function(...) cat(sprintf(...), "\n", sep = "")
  hdr <- function(x) cat("\n\n==================== ", x, " ====================\n", sep = "")
  r3  <- function(x) round(x, 3)
  co  <- readRDS(file.path(RES_DIR, "cohorts.rds"))$cohorts
  EST <- c("ESTIMATE_immune", "ESTIMATE_stromal"); M3 <- c("prolif", EST); M4 <- c(M3, "MCP_T_cell")
  MIN_MEDIAN <- 1
  excl <- unique(c(PROLIF_GENES, CHECKPOINTS))

  # ---------------- expression matrices (same handling as scripts 07/08/11) ----------------
  collapse_symbols <- function(mat, sym) {
    ok <- !is.na(sym) & sym != ""; mat <- mat[ok, , drop = FALSE]; sym <- sym[ok]
    o <- order(rowMeans(mat), decreasing = TRUE); mat <- mat[o, , drop = FALSE]; sym <- sym[o]
    keep <- !duplicated(sym); mat <- mat[keep, , drop = FALSE]; rownames(mat) <- sym[keep]; mat
  }
  tpm <- readRDS(file.path(DATA_DIR, "tcga_tpm.rds")); gi <- readRDS(file.path(DATA_DIR, "tcga_gene_info.rds"))
  mats <- list(TCGA = collapse_symbols(as.matrix(tpm[, co$TCGA$sample, drop = FALSE]), gi$gene_name[match(rownames(tpm), rownames(gi))]))
  rm(tpm); invisible(gc())
  for (id in c("325", "693")) {
    ex <- readRDS(file.path(DATA_DIR, paste0("cgga", id, "_expr_gbm.rds")))
    m <- as.matrix(ex[, -1]); storage.mode(m) <- "numeric"
    mats[[paste0("CGGA_", id)]] <- collapse_symbols(m, as.character(ex[[1]]))
  }
  cc_genes <- tryCatch({
    h <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
                  error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
    setdiff(unique(h$gene_symbol[h$gs_name %in% c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")]), c("EZH2", CHECKPOINTS))
  }, error = function(e) NULL)

  # ---------------- per-cohort preparation ----------------
  resid_cols <- function(Y, Z) { Z1 <- cbind(1, Z); Y - Z1 %*% solve(crossprod(Z1), crossprod(Z1, Y)) }
  cor_mat <- function(A, B) { A <- scale(A); B <- scale(B); crossprod(A, B) / (nrow(A) - 1) }
  pool <- function(rmat, n, k) {                    # rmat: genes x cohorts
    z <- atanh(pmin(pmax(rmat, -0.999999), 0.999999)); w <- n - 3 - k
    zp <- as.vector(z %*% w) / sum(w); se <- 1 / sqrt(sum(w))
    data.frame(r = tanh(zp), p = 2 * pnorm(-abs(zp / se)))
  }
  prep <- list()
  for (k in names(co)) {
    d <- co[[k]]; m <- mats[[k]][, d$sample, drop = FALSE]
    lg <- log2(m + 1); lg <- lg[apply(lg, 1, var) > 0, , drop = FALSE]
    med <- apply(m[rownames(lg), , drop = FALSE], 1, median)
    R  <- apply(lg[names(med)[med >= MIN_MEDIAN], , drop = FALSE], 1, rank)
    R  <- cbind(R[, setdiff(colnames(R), "EZH2"), drop = FALSE], EZH2 = rank(d$EZH2))
    Tk <- sapply(CHECKPOINTS, function(t) rank(d[[t]]))
    S <- ng <- Zcc <- NULL
    if (!is.null(cc_genes)) { g <- intersect(cc_genes, rownames(lg)); Zcc <- scale(t(lg[g, , drop = FALSE])); S <- rowSums(Zcc); ng <- length(g) }
    Xpc <- scale(R[, setdiff(colnames(R), c("EZH2", excl)), drop = FALSE])
    Xpc <- Xpc[, colSums(is.na(Xpc)) == 0, drop = FALSE]
    pcs <- prcomp(Xpc, center = FALSE, scale. = FALSE, rank. = 10)
    ve  <- pcs$sdev^2 / sum(pcs$sdev^2)
    prep[[k]] <- list(d = d, R = R, Tk = Tk, Zcc = Zcc, S = S, ng = ng, PC = pcs$x[, 1:10], ve = ve[1:10],
                      rho = drop(cor(R, rank(d$prolif))),
                      rho_imm = drop(cor(R, rank(d$ESTIMATE_immune))),
                      rho_str = drop(cor(R, rank(d$ESTIMATE_stromal))),
                      rho_T   = drop(cor(R, rank(d$MCP_T_cell))))
  }
  n_c <- sapply(prep, function(p) nrow(p$d))
  genes <- Reduce(intersect, lapply(prep, function(p) colnames(p$R)))
  getm <- function(f) { x <- sapply(prep, function(p) p[[f]][genes]); rownames(x) <- genes; x }
  rho_mean <- rowMeans(getm("rho"))
  strict_non_immune <- apply(abs(getm("rho_imm")) < 0.1 & abs(getm("rho_str")) < 0.1 & abs(getm("rho_T")) < 0.1, 1, all)
  say("Expressed genes common to all cohorts: %d (strictly non-immune: %d)", length(genes), sum(strict_non_immune))
  for (k in names(prep)) say("%-9s variance explained by PC1-PC10: %s (cumulative %.1f%%)", k,
                             paste(sprintf("%.1f", 100 * prep[[k]]$ve), collapse = " "), 100 * sum(prep[[k]]$ve))
  for (k in names(prep)) say("%-9s |rho| of PC1-PC5 with prolif / ESTIMATE immune / MCP T cells: %s", k,
    paste(sapply(1:5, function(j) sprintf("PC%d %.2f/%.2f/%.2f", j,
      abs(cor(prep[[k]]$PC[, j], prep[[k]]$d$prolif, method = "spearman")),
      abs(cor(prep[[k]]$PC[, j], prep[[k]]$d$ESTIMATE_immune, method = "spearman")),
      abs(cor(prep[[k]]$PC[, j], prep[[k]]$d$MCP_T_cell, method = "spearman")))), collapse = "; "))

  # genes x checkpoints x cohorts partial correlations
  partial_all <- function(gs, Zvars, cc = FALSE, npc = 0) {
    arr <- array(NA_real_, c(length(gs), length(CHECKPOINTS), length(co)), dimnames = list(gs, CHECKPOINTS, names(co)))
    for (k in names(co)) {
      P <- prep[[k]]; d <- P$d
      Zr <- if (length(Zvars)) sapply(Zvars, function(z) rank(d[[z]])) else matrix(0, nrow(d), 0)
      if (npc > 0) Zr <- cbind(Zr, apply(P$PC[, seq_len(npc), drop = FALSE], 2, rank))
      if (cc) Zr <- cbind(Zr, cc = rank(P$S / P$ng))
      arr[, , k] <- cor_mat(resid_cols(P$R[, gs, drop = FALSE], Zr), resid_cols(P$Tk, Zr))
      if (cc) for (g in intersect(gs, colnames(P$Zcc))) {
        Zg <- cbind(Zr[, setdiff(colnames(Zr), "cc"), drop = FALSE], cc = rank((P$S - P$Zcc[, g]) / (P$ng - 1)))
        arr[g, , k] <- cor_mat(resid_cols(P$R[, g, drop = FALSE], Zg), resid_cols(P$Tk, Zg))
      }
    }
    arr
  }
  chk <- partial_all("EZH2", M3)
  for (k in names(co)) stopifnot(abs(chk["EZH2", "LAG3", k] - pspear("EZH2", "LAG3", M3, co[[k]])[["r"]]) < 1e-8)
  say("Vectorised estimates match pspear() for EZH2 (|difference| < 1e-8).")

  MODELS <- list(M0 = list(Z = character(0)), M3 = list(Z = M3))
  if (!is.null(cc_genes)) MODELS$M3h <- list(Z = M3, cc = TRUE)
  MODELS$M4 <- list(Z = M4); MODELS$`M4+PC5` <- list(Z = M4, npc = 5); MODELS$`M4+PC10` <- list(Z = M4, npc = 10)
  pooled <- list()
  for (mn in names(MODELS)) {
    s <- MODELS[[mn]]; npc <- if (is.null(s$npc)) 0 else s$npc
    arr <- partial_all(genes, s$Z, isTRUE(s$cc), npc)
    kk <- length(s$Z) + isTRUE(s$cc) + npc
    pooled[[mn]] <- lapply(CHECKPOINTS, function(t) { x <- pool(arr[, t, ], n_c, kk); rownames(x) <- genes; x })
    names(pooled[[mn]]) <- CHECKPOINTS
  }
  ref <- genes[!(genes %in% c("EZH2", excl))]                 # genome-wide reference genes
  ref_hi <- ref[rho_mean[ref] >= 0.75]                        # strongly proliferation-coupled
  ref_ni <- ref[strict_non_immune[ref]]                       # strictly non-immune
  say("Reference genes: all %d | coupling >= 0.75: %d | strictly non-immune: %d", length(ref), length(ref_hi), length(ref_ni))

  # =====================================================================
  # E1. EZH2 against the genome-wide null
  # =====================================================================
  hdr("E1. EZH2 rank among all expressed genes (genome-wide null)")
  rank_row <- function(x, set, e) c(median = median(x[set, "r"]), rank = 1 + sum(x[set, "r"] > e),
                                    n = length(set), emp_p = (1 + sum(x[set, "r"] >= e)) / (1 + length(set)))
  rowsE1 <- list()
  for (mn in names(MODELS)) for (t in CHECKPOINTS) {
    x <- pooled[[mn]][[t]]; e <- x["EZH2", "r"]
    a <- rank_row(x, ref, e); h <- rank_row(x, ref_hi, e); q <- rank_row(x, ref_ni, e)
    rowsE1[[length(rowsE1) + 1]] <- data.frame(model = mn, target = t, EZH2_r = e, EZH2_p = x["EZH2", "p"],
      genome_median = a[["median"]], genome_rank = a[["rank"]], genome_n = a[["n"]], genome_emp_p = a[["emp_p"]],
      genome_percentile = 100 * (1 - (a[["rank"]] - 1) / a[["n"]]),
      coupled_median = h[["median"]], coupled_rank = h[["rank"]], coupled_n = h[["n"]], coupled_emp_p = h[["emp_p"]],
      nonimmune_median = q[["median"]], nonimmune_emp_p = q[["emp_p"]])
  }
  E1 <- do.call(rbind, rowsE1)
  write.csv(E1, file.path(RES_DIR, "Table_S10a_EZH2_genomewide_rank.csv"), row.names = FALSE)
  for (t in c("LAG3", "PDCD1", "CD276")) {
    cat("--- ", t, " ---\n", sep = "")
    x <- subset(E1, target == t)
    print(data.frame(model = x$model, EZH2_r = r3(x$EZH2_r), P = fmt_p(x$EZH2_p), genome_median = r3(x$genome_median),
                     rank = paste0(x$genome_rank, "/", x$genome_n), emp_p = r3(x$genome_emp_p),
                     coupled_median = r3(x$coupled_median), coupled_rank = paste0(x$coupled_rank, "/", x$coupled_n),
                     coupled_p = r3(x$coupled_emp_p), nonimm_median = r3(x$nonimmune_median), nonimm_p = r3(x$nonimmune_emp_p)),
          row.names = FALSE)
  }
  cat("--- Genome-wide empirical P for EZH2, all checkpoints (rows = model) ---\n")
  wp <- reshape(E1[, c("model", "target", "genome_emp_p")], idvar = "model", timevar = "target", direction = "wide")
  names(wp) <- sub("genome_emp_p.", "", names(wp)); print(wp, row.names = FALSE, digits = 2)
  pM3 <- subset(E1, model == "M3"); pM3$BH <- p.adjust(pM3$genome_emp_p, "BH")
  say("M3, 8 checkpoints: smallest genome-wide empirical P = %.3f (%s); after BH over 8 checkpoints the smallest q = %.3f",
      min(pM3$genome_emp_p), pM3$target[which.min(pM3$genome_emp_p)], min(pM3$BH))

  # =====================================================================
  # E2/E3. Background skew: positive vs negative nominal associations
  # =====================================================================
  hdr("E2/E3. Background: share of genes with positive vs negative nominal associations")
  rowsE2 <- list()
  for (mn in names(MODELS)) for (t in CHECKPOINTS) for (sub in c("all genes", "strictly non-immune")) {
    x <- pooled[[mn]][[t]]; g <- if (sub == "all genes") ref else ref_ni
    rowsE2[[length(rowsE2) + 1]] <- data.frame(model = mn, target = t, gene_set = sub, n_genes = length(g),
      median_r = median(x[g, "r"]), frac_pos = mean(x[g, "r"] > 0 & x[g, "p"] < 0.05),
      frac_neg = mean(x[g, "r"] < 0 & x[g, "p"] < 0.05),
      frac_pos_BH = mean(x[g, "r"] > 0 & p.adjust(x$p[match(ref, rownames(x))], "BH")[match(g, ref)] < 0.05))
  }
  E2 <- do.call(rbind, rowsE2)
  write.csv(E2, file.path(RES_DIR, "Table_S10b_background_skew.csv"), row.names = FALSE)
  for (sub in c("all genes", "strictly non-immune")) {
    cat("--- ", sub, ": median r | positive share | negative share (nominal P < 0.05; 0.025 each expected under no association) ---\n", sep = "")
    x <- subset(E2, gene_set == sub)
    x$cell <- sprintf("%+.2f | %.2f | %.2f", x$median_r, x$frac_pos, x$frac_neg)
    w <- reshape(x[, c("target", "model", "cell")], idvar = "target", timevar = "model", direction = "wide")
    names(w) <- sub("cell.", "", names(w)); print(w, row.names = FALSE)
  }

  # =====================================================================
  # E4 summary + figure
  # =====================================================================
  hdr("E4. Summary for LAG3 and PDCD1")
  for (t in c("LAG3", "PDCD1")) for (mn in names(MODELS)) {
    a <- subset(E1, model == mn & target == t); b <- subset(E2, model == mn & target == t & gene_set == "all genes")
    say("%-6s %-8s EZH2 r = %6.3f (P %s) | genome median %6.3f | genome rank %5d/%d, emp P %.3f | background pos/neg %.2f/%.2f",
        t, mn, a$EZH2_r, fmt_p(a$EZH2_p), a$genome_median, a$genome_rank, a$genome_n, a$genome_emp_p, b$frac_pos, b$frac_neg)
  }

  mod_lab <- c(M0 = "M0 unadjusted", M3 = "M3 prolif + ESTIMATE", M3h = "M3h + Hallmark", M4 = "M4 + T cells",
               `M4+PC5` = "M4 + 5 PCs", `M4+PC10` = "M4 + 10 PCs")
  show_m <- intersect(c("M3", "M4", "M4+PC10"), names(MODELS))
  hd <- do.call(rbind, lapply(show_m, function(mn) do.call(rbind, lapply(c("LAG3", "PDCD1"), function(t)
    data.frame(model = mod_lab[[mn]], target = t, r = pooled[[mn]][[t]][ref, "r"])))))
  hd$model <- factor(hd$model, levels = mod_lab[show_m])
  ez <- subset(E1, model %in% show_m & target %in% c("LAG3", "PDCD1"))
  ez$model <- factor(mod_lab[ez$model], levels = mod_lab[show_m])
  ez$lab <- sprintf("EZH2 %.2f\nrank %d/%d", ez$EZH2_r, ez$genome_rank, ez$genome_n)
  pa <- ggplot(hd, aes(r)) + geom_histogram(bins = 60, fill = "grey70", colour = NA) +
    geom_vline(xintercept = 0, linetype = 3, colour = "grey40") +
    geom_vline(data = ez, aes(xintercept = EZH2_r), colour = "#C0392B", linewidth = 0.8) +
    geom_label(data = ez, aes(x = EZH2_r, y = Inf, label = lab), inherit.aes = FALSE, vjust = 1.1, hjust = -0.05,
               size = 2.7, label.size = 0, fill = "white", colour = "#C0392B") +
    facet_grid(target ~ model) +
    labs(x = "Pooled partial Spearman r of each expressed gene with the checkpoint", y = "Genes",
         title = "a  Genome-wide null: every expressed gene tested like EZH2") +
    theme_pub + theme(strip.text.y = element_text(face = "bold.italic"))
  sk <- subset(E2, gene_set == "all genes" & model %in% show_m)
  sk <- rbind(data.frame(sk[, c("model", "target")], dir = "positive", share = sk$frac_pos),
              data.frame(sk[, c("model", "target")], dir = "negative", share = sk$frac_neg))
  sk$model <- factor(mod_lab[sk$model], levels = mod_lab[show_m])
  pb <- ggplot(sk, aes(target, share, colour = dir, shape = dir)) +
    geom_hline(yintercept = 0.025, linetype = 3, colour = "grey50") + geom_point(size = 2.2) +
    facet_wrap(~model, nrow = 1) + scale_colour_manual(values = c(positive = "#C0392B", negative = "#2C3E50")) +
    labs(x = NULL, y = "Share of genes, P < 0.05", colour = "Direction", shape = "Direction",
         title = "b  Background: share of all genes with a nominally significant association",
         subtitle = "Dotted line: 0.025 per direction, expected if genes were unrelated to the checkpoint") +
    theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), legend.position = "bottom")
  save_fig(pa / pb + plot_layout(heights = c(1.15, 1)), "Figure4_genomewide_null", 11, 9)
})

# =====================================================================
# 4c Positive controls, spike-in power, audit of published claims
# =====================================================================
step("4c Positive controls, spike-in power, audit of published claims")
local({
  set.seed(2026)
  suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
  say <- function(...) cat(sprintf(...), "\n", sep = "")
  hdr <- function(x) cat("\n\n==================== ", x, " ====================\n", sep = "")
  r3  <- function(x) round(x, 3)
  co  <- readRDS(file.path(RES_DIR, "cohorts.rds"))$cohorts
  EST <- c("ESTIMATE_immune", "ESTIMATE_stromal"); M3 <- c("prolif", EST); M4 <- c(M3, "MCP_T_cell")
  MIN_MEDIAN <- 1
  excl <- unique(c(PROLIF_GENES, CHECKPOINTS))

  # ---------------- expression matrices (same handling as scripts 07/08/11) ----------------
  collapse_symbols <- function(mat, sym) {
    ok <- !is.na(sym) & sym != ""; mat <- mat[ok, , drop = FALSE]; sym <- sym[ok]
    o <- order(rowMeans(mat), decreasing = TRUE); mat <- mat[o, , drop = FALSE]; sym <- sym[o]
    keep <- !duplicated(sym); mat <- mat[keep, , drop = FALSE]; rownames(mat) <- sym[keep]; mat
  }
  tpm <- readRDS(file.path(DATA_DIR, "tcga_tpm.rds")); gi <- readRDS(file.path(DATA_DIR, "tcga_gene_info.rds"))
  mats <- list(TCGA = collapse_symbols(as.matrix(tpm[, co$TCGA$sample, drop = FALSE]), gi$gene_name[match(rownames(tpm), rownames(gi))]))
  rm(tpm); invisible(gc())
  for (id in c("325", "693")) {
    ex <- readRDS(file.path(DATA_DIR, paste0("cgga", id, "_expr_gbm.rds")))
    m <- as.matrix(ex[, -1]); storage.mode(m) <- "numeric"
    mats[[paste0("CGGA_", id)]] <- collapse_symbols(m, as.character(ex[[1]]))
  }
  cc_genes <- tryCatch({
    h <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
                  error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
    setdiff(unique(h$gene_symbol[h$gs_name %in% c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")]), c("EZH2", CHECKPOINTS))
  }, error = function(e) NULL)

  # ---------------- per-cohort preparation ----------------
  resid_cols <- function(Y, Z) { Z1 <- cbind(1, Z); Y - Z1 %*% solve(crossprod(Z1), crossprod(Z1, Y)) }
  cor_mat <- function(A, B) { A <- scale(A); B <- scale(B); crossprod(A, B) / (nrow(A) - 1) }
  pool <- function(rmat, n, k) {                    # rmat: genes x cohorts
    z <- atanh(pmin(pmax(rmat, -0.999999), 0.999999)); w <- n - 3 - k
    zp <- as.vector(z %*% w) / sum(w); se <- 1 / sqrt(sum(w))
    data.frame(r = tanh(zp), p = 2 * pnorm(-abs(zp / se)))
  }
  prep <- list()
  for (k in names(co)) {
    d <- co[[k]]; m <- mats[[k]][, d$sample, drop = FALSE]
    lg <- log2(m + 1); lg <- lg[apply(lg, 1, var) > 0, , drop = FALSE]
    med <- apply(m[rownames(lg), , drop = FALSE], 1, median)
    R  <- apply(lg[names(med)[med >= MIN_MEDIAN], , drop = FALSE], 1, rank)
    R  <- cbind(R[, setdiff(colnames(R), "EZH2"), drop = FALSE], EZH2 = rank(d$EZH2))
    Tk <- sapply(CHECKPOINTS, function(t) rank(d[[t]]))
    S <- ng <- Zcc <- NULL
    if (!is.null(cc_genes)) { g <- intersect(cc_genes, rownames(lg)); Zcc <- scale(t(lg[g, , drop = FALSE])); S <- rowSums(Zcc); ng <- length(g) }
    Xpc <- scale(R[, setdiff(colnames(R), c("EZH2", excl)), drop = FALSE])
    Xpc <- Xpc[, colSums(is.na(Xpc)) == 0, drop = FALSE]
    pcs <- prcomp(Xpc, center = FALSE, scale. = FALSE, rank. = 10)
    ve  <- pcs$sdev^2 / sum(pcs$sdev^2)
    prep[[k]] <- list(d = d, R = R, Tk = Tk, Zcc = Zcc, S = S, ng = ng, PC = pcs$x[, 1:10], ve = ve[1:10],
                      rho = drop(cor(R, rank(d$prolif))),
                      rho_imm = drop(cor(R, rank(d$ESTIMATE_immune))),
                      rho_str = drop(cor(R, rank(d$ESTIMATE_stromal))),
                      rho_T   = drop(cor(R, rank(d$MCP_T_cell))))
  }
  n_c <- sapply(prep, function(p) nrow(p$d))
  genes <- Reduce(intersect, lapply(prep, function(p) colnames(p$R)))
  getm <- function(f) { x <- sapply(prep, function(p) p[[f]][genes]); rownames(x) <- genes; x }
  rho_mean <- rowMeans(getm("rho"))
  strict_non_immune <- apply(abs(getm("rho_imm")) < 0.1 & abs(getm("rho_str")) < 0.1 & abs(getm("rho_T")) < 0.1, 1, all)
  say("Expressed genes common to all cohorts: %d (strictly non-immune: %d)", length(genes), sum(strict_non_immune))
  for (k in names(prep)) say("%-9s variance explained by PC1-PC10: %s (cumulative %.1f%%)", k,
                             paste(sprintf("%.1f", 100 * prep[[k]]$ve), collapse = " "), 100 * sum(prep[[k]]$ve))
  for (k in names(prep)) say("%-9s |rho| of PC1-PC5 with prolif / ESTIMATE immune / MCP T cells: %s", k,
    paste(sapply(1:5, function(j) sprintf("PC%d %.2f/%.2f/%.2f", j,
      abs(cor(prep[[k]]$PC[, j], prep[[k]]$d$prolif, method = "spearman")),
      abs(cor(prep[[k]]$PC[, j], prep[[k]]$d$ESTIMATE_immune, method = "spearman")),
      abs(cor(prep[[k]]$PC[, j], prep[[k]]$d$MCP_T_cell, method = "spearman")))), collapse = "; "))

  # genes x checkpoints x cohorts partial correlations
  partial_all <- function(gs, Zvars, cc = FALSE, npc = 0) {
    arr <- array(NA_real_, c(length(gs), length(CHECKPOINTS), length(co)), dimnames = list(gs, CHECKPOINTS, names(co)))
    for (k in names(co)) {
      P <- prep[[k]]; d <- P$d
      Zr <- if (length(Zvars)) sapply(Zvars, function(z) rank(d[[z]])) else matrix(0, nrow(d), 0)
      if (npc > 0) Zr <- cbind(Zr, apply(P$PC[, seq_len(npc), drop = FALSE], 2, rank))
      if (cc) Zr <- cbind(Zr, cc = rank(P$S / P$ng))
      arr[, , k] <- cor_mat(resid_cols(P$R[, gs, drop = FALSE], Zr), resid_cols(P$Tk, Zr))
      if (cc) for (g in intersect(gs, colnames(P$Zcc))) {
        Zg <- cbind(Zr[, setdiff(colnames(Zr), "cc"), drop = FALSE], cc = rank((P$S - P$Zcc[, g]) / (P$ng - 1)))
        arr[g, , k] <- cor_mat(resid_cols(P$R[, g, drop = FALSE], Zg), resid_cols(P$Tk, Zg))
      }
    }
    arr
  }
  chk <- partial_all("EZH2", M3)
  for (k in names(co)) stopifnot(abs(chk["EZH2", "LAG3", k] - pspear("EZH2", "LAG3", M3, co[[k]])[["r"]]) < 1e-8)
  say("Vectorised estimates match pspear() for EZH2 (|difference| < 1e-8).")


  MODELS <- list(M0 = list(Z = character(0)), M3 = list(Z = M3), M4 = list(Z = M4), `M4+PC10` = list(Z = M4, npc = 10))
  npc_of <- function(s) if (is.null(s$npc)) 0 else s$npc
  kk_of  <- function(s) length(s$Z) + npc_of(s)
  Zmat <- function(k, s) {
    d <- prep[[k]]$d
    Zr <- if (length(s$Z)) sapply(s$Z, function(z) rank(d[[z]])) else matrix(0, nrow(d), 0)
    if (npc_of(s) > 0) Zr <- cbind(Zr, apply(prep[[k]]$PC[, seq_len(npc_of(s)), drop = FALSE], 2, rank))
    Zr
  }
  # genome-wide reference: pooled r and p, genes x checkpoints, per model
  ref <- genes[!(genes %in% c("EZH2", excl))]
  Rg <- Pg <- list()
  for (mn in names(MODELS)) {
    s <- MODELS[[mn]]; arr <- partial_all(genes, s$Z, FALSE, npc_of(s))
    Rg[[mn]] <- Pg[[mn]] <- matrix(NA_real_, length(genes), length(CHECKPOINTS), dimnames = list(genes, CHECKPOINTS))
    for (t in CHECKPOINTS) { x <- pool(arr[, t, ], n_c, kk_of(s)); Rg[[mn]][, t] <- x$r; Pg[[mn]][, t] <- x$p }
  }
  say("Genome-wide reference genes: %d", length(ref))
  coup_P <- rowMeans(getm("rho")); coup_I <- rowMeans(getm("rho_imm"))

  # ---------------- test-gene machinery ----------------
  test_vec <- function(g, k) {
    d <- prep[[k]]$d
    if (g %in% colnames(prep[[k]]$R)) return(prep[[k]]$R[, g])
    if (g %in% CHECKPOINTS) return(rank(d[[g]]))
    m <- mats[[k]]; if (!g %in% rownames(m)) return(NULL)
    rank(log2(m[g, d$sample] + 1))
  }
  test_median <- function(g, k) {
    d <- prep[[k]]$d
    if (g == "EZH2" || g %in% CHECKPOINTS) return(median(d[[g]]))
    m <- mats[[k]]; if (!g %in% rownames(m)) return(NA_real_); median(m[g, d$sample])
  }
  test_pool <- function(xl, s) {             # xl: list of rank vectors per cohort -> data.frame per checkpoint
    r <- sapply(names(co), function(k) { Zr <- Zmat(k, s); drop(cor_mat(resid_cols(matrix(xl[[k]]), Zr), resid_cols(prep[[k]]$Tk, Zr))) })
    if (is.null(dim(r))) r <- matrix(r, nrow = 1)
    x <- pool(r, n_c, kk_of(s)); rownames(x) <- CHECKPOINTS; x
  }
  rank_against <- function(e, rr) c(median = median(rr), rank = 1 + sum(rr > e), n = length(rr), emp_p = (1 + sum(rr >= e)) / (1 + length(rr)))
  nearest <- function(cp, ci, pool_genes, n = 100) {
    dd <- sqrt((coup_P[pool_genes] - cp)^2 + (coup_I[pool_genes] - ci)^2)
    pool_genes[order(dd)][seq_len(min(n, length(pool_genes)))]
  }

  POS <- c("CD8A", "CD8B", "GZMK", "GZMB", "PRF1", "NKG7", "CXCL13", "TOX", "IFNG", "CD2")
  AUDIT <- data.frame(
    gene   = c("EZH2", "MELK", "CENPA", "FBP1", "WTAP", "ZC3H13", "IFI44", "FCGR3A"),
    source = c("Chen 2021 (PMID 33277782)", "Yang 2022 (PMID 36353126)", "Ling 2024 (PMID 39620895)", "Sun 2023 (PMID 37507483)",
               "Gao 2025 (PMID 39910141)", "Gao 2025 (PMID 39910141)", "Jin 2025 (PMID 41327009)", "Sun 2022, LGG (PMID 39280892)"),
    stringsAsFactors = FALSE)
  TESTS <- rbind(data.frame(group = "positive control (T-cell gene)", gene = POS, source = "", stringsAsFactors = FALSE),
                 data.frame(group = "positive control (checkpoint)", gene = CHECKPOINTS, source = "", stringsAsFactors = FALSE),
                 data.frame(group = "published claim", gene = AUDIT$gene, source = AUDIT$source, stringsAsFactors = FALSE))
  rowsF <- list(); skipped <- character(0)
  for (i in seq_len(nrow(TESTS))) {
    g <- TESTS$gene[i]
    xl <- lapply(names(co), function(k) test_vec(g, k)); names(xl) <- names(co)
    if (any(vapply(xl, is.null, logical(1))) || any(vapply(xl, function(v) var(v) == 0, logical(1)))) { skipped <- c(skipped, g); next }
    cp <- mean(sapply(names(co), function(k) cor(xl[[k]], rank(prep[[k]]$d$prolif))))
    ci <- mean(sapply(names(co), function(k) cor(xl[[k]], rank(prep[[k]]$d$ESTIMATE_immune))))
    ref_g <- setdiff(ref, g); ctrl <- nearest(cp, ci, ref_g)
    med <- sapply(names(co), function(k) test_median(g, k))
    for (mn in names(MODELS)) {
      x <- test_pool(xl, MODELS[[mn]])
      for (t in setdiff(CHECKPOINTS, g)) {
        e <- x[t, "r"]; a <- rank_against(e, Rg[[mn]][ref_g, t]); h <- rank_against(e, Rg[[mn]][ctrl, t])
        rowsF[[length(rowsF) + 1]] <- data.frame(group = TESTS$group[i], gene = g, source = TESTS$source[i],
          median_expr = paste(sprintf("%.2f", med), collapse = "/"), coupling_prolif = cp, coupling_immune = ci,
          model = mn, target = t, r = e, p = x[t, "p"],
          genome_median = a[["median"]], genome_rank = a[["rank"]], genome_n = a[["n"]], genome_emp_p = a[["emp_p"]],
          genome_percentile = 100 * (1 - (a[["rank"]] - 1) / a[["n"]]),
          matched_median = h[["median"]], matched_rank = h[["rank"]], matched_n = h[["n"]], matched_emp_p = h[["emp_p"]])
      }
    }
  }
  F <- do.call(rbind, rowsF)
  F$standard_positive <- F$r > 0 & F$p < 0.05
  F$genome_specific <- F$genome_emp_p < 0.05
  F$matched_specific <- F$matched_emp_p < 0.05
  write.csv(F, file.path(RES_DIR, "Table_S11a_positive_controls_and_audit.csv"), row.names = FALSE)
  if (length(skipped)) say("Not tested (absent or constant in a cohort): %s", paste(skipped, collapse = ", "))

  show <- function(df) print(data.frame(gene = df$gene, model = df$model, r = r3(df$r), P = fmt_p(df$p),
    genome = paste0(df$genome_rank, "/", df$genome_n), genome_p = r3(df$genome_emp_p), pct = round(df$genome_percentile, 1),
    matched = paste0(df$matched_rank, "/", df$matched_n), matched_p = r3(df$matched_emp_p)), row.names = FALSE)

  hdr("F1. Positive controls")
  cat("Median expression per cohort (TCGA TPM / CGGA_325 FPKM / CGGA_693 FPKM) and coupling:\n")
  u <- unique(F[F$group != "published claim", c("gene", "median_expr", "coupling_prolif", "coupling_immune")])
  print(transform(u, coupling_prolif = r3(coupling_prolif), coupling_immune = r3(coupling_immune)), row.names = FALSE)
  for (t in c("LAG3", "PDCD1")) {
    cat("--- target ", t, " (T-cell genes and other checkpoints; models M3 and M4) ---\n", sep = "")
    show(subset(F, target == t & group != "published claim" & model %in% c("M3", "M4")))
  }
  pc <- subset(F, group == "positive control (T-cell gene)" & target %in% c("LAG3", "PDCD1"))
  for (mn in names(MODELS)) { x <- subset(pc, model == mn)
    say("%-8s T-cell positive controls vs LAG3/PDCD1: %d of %d pairs genome-specific (emp P < 0.05); %d of %d matched-specific; median genome percentile %.1f",
        mn, sum(x$genome_specific), nrow(x), sum(x$matched_specific), nrow(x), median(x$genome_percentile)) }
  ck <- subset(F, group == "positive control (checkpoint)")
  for (mn in names(MODELS)) { x <- subset(ck, model == mn)
    say("%-8s checkpoint-checkpoint pairs: %d of %d genome-specific; %d of %d matched-specific (CTLA4 pairs included; CTLA4 is an MCP T-cell marker)",
        mn, sum(x$genome_specific), nrow(x), sum(x$matched_specific), nrow(x)) }

  hdr("F2. Audit of published glioma gene-immune claims")
  u <- unique(F[F$group == "published claim", c("gene", "source", "median_expr", "coupling_prolif", "coupling_immune")])
  print(transform(u, coupling_prolif = r3(coupling_prolif), coupling_immune = r3(coupling_immune)), row.names = FALSE)
  au <- subset(F, group == "published claim")
  sumA <- do.call(rbind, lapply(split(au, list(au$gene, au$model), drop = TRUE), function(x) data.frame(gene = x$gene[1], model = x$model[1],
    standard_positive = sum(x$standard_positive), genome_specific = sum(x$standard_positive & x$genome_specific),
    matched_specific = sum(x$standard_positive & x$matched_specific),
    best_target = x$target[which.min(x$matched_emp_p)], best_matched_p = min(x$matched_emp_p))))
  sumA <- sumA[order(match(sumA$gene, AUDIT$gene), match(sumA$model, names(MODELS))), ]
  write.csv(sumA, file.path(RES_DIR, "Table_S11b_audit_summary.csv"), row.names = FALSE)
  cat("--- per gene and model: checkpoints (of 8) passing the standard test, and of those, specific vs genome-wide null / matched controls ---\n")
  print(transform(sumA, best_matched_p = r3(best_matched_p)), row.names = FALSE)
  for (mn in c("M3", "M4")) { x <- subset(au, model == mn)
    say("%s: %d of %d audit gene-checkpoint pairs pass the standard test; %d genome-specific; %d matched-specific; %d both",
        mn, sum(x$standard_positive), nrow(x), sum(x$standard_positive & x$genome_specific),
        sum(x$standard_positive & x$matched_specific), sum(x$standard_positive & x$genome_specific & x$matched_specific)) }
  for (t in c("LAG3", "PDCD1")) { cat("--- target ", t, " (published-claim genes; M3 and M4) ---\n", sep = ""); show(subset(au, target == t & model %in% c("M3", "M4"))) }

  hdr("F3. Spike-in power: smallest specific association detected by the genome-wide null")
  rowsS <- list()
  cpE <- mean(sapply(names(co), function(k) cor(prep[[k]]$R[, "EZH2"], rank(prep[[k]]$d$prolif))))
  ciE <- mean(sapply(names(co), function(k) cor(prep[[k]]$R[, "EZH2"], rank(prep[[k]]$d$ESTIMATE_immune))))
  ctrlE <- nearest(cpE, ciE, ref)
  for (mn in c("M3", "M4")) for (t in c("LAG3", "PDCD1")) for (dl in seq(0, 1.5, by = 0.025)) {
    s <- MODELS[[mn]]
    xl <- lapply(names(co), function(k) { Zr <- Zmat(k, s); P <- prep[[k]]
      spec <- drop(resid_cols(matrix(P$Tk[, t]), Zr))
      rank(drop(scale(P$R[, "EZH2"])) + dl * drop(scale(spec))) })
    names(xl) <- names(co)
    x <- test_pool(xl, s); e <- x[t, "r"]
    a <- rank_against(e, Rg[[mn]][ref, t]); h <- rank_against(e, Rg[[mn]][ctrlE, t])
    rowsS[[length(rowsS) + 1]] <- data.frame(model = mn, target = t, delta = dl, r = e, p = x[t, "p"],
      genome_median = a[["median"]], genome_emp_p = a[["emp_p"]], matched_emp_p = h[["emp_p"]])
  }
  S <- do.call(rbind, rowsS)
  write.csv(S, file.path(RES_DIR, "Table_S11c_spike_in_power.csv"), row.names = FALSE)
  for (mn in c("M3", "M4")) for (t in c("LAG3", "PDCD1")) {
    x <- subset(S, model == mn & target == t)
    f05 <- x[x$genome_emp_p < 0.05, ]; f01 <- x[x$genome_emp_p < 0.01, ]; m05 <- x[x$matched_emp_p < 0.05, ]
    say("%-3s %-6s EZH2 r without spike %.3f (genome median %.3f) | genome P < 0.05 from r = %s | < 0.01 from r = %s | matched P < 0.05 from r = %s",
        mn, t, x$r[1], x$genome_median[1],
        if (nrow(f05)) sprintf("%.3f", f05$r[1]) else "not reached", if (nrow(f01)) sprintf("%.3f", f01$r[1]) else "not reached",
        if (nrow(m05)) sprintf("%.3f", m05$r[1]) else "not reached")
  }

  # ---------------- figure ----------------
  lab_g <- function(df) factor(df$gene, levels = unique(c(POS, CHECKPOINTS, AUDIT$gene)))
  fa <- subset(F, target %in% c("LAG3", "PDCD1") & model %in% c("M3", "M4"))
  fa$gene_f <- lab_g(fa)
  fa$grp <- factor(ifelse(fa$group == "published claim", "Published claims",
                   ifelse(fa$group == "positive control (checkpoint)", "Other checkpoints", "T-cell genes")),
                   levels = c("T-cell genes", "Other checkpoints", "Published claims"))
  fa$model <- factor(fa$model, levels = c("M3", "M4"), labels = c("M3 prolif + ESTIMATE", "M4 = M3 + T cells"))
  pa <- ggplot(fa, aes(gene_f, genome_percentile, colour = model, shape = model)) +
    geom_hline(yintercept = 95, linetype = 2, colour = "grey50") +
    geom_point(size = 2, position = position_dodge(width = 0.5)) +
    facet_grid(target ~ grp, scales = "free_x", space = "free_x") +
    scale_colour_manual(values = c("#C0392B", "#2C3E50")) + scale_y_continuous(limits = c(0, 100)) +
    labs(x = NULL, y = "Genome-wide percentile", colour = NULL, shape = NULL,
         title = "a  Positive controls and published claims against the genome-wide null",
         subtitle = "Dashed line: 95th percentile (empirical P = 0.05)") +
    theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), legend.position = "bottom",
                      strip.text.y = element_text(face = "bold.italic"))
  floor_p <- 1 / (1 + length(ref))
  Sp <- do.call(rbind, lapply(split(S, list(S$model, S$target), drop = TRUE), function(x) {
    hit <- which(x$genome_emp_p <= floor_p); x[seq_len(if (length(hit)) min(nrow(x), hit[1] + 2) else nrow(x)), ] }))
  pb <- ggplot(Sp, aes(r, genome_emp_p, colour = model, linetype = target)) +
    geom_hline(yintercept = 0.05, linetype = 3, colour = "grey50") + geom_line() + geom_point(size = 1.2) +
    scale_y_log10() + scale_colour_manual(values = c(M3 = "#C0392B", M4 = "#2C3E50")) +
    labs(x = "Pooled partial r of EZH2 after adding a known target-specific signal", y = "Genome-wide empirical P",
         colour = "Model", linetype = "Target", title = "b  Spike-in: how large a specific association the null detects") +
    theme_pub + theme(legend.position = "bottom")
  save_fig(pa / pb + plot_layout(heights = c(1.4, 1)), "Figure5_positive_controls_audit", 11, 10)
})

# =====================================================================
# Copy the figures under their article numbers (figures/article/)
# =====================================================================
step("Article figure numbering")
local({
  map <- c(Fig1 = "Figure3_checkpoints_partial_correlation", Fig2 = "Figure2_specificity_benchmark",
           Fig3 = "Figure4_genomewide_null", Fig4 = "FigureS6_calibration", Fig5 = "Figure5_positive_controls_audit",
           FigS1 = "Figure1_expression_PRC2_proliferation", FigS2 = "Figure2_survival", FigS3 = "FigureS1_immune_landscape",
           FigS4 = "FigureS4_forest_checkpoints", FigS5 = "Figure4_single_cell", FigS6 = "Figure3_simulation")
  out <- file.path(FIG_DIR, "article"); dir.create(out, showWarnings = FALSE)
  for (k in names(map)) for (ext in c("pdf", "png")) {
    src <- file.path(FIG_DIR, paste0(map[[k]], ".", ext))
    if (file.exists(src)) file.copy(src, file.path(out, paste0(k, ".", ext)), overwrite = TRUE)
    else cat("  not found (run the earlier steps):", src, "\n")
  }
  cat("  Article figures written to", out, "\n")
})

cat("\n\u2713 04_genomewide_controls_audit.R complete.\n")
}, finally = end_log())
