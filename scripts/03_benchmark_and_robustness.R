# =====================================================================
# 03_benchmark_and_robustness.R
# Step 3: matched-control benchmark, simulation, robustness and survival sensitivity, figures
# Part of: Covariate-adjusted co-expression does not establish gene-immune
# checkpoint specificity in glioblastoma (EZH2 case study).
#
# HOW TO RUN: set the working directory to the PROJECT folder (the folder
# that contains share/), then source this script from wherever it is, e.g.
#   setwd("D:/Bioinformatics/EZH2_GBM_Bioinfo_2026")
#   source("scripts/03_benchmark_and_robustness.R")
# Outputs go to results/ and figures/ in the project folder; a log is written
# to results/log_03_benchmark_and_robustness.txt.
# =====================================================================
.this <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/"), error = function(e) NA_character_)
source(file.path(if (is.na(.this)) getwd() else dirname(.this), "00_setup.R"))
check_project(c("results/cohorts.rds", "results/Table3a_checkpoint_partial_correlations_by_cohort.csv"))
start_log("log_03_benchmark_and_robustness.txt")
tryCatch({

# =====================================================================
# 3a Measurement-error simulation
# =====================================================================
step("3a Measurement-error simulation")
local({
  set.seed(2026)
  set.seed(2026)
  N_COH   <- c(TCGA = 160, CGGA_325 = 74, CGGA_693 = 109)
  RHO_OBS <- 0.90
  NREP    <- 2000
  say <- function(...) cat(sprintf(...), "\n", sep = "")

  # ---- optional: reliability of the real 10-gene proliferation score ----
  PROLIF_GENES <- c("MKI67","TOP2A","CCNB1","CDK1","BUB1","CENPF","AURKA","PCNA","MCM2","TPX2")
  lamP_data <- NA
  if (file.exists("results/cohorts.rds")) {
    co <- readRDS("results/cohorts.rds")$cohorts
    alpha <- sapply(co, function(d) {                # Cronbach alpha of the z-scored log2 genes
      X <- scale(log2(d[, PROLIF_GENES] + 1)); k <- ncol(X)
      k / (k - 1) * (1 - sum(apply(X, 2, var)) / var(rowSums(X)))
    })
    say("Cronbach alpha of the 10-gene score: %s", paste(sprintf("%s %.3f", names(alpha), alpha), collapse = " | "))
    lamP_data <- sqrt(mean(alpha))
    say("=> approximate correlation of the proxy with the latent factor: lamP = sqrt(alpha) = %.3f", lamP_data)
    say("   (approximate: alpha assumes the 10 genes share ONE factor; real cell-cycle genes also share other variation)\n")
  }

  # ---- fast rank-based partial correlation of x and y given z ----
  pcor_rank <- function(x, y, z) {
    Z <- cbind(1, rank(z))
    ex <- lm.fit(Z, rank(x))$residuals; ey <- lm.fit(Z, rank(y))$residuals
    cor(ex, ey)
  }
  pool_p <- function(r, n, k = 1) {                  # fixed-effect Fisher z, weights n - 3 - k
    w <- n - 3 - k; z <- atanh(r); zp <- sum(w * z) / sum(w); se <- 1 / sqrt(sum(w))
    c(r = tanh(zp), p = 2 * pnorm(-abs(zp / se)))
  }
  sim_cell <- function(lamP, gamma, beta = 0) {
    lamE <- min(0.999, RHO_OBS / lamP)
    out <- replicate(NREP, {
      r <- sapply(N_COH, function(n) {
        P <- rnorm(n)
        prox <- lamP * P + sqrt(1 - lamP^2) * rnorm(n)
        E    <- lamE * P + sqrt(1 - lamE^2) * rnorm(n)
        L    <- gamma * P + beta * E + rnorm(n)
        pcor_rank(E, L, prox)
      })
      pool_p(r, N_COH)
    })
    c(mean_pooled_r = mean(out["r", ]), reject_rate = mean(out["p", ] < 0.05), lamE = lamE)
  }
  analytic <- function(lamP, gamma) {                # Pearson analogue, for reference
    lamE <- min(0.999, RHO_OBS / lamP)
    gamma * lamE * (1 - lamP^2) / sqrt((1 - lamE^2 * lamP^2) * (1 - gamma^2 * lamP^2))
  }

  lamP_grid <- sort(unique(c(0.99, 0.97, 0.95, 0.93, 0.90, if (!is.na(lamP_data)) round(lamP_data, 3))), decreasing = TRUE)
  grid <- rbind(
    cbind(world = "A_confounding_only", expand.grid(lamP = lamP_grid, gamma = c(0.3, 0.5), beta = 0)),
    cbind(world = "B_true_direct_effect", expand.grid(lamP = lamP_grid, gamma = 0.3, beta = 0.3)))
  res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]; s <- sim_cell(g$lamP, g$gamma, g$beta)
    data.frame(g, mean_pooled_r = s[["mean_pooled_r"]], reject_rate = s[["reject_rate"]], lamE = s[["lamE"]],
               analytic_r_world_A = if (g$beta == 0) analytic(g$lamP, g$gamma) else NA, row.names = NULL)
  }))
  dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)
  write.csv(res, "results/Table_S6c_simulation_measurement_error.csv", row.names = FALSE)
  cat("==== Simulation: pooled adjusted EZH2-LAG3 partial r when EZH2 has NO direct effect (world A) ====\n")
  print(transform(res, mean_pooled_r = round(mean_pooled_r, 3), reject_rate = round(reject_rate, 3),
                  lamE = round(lamE, 3), analytic_r_world_A = round(analytic_r_world_A, 3)), row.names = FALSE)
  cat("\nHow to read: reject_rate in world A is the false-positive rate of the paper's test (nominal 0.05).\n",
      "Compare mean_pooled_r in world A with the observed pooled partial r (LAG3 0.29; PDCD1 0.24).\n",
      "If world A reaches that size at a plausible lamP, the adjusted association alone cannot prove specificity;\n",
      "this is why script 08 (negative-control genes) is needed.\n", sep = "")

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
    A <- subset(res, world == "A_confounding_only")
    p <- ggplot(A, aes(lamP, mean_pooled_r, colour = factor(gamma))) +
      geom_hline(yintercept = 0.29, linetype = 2, colour = "#C0392B") +
      annotate("text", x = min(A$lamP), y = 0.30, label = "observed pooled r (LAG3) = 0.29", hjust = 0, vjust = -0.4, size = 3, colour = "#C0392B") +
      geom_line() + geom_point() +
      scale_colour_manual(values = c(`0.3` = "#2C3E50", `0.5` = "#E67E22")) +
      labs(x = "Correlation of the proliferation proxy with latent proliferation (lamP)",
           y = "Mean pooled adjusted EZH2-LAG3 partial r\n(no direct effect simulated)", colour = "LAG3-proliferation\nloading (gamma)",
           title = "Spurious adjusted association left by a noisy proliferation proxy",
           subtitle = sprintf("cor(EZH2, proxy) fixed at %.2f; n = %s; %d replicates per point", RHO_OBS, paste(N_COH, collapse = "/"), NREP)) +
      theme_classic(base_size = 11)
    ggsave("figures/FigureS3_simulation.pdf", p, width = 7, height = 4.5)
    ggsave("figures/FigureS3_simulation.png", p, width = 7, height = 4.5, dpi = 300)
  }
  cat("\n✓ Table S6c saved\n")
})

# =====================================================================
# 3b Matched negative-control benchmark
# =====================================================================
step("3b Matched negative-control benchmark")
local({
  set.seed(2026)
  RHO_MIN    <- 0.75     # coupling to proliferation required of a comparator gene (EZH2: 0.88-0.91)
  MIN_MEDIAN <- 1        # comparator genes must have median TPM/FPKM >= 1 in every cohort
  MATCH_TOL  <- 0.05     # "matched" comparators: mean rho within +/- this of EZH2's mean rho
  TARGETS    <- c("LAG3", "PDCD1", "CD276")
  EST        <- c("ESTIMATE_immune", "ESTIMATE_stromal")
  say <- function(...) cat(sprintf(...), "\n", sep = "")

  co <- readRDS(file.path(RES_DIR, "cohorts.rds"))$cohorts

  # ---------------- expression matrices (same handling as script 07) ----------------
  collapse_symbols <- function(mat, sym) {
    ok <- !is.na(sym) & sym != ""; mat <- mat[ok, , drop = FALSE]; sym <- sym[ok]
    o <- order(rowMeans(mat), decreasing = TRUE); mat <- mat[o, , drop = FALSE]; sym <- sym[o]
    keep <- !duplicated(sym); mat <- mat[keep, , drop = FALSE]; rownames(mat) <- sym[keep]; mat
  }
  tpm <- readRDS(file.path(DATA_DIR, "tcga_tpm.rds")); gi <- readRDS(file.path(DATA_DIR, "tcga_gene_info.rds"))
  mats <- list(TCGA = collapse_symbols(as.matrix(tpm[, co$TCGA$sample, drop = FALSE]),
                                       gi$gene_name[match(rownames(tpm), rownames(gi))]))
  rm(tpm); invisible(gc())
  for (id in c("325", "693")) {
    ex <- readRDS(file.path(DATA_DIR, paste0("cgga", id, "_expr_gbm.rds")))
    m <- as.matrix(ex[, -1]); storage.mode(m) <- "numeric"
    mats[[paste0("CGGA_", id)]] <- collapse_symbols(m, as.character(ex[[1]]))
  }

  # ---------------- Hallmark cell-cycle genes (as in script 07) ----------------
  cc_genes <- tryCatch({
    h <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
                  error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
    setdiff(unique(h$gene_symbol[h$gs_name %in% c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")]),
            c("EZH2", CHECKPOINTS))
  }, error = function(e) { message("msigdbr unavailable -> model M3h skipped: ", conditionMessage(e)); NULL })

  # ---------------- vectorised partial Spearman ----------------
  resid_cols <- function(Y, Z) {                      # residuals of every column of Y on [1, Z]
    Z1 <- cbind(1, Z)
    Y - Z1 %*% solve(crossprod(Z1), crossprod(Z1, Y))
  }
  cor_mat <- function(A, B) { A <- scale(A); B <- scale(B); crossprod(A, B) / (nrow(A) - 1) }
  partial_block <- function(Rk, Tk, Zr) {            # Rk: ranks (samples x genes); Tk: ranks (samples x targets)
    if (is.null(Zr)) Zr <- matrix(0, nrow(Rk), 0)
    cor_mat(resid_cols(Rk, Zr), resid_cols(Tk, Zr))
  }

  # ---------------- per-cohort preparation ----------------
  need <- c("EZH2", TARGETS, "prolif", EST, "MCP_T_cell")
  prep <- list()
  for (k in names(co)) {
    d <- co[[k]]; m <- mats[[k]][, d$sample, drop = FALSE]
    stopifnot(all(need %in% colnames(d)), !anyNA(d[, need]))
    lg <- log2(m + 1); lg <- lg[apply(lg, 1, var) > 0, , drop = FALSE]
    med <- apply(m[rownames(lg), , drop = FALSE], 1, median)
    Zcc <- S <- ng <- NULL
    if (!is.null(cc_genes)) {                        # same definition as cc_meanz in script 07
      g <- intersect(cc_genes, rownames(lg))
      Zcc <- scale(t(lg[g, , drop = FALSE]))          # samples x genes, z-scored per gene
      S <- rowSums(Zcc); ng <- length(g)
      stopifnot(ng > 1, length(S) == ncol(lg))
    }
    keep <- names(med)[med >= MIN_MEDIAN]
    Rk <- apply(lg[keep, , drop = FALSE], 1, rank)   # samples x genes
    Rk <- cbind(Rk[, setdiff(colnames(Rk), "EZH2"), drop = FALSE], EZH2 = rank(d$EZH2))
    Tk <- sapply(TARGETS, function(t) rank(d[[t]]))
    prep[[k]] <- list(d = d, Rk = Rk, Tk = Tk, Zcc = Zcc, S = S, ng = ng,
                      rho = drop(cor(Rk, rank(d$prolif))))
  }
  n_c <- sapply(prep, function(p) nrow(p$d))

  # ---------------- comparator set ----------------
  genes_common <- Reduce(intersect, lapply(prep, function(p) colnames(p$Rk)))
  rho_mat <- sapply(prep, function(p) p$rho[genes_common]); rownames(rho_mat) <- genes_common
  excl <- unique(c("EZH2", PROLIF_GENES, CHECKPOINTS))
  cand <- genes_common[apply(rho_mat >= RHO_MIN, 1, all) & !(genes_common %in% excl)]
  say("Comparator genes (rho with proliferation >= %.2f in all cohorts, median expression >= %g): %d", RHO_MIN, MIN_MEDIAN, length(cand))
  say("EZH2 rho with proliferation: %s", paste(sprintf("%s %.2f", names(co), rho_mat["EZH2", ]), collapse = " | "))
  if (length(cand) < 30) warning("Fewer than 30 comparator genes: lower RHO_MIN (e.g. 0.70) for a usable null distribution.")
  genes_run <- c("EZH2", cand)

  # ---------------- models ----------------
  MODELS <- list(
    M0_raw                       = list(Z = character(0)),
    M3_prolif_ESTIMATE           = list(Z = c("prolif", EST)),
    M4_prolif_ESTIMATE_Tcells    = list(Z = c("prolif", EST, "MCP_T_cell")),
    M3h_prolif_Hallmark_ESTIMATE = list(Z = c("prolif", EST), cc = TRUE))
  if (is.null(cc_genes)) MODELS$M3h_prolif_Hallmark_ESTIMATE <- NULL

  run_model <- function(spec) {
    r_arr <- array(NA_real_, c(length(genes_run), length(TARGETS), length(co)),
                   dimnames = list(genes_run, TARGETS, names(co)))
    for (k in names(co)) {
      P <- prep[[k]]; d <- P$d
      Zr <- if (length(spec$Z)) sapply(spec$Z, function(z) rank(d[[z]])) else NULL
      if (isTRUE(spec$cc)) { stopifnot(!is.null(P$ng)); Zr <- cbind(Zr, cc = rank(P$S / P$ng)); stopifnot("cc" %in% colnames(Zr)) }
      r_arr[, , k] <- partial_block(P$Rk[, genes_run, drop = FALSE], P$Tk, Zr)
      if (isTRUE(spec$cc))                           # comparator genes inside the signature: remove them from it
        for (g in intersect(genes_run, colnames(P$Zcc))) {
          Zg <- cbind(Zr[, setdiff(colnames(Zr), "cc"), drop = FALSE], cc = rank((P$S - P$Zcc[, g]) / (P$ng - 1)))
          r_arr[g, , k] <- partial_block(P$Rk[, g, drop = FALSE], P$Tk, Zg)
        }
    }
    r_arr
  }
  pool <- function(r3, n, k) {                       # Fisher-z fixed effect, weights n - 3 - k (as meta_r in script 00)
    z <- atanh(pmin(pmax(r3, -0.999999), 0.999999)); w <- n - 3 - k
    zp <- as.vector(z %*% w) / sum(w); se <- 1 / sqrt(sum(w))
    Q  <- rowSums(sweep(sweep(z, 1, zp, "-")^2, 2, w, "*")); df <- ncol(r3) - 1
    data.frame(r_pooled = tanh(zp), p_pooled = 2 * pnorm(-abs(zp / se)),
               I2 = ifelse(Q > 0, pmax(0, (Q - df) / Q), 0), row.names = NULL)
  }

  rho_mean <- rowMeans(rho_mat[genes_run, , drop = FALSE])
  rho_ez   <- rho_mean[["EZH2"]]
  long <- list()
  for (mn in names(MODELS)) {
    spec <- MODELS[[mn]]; r_arr <- run_model(spec)
    k_cov <- length(spec$Z) + as.integer(isTRUE(spec$cc))
    # internal validation: vectorised EZH2 estimates must equal the published pspear() implementation
    if (!isTRUE(spec$cc)) for (k in names(co)) for (t in TARGETS) {
      ref <- pspear("EZH2", t, spec$Z, co[[k]])[["r"]]
      stopifnot(abs(ref - r_arr["EZH2", t, k]) < 1e-8)
    }
    for (t in TARGETS) {
      pl <- pool(r_arr[, t, ], n_c, k_cov)
      long[[length(long) + 1]] <- data.frame(model = mn, target = t, gene = genes_run,
        setNames(as.data.frame(r_arr[, t, ]), paste0("r_", names(co))), pl,
        rho_prolif_mean = rho_mean, matched = abs(rho_mean - rho_ez) <= MATCH_TOL, row.names = NULL)
    }
  }
  long <- do.call(rbind, long)
  write.csv(long, file.path(RES_DIR, "Table_S6a_specificity_by_gene.csv"), row.names = FALSE)

  # ---------------- summary: where does EZH2 sit in the empirical null? ----------------
  summ <- do.call(rbind, lapply(split(long, list(long$model, long$target), drop = TRUE), function(x) {
    e <- x[x$gene == "EZH2", ]; cnd <- x[x$gene != "EZH2", ]
    one <- function(s, label) {
      N <- nrow(s); if (N == 0) return(NULL)
      n_ge <- sum(s$r_pooled >= e$r_pooled)
      data.frame(model = e$model, target = e$target, comparator_set = label, N_comparators = N,
        EZH2_r_pooled = e$r_pooled, EZH2_p = e$p_pooled, EZH2_I2 = e$I2,
        null_median_r = median(s$r_pooled), null_p25 = unname(quantile(s$r_pooled, .25)), null_p75 = unname(quantile(s$r_pooled, .75)),
        frac_nominal_p05 = mean(s$p_pooled < 0.05), frac_BH_q05 = mean(p.adjust(s$p_pooled, "BH") < 0.05),
        n_comparators_ge_EZH2 = n_ge, empirical_p_onesided = (1 + n_ge) / (1 + N), row.names = NULL)
    }
    rbind(one(cnd, "all_coupled"), one(cnd[cnd$matched, ], "rho_matched"))
  }))
  write.csv(summ, file.path(RES_DIR, "Table_S6b_specificity_summary.csv"), row.names = FALSE)

  cat("\n==== EZH2 versus proliferation-coupled comparator genes (pooled partial Spearman r) ====\n")
  print(transform(summ, EZH2_r_pooled = round(EZH2_r_pooled, 3), EZH2_p = fmt_p(EZH2_p), EZH2_I2 = round(EZH2_I2, 2),
                  null_median_r = round(null_median_r, 3), null_p25 = round(null_p25, 3), null_p75 = round(null_p75, 3),
                  frac_nominal_p05 = round(frac_nominal_p05, 2), frac_BH_q05 = round(frac_BH_q05, 2),
                  empirical_p_onesided = round(empirical_p_onesided, 4)), row.names = FALSE)
  cat("\nHow to read: null_median_r far from 0 means adjustment leaves systematic residual confounding\n",
      "(the true null is NOT zero). A high frac_nominal_p05 means many proliferation-coupled genes pass the\n",
      "same test as EZH2. EZH2 is specific only if empirical_p_onesided is small in LAG3 under M3, M3h and M4.\n", sep = "")

  # ---------------- figure ----------------
  pl <- subset(long, gene != "EZH2" & model != "M0_raw")
  ez <- subset(long, gene == "EZH2" & model != "M0_raw")
  pl$model <- factor(pl$model); ez$model <- factor(ez$model, levels = levels(pl$model))
  fig <- ggplot(pl, aes(r_pooled)) +
    geom_histogram(bins = 30, fill = "grey70", colour = "white") +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    geom_vline(data = ez, aes(xintercept = r_pooled), colour = "#C0392B", linewidth = 0.9) +
    facet_grid(target ~ model, scales = "free_y") +
    labs(x = "Pooled partial Spearman r with target gene (3 cohorts)", y = "Comparator genes",
         title = "EZH2 (red) within the null distribution of proliferation-coupled genes",
         subtitle = sprintf("Comparators: rho with proliferation >= %.2f in all cohorts (n = %d genes); identical adjustment for each", RHO_MIN, length(cand))) +
    theme_pub + theme(strip.text.x = element_text(size = 7))
  save_fig(fig, "FigureS2_specificity_benchmark", 11, 7)
  cat("\n✓ Table S6a/S6b and Figure S2 saved\n")
})

# =====================================================================
# 3c Survival, robustness, pooling and cohort sensitivity analyses
# =====================================================================
step("3c Survival, robustness, pooling and cohort sensitivity analyses")
local({
  set.seed(2026)
  say  <- function(...) cat(sprintf(...), "\n", sep = "")
  hdr  <- function(x) cat("\n\n==================== ", x, " ====================\n", sep = "")
  r3   <- function(x) round(x, 3)
  num  <- function(x) suppressWarnings(as.numeric(as.character(x)))
  co   <- readRDS(file.path(RES_DIR, "cohorts.rds"))$cohorts
  EST  <- c("ESTIMATE_immune", "ESTIMATE_stromal")
  M3   <- c("prolif", EST)
  MODELS <- list(M0_raw = character(0), M1_proliferation = "prolif", M2_ESTIMATE = EST,
                 M3_prolif_ESTIMATE = M3, M4_prolif_ESTIMATE_Tcells = c(M3, "MCP_T_cell"))

  # =====================================================================
  # A. SURVIVAL
  # =====================================================================
  hdr("A. Survival sensitivity analyses")
  fit_cox <- function(d, covs, label, cohort) {
    f  <- as.formula(paste("Surv(time, status) ~", paste(c("EZH2_z", covs), collapse = " + ")))
    dd <- d[complete.cases(d[, c("time", "status", "EZH2_z", covs)]), ]
    cx <- coxph(f, data = dd); s <- summary(cx)
    zt <- tryCatch(cox.zph(cx)$table, error = function(e) NULL)
    out <- data.frame(cohort = cohort, analysis = label, term = rownames(s$conf.int), n = s$n, events = s$nevent,
                      HR = s$conf.int[, 1], lo = s$conf.int[, 3], hi = s$conf.int[, 4], p = s$coefficients[, 5],
                      row.names = NULL)
    out$PH_p <- if (is.null(zt)) NA else zt[match(sub("MGMTUnmethylated", "MGMT", out$term), rownames(zt)), "p"]
    out$PH_global_p <- if (is.null(zt)) NA else zt["GLOBAL", "p"]
    out
  }
  surv_rows <- list()
  for (k in names(co)) {
    d <- co[[k]]
    d$time   <- d$OS_months; d$status <- num(d$event)
    d <- d[!is.na(d$time) & !is.na(d$status) & d$time > 0, ]
    d$EZH2_z <- as.numeric(scale(d$log2EZH2)); d$age10 <- num(d$age) / 10
    if (k != "TCGA") { d$radio <- num(d$radio); d$chemo <- num(d$chemo) }
    harm <- c("age10", "MGMT")
    cc_h <- complete.cases(d[, c("time", "status", "EZH2_z", harm)])
    surv_rows[[length(surv_rows) + 1]] <- fit_cox(d, character(0), "Univariable, all patients", k)
    surv_rows[[length(surv_rows) + 1]] <- fit_cox(d[cc_h, ], character(0), "Univariable, harmonised complete-case sample", k)
    surv_rows[[length(surv_rows) + 1]] <- fit_cox(d, harm, "Harmonised: age + MGMT", k)
    if (k == "TCGA") {
      surv_rows[[length(surv_rows) + 1]] <- fit_cox(d, c("age10", "KPS", "MGMT"), "Original: age + KPS + MGMT", k)
    } else {
      orig <- c("age10", "MGMT", "radio", "chemo")
      surv_rows[[length(surv_rows) + 1]] <- fit_cox(d, orig, "Original: age + MGMT + radio + chemo", k)
      # 3-month landmark: only patients alive at 3 months, time counted from the landmark
      lm3 <- d[d$time > 3, ]; lm3$time <- lm3$time - 3
      surv_rows[[length(surv_rows) + 1]] <- fit_cox(lm3, orig, "Landmark 3 mo: age + MGMT + radio + chemo", k)
      surv_rows[[length(surv_rows) + 1]] <- fit_cox(lm3, harm, "Landmark 3 mo: age + MGMT", k)
      say("  %s: treatment coding radio %s | chemo %s (values: %s / %s)", k,
          paste(names(table(d$radio, useNA = "ifany")), table(d$radio, useNA = "ifany"), sep = "=", collapse = " "),
          paste(names(table(d$chemo, useNA = "ifany")), table(d$chemo, useNA = "ifany"), sep = "=", collapse = " "),
          paste(sort(unique(d$radio)), collapse = ","), paste(sort(unique(d$chemo)), collapse = ","))
      say("  %s: deaths within 3 months: %d of %d (excluded from landmark analysis)", k, sum(d$time <= 3 & d$status == 1), nrow(d))
    }
    # alternative cut points and continuous time interaction (harmonised covariates)
    dh <- d[cc_h, ]
    for (cut in c(12, 18, 24)) {
      ds <- survSplit(Surv(time, status) ~ ., data = dh, cut = cut, episode = "period")
      s3 <- summary(coxph(Surv(tstart, time, status) ~ EZH2_z:strata(period) + age10 + MGMT, data = ds))
      ii <- grep("EZH2_z", rownames(s3$coefficients))
      surv_rows[[length(surv_rows) + 1]] <- data.frame(cohort = k,
        analysis = paste0("Harmonised, EZH2 split at ", cut, " mo"),
        term = ifelse(grepl("=1", rownames(s3$coefficients)[ii]), paste0("EZH2_z (0-", cut, " mo)"), paste0("EZH2_z (>", cut, " mo)")),
        n = s3$n, events = s3$nevent, HR = s3$conf.int[ii, 1], lo = s3$conf.int[ii, 3], hi = s3$conf.int[ii, 4],
        p = s3$coefficients[ii, 5], PH_p = NA, PH_global_p = NA, row.names = NULL)
    }
    tt_fit <- coxph(Surv(time, status) ~ EZH2_z + tt(EZH2_z) + age10 + MGMT, data = dh,
                    tt = function(x, t, ...) x * log(t + 1))
    st <- summary(tt_fit)$coefficients
    surv_rows[[length(surv_rows) + 1]] <- data.frame(cohort = k, analysis = "Harmonised, EZH2 x log(time+1)",
      term = c("EZH2_z", "EZH2_z x log(time+1)"), n = tt_fit$n, events = tt_fit$nevent,
      HR = exp(st[c("EZH2_z", "tt(EZH2_z)"), "coef"]),
      lo = exp(st[c("EZH2_z", "tt(EZH2_z)"), "coef"] - 1.96 * st[c("EZH2_z", "tt(EZH2_z)"), "se(coef)"]),
      hi = exp(st[c("EZH2_z", "tt(EZH2_z)"), "coef"] + 1.96 * st[c("EZH2_z", "tt(EZH2_z)"), "se(coef)"]),
      p = st[c("EZH2_z", "tt(EZH2_z)"), "Pr(>|z|)"], PH_p = NA, PH_global_p = NA, row.names = NULL)
  }
  surv <- do.call(rbind, surv_rows)

  # pooled EZH2 HR (inverse-variance fixed effect) per analysis available in all 3 cohorts
  ez <- subset(surv, term == "EZH2_z")
  pool_hr <- do.call(rbind, lapply(split(ez, ez$analysis), function(x) {
    if (nrow(x) < 2) return(NULL)
    b <- log(x$HR); se <- (log(x$hi) - log(x$lo)) / (2 * 1.96); w <- 1 / se^2
    bp <- sum(w * b) / sum(w); sp <- 1 / sqrt(sum(w))
    data.frame(cohort = paste0("Pooled (", nrow(x), " cohorts)"), analysis = x$analysis[1], term = "EZH2_z",
               n = sum(x$n), events = sum(x$events), HR = exp(bp), lo = exp(bp - 1.96 * sp), hi = exp(bp + 1.96 * sp),
               p = 2 * pnorm(-abs(bp / sp)), PH_p = NA, PH_global_p = NA)
  }))
  surv <- rbind(surv, pool_hr)
  write.csv(surv, file.path(RES_DIR, "Table_S7a_survival_sensitivity.csv"), row.names = FALSE)
  cat("\n--- EZH2 hazard ratio per SD, all sensitivity models ---\n")
  print(transform(subset(surv, grepl("EZH2", term)), HR = r3(HR), lo = r3(lo), hi = r3(hi), p = fmt_p(p),
                  PH_p = r3(PH_p), PH_global_p = r3(PH_global_p)), row.names = FALSE)
  cat("\n--- Treatment terms (original vs 3-month landmark) ---\n")
  print(transform(subset(surv, term %in% c("radio", "chemo")), HR = r3(HR), lo = r3(lo), hi = r3(hi), p = fmt_p(p)), row.names = FALSE)
  cat("\n--- Proportional-hazards violations (p < 0.05), all models ---\n")
  print(transform(subset(surv, !is.na(PH_p) & PH_p < 0.05), PH_p = r3(PH_p)), row.names = FALSE)

  # power: minimum detectable HR per SD (two-sided alpha 0.05, 80% power; Schoenfeld)
  mdhr <- function(D) exp((qnorm(0.975) + qnorm(0.80)) / sqrt(D))
  pw <- subset(ez, analysis %in% c("Univariable, all patients", "Harmonised: age + MGMT"))
  pw <- rbind(pw[, c("cohort", "analysis", "events")],
              data.frame(cohort = "All three", analysis = unique(pw$analysis),
                         events = tapply(pw$events, pw$analysis, sum)[unique(pw$analysis)]))
  pw$min_detectable_HR_per_SD <- r3(mdhr(pw$events))
  pw$min_detectable_protective_HR <- r3(1 / mdhr(pw$events))
  cat("\n--- Power: minimum HR per SD detectable with 80% power (alpha 0.05, two-sided) ---\n")
  print(pw, row.names = FALSE)
  write.csv(pw, file.path(RES_DIR, "Table_S7b_survival_power.csv"), row.names = FALSE)

  # =====================================================================
  # B. LAG3 / PDCD1 ROBUSTNESS
  # =====================================================================
  hdr("B. LAG3 and PDCD1 robustness")
  TG <- c("LAG3", "PDCD1")
  det <- do.call(rbind, lapply(names(co), function(k) do.call(rbind, lapply(TG, function(g) {
    x <- co[[k]][[g]]
    data.frame(cohort = k, gene = g, n = length(x), median = median(x), IQR_low = unname(quantile(x, .25)),
               IQR_high = unname(quantile(x, .75)), pct_zero = 100 * mean(x == 0), pct_gt_0.1 = 100 * mean(x > 0.1),
               pct_ge_1 = 100 * mean(x >= 1), n_distinct = length(unique(x)), max = max(x))
  }))))
  cat("--- Detectability (TPM in TCGA, FPKM in CGGA) ---\n"); print(transform(det, median = r3(median), IQR_low = r3(IQR_low),
    IQR_high = r3(IQR_high), pct_zero = round(pct_zero, 1), pct_gt_0.1 = round(pct_gt_0.1, 1), pct_ge_1 = round(pct_ge_1, 1), max = r3(max)), row.names = FALSE)

  rob_rows <- list()
  for (g in TG) for (mn in c("M3_prolif_ESTIMATE", "M4_prolif_ESTIMATE_Tcells")) {
    Z <- MODELS[[mn]]
    variants <- list(
      "All samples"                 = function(d) d,
      "Expression > 0 only"         = function(d) d[d[[g]] > 0, ],
      "Top 5% of target excluded"   = function(d) d[d[[g]] <= quantile(d[[g]], 0.95), ],
      "Top 5% of EZH2 excluded"     = function(d) d[d$EZH2 <= quantile(d$EZH2, 0.95), ])
    for (vn in names(variants)) {
      per <- lapply(names(co), function(k) pspear("EZH2", g, Z, variants[[vn]](co[[k]])))
      r <- sapply(per, `[[`, "r"); n <- sapply(per, `[[`, "n"); m <- meta_r(r, n, length(Z))
      rob_rows[[length(rob_rows) + 1]] <- data.frame(gene = g, model = mn, variant = vn,
        n_TCGA = n[1], n_CGGA_325 = n[2], n_CGGA_693 = n[3], r_TCGA = r[1], r_CGGA_325 = r[2], r_CGGA_693 = r[3],
        r_pooled = m[["r"]], lo = m[["lo"]], hi = m[["hi"]], p_pooled = m[["p"]], I2 = m[["I2"]],
        LOO_min_r = NA, LOO_max_r = NA, LOO_min_pooled = NA, LOO_max_pooled = NA)
    }
    # leave-one-out: drop each sample once (per cohort), recompute cohort r and pooled r
    base <- sapply(names(co), function(k) pspear("EZH2", g, Z, co[[k]])[["r"]])
    n0   <- sapply(names(co), function(k) nrow(na.omit(co[[k]][, c("EZH2", g, Z)])))
    loo_r <- list(); loo_pool <- c()
    for (j in seq_along(co)) {
      d <- co[[j]]
      rj <- vapply(seq_len(nrow(d)), function(i) pspear("EZH2", g, Z, d[-i, ])[["r"]], numeric(1))
      loo_r[[j]] <- rj
      for (x in rj) { rr <- base; nn <- n0; rr[j] <- x; nn[j] <- nn[j] - 1; loo_pool <- c(loo_pool, meta_r(rr, nn, length(Z))[["r"]]) }
    }
    i <- length(rob_rows) - 3
    rob_rows[[i]]$LOO_min_r <- min(unlist(loo_r)); rob_rows[[i]]$LOO_max_r <- max(unlist(loo_r))
    rob_rows[[i]]$LOO_min_pooled <- min(loo_pool); rob_rows[[i]]$LOO_max_pooled <- max(loo_pool)
  }
  rob <- do.call(rbind, rob_rows)
  write.csv(det, file.path(RES_DIR, "Table_S7c_checkpoint_detectability.csv"), row.names = FALSE)
  write.csv(rob, file.path(RES_DIR, "Table_S7d_checkpoint_robustness.csv"), row.names = FALSE)
  cat("\n--- Robustness of EZH2-LAG3 / EZH2-PDCD1 (pooled partial Spearman r) ---\n")
  print(transform(rob[, c("gene", "model", "variant", "r_TCGA", "r_CGGA_325", "r_CGGA_693", "r_pooled", "lo", "hi", "p_pooled", "I2",
                          "LOO_min_pooled", "LOO_max_pooled")],
                  r_TCGA = r3(r_TCGA), r_CGGA_325 = r3(r_CGGA_325), r_CGGA_693 = r3(r_CGGA_693), r_pooled = r3(r_pooled),
                  lo = r3(lo), hi = r3(hi), p_pooled = fmt_p(p_pooled), I2 = round(I2, 2),
                  LOO_min_pooled = r3(LOO_min_pooled), LOO_max_pooled = r3(LOO_max_pooled)), row.names = FALSE)

  # =====================================================================
  # C. META-ANALYSIS SENSITIVITY
  # =====================================================================
  hdr("C. Meta-analysis sensitivity (all 8 checkpoints x 5 models)")
  meta_all <- function(r, n, k) {
    # (1) Fisher z, var 1/(n-3-k)            [main analysis]
    # (2) Fisher z, var 1.06/(n-3-k)         [Fieller et al. 1957 variance for Spearman]
    # (3) DerSimonian-Laird random effects on (1); (4) Hartung-Knapp-Sidik-Jonkman CI
    # (5) partial r scale, var (1-r^2)^2/(n-2-k), fixed effect [van Aert & Goos 2023]
    z <- atanh(r); out <- list()
    fe <- function(y, v) { w <- 1 / v; m <- sum(w * y) / sum(w); se <- 1 / sqrt(sum(w)); c(m, se, sum(w * (y - m)^2)) }
    a <- fe(z, 1 / (n - 3 - k)); df <- length(r) - 1
    out$FE      <- c(tanh(a[1]), tanh(a[1] - 1.96 * a[2]), tanh(a[1] + 1.96 * a[2]), 2 * pnorm(-abs(a[1] / a[2])))
    b <- fe(z, 1.06 / (n - 3 - k))
    out$Fieller <- c(tanh(b[1]), tanh(b[1] - 1.96 * b[2]), tanh(b[1] + 1.96 * b[2]), 2 * pnorm(-abs(b[1] / b[2])))
    v <- 1 / (n - 3 - k); w <- 1 / v
    tau2 <- max(0, (a[3] - df) / (sum(w) - sum(w^2) / sum(w)))
    wr <- 1 / (v + tau2); mr <- sum(wr * z) / sum(wr); ser <- 1 / sqrt(sum(wr))
    out$DL <- c(tanh(mr), tanh(mr - 1.96 * ser), tanh(mr + 1.96 * ser), 2 * pnorm(-abs(mr / ser)))
    q_hk <- sum(wr * (z - mr)^2) / df; se_hk <- sqrt(q_hk / sum(wr)); tq <- qt(0.975, df)
    out$HK <- c(tanh(mr), tanh(mr - tq * se_hk), tanh(mr + tq * se_hk), 2 * pt(-abs(mr / se_hk), df))
    vr <- (1 - r^2)^2 / (n - 2 - k); cc <- fe(r, vr)
    out$partial_r_scale <- c(cc[1], cc[1] - 1.96 * cc[2], cc[1] + 1.96 * cc[2], 2 * pnorm(-abs(cc[1] / cc[2])))
    c(unlist(lapply(names(out), function(nm) setNames(out[[nm]], paste0(nm, c("_r", "_lo", "_hi", "_p"))))), tau2 = tau2)
  }
  mrows <- list(); loo_rows <- list(); cohort_rows <- list()
  for (mn in names(MODELS)) for (g in CHECKPOINTS) {
    Z <- MODELS[[mn]]
    per <- lapply(names(co), function(k) pspear("EZH2", g, Z, co[[k]]))
    r <- sapply(per, `[[`, "r"); n <- sapply(per, `[[`, "n")
    mrows[[length(mrows) + 1]] <- data.frame(model = mn, gene = g, t(meta_all(r, n, length(Z))), check.names = FALSE)
    cohort_rows[[length(cohort_rows) + 1]] <- data.frame(model = mn, gene = g, cohort = names(co), n = n, r = r,
      lo = tanh(atanh(r) - 1.96 / sqrt(n - 3 - length(Z))), hi = tanh(atanh(r) + 1.96 / sqrt(n - 3 - length(Z))))
    for (j in seq_along(co)) { m <- meta_r(r[-j], n[-j], length(Z))
      loo_rows[[length(loo_rows) + 1]] <- data.frame(model = mn, gene = g, cohort_left_out = names(co)[j],
        r_pooled = m[["r"]], lo = m[["lo"]], hi = m[["hi"]], p = m[["p"]], I2 = m[["I2"]]) }
  }
  msens <- do.call(rbind, mrows); loo <- do.call(rbind, loo_rows); coh <- do.call(rbind, cohort_rows)
  write.csv(msens, file.path(RES_DIR, "Table_S7e_meta_sensitivity.csv"), row.names = FALSE)
  write.csv(loo,   file.path(RES_DIR, "Table_S7f_leave_one_cohort_out.csv"), row.names = FALSE)
  show <- subset(msens, gene %in% c("LAG3", "PDCD1", "CD276") & model %in% c("M3_prolif_ESTIMATE", "M4_prolif_ESTIMATE_Tcells"))
  cat("--- LAG3 / PDCD1 / CD276 under alternative pooling methods (r [95% CI]; P) ---\n")
  for (i in seq_len(nrow(show))) { x <- show[i, ]
    f <- function(p) sprintf("%.2f [%.2f, %.2f]; P=%s", x[[paste0(p, "_r")]], x[[paste0(p, "_lo")]], x[[paste0(p, "_hi")]], fmt_p(x[[paste0(p, "_p")]]))
    say("%-6s %-26s FE %s | Fieller %s | DL %s | HK %s | partial-r %s", x$gene, x$model, f("FE"), f("Fieller"), f("DL"), f("HK"), f("partial_r_scale")) }
  cat("\n--- Leave-one-cohort-out (LAG3, PDCD1) ---\n")
  print(transform(subset(loo, gene %in% c("LAG3", "PDCD1") & model %in% c("M3_prolif_ESTIMATE", "M4_prolif_ESTIMATE_Tcells")),
                  r_pooled = r3(r_pooled), lo = r3(lo), hi = r3(hi), p = fmt_p(p), I2 = round(I2, 2)), row.names = FALSE)

  # forest plot: LAG3, PDCD1, CD276 across the five models
  lab_m <- c(M0_raw = "Unadjusted", M1_proliferation = "+ Proliferation", M2_ESTIMATE = "+ ESTIMATE",
             M3_prolif_ESTIMATE = "+ Prolif. + ESTIMATE", M4_prolif_ESTIMATE_Tcells = "+ Prolif. + ESTIMATE + T cells")
  fp <- rbind(transform(subset(coh, gene %in% c("LAG3", "PDCD1", "CD276"))[, c("model", "gene", "cohort", "r", "lo", "hi")], type = "Cohort"),
              with(subset(msens, gene %in% c("LAG3", "PDCD1", "CD276")),
                   data.frame(model, gene, cohort = "Pooled (FE)", r = FE_r, lo = FE_lo, hi = FE_hi, type = "Pooled")),
              with(subset(msens, gene %in% c("LAG3", "PDCD1", "CD276")),
                   data.frame(model, gene, cohort = "Pooled (RE, HK)", r = HK_r, lo = HK_lo, hi = HK_hi, type = "Pooled")))
  fp$model  <- factor(lab_m[fp$model], levels = rev(lab_m))
  fp$cohort <- factor(fp$cohort, levels = c(names(co), "Pooled (FE)", "Pooled (RE, HK)"))
  fp$gene   <- factor(fp$gene, levels = c("LAG3", "PDCD1", "CD276"))
  pfp <- ggplot(fp, aes(r, model, colour = cohort, shape = type)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = pmax(lo, -1), xmax = pmin(hi, 1)), position = position_dodge(width = 0.75), size = 0.3) +
    facet_wrap(~gene) +
    scale_colour_manual(values = c(TCGA = "#2C3E50", CGGA_325 = "#E67E22", CGGA_693 = "#16A085",
                                   `Pooled (FE)` = "#C0392B", `Pooled (RE, HK)` = "#8E44AD")) +
    scale_shape_manual(values = c(Cohort = 16, Pooled = 18), guide = "none") +
    labs(x = "Partial Spearman r with EZH2 (95% CI)", y = NULL, colour = NULL,
         title = "EZH2-checkpoint associations by cohort and adjustment model",
         subtitle = "FE: fixed effect (Fisher z); RE, HK: DerSimonian-Laird with Hartung-Knapp CI") +
    theme_pub + theme(legend.position = "bottom")
  save_fig(pfp, "FigureS4_forest_checkpoints", 11, 6)

  # =====================================================================
  # D. COHORT DEFINITION
  # =====================================================================
  hdr("D. Cohort definition and overlap")
  ## D1. CGGA_325 vs CGGA_693 duplicates
  c325 <- readRDS(file.path(DATA_DIR, "cgga325_clinical.rds")); c693 <- readRDS(file.path(DATA_DIR, "cgga693_clinical.rds"))
  for (cc in c("c325", "c693")) { x <- get(cc); for (v in names(x)) if (is.character(x[[v]])) x[[v]] <- trimws(x[[v]]); assign(cc, x) }
  a <- c325[c325$CGGA_ID %in% co$CGGA_325$sample, ]; b <- c693[c693$CGGA_ID %in% co$CGGA_693$sample, ]
  key <- function(x) paste(num(x$Age), tolower(x$Gender), num(x$OS), num(x$Censor))
  clin_match <- intersect(key(a), key(b))
  say("Shared CGGA_IDs between the two analysed cohorts: %d", length(intersect(a$CGGA_ID, b$CGGA_ID)))
  say("Patients with identical age + sex + OS + censoring status in both analysed cohorts: %d", length(clin_match))
  e325 <- readRDS(file.path(DATA_DIR, "cgga325_expr_gbm.rds")); e693 <- readRDS(file.path(DATA_DIR, "cgga693_expr_gbm.rds"))
  to_mat <- function(ex, ids) { m <- as.matrix(ex[, ids, drop = FALSE]); storage.mode(m) <- "numeric"; rownames(m) <- as.character(ex[[1]])
    m[!duplicated(rownames(m)) & !is.na(rownames(m)), , drop = FALSE] }
  m1 <- to_mat(e325, co$CGGA_325$sample); m2 <- to_mat(e693, co$CGGA_693$sample)
  gs <- intersect(rownames(m1), rownames(m2))
  l1 <- log2(m1[gs, ] + 1); l2 <- log2(m2[gs, ] + 1)
  v  <- apply(cbind(l1, l2), 1, var); top <- names(sort(v, decreasing = TRUE))[seq_len(min(2000, length(v)))]
  xc <- cor(l1[top, ], l2[top, ], method = "spearman")
  w1 <- cor(l1[top, ], method = "spearman"); diag(w1) <- NA
  best <- apply(xc, 1, max)
  say("Expression check (%d most variable shared genes, Spearman): max between-cohort sample correlation = %.3f;",
      length(top), max(xc))
  say("  median of each CGGA_325 sample's best match in CGGA_693 = %.3f; within CGGA_325 best-match median = %.3f",
      median(best), median(apply(w1, 1, max, na.rm = TRUE)))
  say("  pairs with between-cohort r > 0.98 (possible duplicates): %d", sum(xc > 0.98))
  ov <- data.frame(check = c("shared CGGA_ID", "identical age/sex/OS/censor", "expression pairs r>0.98", "max between-cohort r",
                             "median best-match r (325 vs 693)", "median best-match r (within 325)"),
                   value = c(length(intersect(a$CGGA_ID, b$CGGA_ID)), length(clin_match), sum(xc > 0.98), r3(max(xc)),
                             r3(median(best)), r3(median(apply(w1, 1, max, na.rm = TRUE)))))
  write.csv(ov, file.path(RES_DIR, "Table_S7g_CGGA_overlap_check.csv"), row.names = FALSE)

  ## D2. TCGA molecular annotation relevant to WHO CNS5
  cl <- readRDS(file.path(DATA_DIR, "tcga_cohort_clinical.rds"))
  cat("\nTCGA curated annotation columns that may carry CNS5 features:\n")
  pat <- "TERT|EGFR|Chr.7|Chr.10|7.gain|10.loss|IDH|Grade|Histology|codel|ATRX|Mutation.Count|TMB|Supervised|Original.Subtype"
  mol_cols <- grep(pat, colnames(cl), value = TRUE, ignore.case = TRUE)
  print(mol_cols)
  for (v in mol_cols) { x <- cl[[v]]
    if (is.numeric(x) && length(unique(x)) > 10) say("  %-45s numeric: median %.2f, missing %d", v, median(x, na.rm = TRUE), sum(is.na(x)))
    else { t <- table(x, useNA = "ifany"); say("  %-45s %s", v, paste(names(t), t, sep = "=", collapse = " | ")) } }

  ## D3. Subtype (and TMB) adjustment in TCGA
  pspear_f <- function(x, y, Z, F = character(0), data) {      # ranks for x, y and numeric Z; factors in F entered as dummies
    dd <- na.omit(data[, c(x, y, Z, F), drop = FALSE]); n <- nrow(dd)
    R  <- as.data.frame(lapply(dd[, c(x, y, Z), drop = FALSE], rank))
    X  <- model.matrix(~ ., data = cbind(R[, Z, drop = FALSE], dd[, F, drop = FALSE]))
    k  <- ncol(X) - 1
    rr <- cor(lm.fit(X, R[[x]])$residuals, lm.fit(X, R[[y]])$residuals)
    tt <- rr * sqrt((n - 2 - k) / (1 - rr^2)); c(n = n, k = k, r = rr, p = 2 * pt(-abs(tt), n - 2 - k))
  }
  tc <- co$TCGA
  tmb_col <- grep("Mutation.Count|TMB", colnames(cl), value = TRUE, ignore.case = TRUE)[1]
  if (!is.na(tmb_col)) tc$TMB <- num(cl[match(tc$sample, rownames(cl)), tmb_col])
  sub_rows <- list()
  for (g in c("LAG3", "PDCD1", "CD276")) {
    specs <- list("M3 (same patients with subtype)" = list(Z = M3, F = character(0), d = tc[!is.na(tc$subtype), ]),
                  "M3 + subtype"                    = list(Z = M3, F = "subtype", d = tc),
                  "M4 + subtype"                    = list(Z = c(M3, "MCP_T_cell"), F = "subtype", d = tc))
    if (!is.na(tmb_col)) specs[["M3 + subtype + mutation count"]] <- list(Z = c(M3, "TMB"), F = "subtype", d = tc)
    for (sn in names(specs)) { s <- specs[[sn]]; v <- pspear_f("EZH2", g, s$Z, s$F, s$d)
      sub_rows[[length(sub_rows) + 1]] <- data.frame(gene = g, model = sn, n = v[["n"]], k = v[["k"]], r = v[["r"]], p = v[["p"]]) }
  }
  subadj <- do.call(rbind, sub_rows)
  write.csv(subadj, file.path(RES_DIR, "Table_S7h_TCGA_subtype_adjustment.csv"), row.names = FALSE)
  cat(sprintf("\n--- TCGA: adjustment for transcriptional subtype%s ---\n", if (!is.na(tmb_col)) paste0(" and ", tmb_col) else " (no mutation-count column found)"))
  print(transform(subadj, r = r3(r), p = fmt_p(p)), row.names = FALSE)

  # =====================================================================
  # E. GENE OVERLAP WITH ADJUSTMENT SIGNATURES
  # =====================================================================
  hdr("E. Overlap of tested genes with ESTIMATE / MCP-counter signatures")
  tested <- unique(c("EZH2", CHECKPOINTS))
  est_genes <- tryCatch({ e <- new.env(); utils::data("SI_geneset", package = "estimate", envir = e)
    list(ESTIMATE_stromal = as.character(unlist(e$SI_geneset[1, -1])), ESTIMATE_immune = as.character(unlist(e$SI_geneset[2, -1]))) },
    error = function(e) { message("estimate gene sets not available: ", conditionMessage(e)); NULL })
  mcp_genes <- tryCatch({
    f <- system.file("extdata", "mcp_counter_genes.txt", package = "immunedeconv")
    if (!nzchar(f)) f <- list.files(system.file(package = "immunedeconv"), pattern = "mcp.*gene", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)[1]
    if (is.na(f) || !nzchar(f)) stop("MCP-counter marker file not found in immunedeconv")
    g <- read.delim(f, check.names = FALSE); split(as.character(g[[grep("HUGO|symbol", names(g), ignore.case = TRUE)[1]]]),
                                                  g[[grep("Cell population|population", names(g), ignore.case = TRUE)[1]]]) },
    error = function(e) { message("MCP-counter markers not read: ", conditionMessage(e)); NULL })
  sigs <- c(est_genes, mcp_genes, list(Proliferation_10gene = PROLIF_GENES))
  ovl <- do.call(rbind, lapply(names(sigs), function(s) data.frame(signature = s, size = length(unique(sigs[[s]])),
    tested_genes_in_signature = paste(intersect(tested, sigs[[s]]), collapse = ", "))))
  ovl$tested_genes_in_signature[ovl$tested_genes_in_signature == ""] <- "none"
  print(ovl, row.names = FALSE)
  write.csv(ovl, file.path(RES_DIR, "Table_S7i_signature_gene_overlap.csv"), row.names = FALSE)
})

# =====================================================================
# 3d Figures for the benchmark, simulation and forest plot
# =====================================================================
step("3d Figures for the benchmark, simulation and forest plot")
local({
  set.seed(2026)
  suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(ggrepel) })
  RES <- "results"; FIG <- "figures"; dir.create(FIG, showWarnings = FALSE)
  theme_pub <- theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11), strip.background = element_blank(),
          strip.text = element_text(face = "bold"), plot.title.position = "plot")
  save_fig <- function(p, name, w, h) {
    ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h)
    ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
    ggsave(file.path(FIG, paste0(name, ".tiff")), p, width = w, height = h, dpi = 300, compression = "lzw")
  }
  MODEL_LAB <- c(M3_prolif_ESTIMATE = "Proliferation + ESTIMATE",
                 M3h_prolif_Hallmark_ESTIMATE = "+ Hallmark cell-cycle score",
                 M4_prolif_ESTIMATE_Tcells = "+ MCP-counter T cells")

  # ---------------- Figure 2: negative-control benchmark ----------------
  a <- read.csv(file.path(RES, "Table_S6a_specificity_by_gene.csv"))
  s <- read.csv(file.path(RES, "Table_S6b_specificity_summary.csv"))
  a <- subset(a, model %in% names(MODEL_LAB)); sm <- subset(s, model %in% names(MODEL_LAB) & comparator_set == "rho_matched")
  s <- subset(s, model %in% names(MODEL_LAB) & comparator_set == "all_coupled")
  sm$model <- factor(MODEL_LAB[sm$model], levels = MODEL_LAB); sm$target <- factor(sm$target, levels = c("LAG3", "PDCD1", "CD276"))
  s$model <- MODEL_LAB[s$model]
  a$model  <- factor(MODEL_LAB[a$model], levels = MODEL_LAB)
  a$target <- factor(a$target, levels = c("LAG3", "PDCD1", "CD276"))
  s$model <- factor(s$model, levels = MODEL_LAB); s$target <- factor(s$target, levels = c("LAG3", "PDCD1", "CD276"))
  s <- merge(s, sm[, c("model", "target", "n_comparators_ge_EZH2", "N_comparators", "empirical_p_onesided", "null_median_r")], by = c("model", "target"), suffixes = c("", "_m"))
  s$label <- sprintf("Matched: rank %d/%d, P = %.2f\nAll: rank %d/%d, P = %.2f", s$n_comparators_ge_EZH2_m + 1, s$N_comparators_m + 1, s$empirical_p_onesided_m,
                     s$n_comparators_ge_EZH2 + 1, s$N_comparators + 1, s$empirical_p_onesided)
  cmp <- subset(a, gene != "EZH2"); ez <- subset(a, gene == "EZH2")
  cmp$set <- factor(ifelse(cmp$matched, "Matched to EZH2 (n = 80)", "Other coupled genes (n = 58)"), levels = c("Matched to EZH2 (n = 80)", "Other coupled genes (n = 58)"))
  pA <- ggplot(cmp, aes(r_pooled, fill = set)) +
    geom_histogram(binwidth = 0.025, colour = "white", boundary = 0) +
    scale_fill_manual(values = c("grey45", "grey80"), name = NULL) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
    geom_vline(data = s, aes(xintercept = null_median_r_m), linetype = 3, colour = "#2C3E50", linewidth = 0.6) +
    geom_vline(data = ez, aes(xintercept = r_pooled), colour = "#C0392B", linewidth = 0.9) +
    geom_label(data = s, aes(x = Inf, y = Inf, label = label), hjust = 1, vjust = 1, size = 2.4, lineheight = 0.9, fill = "white", alpha = 0.9, label.size = 0, inherit.aes = FALSE) +
    facet_grid(target ~ model) +
    labs(x = "Pooled partial Spearman r with the target gene (3 cohorts)", y = "Comparator genes (n = 138)",
         title = "a  EZH2 (red line) among 138 proliferation-coupled comparator genes (dotted line: median of matched comparators)") +
    theme_pub + theme(strip.text.y = element_text(face = "bold.italic"), legend.position = "top")

  w <- reshape(subset(a, model == MODEL_LAB[["M3_prolif_ESTIMATE"]] & target %in% c("LAG3", "PDCD1"))[, c("gene", "target", "r_pooled")],
               idvar = "gene", timevar = "target", direction = "wide")
  names(w) <- c("gene", "LAG3", "PDCD1"); w$is_ez <- w$gene == "EZH2"
  e <- w[w$is_ez, ]
  w$lab <- ifelse(w$is_ez | w$LAG3 >= sort(w$LAG3, TRUE)[6] | w$PDCD1 >= sort(w$PDCD1, TRUE)[5], w$gene, "")
  pB <- ggplot(w, aes(LAG3, PDCD1)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey60") + geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    geom_hline(yintercept = e$PDCD1, colour = "#C0392B", linewidth = 0.3, linetype = 3) +
    geom_vline(xintercept = e$LAG3, colour = "#C0392B", linewidth = 0.3, linetype = 3) +
    geom_point(data = subset(w, !is_ez), colour = "grey45", size = 1.3, alpha = 0.8) +
    geom_point(data = e, colour = "#C0392B", size = 2.8) +
    geom_text_repel(aes(label = lab, colour = is_ez), size = 2.6, fontface = "italic", max.overlaps = 30,
                    min.segment.length = 0, segment.colour = "grey60", seed = 1) +
    scale_colour_manual(values = c(`FALSE` = "grey20", `TRUE` = "#C0392B"), guide = "none") +
    labs(x = "Partial r with LAG3", y = "Partial r with PDCD1",
         title = "b  Residual associations (proliferation + ESTIMATE)") + theme_pub

  rho <- subset(a, model == MODEL_LAB[["M3_prolif_ESTIMATE"]] & target == "LAG3")
  pC <- ggplot(rho, aes(rho_prolif_mean, r_pooled)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    geom_smooth(data = subset(rho, gene != "EZH2"), method = "lm", formula = y ~ x, se = TRUE, colour = "#2C3E50", linewidth = 0.6) +
    geom_point(data = subset(rho, gene != "EZH2"), colour = "grey45", size = 1.3, alpha = 0.8) +
    geom_point(data = subset(rho, gene == "EZH2"), colour = "#C0392B", size = 2.8) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.4, size = 2.8,
             label = sprintf("Spearman rho = %.2f", cor(rho$rho_prolif_mean[rho$gene != "EZH2"], rho$r_pooled[rho$gene != "EZH2"], method = "spearman"))) +
    labs(x = "Mean correlation with the 10-gene proliferation score", y = "Partial r with LAG3",
         title = "c  Residual LAG3 association vs coupling") + theme_pub
  fig4 <- pA / (pB | pC) + plot_layout(heights = c(1.55, 1))
  save_fig(fig4, "Figure2_specificity_benchmark", 10, 10.5)

  # ---------------- Figure 3: measurement-error simulation ----------------
  sim <- read.csv(file.path(RES, "Table_S6c_simulation_measurement_error.csv"))
  A <- subset(sim, world == "A_confounding_only"); A$gamma <- factor(A$gamma)
  null_med <- s$null_median_r_m[s$target == "LAG3" & s$model == MODEL_LAB[["M3_prolif_ESTIMATE"]]]
  lam_alpha <- A$lamP[which.min(abs(A$lamP - 0.989))]
  p1 <- ggplot(A, aes(lamP, mean_pooled_r, colour = gamma)) +
    geom_hline(yintercept = 0.292, linetype = 2, colour = "#C0392B") +
    geom_hline(yintercept = null_med, linetype = 3, colour = "#2C3E50") +
    geom_vline(xintercept = lam_alpha, linetype = 2, colour = "grey55") +
    annotate("text", x = 0.9, y = 0.292, label = "observed EZH2-LAG3, r = 0.29", hjust = 0, vjust = -0.5, size = 2.7, colour = "#C0392B") +
    annotate("text", x = 0.9, y = null_med, label = sprintf("benchmark median, matched comparators, r = %.2f", null_med), hjust = 0, vjust = -0.5, size = 2.7, colour = "#2C3E50") +
    annotate("text", x = lam_alpha, y = 0.005, label = "sqrt(Cronbach alpha)", hjust = 1.05, size = 2.6, colour = "grey40") +
    geom_line() + geom_point(size = 1.6) +
    scale_colour_manual(values = c(`0.3` = "#2C3E50", `0.5` = "#E67E22")) + scale_y_continuous(limits = c(0, 0.33)) +
    labs(x = "Correlation of proxy with latent proliferation", y = "Mean pooled partial r\n(EZH2 has no direct effect)",
         colour = "LAG3 loading on\nproliferation", title = "a  Spurious adjusted association") + theme_pub
  p2 <- ggplot(A, aes(lamP, reject_rate, colour = gamma)) +
    geom_hline(yintercept = 0.05, linetype = 2, colour = "grey45") + geom_line() + geom_point(size = 1.6) +
    scale_colour_manual(values = c(`0.3` = "#2C3E50", `0.5` = "#E67E22")) + scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Correlation of proxy with latent proliferation", y = "False-positive rate (nominal 0.05)",
         colour = "LAG3 loading on\nproliferation", title = "b  False-positive rate of the pooled test") + theme_pub
  save_fig((p1 | p2) + plot_layout(guides = "collect"), "Figure3_simulation", 10, 4)

  # ---------------- Fig. S4: forest plot ----------------
  t3 <- read.csv(file.path(RES, "Table3a_checkpoint_partial_correlations_by_cohort.csv"))
  ms <- read.csv(file.path(RES, "Table_S7e_meta_sensitivity.csv"))
  K  <- c(M0_raw = 0, M1_proliferation = 1, M2_ESTIMATE = 2, M3_prolif_ESTIMATE = 3, M4_prolif_ESTIMATE_Tcells = 4)
  lab_m <- c(M0_raw = "Unadjusted", M1_proliferation = "+ Proliferation", M2_ESTIMATE = "+ ESTIMATE",
             M3_prolif_ESTIMATE = "+ Prolif. + ESTIMATE", M4_prolif_ESTIMATE_Tcells = "+ Prolif. + ESTIMATE + T cells")
  t3 <- subset(t3, model %in% names(K) & gene %in% c("LAG3", "PDCD1", "CD276"))
  se <- 1 / sqrt(t3$n - 3 - K[t3$model])
  fp <- rbind(data.frame(model = t3$model, gene = t3$gene, est = t3$cohort, r = t3$r,
                         lo = tanh(atanh(t3$r) - 1.96 * se), hi = tanh(atanh(t3$r) + 1.96 * se)),
              with(subset(ms, gene %in% c("LAG3", "PDCD1", "CD276")),
                   rbind(data.frame(model, gene, est = "Pooled, fixed effect", r = FE_r, lo = FE_lo, hi = FE_hi),
                         data.frame(model, gene, est = "Pooled, random effects (HK)", r = HK_r, lo = HK_lo, hi = HK_hi))))
  fp$model <- factor(lab_m[fp$model], levels = rev(lab_m))
  fp$est   <- factor(fp$est, levels = c("TCGA", "CGGA_325", "CGGA_693", "Pooled, fixed effect", "Pooled, random effects (HK)"))
  fp$gene  <- factor(fp$gene, levels = c("LAG3", "PDCD1", "CD276"))
  pf <- ggplot(fp, aes(r, model, colour = est)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = pmax(lo, -1), xmax = pmin(hi, 1), shape = grepl("Pooled", est)),
                    position = position_dodge(width = 0.75), size = 0.25) +
    facet_wrap(~gene) + scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
    scale_colour_manual(values = c(TCGA = "#2C3E50", CGGA_325 = "#E67E22", CGGA_693 = "#16A085",
                                   `Pooled, fixed effect` = "#C0392B", `Pooled, random effects (HK)` = "#8E44AD")) +
    labs(x = "Partial Spearman r with EZH2 (95% CI)", y = NULL, colour = NULL) +
    theme_pub + theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic"))
  save_fig(pf, "FigureS4_forest_checkpoints", 10, 5.5)
  cat("✓ Figure 2, Figure 3 and Fig. S4 saved in figures/\n")
})

cat("\n\u2713 03_benchmark_and_robustness.R complete.\n")
}, finally = end_log())
