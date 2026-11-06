# =====================================================================
# 02_main_analyses.R
# Step 2: standard analyses (expression, survival, checkpoints, microenvironment, single cell, Hallmark sensitivity)
# Part of: Covariate-adjusted co-expression does not establish gene-immune
# checkpoint specificity in glioblastoma (EZH2 case study).
#
# HOW TO RUN: set the working directory to the PROJECT folder (the folder
# that contains share/), then source this script from wherever it is, e.g.
#   setwd("D:/Bioinformatics/EZH2_GBM_Bioinfo_2026")
#   source("scripts/02_main_analyses.R")
# Outputs go to results/ and figures/ in the project folder; a log is written
# to results/log_02_main_analyses.txt.
# =====================================================================
.this <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/"), error = function(e) NA_character_)
source(file.path(if (is.na(.this)) getwd() else dirname(.this), "00_setup.R"))
check_project(c("results/cohorts.rds"))
start_log("log_02_main_analyses.txt")
tryCatch({

# =====================================================================
# 2a EZH2 expression, PRC2 and proliferation
# =====================================================================
step("2a EZH2 expression, PRC2 and proliferation")
local({
  set.seed(2026)
  obj <- readRDS(file.path(RES_DIR, "cohorts.rds")); L <- obj$cohorts; N <- obj$tcga_normals
  tcga <- L$TCGA

  # ---- A. Tumor vs normal (TCGA) ----
  tn <- rbind(data.frame(group = "GBM (IDH-wt)", EZH2 = tcga$EZH2),
              data.frame(group = "Normal brain", EZH2 = N$EZH2))
  tn$group <- factor(tn$group, levels = c("Normal brain", "GBM (IDH-wt)"))
  w_tn <- wilcox.test(EZH2 ~ group, data = tn)
  fc   <- median(tcga$EZH2) / median(N$EZH2)
  cat(sprintf("Tumor vs normal: median %.2f vs %.2f TPM, fold = %.1f, Wilcoxon p = %s\n",
              median(tcga$EZH2), median(N$EZH2), fc, fmt_p(w_tn$p.value)))
  pA <- ggplot(tn, aes(group, log2(EZH2 + 1), fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) +
    scale_fill_manual(values = c("grey70", "#C0392B"), guide = "none") +
    annotate("text", x = 1.5, y = max(log2(tn$EZH2 + 1)) * 1.08,
             label = paste0("Wilcoxon p = ", fmt_p(w_tn$p.value)), size = 3.2) +
    labs(x = NULL, y = "EZH2 log2(TPM+1)", title = "A  TCGA: tumor vs normal") + theme_pub

  # ---- B. PRC2 + proliferation correlations, all cohorts ----
  rows <- list()
  for (k in names(L)) for (g in c(PRC2_GENES, "prolif")) {
    s <- pspear("EZH2", g, character(0), L[[k]])
    rows[[length(rows) + 1]] <- data.frame(cohort = k, variable = g, n = s[["n"]], rho = s[["r"]], p = s[["p"]])
  }
  cor_tab <- do.call(rbind, rows); cor_tab$q <- p.adjust(cor_tab$p, "BH")
  print(transform(cor_tab, rho = round(rho, 3), p = fmt_p(p), q = fmt_p(q)))

  bdat <- do.call(rbind, lapply(L, function(d) d[, c("cohort", "log2EZH2", "prolif")])); bdat$cohort <- factor(bdat$cohort, levels = names(L))
  cor_tab$cohort <- factor(cor_tab$cohort, levels = names(L))
  pB <- ggplot(bdat,
               aes(prolif, log2EZH2)) +
    geom_point(size = 0.8, alpha = 0.6) + geom_smooth(method = "lm", se = TRUE, colour = "#C0392B", linewidth = 0.7) +
    facet_wrap(~cohort, scales = "free") +
    geom_text(data = subset(cor_tab, variable == "prolif"),
              aes(x = -Inf, y = Inf, label = sprintf("rho = %.2f\np %s", rho, ifelse(p < 1e-4, "< 0.0001", paste("=", fmt_p(p))))),
              hjust = -0.1, vjust = 1.2, size = 3, inherit.aes = FALSE) +
    labs(x = "Proliferation score (10 genes)", y = "EZH2 log2(expression+1)",
         title = "C  EZH2 tracks tumor proliferation") + theme_pub

  pC <- ggplot(subset(cor_tab, variable != "prolif"), aes(variable, cohort, fill = rho)) +
    geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.2f", rho)), size = 3.2) +
    scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C0392B", limits = c(-1, 1)) +
    labs(x = NULL, y = NULL, fill = "Spearman\nrho", title = "B  EZH2 vs PRC2 genes") + theme_pub

  # ---- D. Clinical correlates (TCGA) ----
  kw  <- kruskal.test(EZH2 ~ subtype, data = tcga)
  pw  <- pairwise.wilcox.test(tcga$EZH2, tcga$subtype, p.adjust.method = "BH")
  w_m <- wilcox.test(EZH2 ~ MGMT, data = tcga); w_s <- wilcox.test(EZH2 ~ sex, data = tcga)
  cat(sprintf("Subtype KW p = %s | MGMT p = %s | sex p = %s\n", fmt_p(kw$p.value), fmt_p(w_m$p.value), fmt_p(w_s$p.value)))
  print(pw$p.value)
  pD <- ggplot(subset(tcga, !is.na(subtype)), aes(subtype, log2EZH2, fill = subtype)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) + geom_jitter(width = 0.15, size = 0.7, alpha = 0.6) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(x = NULL, y = "EZH2 log2(TPM+1)", title = paste0("D  TCGA subtype (Kruskal-Wallis p = ", fmt_p(kw$p.value), ")")) + theme_pub

  fig1 <- (pA | pC) / pB / pD + plot_layout(heights = c(1, 1, 1))
  save_fig(fig1, "Figure1_expression_PRC2_proliferation", 10, 11)

  clin_tab <- data.frame(test = c("Tumor vs normal (TCGA)", "MGMT unmeth vs meth (TCGA)", "Sex (TCGA)", "Subtype KW (TCGA)"),
                         p = c(w_tn$p.value, w_m$p.value, w_s$p.value, kw$p.value))
  write.csv(cor_tab, file.path(RES_DIR, "Table_S1a_EZH2_PRC2_proliferation_correlations.csv"), row.names = FALSE)
  write.csv(clin_tab, file.path(RES_DIR, "Table_S1b_EZH2_clinical_tests.csv"), row.names = FALSE)
  cat("✓ Figure 1 and Table S1 saved\n")
})

# =====================================================================
# 2b Survival
# =====================================================================
step("2b Survival")
local({
  set.seed(2026)
  L <- readRDS(file.path(RES_DIR, "cohorts.rds"))$cohorts

  km_plots <- list(); cox_rows <- list(); ph_rows <- list(); sens_rows <- list()
  for (k in names(L)) {
    d <- subset(L[[k]], !is.na(OS_months) & !is.na(event))
    d$EZH2_z <- as.numeric(scale(d$log2EZH2))
    d$age10  <- d$age / 10

    fit <- survfit(Surv(OS_months, event) ~ EZH2_group, data = d)
    lr  <- survdiff(Surv(OS_months, event) ~ EZH2_group, data = d)
    p_lr <- pchisq(lr$chisq, 1, lower.tail = FALSE)
    fu  <- summary(survfit(Surv(OS_months, 1 - event) ~ 1, data = d))$table["median"]
    med <- summary(fit)$table[, "median"]
    cat(sprintf("\n%s: n=%d, events=%d, median follow-up (reverse KM) = %.1f mo; median OS Low/High = %.1f/%.1f; log-rank p = %s\n",
                k, nrow(d), sum(d$event), fu, med[1], med[2], fmt_p(p_lr)))

    g <- ggsurvplot(fit, data = d, palette = unname(COL_LOWHIGH), legend.labs = c("EZH2-low", "EZH2-high"),
                    legend.title = "", xlab = "Months", ylab = "Overall survival", risk.table = FALSE,
                    censor.size = 2, ggtheme = theme_pub)
    km_plots[[k]] <- g$plot + labs(title = sprintf("%s  %s (n = %d)", LETTERS[match(k, names(L))], k, nrow(d))) +
      annotate("text", x = max(d$OS_months) * 0.6, y = 0.9, size = 3.2, label = paste0("log-rank p = ", fmt_p(p_lr)))

    covs <- if (k == "TCGA") c("age10", "KPS", "MGMT") else c("age10", "MGMT", "radio", "chemo")
    for (model in c("Univariable", "Adjusted")) {
      f  <- as.formula(paste("Surv(OS_months, event) ~ EZH2_z", if (model == "Adjusted") paste("+", paste(covs, collapse = " + ")) else ""))
      cx <- coxph(f, data = d); s <- summary(cx)
      zt <- cox.zph(cx)$table; ph <- zt["GLOBAL", "p"]
      ph_rows[[length(ph_rows) + 1]] <- data.frame(cohort = k, model = model, term = rownames(zt), chisq = zt[, 1], p = zt[, 3])
      for (v in rownames(s$coefficients)) cox_rows[[length(cox_rows) + 1]] <- data.frame(
        cohort = k, model = model, n = s$n, events = s$nevent, term = v,
        HR = s$conf.int[v, 1], lo = s$conf.int[v, 3], hi = s$conf.int[v, 4], p = s$coefficients[v, 5],
        C_index = s$concordance[1], PH_global_p = ph)
    }
    # ---- sensitivity analyses for proportional-hazards violations ----
    if (k == "TCGA") {   # age violates PH -> stratify baseline hazard by age tertile
      d$age_strat <- cut(d$age, quantile(d$age, c(0, 1/3, 2/3, 1), na.rm = TRUE), include.lowest = TRUE)
      s2 <- summary(coxph(Surv(OS_months, event) ~ EZH2_z + KPS + MGMT + strata(age_strat), data = d))
      sens_rows[[length(sens_rows) + 1]] <- data.frame(cohort = k, analysis = "Adjusted, stratified by age tertile",
        HR = s2$conf.int["EZH2_z", 1], lo = s2$conf.int["EZH2_z", 3], hi = s2$conf.int["EZH2_z", 4], p = s2$coefficients["EZH2_z", 5])
    }
    # time-split EZH2 effect (0-24 months vs >24 months), all cohorts
    ds <- survSplit(Surv(OS_months, event) ~ ., data = d, cut = 24, episode = "period")
    s3 <- summary(coxph(as.formula(paste("Surv(tstart, OS_months, event) ~ EZH2_z:strata(period) +", paste(covs, collapse = " + "))), data = ds))
    for (i in grep("EZH2_z", rownames(s3$coefficients))) sens_rows[[length(sens_rows) + 1]] <- data.frame(
      cohort = k, analysis = paste0("Adjusted, EZH2 effect ", ifelse(grepl("=1", rownames(s3$coefficients)[i]), "0-24 mo", ">24 mo")),
      HR = s3$conf.int[i, 1], lo = s3$conf.int[i, 3], hi = s3$conf.int[i, 4], p = s3$coefficients[i, 5])
  }
  cox_tab <- do.call(rbind, cox_rows)
  print(transform(cox_tab, HR = round(HR, 2), lo = round(lo, 2), hi = round(hi, 2), p = fmt_p(p), C_index = round(C_index, 3), PH_global_p = round(PH_global_p, 3)))
  write.csv(cox_tab, file.path(RES_DIR, "Table2_Cox_models.csv"), row.names = FALSE)
  ph_tab <- do.call(rbind, ph_rows); write.csv(ph_tab, file.path(RES_DIR, "Table_S2a_PH_assumption.csv"), row.names = FALSE)
  sens_tab <- do.call(rbind, sens_rows); write.csv(sens_tab, file.path(RES_DIR, "Table_S2b_Cox_sensitivity.csv"), row.names = FALSE)
  cat("\nPH violations (p<0.05):\n"); print(subset(ph_tab, p < 0.05 & term != "GLOBAL"))
  cat("\nSensitivity analyses:\n"); print(transform(sens_tab, HR = round(HR, 2), lo = round(lo, 2), hi = round(hi, 2), p = fmt_p(p)))

  # ---- pooled HR for EZH2 (inverse-variance, log scale) ----
  ez <- subset(cox_tab, term == "EZH2_z")
  pool <- do.call(rbind, lapply(split(ez, ez$model), function(x) {
    b <- log(x$HR); se <- (log(x$hi) - log(x$lo)) / (2 * 1.96); w <- 1 / se^2
    bp <- sum(w * b) / sum(w); sp <- 1 / sqrt(sum(w))
    data.frame(cohort = "Pooled", model = x$model[1], HR = exp(bp), lo = exp(bp - 1.96 * sp), hi = exp(bp + 1.96 * sp),
               p = 2 * pnorm(-abs(bp / sp)))
  }))
  print(pool)
  forest <- rbind(ez[, c("cohort", "model", "HR", "lo", "hi", "p")], pool)
  forest$cohort <- factor(forest$cohort, levels = rev(c(names(L), "Pooled")))
  pF <- ggplot(forest, aes(HR, cohort, colour = model)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5), size = 0.35) +
    scale_colour_manual(values = c(Adjusted = "#C0392B", Univariable = "grey40")) +
    scale_x_log10() + labs(x = "HR per 1 SD of log2 EZH2 (95% CI)", y = NULL, colour = NULL,
                           title = "D  EZH2 and overall survival") + theme_pub + theme(legend.position = "bottom")

  fig2 <- (km_plots[[1]] | km_plots[[2]]) / (km_plots[[3]] | pF)
  save_fig(fig2, "Figure2_survival", 10, 8)
  cat("✓ Figure 2 and Table 2 saved\n")
})

# =====================================================================
# 2c Checkpoint partial correlations
# =====================================================================
step("2c Checkpoint partial correlations")
local({
  set.seed(2026)
  L <- readRDS(file.path(RES_DIR, "cohorts.rds"))$cohorts

  MODELS <- list(
    M0_raw                    = character(0),
    M1_proliferation          = "prolif",
    M2_ESTIMATE               = c("ESTIMATE_immune", "ESTIMATE_stromal"),
    M3_prolif_ESTIMATE        = c("prolif", "ESTIMATE_immune", "ESTIMATE_stromal"),
    M4_prolif_ESTIMATE_Tcells = c("prolif", "ESTIMATE_immune", "ESTIMATE_stromal", "MCP_T_cell"))

  rows <- list()
  for (m in names(MODELS)) for (k in names(L)) for (g in CHECKPOINTS) {
    s <- pspear("EZH2", g, MODELS[[m]], L[[k]])
    rows[[length(rows) + 1]] <- data.frame(model = m, cohort = k, gene = g, n = s[["n"]], r = s[["r"]], p = s[["p"]])
  }
  tab <- do.call(rbind, rows)
  tab$q <- ave(tab$p, tab$model, tab$cohort, FUN = function(p) p.adjust(p, "BH"))

  # TCGA-only DNA purity model
  tp <- do.call(rbind, lapply(CHECKPOINTS, function(g) {
    s <- pspear("EZH2", g, c(MODELS$M3_prolif_ESTIMATE, "ABSOLUTE_purity"), L$TCGA)
    data.frame(model = "M3 + ABSOLUTE purity (TCGA)", cohort = "TCGA", gene = g, n = s[["n"]], r = s[["r"]], p = s[["p"]]) }))
  tp$q <- p.adjust(tp$p, "BH")

  # meta-analysis per gene x model
  meta_tab <- do.call(rbind, lapply(split(tab, list(tab$model, tab$gene)), function(x) {
    m <- meta_r(x$r, x$n, k = length(MODELS[[x$model[1]]]))
    data.frame(model = x$model[1], gene = x$gene[1], r_pooled = m[["r"]], lo = m[["lo"]], hi = m[["hi"]],
               p_pooled = m[["p"]], I2 = m[["I2"]], Q_p = m[["Q_p"]],
               n_cohorts_sig = sum(x$q < 0.05 & sign(x$r) == sign(m[["r"]])),
               same_direction_all = length(unique(sign(x$r))) == 1)
  }))
  meta_tab$q_pooled <- ave(meta_tab$p_pooled, meta_tab$model, FUN = function(p) p.adjust(p, "BH"))
  meta_tab <- meta_tab[order(meta_tab$model, -meta_tab$r_pooled), ]

  write.csv(rbind(tab, tp), file.path(RES_DIR, "Table3a_checkpoint_partial_correlations_by_cohort.csv"), row.names = FALSE)
  write.csv(meta_tab, file.path(RES_DIR, "Table3b_checkpoint_meta_analysis.csv"), row.names = FALSE)

  cat("\n=== Pooled partial correlations (EZH2 vs checkpoint) ===\n")
  print(transform(meta_tab[, c("model", "gene", "r_pooled", "lo", "hi", "p_pooled", "q_pooled", "I2", "n_cohorts_sig")],
                  r_pooled = round(r_pooled, 3), lo = round(lo, 3), hi = round(hi, 3), I2 = round(I2, 2),
                  p_pooled = fmt_p(p_pooled), q_pooled = fmt_p(q_pooled)), row.names = FALSE)
  cat("\n=== LAG3 by cohort ===\n")
  print(transform(subset(rbind(tab, tp), gene == "LAG3"), r = round(r, 3), p = fmt_p(p), q = fmt_p(q)), row.names = FALSE)

  # ---- supporting analyses reported in the text ----
  supp <- do.call(rbind, lapply(names(L), function(k) {
    a1 <- pspear("LAG3", "MCP_T_cell", character(0), L[[k]])        # LAG3 tracks T-cell abundance
    a2 <- pspear("prolif", "LAG3", "EZH2", L[[k]])                  # proliferation beyond EZH2
    a3 <- pspear("EZH2", "LAG3", "prolif", L[[k]])                  # EZH2 beyond proliferation
    data.frame(cohort = k, LAG3_vs_Tcells_rho = a1[["r"]], LAG3_vs_Tcells_p = a1[["p"]],
               prolif_LAG3_given_EZH2_r = a2[["r"]], prolif_LAG3_given_EZH2_p = a2[["p"]],
               EZH2_LAG3_given_prolif_r = a3[["r"]], EZH2_LAG3_given_prolif_p = a3[["p"]]) }))
  write.csv(supp, file.path(RES_DIR, "Table3c_LAG3_supporting_analyses.csv"), row.names = FALSE)
  cat("\n=== LAG3 supporting analyses ===\n"); print(transform(supp, LAG3_vs_Tcells_rho = round(LAG3_vs_Tcells_rho, 3),
    prolif_LAG3_given_EZH2_r = round(prolif_LAG3_given_EZH2_r, 3), EZH2_LAG3_given_prolif_r = round(EZH2_LAG3_given_prolif_r, 3)))

  # ------------------------------ Figure 3 ------------------------------
  lab_models <- c(M0_raw = "Unadjusted", M1_proliferation = "+ Proliferation", M2_ESTIMATE = "+ Immune/stromal (ESTIMATE)",
                  M3_prolif_ESTIMATE = "+ Proliferation + ESTIMATE", M4_prolif_ESTIMATE_Tcells = "+ Prolif. + ESTIMATE + T cells")
  # A: pooled r per gene, unadjusted vs fully adjusted (M3)
  fa <- subset(meta_tab, model %in% c("M0_raw", "M3_prolif_ESTIMATE"))
  fa$model <- factor(lab_models[fa$model], levels = lab_models[c("M0_raw", "M3_prolif_ESTIMATE")])
  fa$gene  <- factor(fa$gene, levels = rev(subset(meta_tab, model == "M3_prolif_ESTIMATE")$gene))
  fa$sig   <- ifelse(fa$q_pooled < 0.05 & fa$I2 < 0.5, "q<0.05 & I2<50%", "not robust")
  pA <- ggplot(fa, aes(r_pooled, gene, colour = model, shape = sig)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), size = 0.35) +
    scale_colour_manual(values = c("grey55", "#C0392B")) + scale_shape_manual(values = c(`not robust` = 1, `q<0.05 & I2<50%` = 16)) +
    labs(x = "Pooled partial Spearman r with EZH2 (95% CI)", y = NULL, colour = NULL, shape = NULL,
         title = "A  EZH2-checkpoint associations, 3 cohorts pooled") + theme_pub + theme(legend.position = "bottom", legend.box = "vertical")

  # B: LAG3 and CD276 per cohort across models
  fb <- subset(tab, gene %in% c("LAG3", "CD276")); fb$cohort <- factor(fb$cohort, levels = names(L))
  fb$model <- factor(lab_models[fb$model], levels = rev(lab_models))
  pB <- ggplot(fb, aes(r, model, colour = cohort)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_point(size = 2.2, position = position_dodge(width = 0.5)) + facet_wrap(~gene) +
    scale_colour_manual(values = c(TCGA = "#2C3E50", CGGA_325 = "#E67E22", CGGA_693 = "#16A085")) +
    labs(x = "Partial Spearman r with EZH2", y = NULL, colour = NULL, title = "B  Effect of confounder adjustment") +
    theme_pub + theme(legend.position = "bottom")

  # C: residual plots (LAG3 vs EZH2 after M3 adjustment)
  resid_df <- do.call(rbind, lapply(names(L), function(k) {
    dd <- na.omit(L[[k]][, c("EZH2", "LAG3", MODELS$M3_prolif_ESTIMATE)]); r <- as.data.frame(lapply(dd, rank))
    Z <- r[, MODELS$M3_prolif_ESTIMATE]
    data.frame(cohort = k, EZH2_res = resid(lm(r$EZH2 ~ ., data = Z)), LAG3_res = resid(lm(r$LAG3 ~ ., data = Z)))
  }))
  lab_c <- subset(tab, gene == "LAG3" & model == "M3_prolif_ESTIMATE")
  resid_df$cohort <- factor(resid_df$cohort, levels = names(L)); lab_c$cohort <- factor(lab_c$cohort, levels = names(L))
  pC <- ggplot(resid_df, aes(EZH2_res, LAG3_res)) + geom_point(size = 0.8, alpha = 0.6) +
    geom_smooth(method = "lm", colour = "#C0392B", linewidth = 0.7) + facet_wrap(~cohort, scales = "free") +
    geom_text(data = lab_c, aes(x = -Inf, y = Inf, label = sprintf("partial r = %.2f\np %s", r, ifelse(p < 1e-4, "< 0.0001", paste("=", fmt_p(p))))),
              hjust = -0.1, vjust = 1.2, size = 3, inherit.aes = FALSE) +
    labs(x = "EZH2 (rank residual)", y = "LAG3 (rank residual)",
         title = "C  EZH2-LAG3 after adjusting for proliferation, immune and stromal content") + theme_pub

  fig3 <- (pA | pB) / pC + plot_layout(heights = c(1.3, 1))
  save_fig(fig3, "Figure3_checkpoints_partial_correlation", 12, 10)
  cat("✓ Figure 3 and Table 3 saved\n")
})

# =====================================================================
# 2d Immune microenvironment
# =====================================================================
step("2d Immune microenvironment")
local({
  set.seed(2026)
  L <- readRDS(file.path(RES_DIR, "cohorts.rds"))$cohorts

  CELL_VARS <- c(ESTIMATE_immune = "ESTIMATE immune score", ESTIMATE_stromal = "ESTIMATE stromal score",
    MCP_T_cell = "T cells", MCP_T_cell_CD8_ = "CD8 T cells", MCP_cytotoxicity_score = "Cytotoxic lymphocytes",
    MCP_NK_cell = "NK cells", MCP_B_cell = "B lineage", MCP_Monocyte = "Monocytic lineage",
    MCP_Myeloid_dendritic_cell = "Myeloid dendritic cells", MCP_Neutrophil = "Neutrophils",
    MCP_Endothelial_cell = "Endothelial cells", MCP_Cancer_associated_fibroblast = "Fibroblasts")
  ADJ <- list(Unadjusted = character(0), `Adjusted for proliferation` = "prolif")

  rows <- list()
  for (a in names(ADJ)) for (k in names(L)) for (v in names(CELL_VARS)) {
    s <- pspear("EZH2", v, ADJ[[a]], L[[k]])
    rows[[length(rows) + 1]] <- data.frame(adjustment = a, cohort = k, variable = v, label = CELL_VARS[[v]],
                                           n = s[["n"]], r = s[["r"]], p = s[["p"]])
  }
  tab <- do.call(rbind, rows)
  tab$q <- ave(tab$p, tab$adjustment, tab$cohort, FUN = function(p) p.adjust(p, "BH"))
  meta <- do.call(rbind, lapply(split(tab, list(tab$adjustment, tab$variable)), function(x) {
    m <- meta_r(x$r, x$n, length(ADJ[[x$adjustment[1]]]))
    data.frame(adjustment = x$adjustment[1], variable = x$variable[1], label = x$label[1],
               TCGA = x$r[x$cohort == "TCGA"], CGGA_325 = x$r[x$cohort == "CGGA_325"], CGGA_693 = x$r[x$cohort == "CGGA_693"],
               r_pooled = m[["r"]], p_pooled = m[["p"]], I2 = m[["I2"]], same_direction = length(unique(sign(x$r))) == 1)
  }))
  meta$q_pooled <- ave(meta$p_pooled, meta$adjustment, FUN = function(p) p.adjust(p, "BH"))
  meta$robust <- meta$q_pooled < 0.05 & meta$I2 < 0.5 & meta$same_direction
  write.csv(tab,  file.path(RES_DIR, "Table_S3a_immune_by_cohort.csv"), row.names = FALSE)
  write.csv(meta, file.path(RES_DIR, "Table_S3b_immune_meta.csv"), row.names = FALSE)
  print(transform(meta[order(meta$adjustment, meta$r_pooled), c("adjustment", "label", "TCGA", "CGGA_325", "CGGA_693", "r_pooled", "q_pooled", "I2", "robust")],
                  TCGA = round(TCGA, 2), CGGA_325 = round(CGGA_325, 2), CGGA_693 = round(CGGA_693, 2),
                  r_pooled = round(r_pooled, 2), q_pooled = fmt_p(q_pooled), I2 = round(I2, 2)), row.names = FALSE)

  tab$label <- factor(tab$label, levels = rev(unname(CELL_VARS)))
  tab$star  <- ifelse(tab$q < 0.05, "*", "")
  pS <- ggplot(tab, aes(cohort, label, fill = r)) + geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f%s", r, star)), size = 2.8) +
    facet_wrap(~adjustment) +
    scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C0392B", limits = c(-0.7, 0.7)) +
    labs(x = NULL, y = NULL, fill = "Spearman r\nwith EZH2",
         title = "EZH2 and the tumor microenvironment (ESTIMATE / MCP-counter)",
         caption = "* BH q < 0.05 within cohort") + theme_pub + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig(pS, "FigureS1_immune_landscape", 9, 6)
  cat("✓ Figure S1 and Table S3 saved\n")
})

# =====================================================================
# 2e Single-cell expression (TISCH2)
# =====================================================================
step("2e Single-cell expression (TISCH2)")
local({
  set.seed(2026)
  read_tisch <- function(g) {
    x <- read.csv(file.path(DATA_DIR, paste0("TISCH_", g, "_heatmap.csv")), check.names = FALSE, fileEncoding = "UTF-8-BOM")
    colnames(x) <- c("cell_type", "dataset", "value"); x$gene <- g; x
  }
  sc <- do.call(rbind, lapply(c("EZH2", "LAG3", "PDCD1"), read_tisch))
  DATASETS <- c("Glioma_GSE131928_Smartseq2" = "Neftel Smart-seq2")
  CELLS <- c("CD8Tex" = "Exhausted CD8 T", "Mono/Macro" = "Monocyte/macrophage",
             "AC-like Malignant" = "Malignant AC-like", "MES-like Malignant" = "Malignant MES-like",
             "NPC-like Malignant" = "Malignant NPC-like", "OPC-like Malignant" = "Malignant OPC-like")
  sc <- subset(sc, dataset %in% names(DATASETS) & cell_type %in% names(CELLS))
  sc$dataset   <- factor(DATASETS[sc$dataset], levels = DATASETS)
  sc$cell_type <- factor(CELLS[sc$cell_type], levels = CELLS)
  sc$gene      <- factor(sc$gene, levels = c("PDCD1", "LAG3", "EZH2"))
  write.csv(sc, file.path(RES_DIR, "Table_S4_single_cell_TISCH2.csv"), row.names = FALSE)
  print(reshape(sc[, c("dataset", "gene", "cell_type", "value")], idvar = c("dataset", "gene"), timevar = "cell_type", direction = "wide"))

  p4 <- ggplot(sc, aes(cell_type, gene, fill = value)) + geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3.2) +
    scale_fill_gradient(low = "grey97", high = "#B03A2E", name = "Mean\nlog(TPM/10+1)") +
    labs(x = NULL, y = NULL, title = "Cell-type expression in adult GBM single-cell data (GSE131928, Smart-seq2)") +
    theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_text(face = "italic"))
  save_fig(p4, "Figure4_single_cell", 7.5, 3.4)
  cat("✓ Figure 4 saved\n")
})

# =====================================================================
# 2f Hallmark proliferation sensitivity
# =====================================================================
step("2f Hallmark proliferation sensitivity")
local({
  set.seed(2026)
  suppressPackageStartupMessages(library(data.table))
  if (!requireNamespace("GSVA", quietly = TRUE)) stop("Package GSVA needed: install.packages('BiocManager'); BiocManager::install('GSVA')")
  say <- function(...) cat(sprintf(...), "\n", sep = "")
  fmt_p <- function(p) ifelse(p < 1e-4, "<0.0001", formatC(p, format = "g", digits = 2))

  # ---- helpers (same partial Spearman and Fisher-z pooling as the main analysis; + random effects) ----
  pspear <- function(x, y, Z = character(0), data) {
    dd <- na.omit(data[, c(x, y, Z), drop = FALSE]); r <- as.data.frame(lapply(dd, rank))
    n <- nrow(dd); k <- length(Z)
    rr <- if (k == 0) cor(r[[x]], r[[y]]) else
      cor(resid(lm(r[[x]] ~ ., data = r[, Z, drop = FALSE])), resid(lm(r[[y]] ~ ., data = r[, Z, drop = FALSE])))
    tt <- rr * sqrt((n - 2 - k) / (1 - rr^2))
    c(n = n, r = rr, p = 2 * pt(-abs(tt), n - 2 - k))
  }
  meta_r <- function(r, n, k = 0) {
    z <- atanh(r); w <- n - 3 - k; zp <- sum(w * z) / sum(w); se <- 1 / sqrt(sum(w))
    Q <- sum(w * (z - zp)^2); df <- length(r) - 1
    tau2 <- if (df > 0) max(0, (Q - df) / (sum(w) - sum(w^2) / sum(w))) else 0
    wr <- 1 / (1 / w + tau2); zr <- sum(wr * z) / sum(wr); ser <- 1 / sqrt(sum(wr))
    c(r = tanh(zp), lo = tanh(zp - 1.96 * se), hi = tanh(zp + 1.96 * se), p = 2 * pnorm(-abs(zp / se)),
      I2 = ifelse(Q > 0, max(0, (Q - df) / Q), 0),
      r_RE = tanh(zr), lo_RE = tanh(zr - 1.96 * ser), hi_RE = tanh(zr + 1.96 * ser))
  }
  collapse_symbols <- function(mat, sym) {
    ok <- !is.na(sym) & sym != ""; mat <- mat[ok, , drop = FALSE]; sym <- sym[ok]
    o <- order(rowMeans(mat), decreasing = TRUE); mat <- mat[o, , drop = FALSE]; sym <- sym[o]
    keep <- !duplicated(sym); mat <- mat[keep, , drop = FALSE]; rownames(mat) <- sym[keep]; mat
  }

  # ---- data ----
  co <- readRDS(file.path("results", "cohorts.rds"))$cohorts
  tpm <- readRDS(file.path("share", "tcga_tpm.rds")); gi <- readRDS(file.path("share", "tcga_gene_info.rds"))
  mats <- list(TCGA = collapse_symbols(as.matrix(tpm), gi$gene_name[match(rownames(tpm), rownames(gi))]))
  rm(tpm); invisible(gc())
  for (id in c("325", "693")) {
    ex <- readRDS(file.path("share", paste0("cgga", id, "_expr_gbm.rds")))
    m <- as.matrix(ex[, -1]); storage.mode(m) <- "numeric"
    mats[[paste0("CGGA_", id)]] <- collapse_symbols(m, as.character(ex[[1]]))
  }

  # ---- published cell-cycle signature (EZH2 removed) ----
  if (!exists("get_hallmark")) get_hallmark <- function() {
    if (!requireNamespace("msigdbr", quietly = TRUE)) stop("Package msigdbr needed: install.packages('msigdbr')")
    h <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
                  error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
    unique(h$gene_symbol[h$gs_name %in% c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")])
  }
  cc_genes <- setdiff(get_hallmark(), c("EZH2", "LAG3", "PDCD1", "BTLA", "CD276"))
  say("Cell-cycle signature: Hallmark E2F targets + G2M checkpoint, %d genes after removing EZH2 (and target genes)", length(cc_genes))

  for (k in names(co)) {
    d <- co[[k]]; m <- mats[[k]][, d$sample]; lm2 <- log2(m + 1); lm2 <- lm2[apply(lm2, 1, var) > 0, ]
    g <- intersect(cc_genes, rownames(lm2))
    es <- suppressMessages(GSVA::gsva(GSVA::ssgseaParam(lm2, list(cc = g)), verbose = FALSE))
    d$cc_ssgsea <- as.numeric(es["cc", ])
    d$cc_meanz  <- colMeans(t(scale(t(lm2[g, , drop = FALSE]))))
    co[[k]] <- d
    say("  [%s] %d/%d signature genes found | rho(EZH2, cell-cycle ssGSEA) = %.2f | rho(EZH2, 10-gene score) = %.2f | rho(10-gene, cell-cycle) = %.2f",
        k, length(g), length(cc_genes), cor(d$EZH2, d$cc_ssgsea, method = "spearman"),
        cor(d$EZH2, d$prolif, method = "spearman"), cor(d$prolif, d$cc_ssgsea, method = "spearman"))
  }

  # ---- models ----
  EST <- c("ESTIMATE_immune", "ESTIMATE_stromal")
  MODELS <- list(
    "Original: 10-gene + ESTIMATE"                 = c("prolif", EST),
    "Hallmark cell-cycle (ssGSEA) + ESTIMATE"      = c("cc_ssgsea", EST),
    "Hallmark cell-cycle (mean z) + ESTIMATE"      = c("cc_meanz", EST),
    "10-gene + Hallmark (ssGSEA) + ESTIMATE"       = c("prolif", "cc_ssgsea", EST),
    "MKI67 alone + ESTIMATE"                       = c("MKI67", EST),
    "10-gene + Hallmark + ESTIMATE + T cells"      = c("prolif", "cc_ssgsea", EST, "MCP_T_cell"))
  GENES <- c("LAG3", "PDCD1", "BTLA", "CD276")

  rows <- list()
  for (g in GENES) for (mn in names(MODELS)) {
    Z <- MODELS[[mn]]
    per <- rbindlist(lapply(names(co), function(k) { v <- pspear("EZH2", g, Z, co[[k]])
      data.table(cohort = k, n = v[["n"]], r = v[["r"]], p = v[["p"]]) }))
    m <- meta_r(per$r, per$n, length(Z))
    rows[[length(rows) + 1]] <- data.table(gene = g, model = mn,
      r_TCGA = per$r[1], r_CGGA_325 = per$r[2], r_CGGA_693 = per$r[3],
      p_TCGA = per$p[1], p_CGGA_325 = per$p[2], p_CGGA_693 = per$p[3],
      r_pooled = m[["r"]], lo = m[["lo"]], hi = m[["hi"]], p_pooled = m[["p"]], I2 = m[["I2"]],
      r_RE = m[["r_RE"]], lo_RE = m[["lo_RE"]], hi_RE = m[["hi_RE"]],
      same_direction = length(unique(sign(per$r))) == 1)
  }
  res <- rbindlist(rows)
  res[, q_pooled := p.adjust(p_pooled, "BH"), by = model]

  cat("\n---- EZH2 vs checkpoints: partial Spearman r (per cohort: TCGA CGGA_325 CGGA_693) ----\n")
  for (i in seq_len(nrow(res))) { t <- res[i]
    say("%-6s %-42s | %5.2f %5.2f %5.2f | pooled r = %5.2f [%5.2f, %5.2f] q = %-7s I2 = %3.0f%% | RE r = %5.2f [%5.2f, %5.2f]",
        t$gene, t$model, t$r_TCGA, t$r_CGGA_325, t$r_CGGA_693, t$r_pooled, t$lo, t$hi, fmt_p(t$q_pooled),
        100 * t$I2, t$r_RE, t$lo_RE, t$hi_RE) }

  # reverse check: does the Hallmark score itself track LAG3 once EZH2 is accounted for?
  cat("\n---- Reverse check: Hallmark cell-cycle score vs LAG3, adjusted for EZH2 + ESTIMATE ----\n")
  for (k in names(co)) { v <- pspear("cc_ssgsea", "LAG3", c("EZH2", EST), co[[k]])
    say("  %-9s r = %5.2f  P = %s", k, v[["r"]], fmt_p(v[["p"]])) }

  fwrite(res, file.path("results", "Table_S5_proliferation_sensitivity.csv"))
  cat("\nOK results/Table_S5_proliferation_sensitivity.csv saved\n")
})

cat("\n\u2713 02_main_analyses.R complete.\n")
}, finally = end_log())
