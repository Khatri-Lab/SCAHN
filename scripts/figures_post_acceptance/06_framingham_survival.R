# ==============================================================================
# Framingham survival, lab correlations, disease effect sizes and GSE184050
#
# Cox proportional hazards models of the neutrophil signature scores against
# 5-, 7- and 10-year mortality in the Framingham cohort, Spearman correlation of
# the scores with plasma lab measures, effect sizes for disease,
# and a longitudinal T2D validation in GSE184050.
#
# Requires 00_setup.R (packages, figure settings, palettes) and 01_load_data.R
# (gene_sets, genes.* signature vectors, dir.data / dir.results / dir.fig).
# Also uses ggbio and ggfortify, loaded below.
#
# Reads, on top of what 01_load_data.R loads:
#   results/Framingham/FH.score_pheno.csv
#   results/Framingham/FH.pheno.csv
#   results/Framingham/FH.diseaseES.csv
#   data/GSE184050_exprlist.rds
#
# Writes to figures/post_acceptance/:
#   Figure6.pdf   Figure 6
#   Fig.e8.pdf    Extended figure 8 (stratified correlations, Cox model tables)
# ==============================================================================

SCAHN_SCRIPTS = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/scripts/figures_post_acceptance"
source(file.path(SCAHN_SCRIPTS, "00_setup.R"))
source(file.path(SCAHN_SCRIPTS, "01_load_data.R"))


# -- Framingham functions for survival -----------------------------------------

library(ggbio)
library(ggfortify)

print_cox_model_summary = function(cox_model, test = c("wald", "lrt"), data = NULL) {
  test = match.arg(test)
  summary_cox_model = summary(cox_model)
  summary_cox_model_data_frame =
    as.data.frame(summary_cox_model$coefficients)[, c(2, 5)]
  colnames(summary_cox_model_data_frame) = c("hazard ratio", "pval")
  summary_cox_model_data_frame$`hazard ratio` = paste0(
    round(summary_cox_model_data_frame$`hazard ratio`, 2), " (",
    round(summary_cox_model$conf.int[, 3], 2), ", ", round(summary_cox_model$conf.int[, 4], 2), ")")

  if (test == "lrt") {
    if (is.null(data)) stop("'data' must be provided when test = 'lrt'")
    n_coefs = nrow(summary_cox_model$coefficients)
    lrt_pvals = numeric(n_coefs)
    coef_names = rownames(summary_cox_model$coefficients)
    model_terms = attr(terms(cox_model), "term.labels")

    for (i in seq_len(n_coefs)) {
      matched_term = NULL
      for (tm in model_terms) {
        if (startsWith(coef_names[i], tm) || coef_names[i] == tm) {
          matched_term = tm
          break
        }
      }
      if (is.null(matched_term)) {
        lrt_pvals[i] = summary_cox_model$coefficients[i, 5]
        next
      }
      reduced_terms = setdiff(model_terms, matched_term)
      if (length(reduced_terms) == 0L) {
        null_ll = cox_model$loglik[1]
      } else {
        resp_str = deparse(formula(cox_model)[[2]])
        fml_reduced = as.formula(paste0(resp_str, " ~ ", paste(reduced_terms, collapse = " + ")))
        fit_reduced = coxph(fml_reduced, data = data)
        null_ll = fit_reduced$loglik[2]
      }
      lrt_stat = 2 * (cox_model$loglik[2] - null_ll)
      lrt_pvals[i] = pchisq(lrt_stat, df = 1, lower.tail = FALSE)
    }
    summary_cox_model_data_frame$pval = lrt_pvals
  }

  summary_cox_model_data_frame$pval = as.character(ifelse(
    as.numeric(summary_cox_model_data_frame$pval) < 0.01,
    format(as.numeric(summary_cox_model_data_frame$pval), digits = 2, scientific = TRUE),
    round(as.numeric(summary_cox_model_data_frame$pval), 3)))
  return(summary_cox_model_data_frame)
}

### full cohort: score + age + sex + BMI + disease + smoking
evaluate_framingham_mortality_models = function(n_years_cutoff, score_pheno_matrix,
                                                score_name, score_name_print) {
  cutoff_days = 365 * n_years_cutoff

  # ---- Censoring at n_years_cutoff ----
  score_pheno_matrix$status =
    as.integer(!is.na(score_pheno_matrix$death_date) & score_pheno_matrix$death_date <= cutoff_days)
  score_pheno_matrix$death_date =
    pmin(score_pheno_matrix$death_date, cutoff_days, na.rm = TRUE) / 365
  score_pheno_matrix$disease_status = as.factor(1 - score_pheno_matrix$healthy)

  cat(sprintf("  Events: %d  |  Censored: %d\n",
              sum(score_pheno_matrix$status), sum(!score_pheno_matrix$status)))

  # ---- Full Cox model (score + covariates) ----
  formula_full = as.formula(sprintf(
    paste0("Surv(death_date, status) ~ %s + age + as.factor(sex) + bmi + ",
           "disease_status + as.factor(regular_smokers)"), score_name))
  formula_no_score = update(formula_full, paste(". ~ . -", score_name))

  cox_full     = coxph(formula_full,     data = score_pheno_matrix)
  cox_no_score = coxph(formula_no_score, data = score_pheno_matrix)
  ph_test      = cox.zph(cox_full)

  cox_full_summary = print_cox_model_summary(cox_full, test = "wald", data = score_pheno_matrix)
  rownames(cox_full_summary) = c(score_name_print, "age", "sex", "BMI", "disease", "smoking")

  lrt_p_all = pchisq(2 * (cox_full$loglik[2] - cox_no_score$loglik[2]), df = 1, lower.tail = FALSE)

  # ---- Covariate-adjusted residual score ----
  lm_adjust = lm(as.formula(sprintf(
    "%s ~ age + as.factor(sex) + bmi + disease_status + as.factor(regular_smokers)",
    score_name)), data = score_pheno_matrix)

  score_pheno_matrix$residuals = NA_real_
  score_pheno_matrix$residuals[as.integer(names(lm_adjust$residuals))] =
    lm_adjust$residuals
  score_pheno_matrix$score_discretized = factor(
    ifelse(is.na(score_pheno_matrix$residuals), NA,
           ifelse(score_pheno_matrix$residuals >= 0, "high", "low")), levels = c("low", "high"))

  # ---- Residual Cox models ----
  cox_resid = coxph(Surv(death_date, status) ~ residuals, data = score_pheno_matrix)
  cox_resid_disc = coxph(Surv(death_date, status) ~ score_discretized, data = score_pheno_matrix)

  cox_resid_summary = print_cox_model_summary(cox_resid, test = "lrt", data = score_pheno_matrix)
  rownames(cox_resid_summary) = paste("adjusted", score_name_print)

  cox_resid_disc_summary = print_cox_model_summary(cox_resid_disc, test = "lrt",
                                                   data = score_pheno_matrix)
  rownames(cox_resid_disc_summary) = paste("adjusted discretized", score_name_print)

  lrt_p_residual = pchisq(2 * diff(cox_resid$loglik),      df = 1, lower.tail = FALSE)
  lrt_p_discretized = pchisq(2 * diff(cox_resid_disc$loglik), df = 1, lower.tail = FALSE)

  # ---- KM curve ----
  km_fit = survfit(Surv(death_date, status) ~ score_discretized, data = score_pheno_matrix)

  list(
    cox_model_all_var_summary              = cox_full_summary,
    lrt_p_all                              = lrt_p_all,
    ph_test                                = ph_test,
    residual_cox_model_summary             = cox_resid_summary,
    lrt_p_residual                         = lrt_p_residual,
    residual_cox_model_discretized_summary = cox_resid_disc_summary,
    lrt_p_discretized = lrt_p_discretized, km_fit = km_fit)
}

### within one sex: same model without the sex term, and Wald p values
evaluate_framingham_mortality_models.sex = function(n_years_cutoff, score_pheno_matrix,
                                                    score_name, score_name_print) {
  cutoff_days = 365 * n_years_cutoff

  # ---- Censoring at n_years_cutoff ----
  score_pheno_matrix$status =
    as.integer(!is.na(score_pheno_matrix$death_date) & score_pheno_matrix$death_date <= cutoff_days)
  score_pheno_matrix$death_date =
    pmin(score_pheno_matrix$death_date, cutoff_days, na.rm = TRUE) / 365
  score_pheno_matrix$disease_status = as.factor(1 - score_pheno_matrix$healthy)

  cat(sprintf("  Events: %d  |  Censored: %d\n",
              sum(score_pheno_matrix$status), sum(!score_pheno_matrix$status)))

  # ---- Full Cox model (score + covariates) ----
  formula_full = as.formula(sprintf(
    "Surv(death_date, status) ~ %s + age +  bmi + disease_status + as.factor(regular_smokers)",
    score_name))
  formula_no_score = update(formula_full, paste(". ~ . -", score_name))

  cox_full     = coxph(formula_full,     data = score_pheno_matrix)
  cox_no_score = coxph(formula_no_score, data = score_pheno_matrix)
  ph_test      = cox.zph(cox_full)

  cox_full_summary = print_cox_model_summary(cox_full, test = "wald", data = score_pheno_matrix)
  rownames(cox_full_summary) = c(score_name_print, "age", "BMI", "disease", "smoking")

  lrt_p_all = pchisq(2 * (cox_full$loglik[2] - cox_no_score$loglik[2]), df = 1, lower.tail = FALSE)

  # ---- Covariate-adjusted residual score ----
  lm_adjust = lm(as.formula(sprintf(
    "%s ~ age + bmi + disease_status + as.factor(regular_smokers)",
    score_name)), data = score_pheno_matrix)

  score_pheno_matrix$residuals = NA_real_
  score_pheno_matrix$residuals[as.integer(names(lm_adjust$residuals))] =
    lm_adjust$residuals
  score_pheno_matrix$score_discretized = factor(
    ifelse(is.na(score_pheno_matrix$residuals), NA,
           ifelse(score_pheno_matrix$residuals >= 0, "high", "low")), levels = c("low", "high"))

  # ---- Residual Cox models ----
  cox_resid = coxph(Surv(death_date, status) ~ residuals, data = score_pheno_matrix)
  cox_resid_disc = coxph(Surv(death_date, status) ~ score_discretized, data = score_pheno_matrix)

  cox_resid_summary = print_cox_model_summary(cox_resid)
  rownames(cox_resid_summary) = paste("adjusted", score_name_print)

  cox_resid_disc_summary = print_cox_model_summary(cox_resid_disc)
  rownames(cox_resid_disc_summary) = paste("adjusted discretized", score_name_print)

  lrt_p_residual = pchisq(2 * diff(cox_resid$loglik),      df = 1, lower.tail = FALSE)
  lrt_p_discretized = pchisq(2 * diff(cox_resid_disc$loglik), df = 1, lower.tail = FALSE)

  # ---- KM curve ----
  km_fit = survfit(Surv(death_date, status) ~ score_discretized, data = score_pheno_matrix)

  list(
    cox_model_all_var_summary              = cox_full_summary,
    lrt_p_all                              = lrt_p_all,
    ph_test                                = ph_test,
    residual_cox_model_summary             = cox_resid_summary,
    lrt_p_residual                         = lrt_p_residual,
    residual_cox_model_discretized_summary = cox_resid_disc_summary,
    lrt_p_discretized = lrt_p_discretized, km_fit = km_fit)
}

### run every score at 5, 7 and 10 years, and flatten to one table
survival_framingham = function(basedata, bysex = FALSE) {
  survival_summary = list()
  survival_dt = list()
  for (year in c(5, 7, 10)) {
    survival_summary[[paste0("year", year)]] = list()
    for (score in c(scores_to_test, scores_to_test_z)) {

      if (!bysex) {
        cox_model_summaries =
          evaluate_framingham_mortality_models(year, basedata, score, score)
      } else {
        cox_model_summaries =
          evaluate_framingham_mortality_models.sex(year, basedata, score, score)
      }
      survival_summary[[paste0("year", year)]][[score]] = cox_model_summaries
    }

    survival_dt[[paste0("year", year)]] = rbindlist(lapply(
      names(survival_summary[[paste0("year", year)]]), function(gene) {
        s = survival_summary[[paste0("year", year)]][[gene]]

        row_all  = s$cox_model_all_var_summary[gene, ]
        row_res  = s$residual_cox_model_summary[1, ]
        row_disc = s$residual_cox_model_discretized_summary[1, ]

        data.table(
          gene = gene, hr_all = row_all[["hazard ratio"]],
          pval_all           = as.numeric(row_all[["pval"]]),
          lrt_p_all          = as.numeric(s$lrt_p_all),
          hr_residual        = row_res[["hazard ratio"]],
          pval_residual      = as.numeric(row_res[["pval"]]),
          lrt_p_residual     = as.numeric(s$lrt_p_residual),
          hr_discretized     = row_disc[["hazard ratio"]],
          pval_discretized   = as.numeric(row_disc[["pval"]]),
          lrt_p_discretized = as.numeric(s$lrt_p_discretized), year = year)
      }))
  }
  survival_dt = do.call(rbind, survival_dt)
  return(list(survival_summary = survival_summary, survival_dt = survival_dt))
}

theme_legend_right_box = theme_legend_right + theme(legend.box = "vertical")

### forest plot of hazard ratios, one row per score, shape per follow-up length
plot_forest_framingham = function(survival_res, hr_col = "hr_all",
                                  pval_col = "lrt_p_all", scores = NULL) {
  dt = copy(survival_res$survival_dt)
  if (!is.null(scores)) dt = dt[dt$gene %in% scores, ]

  dt$gene = gsub("_z", "", dt$gene)

  # FDR per year (across all scores in that year)
  #dt = dt[, fdr := p.adjust(.SD[[pval_col]], method = "BH"), by = year]
  dt = dt[, fdr := p.adjust(.SD[[pval_col]], method = "BH")]

  # Parse "HR (lo, hi)" string
  hr_str = dt[[hr_col]]
  dt$hr  = as.numeric(sub(" \\(.*", "", hr_str))
  dt$lo  = as.numeric(sub(".*\\(([^,]+),.*", "\\1", hr_str))
  dt$hi  = as.numeric(sub(".*,\\s*([^)]+)\\)", "\\1", hr_str))

  dt$sig = ifelse(dt$fdr <= 0.05, "FDR ≤ 0.05", ifelse(dt$fdr <= 0.1, "FDR ≤ 0.1", "ns"))
  dt$sig = factor(dt$sig, levels = c("FDR ≤ 0.05", "FDR ≤ 0.1", "ns"))

  dt$year = factor(dt$year, levels = rev(c(5, 7, 10)))

  # Order genes by 5-year HR
  gene_order = dt[dt$year == 5, ][order(hr), gene]
  missing    = setdiff(unique(dt$gene), gene_order)
  dt$gene    = factor(dt$gene, levels = c(gene_order, missing))

  shape_map  = c("5" = 15, "7" = 16, "10" = 18)   # filled circle, filled circle, diamond
  dodge      = position_dodge(width = 0.6)

  ggplot(dt, aes(x = hr, y = gene, color = sig, shape = year, group = year)) +
    geom_vline(xintercept = 1, linetype = "solid", color = "grey75", linewidth = 0.3,
               alpha = 0.5) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.4, size = 0.3, position = dodge) +
    geom_point(aes(fill = sig), size = 1, position = dodge, stroke = 0.4) +
    scale_x_log10(breaks = seq(0.8, 1.4, 0.2)) +
    scale_y_discrete(labels = cap_first) +
    scale_shape_manual(values = shape_map, name = "Follow-up",
                       labels = c("10" = "10-year mortality", "7" = "7-year mortality",
                                  "5"  = "5-year mortality")) +
    scale_color_manual(values = c("FDR ≤ 0.05" = "#C0392B", "FDR ≤ 0.1" = "#E67E22",
                                  "ns" = "grey50"), name = "Significance") +
    scale_fill_manual(values = c("FDR ≤ 0.05" = "#C0392B", "FDR ≤ 0.1" = "#E67E22",
                                  "ns" = "grey50"), name = "Significance") +
    labs(x = "Hazard ratio", y = NULL) +
    theme_expresso(legend_text_size = textsize, axis_text_size = textsize,
                   plot_background = "white", panel_background = "white") +
    theme(panel.grid.major.y = element_blank()) +
    theme_legend_right_box +
    guides(shape = guide_legend(reverse = TRUE, ncol = 1, order = 1),
           color = guide_legend(ncol = 1, order = 2,
                                override.aes = list(shape = NA, linewidth = 0.6)),
           fill  = "none")
}

### Kaplan-Meier curve for the discretized covariate-adjusted score
plot_km_residual = function(survival_summary, n_year, score, col_survival) {
  res     = survival_summary[[paste0("year", n_year)]][[score]]
  km_fit  = res$km_fit
  cox_row = res$residual_cox_model_discretized_summary[1, ]

  hr_str   = cox_row[["hazard ratio"]]
  pval_num = as.numeric(cox_row[["pval"]])
  pval_fmt = ifelse(pval_num < 0.001, formatC(pval_num, digits = 2, format = "e"),
                    round(pval_num, 3))
  label = sprintf("HR: %s\np = %s", hr_str, pval_fmt)

  autoplot(km_fit, conf.int = TRUE) +
    labs(color = paste0(score, " adjusted\n(discretized)"),
         fill  = paste0(score, " adjusted\n(discretized)")) +
    labs(x = "Time (years)", y = "Survival (%)", title = cap_first(score)) +
    scale_y_continuous(labels = label_percent(suffix = "")) +
    annotate("text", x = 0, y = -Inf, label = label, hjust = 0, vjust = -0.2, lineheight = 0.8,
             size = pt2mm(textsize), color = textcolor) +
    scale_fill_manual(values   = col_survival) +
    scale_colour_manual(values = col_survival) +
    theme_bw() +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
          panel.grid.major.y   = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.border         = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          plot.margin          = margin(0, 0, 0, 0, "line"),
          plot.title = element_text(size = textsize, color = textcolor,
                                    margin = margin(b = 1)),
          axis.text = element_text(size = textsize, color = textcolor),
          axis.title = element_text(size = textsize, color = textcolor),
          axis.title.x = element_text(margin = margin(t = 1)),
          axis.title.y = element_text(margin = margin(r = 1)),
          axis.ticks = element_line(color = linecolor, linewidth = linesize),
          axis.line = element_line(color = linecolor, linewidth = linesize))
}

### the full Cox model table, drawn as a grob in extended figure 8
make_cox_table = function(surv_obj, score, n_year, row_label = NULL) {
  dt = surv_obj$survival_summary[[paste0("year", n_year)]][[score]]$cox_model_all_var_summary
  if (is.null(row_label)) row_label = score
  rownames(dt)[1] = row_label
  rownames(dt) = cap_first(rownames(dt))
  colnames(dt) = replace(cap_first(colnames(dt)), colnames(dt) == "pval", "pval")

  tg = tableGrob(dt, theme = ttheme_default(
      core = list(
        fg_params = list(fontsize = textsize, fontfamily = "Arial",
                         fontface = "plain", hjust = 1, x = 0.99),
        bg_params = list(fill = "white"), padding = unit(c(0.35, 0.3), "line")), colhead = list(
        fg_params = list(fontsize = textsize, fontfamily = "Arial", fontface = "plain"),
        bg_params = list(fill = "white"), padding = unit(c(0.35, 0.3), "line")), rowhead = list(
        fg_params = list(fontsize = textsize, fontfamily = "Arial", fontface = "plain"),
        padding = unit(c(0.3, 0.3), "line"))))

  gtable_add_grob(
    tg, grobs = replicate(nrow(tg), segmentsGrob(
      x0 = unit(0, "npc"), x1 = unit(1, "npc"), y0 = unit(0, "npc"), y1 = unit(0, "npc"),
      gp = gpar(col = linecolor, lwd = linesize * ggplot2::.pt)), simplify = FALSE),
    t = seq_len(nrow(tg)), l = 1, r = ncol(tg))
}


# -- Framingham survival -------------------------------------------------------

col_survival = c("#214DC8", "#FF7F0E")

dichromat_pal("DarkRedtoBlue.12")(12)[c(2, 11)]

### all subjects; age, BMI and the scores are z-scored before fitting
framingham_score_pheno = fread(paste0(
  dir.results, "Framingham/FH.score_pheno.csv"))
framingham_score_pheno$immature_mature_ratio =
  framingham_score_pheno$immature / framingham_score_pheno$mature

base_rows   = c("age", "sex", "BMI", "smoking")
scores_to_test = c("immature", "mature", "degranulating", "antiprotease",
                   "total", "IFN", "immature_mature_ratio")
scores_to_test_z = paste0(scores_to_test, "_z")

framingham_score_pheno[, age := as.numeric(scale(age))]
framingham_score_pheno[, bmi := as.numeric(scale(bmi))]

framingham_score_pheno[, (scores_to_test_z) := lapply(.SD, function(x) as.numeric(scale(x))),
                       .SDcols = scores_to_test]

surv.all = survival_framingham(framingham_score_pheno, bysex = FALSE)

### forest plot
FH.survival.forest = plot_forest_framingham(surv.all, scores = c(scores_to_test_z[1:6]))

surv.all$survival_summary$year5

df = as.data.table(surv.all$survival_dt)[grep("_z", gene)]
df$hr_all.numeric = as.numeric(gsub(" \\(.*", "", df$hr_all))

### the same models within each sex
surv.F = survival_framingham(framingham_score_pheno[sex == "F"], bysex = TRUE)
surv.M = survival_framingham(framingham_score_pheno[sex == "M"], bysex = TRUE)

surv.F$survival_dt$sex = "F"
surv.M$survival_dt$sex = "M"

# -- Framingham correlations ---------------------------------------------------

framingham_score_pheno = fread(paste0(
  dir.results, "Framingham/FH.score_pheno.csv"))
framingham_score_pheno2 = fread(paste0(
  dir.results, "Framingham/FH.pheno.csv"))
framingham_score_pheno2 = framingham_score_pheno2[
  order(match(shareid, framingham_score_pheno$shareid))]
all.equal(framingham_score_pheno2$shareid, framingham_score_pheno$shareid)
framingham_score_pheno2 = cbind(
  framingham_score_pheno2, framingham_score_pheno[, c(scores_to_test[1:6],
                             paste0(scores_to_test[1:6], "_adj")), with = FALSE])

framingham_score_pheno2$disease = as.factor(1 - framingham_score_pheno2$healthy)
framingham_score_pheno2$LDL_plasma_estimated =
  framingham_score_pheno2$total_cholestrol_plasma -
  framingham_score_pheno2$HDL_cholestrol_plasma -
  framingham_score_pheno2$triglycerides_plasma / 5
list_of_lab_scores = c("glucose_plasma", "A1C", "triglycerides_plasma",
                       "albumin_urine", "creatinine_serum", "HDL_cholestrol_plasma")

framingham_score_pheno2$diabetes_yesno = factor(framingham_score_pheno2$diabetes)
framingham_score_pheno2[, diabetes_yesno := ifelse(diabetes %in% c(1), "diabetes", "no diabetes")]
framingham_score_pheno2[, BMI_categ := ifelse(
  bmi > median(framingham_score_pheno2$bmi, na.rm = TRUE), "higher BMI", "lower BMI")]
framingham_score_pheno2[, sex_more := ifelse(sex %in% c("F"), "female", "male")]

### Spearman correlation of each score against the plasma lab measures
Framingham_correlation_ht = function(dtf, ptitle, list_of_gene_scores = scores,
                                     fdr_threshold = 0.05, circlesize = 0.5) {
  list_of_lab_scores = c("glucose_plasma", "A1C", "triglycerides_plasma", "HDL_cholestrol_plasma")

  corr_estimate = matrix(NA_real_, nrow = length(list_of_gene_scores),
                         ncol = length(list_of_lab_scores), dimnames = list(list_of_gene_scores,
                                         list_of_lab_scores))
  corr_p_value  = corr_estimate

  for (gene_score in list_of_gene_scores) {
    for (lab_score in list_of_lab_scores) {
      ct = cor.test(dtf[[lab_score]], dtf[[gene_score]], method = "spearman")
      corr_estimate[gene_score, lab_score] = ct$estimate
      corr_p_value[gene_score, lab_score]  = ct$p.value
    }
  }

  # FDR correction across all tests
  fdr_mat = matrix(p.adjust(as.vector(corr_p_value), method = "BH"),
                   nrow = nrow(corr_p_value), ncol = ncol(corr_p_value),
                   dimnames = dimnames(corr_p_value))

  col_labels = cap_first(c("glucose", "HbA1C", "triglycerides", "HDL"))

  # Circle sizes: large = significant; outline drawn at fdr_threshold boundary
  sz       = 1.2 - fdr_mat
  fdr_clip = fdr_mat; fdr_clip[fdr_clip > fdr_threshold] = 1.2
  sz2      = 1.2 - fdr_clip

  clamp = max(abs(corr_estimate), na.rm = TRUE)
  clamp = min(ceiling(clamp * 10) / 10, 1)   # round up to nearest 0.1, cap at 1
  pal12 = dichromat::colorschemes$DarkRedtoBlue.12
  cf    = circlize::colorRamp2(seq(-clamp, clamp, length.out = 12), pal12)

  Heatmap(
    corr_estimate, rect_gp = gpar(type = "none", fill = "white"),
    col = cf, na_col = "grey", cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_dend = FALSE, show_column_dend = FALSE,
    show_column_names = TRUE, column_labels = col_labels,
    column_names_side = "bottom", column_names_rot = 90,
    column_names_gp = gpar(fontsize = textsize, col = textcolor), show_row_names = TRUE,
    row_labels = cap_first(list_of_gene_scores),
    row_names_side = "left", row_names_gp = gpar(fontsize = textsize, col = textcolor,
                               fontface = "plain"), row_names_max_width = unit(30, "line"),
    column_title = ptitle, column_title_side = "top",
    column_title_gp = gpar(fontsize = textsize, col = textcolor, fontface = "plain"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.rect(x, y, width, height, gp = gpar(col = NA, fill = "white", alpha = 0.2))
      r = corr_estimate[i, j]
      if (is.na(r) || is.na(sz[i, j])) return(invisible(NULL))
      grid.circle(x, y, r = circlesize * sz[i, j] * min(unit.c(width, height)),
                  gp = gpar(fill = cf(r), col = NA, alpha = 0.9))
      grid.circle(x, y, r = circlesize * sz2[i, j] * min(unit.c(width, height)),
                  gp = gpar(fill = NA, col = "grey17", alpha = 1, lwd = 0.6))
    }, show_heatmap_legend = TRUE, heatmap_legend_param = list(
      title = "Correlation", title_gp = gpar(fontsize = textsize, col = textcolor),
      labels_gp = gpar(fontsize = textsize, col = textcolor), legend_height = unit(1, "line"),
      grid_width = unit(0.1, "line"), direction = "vertical",
      title_position = "leftcenter-rot", at = seq(-0.3, 0.3, length = 3)))
}

### heatmaps, whole cohort and stratified

scores = c("immature", "mature", "degranulating", "antiprotease", "total", "IFN")
scores2 = paste0(scores_to_test[1:5], "_adj")

pht.all = Framingham_correlation_ht(framingham_score_pheno2, NULL,
                                    list_of_gene_scores = scores)
pht.all2 = Framingham_correlation_ht(framingham_score_pheno2, "Framingham cohort, 5321 subjects",
                                     list_of_gene_scores = scores2)

pht.females = Framingham_correlation_ht(framingham_score_pheno2[sex %in% "F"],
                                        "2861 females",
                                        list_of_gene_scores = scores)
pht.males = Framingham_correlation_ht(framingham_score_pheno2[sex %in% "M"],
                                      "2460 males", list_of_gene_scores = scores)
pht.nodiabetes = Framingham_correlation_ht(
  framingham_score_pheno2[diabetes == 0], "4846 subjects without diabetes",
  list_of_gene_scores = scores)


# -- Framingham disease --------------------------------------------------------

FH.diseaseES = fread(paste0(dir.results, "Framingham/FH.diseaseES.csv"))

FH.diseaseES = as.data.table(FH.diseaseES)[grepl("adj", score_name)][
  grepl("future|dia", condition)]
FH.diseaseES$fdr = p.adjust(FH.diseaseES$p, method = "BH")

FH.diseaseES$sig = ifelse(FH.diseaseES$fdr <= 0.05, "FDR ≤ 0.05",
                          ifelse(FH.diseaseES$fdr <= 0.1, "FDR ≤ 0.1", "ns"))
FH.diseaseES$sig = factor(FH.diseaseES$sig, levels = c("FDR ≤ 0.05", "FDR ≤ 0.1", "ns"))

shape_map  = c("diabetes" = 15, "future cancer" = 16, "future CVD" = 18)
FH.diseaseES$condition = factor(FH.diseaseES$condition, levels = rev(names(shape_map)))
FH.diseaseES$score_name = gsub("_adj", "", FH.diseaseES$score_name)

# Order by diabetes ES
gene_order = FH.diseaseES[condition == "diabetes"][order(effsize)]$score_name
FH.diseaseES$score_name    = factor(FH.diseaseES$score_name, levels = gene_order)
dodge      = position_dodge(width = 0.6)

FH.disease.forest =
  ggplot(FH.diseaseES, aes(x = effsize, y = score_name, color = sig,
                           shape = condition, group = condition)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey75", linewidth = 0.3,
             alpha = 0.5) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.5, size = 0.3, position = dodge) +
  geom_point(aes(fill = sig), size = 1.1, position = dodge, stroke = 0.3) +
  scale_shape_manual(values = shape_map, name = "condition", labels = cap_first) +
  scale_color_manual(values = c("FDR ≤ 0.05" = "#C0392B", "FDR ≤ 0.1" = "#E67E22",
                                "ns" = "grey50"), name = "Significance") +
  scale_fill_manual(values = c("FDR ≤ 0.05" = "#C0392B", "FDR ≤ 0.1" = "#E67E22",
                                "ns" = "grey50"), name = "Significance") +
  scale_y_discrete(labels = cap_first) +
  labs(x = "Effect size\ndisease vs. healthy", y = NULL) +
  theme_expresso(legend_text_size = textsize, axis_text_size = textsize,
                 plot_background = "white", panel_background = "white") +
  theme(legend.title = element_blank(), panel.grid.major.y = element_blank(),
        axis.title.x = element_text(hjust = 0.5, margin = margin(0, 0.1, 0.4, 0.1,
                                                                 unit = "line"))) +
  theme_legend_right_box +
  guides(shape = guide_legend(reverse = TRUE, ncol = 1, order = 1),
         color = guide_legend(ncol = 1, order = 2,
                              override.aes = list(shape = NA, linewidth = 0.6)),
         fill  = "none")


# -- GSE184050, T2D longitudinal -----------------------------------------------

expr.184050 = readRDS(paste0(dir.data, "GSE184050_exprlist.rds"))
expr.184050$logtpm = expr.184050$logtpm + 1

gene_sets2 = gene_sets
gene_sets2$IFN = c(genes.ifn, "DDX58", "IRF7", "SAMD9L")

# ── Compute gene set scores (gene-sets × samples) ─────────────────────────
scores_to_test_184050 = names(gene_sets2)
score_mat_184050 = do.call(rbind, lapply(scores_to_test_184050, function(gs) {
  genes = intersect(gene_sets2[[gs]], rownames(expr.184050$logtpm))
  apply(expr.184050$logtpm[genes, , drop = FALSE], 2, geom_mean)
}))
rownames(score_mat_184050) = scores_to_test_184050

# samples × scores + pheno
score_pheno_184050 = cbind(expr.184050$pheno, as.data.frame(t(score_mat_184050)))

# ── Compute deltas (follow-up − baseline) per patient ─────────────────────
bl  = score_pheno_184050[score_pheno_184050$timepoint == "baseline",  ]
fu  = score_pheno_184050[score_pheno_184050$timepoint == "follow-up", ]
pts = intersect(bl$patient_id, fu$patient_id)
bl  = bl[match(pts, bl$patient_id), ]
fu  = fu[match(pts, fu$patient_id), ]

delta_mat_184050 = as.matrix(fu[, scores_to_test_184050]) -
                   as.matrix(bl[, scores_to_test_184050])
delta_df_184050 = data.frame(patient_id = pts, disease_state = bl$disease_state,
                              as.data.frame(delta_mat_184050))

# ── Helper: paired effect size (Cohen's d_z) + Wilcoxon signed-rank ───────
.paired_stats_184050 = function(delta_sub, scores2test, condition) {
  do.call(rbind, lapply(scores2test, function(s) {
    d  = delta_sub[[s]]
    n  = length(d)
    dz = mean(d) / sd(d)
    se = sqrt((1 + dz^2 / 2) / n)
    data.frame(score_name = s, condition = condition, effsize = dz,
               effsize_se = se, conf.low = dz - 1.96 * se,
               conf.high = dz + 1.96 * se, p = wilcox.test(d, mu = 0)$p.value,
               n1 = n, n2 = NA_integer_)
  }))
}

# ── 1) Follow-up vs baseline in controls (paired) ─────────────────────────
res_ctrl_184050 = .paired_stats_184050(
  delta_df_184050[delta_df_184050$disease_state == "control", ],
  scores_to_test_184050, "control: follow-up vs. baseline")

# ── 2) Follow-up vs baseline in T2D (paired) ──────────────────────────────
res_t2d_184050 = .paired_stats_184050(
  delta_df_184050[delta_df_184050$disease_state == "T2D", ],
  scores_to_test_184050, "T2D: follow-up vs. baseline")

# ── 3) T2D delta vs control delta (unpaired, Hedges' g) ───────────────────
cls_184050 = as.integer(delta_df_184050$disease_state == "T2D")
hg_184050  = expresso::gene_hedges_g(
  t(as.matrix(delta_df_184050[, scores_to_test_184050])), cls_184050)
se_184050  = sqrt(hg_184050$var)
pvals_184050 = sapply(scores_to_test_184050, function(s)
  wilcox.test(delta_df_184050[[s]][cls_184050 == 1], delta_df_184050[[s]][cls_184050 == 0])$p.value)

res_delta_184050 = data.frame(
  score_name = scores_to_test_184050, condition = "T2D delta vs. control delta",
  effsize = hg_184050$es, effsize_se = se_184050,
  conf.low = hg_184050$es - 1.96 * se_184050, conf.high = hg_184050$es + 1.96 * se_184050,
  p = pvals_184050, n1 = sum(cls_184050 == 1), n2 = sum(cls_184050 == 0))

# ── Combine & FDR ─────────────────────────────────────────────────────────
gse184050_effsizes = rbind(res_ctrl_184050, res_t2d_184050, res_delta_184050)
gse184050_effsizes$fdr = p.adjust(gse184050_effsizes$p, method = "BH")

# ── Heatmap ────────────────────────────────────────────────────────────────
make_cellfun =
  function(val_mat, fdr_mat, cf, dot_scale = 0.35, fdr_cutoff = 0.1) {
    sz_mat  = 1.2 - fdr_mat
    fd2     = fdr_mat; fd2[fd2 > fdr_cutoff] = 1.2
    sz2_mat = 1.2 - fd2          # 0 when FDR > fdr_cutoff → outline invisible

    function(j, i, x, y, width, height, fill) {
      grid.rect(x, y, width, height, gp = gpar(col = NA, fill = "white", alpha = 0.2))
      v = val_mat[i, j]; f = fdr_mat[i, j]
      if (is.na(v) || is.na(f)) return(invisible(NULL))
      r_base = dot_scale * min(unit.c(width, height))
      grid.circle(x, y, r = sz_mat[i, j] * r_base, gp = gpar(fill = cf(v), col = NA, alpha = 0.9))
      grid.circle(x, y, r = sz2_mat[i, j] * r_base,
                  gp = gpar(fill = NA, col = "grey17", alpha = 1, lwd = 0.6))
    }
  }

col_order = c("control: follow-up vs. baseline", "T2D: follow-up vs. baseline",
              "T2D delta vs. control delta")

es_mat_184050 = reshape2::acast(gse184050_effsizes, score_name ~ condition,
                                 value.var = "effsize")[, col_order]
fdr_mat_184050 = reshape2::acast(gse184050_effsizes, score_name ~ condition,
                                 value.var = "fdr")[, col_order]

row_order_184050 = order(es_mat_184050[, "T2D delta vs. control delta"],
                         decreasing = TRUE, na.last = TRUE)
es_mat_184050    = es_mat_184050[row_order_184050, ]
fdr_mat_184050   = fdr_mat_184050[row_order_184050, ]

clamp_184050 = min(ceiling(max(abs(es_mat_184050), na.rm = TRUE) * 10) / 10, 2)
pal12        = dichromat::colorschemes$DarkRedtoBlue.12
cf_184050    = circlize::colorRamp2(
  seq(-clamp_184050, clamp_184050, length.out = 12), pal12)

ht_184050 = Heatmap(
  es_mat_184050, rect_gp = gpar(type = "none"), col = cf_184050, na_col = "grey",
  cluster_rows     = FALSE, cluster_columns  = FALSE,
  show_row_dend = FALSE, show_column_dend = FALSE,
  show_column_names = TRUE, column_names_side = "bottom", column_names_rot = 90,
  column_names_gp = gpar(fontsize = textsize, col = textcolor, fontface = "plain",
                         lineheight = 0.7),
  column_labels = c(
    "control: follow-up vs. baseline" = "Control\nfollow-up vs. baseline",
    "T2D: follow-up vs. baseline"     = "T2D\nfollow-up vs. baseline",
    "T2D delta vs. control delta"     = "T2D delta vs.\ncontrol delta")[
      colnames(es_mat_184050)],
  show_row_names = TRUE, row_names_side = "left",
  row_labels = cap_first(rownames(es_mat_184050)),
  row_names_gp = gpar(fontsize = textsize, col = textcolor, fontface = "plain"),
  row_names_max_width = unit(30, "line"),
  cell_fun          = make_cellfun(es_mat_184050, fdr_mat_184050, cf_184050,
                                   dot_scale = 0.45, fdr_cutoff = 0.1), show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = "Effect size", title_gp = gpar(fontsize = textsize, col = textcolor),
    labels_gp = gpar(fontsize = textsize, col = textcolor), legend_height = unit(1, "line"),
    grid_width = unit(0.1, "line"), direction = "vertical",
    title_position = "leftcenter-rot", at = seq(-0.6, 0.6, length.out = 3)))


# -- Keys for figure 6 and extended figure 8 -----------------------------------

### FDR outline legends, at 0.05 and at 0.1

df.twocircles = data.table(
  x = rep(1, 2), y = 1:2, color = c("grey88", "grey17"), labeltext = c("FDR > 0.05", "FDR ≤ 0.05"))

p.twocircles = ggplot(df.twocircles, aes(x = x, y = y, label = labeltext, color = color)) +
  geom_point(shape = 21, fill = "grey88", size = 2) +
  geom_text(angle = 0, size = pt2mm(textsize), hjust = -0.12, color = textcolor) +
  theme(text = element_text(size = textsize, color = textcolor)) +
  scale_color_identity() +
  theme_void() +
  coord_cartesian(clip = "off")

p.twocircles.t = ggplot(df.twocircles, aes(x = y, y = x, label = labeltext, color = color)) +
  geom_point(shape = 21, fill = "grey88", size = 2, stroke = 0.3) +
  geom_text(angle = 90, size = pt2mm(textsize), hjust = -0.12, color = textcolor) +
  theme(text = element_text(size = textsize, color = textcolor)) +
  scale_color_identity() +
  scale_x_reverse() +
  theme_void() +
  coord_cartesian(clip = "off")

df.twocircles0.1 = data.table(
  x = rep(1, 2), y = 1:2, color = c("grey88", "grey17"), labeltext = c("FDR > 0.1", "FDR ≤ 0.1"))

p.twocircles0.1 = ggplot(df.twocircles0.1, aes(x = x, y = y, label = labeltext, color = color)) +
  geom_point(shape = 21, fill = "grey88", size = 2) +
  geom_text(angle = 0, size = pt2mm(textsize), hjust = -0.135, color = textcolor) +
  theme(text = element_text(size = textsize, color = textcolor)) +
  scale_color_identity() +
  theme_void() +
  coord_cartesian(clip = "off")

p.twocircles0.1.t =
  ggplot(df.twocircles0.1, aes(x = y, y = x, label = labeltext, color = color)) +
  geom_point(shape = 21, fill = "grey88", size = 2, stroke = 0.3) +
  geom_text(angle = 90, size = pt2mm(textsize), hjust = -0.135, color = textcolor) +
  theme(text = element_text(size = textsize, color = textcolor)) +
  scale_color_identity() +
  scale_x_reverse() +
  theme_void() +
  coord_cartesian(clip = "off")

### KM curve colour legend

df.legend_survival = data.table(
  y = c(2, 1), color = c("#214DC8", "#FF7F0E"),
  label = c("< median", "≥ median"))

p.legend_survival = ggplot(data = df.legend_survival) +
  ggplot2::geom_segment(aes(x = 0, xend = 0.3, y = y, yend = y, color = color),
                        linewidth = 0.5) +
  geom_text(aes(x = 0.42, y = y, label = label), hjust = 0,
            size = pt2mm(textsize), color = textcolor) +
  scale_color_identity() +
  xlim(0, 2.6) +
  ylim(0.4, 2.6) +
  theme_void()


# -- Figure 6 ------------------------------------------------------------------

cairo_pdf(file = paste0(dir.fig, "Figure6.pdf"), width = 5.28, height = 5)
pushViewport(viewport(layout = grid.layout(nrow = 93, ncol = 67)))

### lab correlations
pushViewport(viewport(layout.pos.row = 2:25, layout.pos.col = 4:23))
draw(pht.all, merge_legend = FALSE, heatmap_legend_side = "right",
     align_heatmap_legend = "heatmap_top",
     padding = unit(c(0.14, 0.5, 0.5, 0.2), "line"), annotation_legend_side = "top",
     newpage = FALSE, background = "white")
upViewport(1)


### GSE184050
pushViewport(viewport(layout.pos.row = 2:31, layout.pos.col = 38:58))
draw(ht_184050, merge_legend = FALSE, heatmap_legend_side = "right",
     align_heatmap_legend = "heatmap_top",
     padding = unit(c(0, 0.5, 0.5, 0.2), "line"), annotation_legend_side = "top",
     newpage = FALSE, background = "white")
upViewport(1)

### forests
print(FH.disease.forest + theme(plot.margin = margin(0.04, 0.04, 0.04, 0.04, "line")),
      vp = viewport(layout.pos.row = 34:66, layout.pos.col = 2:20))
print(FH.survival.forest + theme(plot.margin = margin(0.04, 0.04, 0.08, 0.04, "line")),
      vp = viewport(layout.pos.row = 34:64, layout.pos.col = 33:51))


print(p.twocircles + theme(plot.margin = margin(0.04, 0.04, 0.04, 0.04, "line")),
      vp = viewport(layout.pos.row = 12:13, layout.pos.col = 19:20))
print(p.twocircles0.1 + theme(plot.margin = margin(0.04, 0.04, 0.04, 0.04, "line")),
      vp = viewport(layout.pos.row = 12:13, layout.pos.col = 54:55))

### KM curves
print(plot_km_residual(surv.all$survival_summary, n_year = 5,
                       score = "degranulating", col_survival = col_survival) +
        theme(plot.margin = margin(0.04, 0.04, 0.04, 0.04, "line")),
      vp = viewport(layout.pos.row = 68:92, layout.pos.col = 2:17))

print(plot_km_residual(surv.all$survival_summary, n_year = 5,
                       score = "antiprotease", col_survival = col_survival) +
        theme(plot.margin = margin(0.04, 0.04, 0.04, 0.04, "line")),
      vp = viewport(layout.pos.row = 68:92, layout.pos.col = 18:33))

print(plot_km_residual(surv.all$survival_summary, n_year = 5,
                       score = "total", col_survival = col_survival) +
        theme(plot.margin = margin(0.04, 0.04, 0.04, 0.04, "line")),
      vp = viewport(layout.pos.row = 68:92, layout.pos.col = 34:49))

print(plot_km_residual(surv.all$survival_summary, n_year = 5,
                       score = "immature", col_survival = col_survival) +
        theme(plot.margin = margin(0.04, 0.04, 0.04, 0.04, "line")),
      vp = viewport(layout.pos.row = 68:92, layout.pos.col = 50:65))

print(p.legend_survival, vp = viewport(layout.pos.row = 68:70, layout.pos.col = 14:19))
print(p.legend_survival, vp = viewport(layout.pos.row = 68:70, layout.pos.col = 30:35))
print(p.legend_survival, vp = viewport(layout.pos.row = 68:70, layout.pos.col = 46:51))
print(p.legend_survival, vp = viewport(layout.pos.row = 68:70, layout.pos.col = 62:67))

### panel labels
grid.text(x = unit(0.02, "npc"), y = unit(0.98, "npc"), label = "a",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", fontfamily = "Arial",
                    col = textcolor))
grid.text(x = unit(0.02, "npc"), y = unit(0.66, "npc"), label = "b",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", fontfamily = "Arial",
                    col = textcolor))
grid.text(x = unit(0.54, "npc"), y = unit(0.98, "npc"), label = "c",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", fontfamily = "Arial",
                    col = textcolor))
grid.text(x = unit(0.48, "npc"), y = unit(0.66, "npc"), label = "d",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", fontfamily = "Arial",
                    col = textcolor))
grid.text(x = unit(0.02, "npc"), y = unit(0.28, "npc"), label = "e",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", fontfamily = "Arial",
                    col = textcolor))

dev.off()


# -- Extended figure 8 ---------------------------------------------------------

cairo_pdf(file = paste0(dir.fig, "Fig.e8.pdf"), width = 5.2, height = 4.6)
pushViewport(viewport(layout = grid.layout(nrow = 83, ncol = 84)))

### stratified lab correlations
pushViewport(viewport(layout.pos.row = 1:36, layout.pos.col = 2:28))
ht_opt(TITLE_PADDING = unit(c(0.5, 0.5), "line"))
draw(pht.females, merge_legend = FALSE, heatmap_legend_side = "right",
     align_heatmap_legend = "heatmap_top",
     annotation_legend_side = "top", newpage = FALSE)
upViewport(1)
ht_opt(TITLE_PADDING = unit(c(0.2, 0.2), "line"))

pushViewport(viewport(layout.pos.row = 1:36, layout.pos.col = 29:55))
ht_opt(TITLE_PADDING = unit(c(0.5, 0.5), "line"))
draw(pht.males, merge_legend = FALSE, heatmap_legend_side = "right",
     align_heatmap_legend = "heatmap_top",
     annotation_legend_side = "top", newpage = FALSE)
upViewport(1)
ht_opt(TITLE_PADDING = unit(c(0.2, 0.2), "line"))

pushViewport(viewport(layout.pos.row = 1:36, layout.pos.col = 56:82))
ht_opt(TITLE_PADDING = unit(c(0.5, 0.5), "line"))
draw(pht.nodiabetes, merge_legend = FALSE, heatmap_legend_side = "right",
     align_heatmap_legend = "heatmap_top",
     annotation_legend_side = "top", newpage = FALSE)
upViewport(1)
ht_opt(TITLE_PADDING = unit(c(0.2, 0.2), "line"))

print(p.twocircles + theme(plot.margin = margin(0.04, 0.04, 0.04, 0.04, "line")),
      vp = viewport(layout.pos.row = 16:17, layout.pos.col = 21:23))
print(p.twocircles + theme(plot.margin = margin(0.04, 0.04, 0.04, 0.04, "line")),
      vp = viewport(layout.pos.row = 16:17, layout.pos.col = 48:50))
print(p.twocircles + theme(plot.margin = margin(0.04, 0.04, 0.04, 0.04, "line")),
      vp = viewport(layout.pos.row = 16:17, layout.pos.col = 75:77))

### Cox proportional hazards model tables, 5-year mortality
pushViewport(viewport(layout.pos.row = 40:60, layout.pos.col = 3:28))
grid.draw(make_cox_table(surv.all, "degranulating_z", n_year = 5, row_label = "degranulating"))
popViewport()
pushViewport(viewport(layout.pos.row = 62:82, layout.pos.col = 3:28))
grid.draw(make_cox_table(surv.all, "antiprotease_z", n_year = 5, row_label = "antiprotease"))
popViewport()

pushViewport(viewport(layout.pos.row = 40:60, layout.pos.col = 30:56))
grid.draw(make_cox_table(surv.all, "total_z", n_year = 5, row_label = "total"))
popViewport()

pushViewport(viewport(layout.pos.row = 62:82, layout.pos.col = 30:55))
grid.draw(make_cox_table(surv.all, "immature_z", n_year = 5, row_label = "immature"))
popViewport()

pushViewport(viewport(layout.pos.row = 40:60, layout.pos.col = 58:82))
grid.draw(make_cox_table(surv.all, "mature_z", n_year = 5, row_label = "mature"))
popViewport()

pushViewport(viewport(layout.pos.row = 62:82, layout.pos.col = 57:82))
grid.draw(make_cox_table(surv.all, "IFN_z", n_year = 5, row_label = "IFN"))
popViewport()

### panel labels
grid.text(x = unit(0.03, "npc"), y = unit(0.98, "npc"), label = "a",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", fontfamily = "Arial",
                    col = textcolor))
grid.text(x = unit(0.35, "npc"), y = unit(0.98, "npc"), label = "b",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", fontfamily = "Arial",
                    col = textcolor))
grid.text(x = unit(0.67, "npc"), y = unit(0.98, "npc"), label = "c",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", fontfamily = "Arial",
                    col = textcolor))

grid.text(x = unit(0.02, "npc"), y = unit(0.52, "npc"), label = "d",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", fontfamily = "Arial",
                    col = textcolor))


dev.off()
