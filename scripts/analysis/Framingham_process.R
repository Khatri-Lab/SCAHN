# ==============================================================================
# Framingham signature scores, covariate adjustment and disease effect sizes
#
# Scores the six neutrophil signatures on the co-normalised Framingham
# expression matrix, regresses out age / sex / BMI / smoking to give the `_adj`
# variants, then computes Hedges' g and a Wilcoxon p value for prevalent and
# incident cancer, CVD and diabetes against the healthy samples.
#
# Two blocks, run in order: the scoring block writes score_pheno.csv, the
# effect-size block reads it back and writes diseaseES.csv. The input is the
# 52-gene subset.
#
# Requires 00_setup.R (data.table, dir.results), 01_load_data.R (genes.*
# signature vectors) and expresso, loaded explicitly below.
#
# Reads:
#   results/Framingham/FH.expr.selectedgenes.rds
#
# Writes:
#   results/Framingham/FH.score_pheno.csv
#     read by 05_sexdifference.R and 06_framingham_survival.R
#   results/Framingham/FH.diseaseES.csv
#     read by 06_framingham_survival.R
# ==============================================================================

SCAHN_SCRIPTS = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/scripts/figures"
source(file.path(SCAHN_SCRIPTS, "00_setup.R"))      # packages, paths
source(file.path(SCAHN_SCRIPTS, "01_load_data.R"))  # genes.* signature vectors

library(expresso)  # get_gene_scores, gene_hedges_g


# -- Covariate adjustment ------------------------------------------------------

# Residuals of the score on age, sex, BMI and smoking status. lm() drops rows
# with any missing covariate, so the residuals are scattered back by name onto a
# full-length NA vector.
compute_adjusted_score = function(score_pheno_matrix, score_name) {
  lm_score_adjust = lm(as.formula(paste0(score_name,
    " ~ age + as.factor(sex) + bmi + as.factor(regular_smokers)")), data = score_pheno_matrix)
  score_pheno_matrix$residuals = NA
  score_pheno_matrix$residuals[as.numeric(names(lm_score_adjust$residuals))] =
    lm_score_adjust$residuals
  return(score_pheno_matrix$residuals)
}


# -- Effect sizes for prevalent and incident disease ---------------------------

# Internal: compute Hedges' g + Wilcoxon for one comparison
.eff_size_one = function(mat, list_of_classes, condition, scores2test) {
  cls       = as.integer(mat$class == list_of_classes[1])
  score_mat = t(as.matrix(mat[, scores2test, drop = FALSE]))

  hg  = gene_hedges_g(score_mat, cls)
  se  = sqrt(hg$var)
  n1  = sum(cls == 1)
  n2  = sum(cls == 0)

  pvals = sapply(scores2test, function(s)
    wilcox.test(mat[[s]][cls == 1], mat[[s]][cls == 0])$p.value)

  data.frame(score_name = scores2test, condition = condition, effsize = hg$es, effsize_se = se,
    conf.low = hg$es - 1.96 * se, conf.high = hg$es + 1.96 * se, p = pvals, n1 = n1, n2 = n2,
    row.names = NULL)
}

# Internal: split a disease group into future/past comparisons vs healthy
.future_past_eff = function(pheno, disease_col, date_col, label_future, label_past, scores2test,
                            window_days = 5 * 365) {
  sub       = pheno[pheno$healthy == 1 | pheno[[disease_col]] == 1, ]
  sub$class = "healthy"
  sub$class[sub[[disease_col]] == 1 & sub[[date_col]] >  0 & sub[[date_col]] <= window_days]  =
    label_future
  sub$class[sub[[disease_col]] == 1 & sub[[date_col]] <  0 & sub[[date_col]] >= -window_days] =
    label_past

  # Drop disease cases outside the window (they stay unlabeled, remove them)
  sub = sub[!(sub[[disease_col]] == 1 & sub$class == "healthy"), ]

  rbind(
    .eff_size_one(sub[sub$class %in% c("healthy", label_future), ],
      c(label_future, "healthy"), label_future, scores2test),
    .eff_size_one(sub[sub$class %in% c("healthy", label_past), ],
      c(label_past, "healthy"), label_past, scores2test)
  )
}

perform_case_control_framingham = function(framingham_score_pheno_matrix, scores2test,
                                           window_years = 5) {
  pheno       = framingham_score_pheno_matrix
  window_days = window_years * 365

  diab       = pheno[pheno$healthy == 1 | pheno$diabetes == 1, ]
  diab$class = ifelse(diab$diabetes == 1, "diabetes", "healthy")

  rbind(
    .future_past_eff(pheno, "cancer", "first_cancer_diagnosis_date", "future cancer",
      "past cancer", scores2test, window_days),
    .future_past_eff(pheno, "cardio", "first_cardio_diagnosis_date", "future CVD",
      "past CVD", scores2test, window_days),
    .eff_size_one(diab, c("diabetes", "healthy"), "diabetes", scores2test)
  )
}


# -- Signature scores and covariate adjustment ---------------------------------

fh.expr = readRDS(paste0(dir.results, "Framingham/FH.expr.selectedgenes.rds"))

framingham_score_pheno = fh.expr$pheno
table(framingham_score_pheno$sex)

dim(fh.expr$expr)   # 52  5321

framingham_score_pheno = as.data.table(framingham_score_pheno)
framingham_score_pheno$immature = get_gene_scores(fh.expr$expr, genes.immature.neu, "")
framingham_score_pheno$mature = get_gene_scores(fh.expr$expr, genes.mature.neu, "")
framingham_score_pheno$degranulating = get_gene_scores(fh.expr$expr, genes.degranulating.neu, "")
framingham_score_pheno$antiprotease = get_gene_scores(fh.expr$expr, genes.antiprotease.neu, "")
framingham_score_pheno$total = get_gene_scores(fh.expr$expr, genes.total.neu, "")
framingham_score_pheno$IFN = get_gene_scores(fh.expr$expr, genes.ifn, "")

rownames(framingham_score_pheno) = NULL

framingham_score_pheno$immature_adj = compute_adjusted_score(framingham_score_pheno, "immature")
framingham_score_pheno$mature_adj = compute_adjusted_score(framingham_score_pheno, "mature")
framingham_score_pheno$degranulating_adj = compute_adjusted_score(framingham_score_pheno,
  "degranulating")
framingham_score_pheno$antiprotease_adj = compute_adjusted_score(framingham_score_pheno,
  "antiprotease")
framingham_score_pheno$total_adj = compute_adjusted_score(framingham_score_pheno, "total")
framingham_score_pheno$IFN_adj = compute_adjusted_score(framingham_score_pheno, "IFN")

fwrite(framingham_score_pheno, paste0(dir.results,
  "Framingham/FH.score_pheno.csv"))


# -- Disease effect sizes ------------------------------------------------------

framingham_score_pheno = fread(paste0(dir.results,
  "Framingham/FH.score_pheno.csv"))
scores_to_test = c("immature", "mature", "degranulating", "antiprotease", "total", "IFN")
framingham_score_pheno = as.data.frame(framingham_score_pheno)
eff_size_p_val_summary = perform_case_control_framingham(framingham_score_pheno,
  c(scores_to_test, paste0(scores_to_test, "_adj")), window_years = 7)

fwrite(eff_size_p_val_summary, paste0(dir.results,
  "Framingham/FH.diseaseES.csv"))
