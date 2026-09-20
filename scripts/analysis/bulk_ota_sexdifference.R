# ==============================================================================
# Sex differences in neutrophil gene expression, ota2021 bulk cohort
#
# Healthy controls only, male vs. female, on the Neu_LDG arm of the ota2021
# bulk expression list. Two independent models over the same cohort: a DESeq2
# Wald test on raw counts with design = ~ sex, and a per-gene logistic
# regression of sex on expression with age as a covariate, which is the one
# 05_sexdifference.R reads.
#
# In the logistic block, expression is VST-transformed then z-scored per gene,
# so log_or reads as "per 1 SD" and is comparable across genes -- the same
# convention as the age-adjusted meta-analysis. log_or > 0 is higher in males,
# matching the sign of the DESeq2 (M vs. F) log2FoldChange. The unadjusted fit
# is merged in alongside so the effect of the age term can be read off directly.
#
# Requires 00_setup.R (data.table, expresso, matrixStats, dir.data /
# dir.results) and DESeq2, loaded explicitly below.
#
# Reads:
#   data/sex/ota2021.clinical_diagnosis_age_sex_v2.txt
#   data/sex/ota2021.exprlist.rds
#
# Writes:
#   results/sex_difference/ota2021.Neu_LDG.sex.DEG.csv
#     read by no script; the unadjusted DESeq2 record
#   results/sex_difference/ota2021.Neu_LDG.sex.ageAdj.logit.DEG.csv
#     read by 05_sexdifference.R
# ==============================================================================

SCAHN_SCRIPTS = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/scripts/figures"
source(file.path(SCAHN_SCRIPTS, "00_setup.R"))  # packages, paths

library(DESeq2)

meta = fread(paste0(dir.data, "sex/ota2021.clinical_diagnosis_age_sex_v2.txt"))
exprobj = readRDS(paste0(dir.data, "sex/ota2021.exprlist.rds"))


# -- DESeq2, male vs. female ---------------------------------------------------

meta.hc = meta[disease == "HC"]
dim(exprobj$Neu_LDG$count)
countmat = round(exprobj$Neu_LDG$count, 0)

# what the gene-name filter drops
rownames(countmat)[grepl("^AC0", rownames(countmat))]
rownames(countmat)[grepl("^AL[0-9]", rownames(countmat))]
countmat = countmat[rownames(countmat)[
  !grepl("^AP[0-9]|^AF[0-9]|^AC[0-9]|^AL[0-9]|-AS1$", rownames(countmat))], ]

meta.hc = meta.hc[id %in% colnames(exprobj$Neu_LDG$count)]
meta.hc = as.data.frame(meta.hc)
rownames(meta.hc) = meta.hc$id
meta.hc$sex = factor(meta.hc$sex, levels = c("F", "M"))
dds = DESeqDataSetFromMatrix(countData = countmat[, meta.hc$id], colData = meta.hc, design = ~ sex)
dds = DESeq(dds)
res = results(dds, contrast = c("sex", "M", "F"))

res = as.data.frame(res)
res$FDR = p.adjust(res$pvalue, method = "fdr")
res$gene = rownames(res)
res = as.data.table(res)

fwrite(res, paste0(dir.results, "sex_difference/ota2021.Neu_LDG.sex.DEG.csv"))


# -- Logistic regression of sex on expression, age-adjusted --------------------

# Model per gene: sex(M=1) ~ gene + age (adjusted) vs. sex(M=1) ~ gene.

# healthy controls with sex + age available
meta.hc = meta[disease == "HC" & sex %in% c("M", "F") & !is.na(age)]
meta.hc = as.data.frame(meta.hc)
rownames(meta.hc) = meta.hc$id
meta.hc$sex01 = ifelse(meta.hc$sex == "M", 1L, 0L)
meta.hc$age   = as.numeric(meta.hc$age)

# VST expression, same gene-name filter as the DESeq2 block
exprmat = exprobj$Neu_LDG$expr
exprmat = exprmat[!grepl("^AP[0-9]|^AF[0-9]|^AC[0-9]|^AL[0-9]|-AS1$", rownames(exprmat)), ]
ids = intersect(colnames(exprmat), rownames(meta.hc))
exprmat = exprmat[, ids]
meta.hc = meta.hc[ids, , drop = FALSE]

# z-score each gene within the dataset -> log_or per 1 SD
exprmat = (exprmat - rowMeans(exprmat, na.rm = TRUE)) / rowSds(exprmat, na.rm = TRUE)
exprmat[!is.finite(exprmat)] = NA

fit.adj   = gene_logistic(exprmat, meta.hc, event_col = "sex01", covariates = "age", test = "wald")
fit.unadj = gene_logistic(exprmat, meta.hc, event_col = "sex01", covariates = NULL,  test = "wald")

res = as.data.table(fit.adj$summary)   # gene, or, log_or, se, p_wald, p_lrt, p_value, n, ...
res[, FDR := p.adjust(p_value, method = "fdr")]
res = merge(
  res, as.data.table(fit.unadj$summary)[, .(gene, log_or.unadj = log_or, p.unadj = p_value)],
  by = "gene", all.x = TRUE)
setorder(res, p_value)

fwrite(res, paste0(dir.results, "sex_difference/ota2021.Neu_LDG.sex.ageAdj.logit.DEG.csv"))
