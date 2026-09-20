# ==============================================================================
# Sex differences in whole-blood and PBMC bulk transcriptomes
#
# Male vs. female across the curated bulk sex cohort, split by tissue into whole
# blood and PBMC and meta-analysed within each. Two independent arms over the
# same datasets: a MetaIntegrator effect-size meta-analysis, and a per-gene
# logistic meta-analysis with age as a covariate, which is the one
# 05_sexdifference.R reads.
#
# In the logistic arm each dataset is fitted on its own, then pooled per tissue
# under a random-effects model. Expression is z-scored per gene within each
# dataset before fitting, so log-OR reads as "per 1 SD" and is comparable across
# genes and datasets. A dataset is skipped when it has no usable age or fewer
# than six samples with both sex and age, and the adjusted and unadjusted fits
# are given the identical sample set so their estimates are comparable.
#
# Both arms restrict to genes present in both tissues, which is what makes the
# WB-vs-PBMC columns a direct side-by-side.
#
# Requires 00_setup.R (data.table, expresso, matrixStats, dir.data /
# dir.results) and MetaIntegrator, called namespace-qualified below.
#
# Reads:
#   data/sex/datasets.sexMeta.rds
#   data/sex/allGeneLocations_dict.rds
#
# Writes:
#   results/sex_difference/dtf.sexMeta.rds
#     read by no script; the unadjusted MetaIntegrator record
#   results/sex_difference/dtf.sexMeta.ageAdj.shared.rds
#     read by 05_sexdifference.R
# ==============================================================================

SCAHN_SCRIPTS = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/scripts/figures"
source(file.path(SCAHN_SCRIPTS, "00_setup.R"))  # packages, paths


# -- Datasets by tissue --------------------------------------------------------

sexdatasets = readRDS(paste0(dir.data, "sex/datasets.sexMeta.rds"))

datasettissuetype = sapply(sexdatasets, function(x) { unique(x$pheno$tissue) })
sexdatasets.WB = sexdatasets[names(which(datasettissuetype == "Whole blood"))]
sexdatasets.PBMC = sexdatasets[names(which(datasettissuetype == "PBMC"))]

# samples of each sex per tissue
lapply(sexdatasets.WB, function(x) { table(x$pheno$sex) })
sum(sapply(sexdatasets.WB, function(x) { sum(x$pheno$sex == "Male") }))
sum(sapply(sexdatasets.WB, function(x) { sum(x$pheno$sex == "Female") }))
lapply(sexdatasets.PBMC, function(x) { table(x$pheno$sex) })
sum(sapply(sexdatasets.PBMC, function(x) { sum(x$pheno$sex == "Male") }))
sum(sapply(sexdatasets.PBMC, function(x) { sum(x$pheno$sex == "Female") }))


# -- MetaIntegrator effect-size meta-analysis ----------------------------------

sexMeta.WB = MetaIntegrator::runMetaAnalysis(list(originalData = sexdatasets.WB),
                                             runLeaveOneOutAnalysis = FALSE)
sexMeta.PBMC = MetaIntegrator::runMetaAnalysis(list(originalData = sexdatasets.PBMC),
                                               runLeaveOneOutAnalysis = FALSE)

dtf.sexMeta.WB = sexMeta.WB$metaAnalysis$pooledResults
dtf.sexMeta.PBMC = sexMeta.PBMC$metaAnalysis$pooledResults

dtf.sexMeta.WB$gene = rownames(dtf.sexMeta.WB)
dtf.sexMeta.PBMC$gene = rownames(dtf.sexMeta.PBMC)

dtf.sexMeta.WB = as.data.table(dtf.sexMeta.WB)
dtf.sexMeta.PBMC = as.data.table(dtf.sexMeta.PBMC)

sharedgenes = intersect(dtf.sexMeta.WB$gene, dtf.sexMeta.PBMC$gene)

dtf.sexMeta.WB = dtf.sexMeta.WB[gene %in% sharedgenes][order(gene)]
dtf.sexMeta.PBMC = dtf.sexMeta.PBMC[gene %in% sharedgenes][order(gene)]
all.equal(dtf.sexMeta.WB$gene, dtf.sexMeta.PBMC$gene)

dtf.sexMeta.WB.tmp = dtf.sexMeta.WB[, c(
  "effectSize", "effectSizeStandardError", "effectSizePval", "effectSizeFDR", "numStudies")]
dtf.sexMeta.PBMC.tmp = dtf.sexMeta.PBMC[, c(
  "effectSize", "effectSizeStandardError", "effectSizePval", "effectSizeFDR", "numStudies")]

colnames(dtf.sexMeta.WB.tmp) = paste0(colnames(dtf.sexMeta.WB.tmp), ".WB")
colnames(dtf.sexMeta.PBMC.tmp) = paste0(colnames(dtf.sexMeta.PBMC.tmp), ".PBMC")

dtf.sexMeta = data.table(gene = dtf.sexMeta.WB$gene)

allGeneLocations_dict = readRDS(paste0(dir.data, "sex/allGeneLocations_dict.rds"))
dtf.sexMeta$chr = allGeneLocations_dict[dtf.sexMeta$gene]
dtf.sexMeta = cbind(dtf.sexMeta, dtf.sexMeta.WB.tmp, dtf.sexMeta.PBMC.tmp)
saveRDS(dtf.sexMeta, paste0(dir.results, "sex_difference/dtf.sexMeta.rds"))


# -- Age-adjusted logistic meta-analysis ---------------------------------------

# Per-dataset logistic fit, restricted to samples with both sex and age available
# (identical sample set for the adjusted and unadjusted runs). Returns NULL for
# datasets without usable age (e.g. GSE18323).
run_one = function(x, covariates) {
  g = x$genes
  g = g[rownames(g) != "", , drop = FALSE]
  ph = as.data.frame(x$pheno)
  ph$sex01 = ifelse(ph$sex == "Male", 1L,
                    ifelse(ph$sex == "Female", 0L, NA_integer_))
  ph$age = suppressWarnings(as.numeric(ph$age))
  keep = !is.na(ph$sex01) & !is.na(ph$age)
  ph = ph[keep, , drop = FALSE]
  g  = g[, rownames(ph), drop = FALSE]
  if (nrow(ph) < 6L || length(unique(ph$sex01)) < 2L) return(NULL)
  # z-score each gene (row) within this dataset -> log-OR is per 1 SD (see header)
  g = (g - rowMeans(g, na.rm = TRUE)) / rowSds(g, na.rm = TRUE)
  g[!is.finite(g)] = NA                                     # constant/degenerate genes -> NA
  gene_logistic(g, ph, event_col = "sex01", covariates = covariates, test = "wald")
}

adj_list = list()
unadj_list = list()
nsamp = data.table()
for (nm in names(sexdatasets)) {
  a = tryCatch(run_one(sexdatasets[[nm]], covariates = "age"), error = function(e) NULL)
  if (is.null(a)) { message("[skip] ", nm, " (no age / too few)"); next }
  u = run_one(sexdatasets[[nm]], covariates = NULL)
  adj_list[[nm]] = a
  unadj_list[[nm]] = u
  nsamp = rbind(nsamp, data.table(dataset = nm, tissue = datasettissuetype[[nm]],
                                  n_used = a$summary$n[1], n_genes = nrow(a$summary)))
}
nsamp

# Pool per tissue
by_tissue = function(lst, tis) {
  lst[names(lst)[sapply(names(lst), function(n) datasettissuetype[[n]] == tis)]]
}

meta_adj_WB   = meta_logistic(by_tissue(adj_list,   "Whole blood"), method = "random")
meta_adj_PBMC = meta_logistic(by_tissue(adj_list,   "PBMC"),        method = "random")
meta_un_WB    = meta_logistic(by_tissue(unadj_list, "Whole blood"), method = "random")
meta_un_PBMC  = meta_logistic(by_tissue(unadj_list, "PBMC"),        method = "random")

# Tidy merged table, analogous to dtf.sexMeta
mk = function(m, tag) {
  s = as.data.table(m$summary)
  s = s[, .(gene = gene_set, pooled_or, pooled_or_lo, pooled_or_hi, pooled_p, FDR, I2, k)]
  setnames(s, setdiff(names(s), "gene"), paste0(setdiff(names(s), "gene"), tag))
  s
}
dtf.sexMeta.ageAdj = Reduce(function(a, b) { merge(a, b, by = "gene", all = TRUE) },
                            list(mk(meta_adj_WB, ".adj.WB"), mk(meta_adj_PBMC, ".adj.PBMC"),
                                 mk(meta_un_WB, ".unadj.WB"), mk(meta_un_PBMC, ".unadj.PBMC")))
dtf.sexMeta.ageAdj[, chr := allGeneLocations_dict[gene]]
setcolorder(dtf.sexMeta.ageAdj, c("gene", "chr"))

# Restrict to genes fitted in both tissues under the adjusted model, mirroring the
# shared-gene intersect used for dtf.sexMeta above.
dtf.sexMeta.ageAdj.shared = dtf.sexMeta.ageAdj[
  is.finite(pooled_or.adj.WB) & is.finite(pooled_or.adj.PBMC)][order(gene)]
saveRDS(dtf.sexMeta.ageAdj.shared,
        paste0(dir.results, "sex_difference/dtf.sexMeta.ageAdj.shared.rds"))
