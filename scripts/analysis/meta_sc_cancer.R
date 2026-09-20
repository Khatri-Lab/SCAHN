# ==============================================================================
# Meta-analysis of tumour vs. adjacent-tissue neutrophil gene expression
#
# Tumour-infiltrating neutrophils compared with neutrophils from adjacent normal
# tissue of the same studies. The per-sample geometric-mean expression of every
# gene is meta-analysed across the three studies that sampled both sides, with
# wu2024 split by cancer type so each contributes its own within-study contrast.
#
# The cohort selection below is the same one used for the proportion analysis in
# 04_infection_cancer.R; edit both if the sample set changes.
#
# Requires 00_setup.R (data.table, expresso, dir.data / dir.results) and
# 01_load_data.R (samplemeta).
#
# Reads:
#   data/samplemeta.cancer.csv
#   results/annotation/SCAHN.subsetproportion.rds
#   results/meta_analysis/dtfgenes.allneu.persample.geommean.rds
#
# Writes:
#   results/meta_analysis/SCAHN.cancer.exprmeta.allgenes.allneu.rds
#     read by 04_infection_cancer.R
# ==============================================================================

SCAHN_SCRIPTS = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/scripts/figures"
source(file.path(SCAHN_SCRIPTS, "00_setup.R"))      # packages, paths
source(file.path(SCAHN_SCRIPTS, "01_load_data.R"))  # samplemeta

# -- Tumour and adjacent-tissue samples ----------------------------------------

dtf.prop = readRDS(paste0(dir.results, "annotation/SCAHN.subsetproportion.rds"))

dtf.prop.cancertissue = dtf.prop[!grepl("peripheral blood", tissue)][
  !grepl("han2020", pid)]
dtf.prop.cancertissue = dtf.prop.cancertissue[grepl("cancer", group)]

samplemeta.cancer = fread(paste0(dir.data, "samplemeta.cancer.csv"))

dtf.prop.cancertissue$stage = samplemeta.cancer[
  match(dtf.prop.cancertissue$sampleid.pid, samplemeta.cancer$sampleid.pid)]$stage
dtf.prop.cancertissue$type = samplemeta.cancer[
  match(dtf.prop.cancertissue$sampleid.pid, samplemeta.cancer$sampleid.pid)]$type
dtf.prop.cancertissue$treatment = samplemeta.cancer[
  match(dtf.prop.cancertissue$sampleid.pid, samplemeta.cancer$sampleid.pid)]$treatment

dtf.prop.cancertissue = dtf.prop.cancertissue[!type %in% "metastasis"][
  n_cells >= 50]
table(dtf.prop.cancertissue$pid, dtf.prop.cancertissue$type)

### cohort size
dtf.prop.cancertissue[pid %in% c("hu2022", "salcher2022", "wu2024")][, .N]

# Only these three studies profiled both tumour and adjacent tissue; all of
# their samples happen to be pre-treatment, so the filter drops nothing.
dtf.prop.cancertissue = dtf.prop.cancertissue[
  pid %in% c("hu2022", "salcher2022", "wu2024")][treatment == "pre-treatment"]


# -- Per-sample expression, tumour vs. adjacent --------------------------------

dtfgenes.mean = readRDS(paste0(
  dir.results, "meta_analysis/dtfgenes.allneu.persample.geommean.rds"))

dtfgenes.mean$group = samplemeta[match(dtfgenes.mean$sampleid.pid, samplemeta$sampleid.pid)]$group
dtfgenes.mean$group.detail = samplemeta[
  match(dtfgenes.mean$sampleid.pid, samplemeta$sampleid.pid)]$group.detail
dtfgenes.mean$tissue = samplemeta[match(dtfgenes.mean$sampleid.pid, samplemeta$sampleid.pid)]$tissue
dtfgenes.mean$source = samplemeta[match(dtfgenes.mean$sampleid.pid, samplemeta$sampleid.pid)]$source

dtfgenes.mean$pid = gsub("[ab]$", "", dtfgenes.mean$pid)

dtf.score.cancer = dtfgenes.mean[
  sampleid.pid %in% dtf.prop.cancertissue$sampleid.pid]
dtf.score.cancer$condition = samplemeta.cancer[
  match(dtf.score.cancer$sampleid.pid, samplemeta.cancer$sampleid.pid)]$type

# every gene in the pseudobulk table, i.e. everything that is not an annotation
features = setdiff(colnames(dtfgenes.mean),
                   c("sampleid.pid", "pid", "group", "tissue", "source", "sex",
                     "age", "condition", "group.detail"))

table(dtf.score.cancer$condition)

dtf.score.cancer[, class := ifelse(condition %in% "cancer", 1, NA)]
dtf.score.cancer[, class := ifelse(condition %in% "adjacent", 0, class)]

# wu2024 pooled three cancer types; split it so each is its own dataset
dtf.score.cancer[, pid := ifelse(pid == "wu2024" & group.detail == "GBC", "wu2024-GBC", pid)]
dtf.score.cancer[, pid := ifelse(pid == "wu2024" & group.detail == "HCC", "wu2024-HCC", pid)]
dtf.score.cancer[, pid := ifelse(pid == "wu2024", "wu2024-other", pid)]


# -- Meta-analysis -------------------------------------------------------------

pids = c("hu2022", "salcher2022", "wu2024-GBC", "wu2024-HCC", "wu2024-other")

metaobj = lapply(pids, function(i) {
  sub = dtf.score.cancer[pid == i]
  meta_dataset(
    expr = t(as.matrix(sub[, features, with = FALSE])), class = sub$class, label = i)
})

meta.expr.cancer = meta_analysis(metaobj, outcome_type = "binary")

saveRDS(meta.expr.cancer, paste0(
  dir.results, "meta_analysis/SCAHN.cancer.exprmeta.allgenes.allneu.rds"))
