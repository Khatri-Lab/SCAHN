# ==============================================================================
# Meta-analysis of infection vs. healthy neutrophil gene expression
#
# Whole-blood samples from the infection studies, in three contrasts against a
# pooled healthy / convalescent baseline: severe COVID-19, bacterial sepsis, and
# non-severe COVID-19 / flu. Each contrast is run twice:
#   1. on per-sample pseudobulk within the granule and IFN subset categories
#   2. on per-sample pseudobulk over all neutrophils pooled
#
# The cohort selection below is the same one used for the proportion analysis in
# 04_infection_cancer.R; edit both if the sample set changes.
#
# Requires 00_setup.R (data.table, expresso, dir.data / dir.results) and
# 01_load_data.R (dtf, samplemeta).
#
# Reads:
#   results/annotation/SCAHN.subsetproportion.rds
#   results/meta_analysis/dtfgenes.subsets_categ.persample.geommean.rds
#   results/meta_analysis/dtfgenes.allneu.persample.geommean.rds
#
# Writes:
#   results/meta_analysis/SCAHN.infect.exprmeta.selectedgenes.subsets.rds
#   results/meta_analysis/SCAHN.infect.exprmeta.selectedgenes.allneu.rds
#     both read by 04_infection_cancer.R
# ==============================================================================

SCAHN_SCRIPTS = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/scripts/figures"
source(file.path(SCAHN_SCRIPTS, "00_setup.R"))      # packages, paths
source(file.path(SCAHN_SCRIPTS, "01_load_data.R"))  # dtf, samplemeta


# -- Per-sample expression within subset categories ----------------------------

dtf.prop = readRDS(paste0(dir.results, "annotation/SCAHN.subsetproportion.rds"))

dtfgenes.mean = readRDS(paste0(
  dir.results, "meta_analysis/dtfgenes.subsets_categ.persample.geommean.rds"))

dtfgenes.mean$group = samplemeta[match(dtfgenes.mean$sampleid.pid, samplemeta$sampleid.pid)]$group
dtfgenes.mean$tissue = samplemeta[match(dtfgenes.mean$sampleid.pid, samplemeta$sampleid.pid)]$tissue
dtfgenes.mean$source = samplemeta[match(dtfgenes.mean$sampleid.pid, samplemeta$sampleid.pid)]$source

### whole-blood samples from the infection studies
dtf.prop.infection =
  dtf.prop[pid %in% unique(dtf[grepl("COVID|sepsis|flu|conv", group)]$pid)][
    grepl("COVID|sepsis|flu|conv|healthy", group)][
    tissue %in% c("peripheral blood")][
    grepl("WB", source)]
dtf.prop.infection = dtf.prop.infection[n_cells >= 50]

dtf.score.infect =
  dtfgenes.mean[sampleid.pid %in% dtf.prop.infection$sampleid.pid][
    grepl("COVID|sepsis|flu|conv|healthy", group)][
    tissue %in% c("peripheral blood")][
    grepl("WB", source)]

# Collapse the clinical groups into the three contrasts plus one baseline. Flu
# and non-severe COVID-19 are pooled; convalescent samples join the healthy
# baseline.
dtf.score.infect[, condition := ifelse(
  group %in% c("COVID-19,nonsevere", "flu", "flu,pregnancy"), "COVID-19/flu,nonsevere",
  ifelse(group %in% c("healthy", "convalescent"), "healthy/convalescent", group))]

length(unique(dtf.prop.infection$sampleid.pid))
length(unique(dtf.score.infect$sampleid.pid))

dtf.score.infect$pid = gsub("[ab]$", "", dtf.score.infect$pid)

# every gene in the pseudobulk table, i.e. everything that is not an annotation
features = setdiff(colnames(dtfgenes.mean), c("sampleid.pid", "pid", "celltype2", "group", "tissue",
                     "source"))

run_meta_pipeline = function(target_condition, target_pids, data, cell_type, feature_list) {

  # subsetting makes a copy, so the caller's data.table is left alone
  sub_data = data[pid %in% target_pids & celltype2 == cell_type]

  # 1 for the target condition, 0 for the baseline, NA for everything else
  sub_data[, class := fcase(
    condition == target_condition, 1, condition == "healthy/convalescent", 0, default = NA)]

  metaobj = lapply(target_pids, function(p) {
    p_data = sub_data[pid == p & !is.na(class)]
    meta_dataset(
      expr = t(as.matrix(p_data[, feature_list, with = FALSE])), class = p_data$class, label = p)
  })

  meta_analysis(metaobj, outcome_type = "binary")
}

subsets = c("granule subsets", "IFN subsets")

# the studies contributing both the condition and the baseline
configs = list(
  "COVID-19,severe" = c("combes2021", "schrepping2020", "sinha2021", "wilk2021"),
  "bacterial sepsis" = c("combes2021", "kaiser2024", "kwok2023", "sinha2021"),
  "COVID-19/flu,nonsevere" = c("combes2021", "schrepping2020", "zhang2023"))

all_results = setNames(lapply(subsets, function(celltypevar) {
  lapply(names(configs), function(cond) {
    run_meta_pipeline(
      target_condition = cond, target_pids = configs[[cond]],
      data = dtf.score.infect, cell_type = celltypevar, feature_list = features)
  }) |> setNames(names(configs))
}), subsets)

saveRDS(all_results, paste0(
  dir.results, "meta_analysis/SCAHN.infect.exprmeta.selectedgenes.subsets.rds"))


# -- Per-sample expression over all neutrophils --------------------------------

dtfgenes.mean = readRDS(paste0(
  dir.results, "meta_analysis/dtfgenes.allneu.persample.geommean.rds"))

dtfgenes.mean$group = samplemeta[match(dtfgenes.mean$sampleid.pid, samplemeta$sampleid.pid)]$group
dtfgenes.mean$tissue = samplemeta[match(dtfgenes.mean$sampleid.pid, samplemeta$sampleid.pid)]$tissue
dtfgenes.mean$source = samplemeta[match(dtfgenes.mean$sampleid.pid, samplemeta$sampleid.pid)]$source

dtf.score.infect =
  dtfgenes.mean[sampleid.pid %in% dtf.prop.infection$sampleid.pid][
    grepl("COVID|sepsis|flu|conv|healthy", group)][
    tissue %in% c("peripheral blood")][
    grepl("WB", source)]

dtf.score.infect[, condition := ifelse(
  group %in% c("COVID-19,nonsevere", "flu", "flu,pregnancy"), "COVID-19/flu,nonsevere",
  ifelse(group %in% c("healthy", "convalescent"), "healthy/convalescent", group))]

dtf.score.infect = dtf.score.infect[
  sampleid.pid %in% unique(dtf.prop.infection$sampleid.pid)]
length(unique(dtf.prop.infection$sampleid.pid))
length(unique(dtf.score.infect$sampleid.pid))
dtf.score.infect$pid = gsub("[ab]$", "", dtf.score.infect$pid)

features = setdiff(colnames(dtfgenes.mean), c("sampleid.pid", "pid", "celltype2", "group", "tissue",
                     "source"))

# Same helper as above without the celltype2 filter: one row per sample, all
# neutrophils pooled.
run_meta_pipeline = function(target_condition, target_pids, data, feature_list) {

  # subsetting makes a copy, so the caller's data.table is left alone
  sub_data = data[pid %in% target_pids]

  # 1 for the target condition, 0 for the baseline, NA for everything else
  sub_data[, class := fcase(
    condition == target_condition, 1, condition == "healthy/convalescent", 0, default = NA)]

  metaobj = lapply(target_pids, function(p) {
    p_data = sub_data[pid == p & !is.na(class)]
    meta_dataset(
      expr = t(as.matrix(p_data[, feature_list, with = FALSE])), class = p_data$class, label = p)
  })

  meta_analysis(metaobj, outcome_type = "binary")
}

configs = list(
  "COVID-19,severe" = c("combes2021", "schrepping2020", "sinha2021", "wilk2021"),
  "bacterial sepsis" = c("combes2021", "kaiser2024", "kwok2023", "sinha2021"),
  "COVID-19/flu,nonsevere" = c("combes2021", "schrepping2020", "zhang2023"))

all_results =
  lapply(names(configs), function(cond) {
    run_meta_pipeline(
      target_condition = cond, target_pids = configs[[cond]],
      data = dtf.score.infect, feature_list = features)
  }) |> setNames(names(configs))

saveRDS(all_results, paste0(
  dir.results, "meta_analysis/SCAHN.infect.exprmeta.selectedgenes.allneu.rds"))
