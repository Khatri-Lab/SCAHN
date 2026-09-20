# ==============================================================================
# Meta-analysis of sex differences in neutrophil gene expression
#
# Healthy whole-blood samples, male vs. female. The per-sample geometric-mean
# expression of every gene is meta-analysed across the five studies that
# contribute enough samples of both sexes.
#
# The cohort selection below is the same one used for the proportion analysis in
# 05_sexdifference.R; edit both if the sample set changes.
#
# Requires 00_setup.R (data.table, expresso, dir.data / dir.results) and
# 01_load_data.R (samplemeta).
#
# Reads:
#   results/annotation/SCAHN.subsetproportion.rds
#   results/meta_analysis/dtfgenes.allneu.persample.geommean.rds
#
# Writes:
#   results/meta_analysis/SCAHN.sex.exprmeta.allgenes.allneu.rds   read by 05_sexdifference.R
# ==============================================================================

SCAHN_SCRIPTS = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/scripts/figures"
source(file.path(SCAHN_SCRIPTS, "00_setup.R"))      # packages, paths
source(file.path(SCAHN_SCRIPTS, "01_load_data.R"))  # samplemeta


# -- Healthy whole-blood samples with a known sex ------------------------------

dtf.prop = readRDS(paste0(dir.results, "annotation/SCAHN.subsetproportion.rds"))

dtf.prop[pid == "wigerblad2022"]$sampleid

dtf.prop.sex.a = dtf.prop[group == "healthy"][
  tissue %in% c("peripheral blood")][
  grepl("WB", source) | grepl("gupta2020|montaldo2022", pid)][
  pid != "wigerblad2022"]
dtf.prop.sex.b = dtf.prop[group == "healthy"][grepl("wigerblad2022", pid)][
  sampleid %in% c("Gran_F1", "Gran_F2", "Gran_F3", "Gran_M1", "Gran_M2",
                  "Gran_M3", "Neut_M0", "WB_oF4", "WB_oF5", "WB_oM4")]
dtf.prop.sex = rbind(dtf.prop.sex.a, dtf.prop.sex.b)

dtf.prop.sex$sex = samplemeta[match(dtf.prop.sex$sampleid.pid, samplemeta$sampleid.pid)]$sex
dtf.prop.sex$age = samplemeta[match(dtf.prop.sex$sampleid.pid, samplemeta$sampleid.pid)]$age
dtf.prop.sex = dtf.prop.sex[!sex %in% c("unknown")]


# -- Per-sample expression, male vs. female ------------------------------------

dtfgenes.mean = readRDS(paste0(
  dir.results, "meta_analysis/dtfgenes.allneu.persample.geommean.rds"))

dtf.score.sex = dtfgenes.mean[sampleid.pid %in% dtf.prop.sex$sampleid.pid]
dtf.score.sex$sex = dtf.prop.sex[match(dtf.score.sex$sampleid.pid, dtf.prop.sex$sampleid.pid)]$sex

dtf.score.sex[, class := ifelse(sex %in% "M", 1, NA)]
dtf.score.sex[, class := ifelse(sex %in% "F", 0, class)]


# -- Meta-analysis -------------------------------------------------------------

# the studies with enough samples of both sexes
pids = c("gupta2020", "kaiser2024", "kwok2023", "wigerblad2022", "sinha2021")

# every gene in the pseudobulk table
features = setdiff(colnames(dtfgenes.mean), c("sampleid.pid", "pid"))

metaobj = lapply(pids, function(i) {
  sub = dtf.score.sex[pid == i]
  meta_dataset(
    expr = t(as.matrix(sub[, features, with = FALSE])), class = sub$class, label = i)
})

meta.expr.sex = meta_analysis(metaobj, outcome_type = "binary")

saveRDS(meta.expr.sex, paste0(dir.results, "meta_analysis/SCAHN.sex.exprmeta.allgenes.allneu.rds"))
