# ==============================================================================
# Table S1: per-study sample, tissue and group summary
#
# One row per study, giving its cell count, its sample count, and the tissues
# and clinical groups it contributes as "value(n)" strings. Cell counts come
# from the per-cell table and everything else from one row per sample, which is
# what the unique() below builds.
#
# Counts are over the atlas as 01_load_data.R leaves it: samples under 30 cells
# are already dropped, and the a/b/c/d suffixes are already stripped from pid,
# so a study's sub-cohorts are pooled into one row here.
#
# The pid values go out raw -- myin2023, reyfman2018, xue2022 -- where the
# figures draw the corrected publication years myin2024, reyfman2019, xue2023.
# That rename is done per script after 01_load_data.R, in 02_atlas.R and
# 04_infection_cancer.R, and this script does not do it.
#
# Requires 00_setup.R (data.table, dir.results) and 01_load_data.R (dtf).
# Reads nothing else: every input arrives through dtf.
#
# Writes:
#   results/tables/SCAHN.tableS1.csv
#     read by no script; the manuscript's supplementary table 1
# ==============================================================================

SCAHN_SCRIPTS = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/scripts/figures"
source(file.path(SCAHN_SCRIPTS, "00_setup.R"))      # packages, paths
source(file.path(SCAHN_SCRIPTS, "01_load_data.R"))  # dtf


# -- One row per sample, and the summary helper --------------------------------

sample_dt = unique(dtf[, .(pid, sampleid.pid, tissue, group)])

# "liver(3); lung(5); other(2)": rarest first, with "other" forced last rather
# than into the place its own count would give it.
make_summary = function(x) {
  counts = sort(table(x))
  nms = names(counts)
  if ("other" %in% nms) {
    nms = c(nms[nms != "other"], "other")
  }
  paste0(nms, "(", counts[nms], ")", collapse = "; ")
}


# -- One row per study ---------------------------------------------------------

summary_dt = dtf[, .(n_cells = .N), by = pid][
  sample_dt[, .(n_samples = .N, tissue = make_summary(tissue),
                group = make_summary(group)), by = pid], on = "pid"]

setorder(summary_dt, pid)
fwrite(summary_dt, paste0(dir.results, "tables/SCAHN.tableS1.csv"))
