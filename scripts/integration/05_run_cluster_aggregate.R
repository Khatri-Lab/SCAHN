# ==============================================================================
# Clustering sweep: aggregate the runs, pick one resolution per pipeline
#
# First of the three clustering-evaluation steps. Collects the 81 clusterings of
# the sweep into one table, records how many clusters each produced, and picks
# the single resolution per pipeline that brings every pipeline closest to a
# common cluster count.
#
# The three steps run in order, each from its own script:
#   1. 05_run_cluster_aggregate.R   this script
#   2. 06a_run_evalcluster.py        concordance of the selection made here
#   3. 07_plot_cluster_eval.R       the two summary figures
#
# Requires data.table only.
#
# Reads:
#   results/integration/clustering/SCAHN.clustering.*.rds
#     written by 04a_run_clusteringlouvain.R
#   results/integration/clustering/SCAHN.clustering.*.csv
#     written by 04b_run_clusteringleiden.py
#
# Writes:
#   results/integration/clustering/clusterings.csv
#     read by 06a_run_evalcluster.py, 06b_run_consensuscluster.py and
#     07_plot_cluster_eval.R
#   results/integration/clustering/clusternumber.csv
#   results/integration/clustering/clusternumber.sub.csv
#     read by 06a_run_evalcluster.py and 07_plot_cluster_eval.R
# ==============================================================================

library(data.table)

dir.results    = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/results/"
dir.clustering = paste0(dir.results, "integration/clustering/")


# -- Collect the sweep ---------------------------------------------------------

# Louvain wrote rds and Leiden wrote csv, but both hold one label per cell in
# the same cell order, so the two sets can be bound side by side. The pattern
# picks up only the per-run outputs, not the aggregates written further down.
files = list.files(dir.clustering, pattern = "SCAHN.clustering")
# 81 = 3 feature sets x 3 PC counts x 3 algorithms x 3 resolutions
length(files) # 81

files.rds = files[grepl("rds", files)]
files.csv = files[grepl("csv", files)]

clustering = list()

for (file in files.rds) {
  pre = gsub("SCAHN.clustering.|.rds", "", file)
  cl = readRDS(paste0(dir.clustering, file))
  # FindClusters returned a factor, and unlist keeps it one, so the labels have
  # to be recovered by value with as.character(): the levels sort lexically
  # ("0", "1", "10", "11", ...), so the underlying integer codes bear no relation
  # to the labels and would not line up with the 0-based Leiden ones.
  clustering[[pre]] = as.character(unlist(cl))
}

for (file in files.csv) {
  pre = gsub("SCAHN.clustering.|.csv", "", file)
  clustering[[pre]] = fread(paste0(dir.clustering, file))
}

df.clustering = as.data.table(do.call(cbind, clustering))
df.clustering[] = lapply(df.clustering, as.integer)
colnames(df.clustering) = names(clustering)

fwrite(df.clustering, paste0(dir.clustering, "clusterings.csv"))


# -- How many clusters each run produced ---------------------------------------

clusternumber = as.data.frame(t(df.clustering[, lapply(.SD, uniqueN)]))
clusternumber$name = rownames(clusternumber)
clusternumber = as.data.table(clusternumber)
clusternumber$pipeline = gsub(".res.*", "", clusternumber$name)
clusternumber$resolution = gsub(".*res", "", clusternumber$name)
colnames(clusternumber)[1] = "number"

fwrite(clusternumber, paste0(dir.clustering, "clusternumber.csv"))


# -- One resolution per pipeline -----------------------------------------------

# Comparing pipelines is only meaningful at a comparable granularity, so each
# pipeline contributes the one resolution that lands nearest a shared target
# cluster count. PC20 is dropped as the low end of the PC sweep.
df = clusternumber[grep("PC20", pipeline, invert = TRUE)]

# Sweep the target over every attainable count, and keep whichever target makes
# the resulting set of cluster numbers tightest.
targets = seq(min(df$number), max(df$number), by = 1)

best_var = Inf
best_target = NA

for (tgt in targets) {
  # One row per pipeline, the one nearest the target. which.min takes the first
  # minimum, so ties go to whichever resolution comes first in the sweep table.
  # order(pipeline) is needed to match the published row order: the sweep table
  # lists the Louvain runs before the Leiden ones, not alphabetically.
  df[, diff := abs(number - tgt)]
  selected = df[df[, .I[which.min(diff)], by = pipeline]$V1][order(pipeline)]

  v = var(selected$number)
  if (v < best_var) {
    best_var = v
    best_target = tgt
    best_selection = copy(selected)
  }
}

cat("Best target:", best_target, "\n")     # 58
cat("Variance:", best_var, "\n")           # 58.58824
cat("Range:", range(best_selection$number), "\n") # 51 78

print(best_selection[, .(pipeline, resolution, number)], nrows = Inf)

fwrite(best_selection, paste0(dir.clustering, "clusternumber.sub.csv"))
