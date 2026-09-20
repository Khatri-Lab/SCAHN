# ==============================================================================
# Clustering sweep: alluvial of the selected clusterings, and cell stability
#
# Last of the three clustering-evaluation steps, and the only one that draws
# anything. The alluvial shows how cells move between the selected clusterings;
# the histogram shows how stable each cell's assignment was across them.
#
# The three steps run in order, each from its own script:
#   1. 05_run_cluster_aggregate.R   the sweep table and the selection
#   2. 06a_run_evalcluster.py        concordance, including the stability scores
#   3. 07_plot_cluster_eval.R       this script
#
# Requires data.table, ggplot2, ggalluvial, grid and expresso.
#
# Reads:
#   results/integration/clustering/clusterings.csv
#   results/integration/clustering/clusternumber.sub.csv
#     both written by 05_run_cluster_aggregate.R
#   results/integration/clustering/concordance_results/cell_stability_scores.csv
#     written by 06a_run_evalcluster.py
#
# Writes:
#   figures/SCAHN.clustering.alluvial.pdf
#   figures/SCAHN.cellstabilityscore.pdf
# ==============================================================================

library(data.table)
library(ggplot2)
library(ggalluvial)
library(grid)
library(expresso)

dir.results    = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/results/"
dir.fig        = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/figures/original/"
dir.clustering = paste0(dir.results, "integration/clustering/")
dir.concord    = paste0(dir.clustering, "concordance_results/")


# -- Helper --------------------------------------------------------------------

# Replace any cluster smaller than N cells with replace_val. Note it coerces
# every column to character, so the labels come back as strings and are
# compared as "-1" below.
replace_small_clusters = function(mat, N, replace_val = -1L) {
  dt = data.table::as.data.table(mat)

  # coerce all columns once; avoids set()
  dt[, (names(dt)) := lapply(.SD, as.character)]
  repv = as.character(replace_val)

  for (col in names(dt)) {
    counts = dt[, .(n = .N), by = .(val = get(col))]
    small  = counts[n < N, val]
    if (length(small)) dt[get(col) %chin% small, (col) := repv]
  }
  dt[]
}


# -- Alluvial of the selected clusterings --------------------------------------

multiclustering = fread(paste0(dir.clustering, "clusterings.csv"))
selected = fread(paste0(dir.clustering, "clusternumber.sub.csv"))
setdiff(selected$name, colnames(multiclustering)) # character(0)
multiclustering = multiclustering[, selected$name, with = FALSE]

# Clusters under 100 cells become -1, and any cell carrying that label in any of
# the selected clusterings is dropped, so the ribbons are not shredded by strata
# too small to read.
multiclustering2 = replace_small_clusters(multiclustering, N = 100, replace_val = -1)
multiclustering2[, sum(Reduce("|", lapply(.SD, function(x) x == "-1")))] # 2194

mask = Reduce("|", lapply(multiclustering2, function(x) x == "-1"))
multiclustering2 = multiclustering2[!mask]

# Axis labels: strip the resolution, and name the two Louvain variants rather
# than numbering them. The .seed alternative is a leftover from the seed sweep
# and no longer matches anything.
colnames(multiclustering2) = gsub(".seed[0-9]+|.res.*$", "", colnames(multiclustering2))
colnames(multiclustering2) = gsub("algo1", "louvain1", colnames(multiclustering2))
colnames(multiclustering2) = gsub("algo2", "louvain2", colnames(multiclustering2))

p = plot_alluvial(multiclustering2, colnames(multiclustering2),
                  top_n_flows = 100000, min_count = 20, stratum_linewidth = 0.1) +
  theme(axis.text.x = element_text(angle = 90), legend.key.height = unit(1, "line"))

pdf(file = paste0(dir.fig, "SCAHN.clustering.alluvial.pdf"), width = 8, height = 5)
pushViewport(viewport(layout = grid.layout(nrow = 100, ncol = 100)))
print(p, vp = viewport(layout.pos.row = 1:100, layout.pos.col = 1:100))
dev.off()


# -- Per-cell stability score --------------------------------------------------

df.cellstability = fread(paste0(dir.concord, "cell_stability_scores.csv"))

meanscore = mean(df.cellstability$stability_score)

p = ggplot(df.cellstability, aes(x = stability_score)) +
  geom_histogram(bins = 100, fill = "#214DC8", color = "white", alpha = 0.8) +
  geom_vline(xintercept = meanscore, color = "#A50F15", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 0.75, y = Inf, vjust = 1.5, hjust = 0,
           label = paste0("mean = ", round(meanscore, 3)), color = "#A50F15", size = 4) +
  theme_expresso() +
  labs(title = "distribution of per-cell clustering stability", x = "cell stability score",
       y = "count")

cairo_pdf(file = paste0(dir.fig, "SCAHN.cellstabilityscore.pdf"), width = 5, height = 3)
pushViewport(viewport(layout = grid.layout(nrow = 100, ncol = 100)))
print(p, vp = viewport(layout.pos.row = 1:100, layout.pos.col = 1:100))
dev.off()
