# ==============================================================================
# Louvain clustering of the shared-nearest-neighbour graph
#
# Partitions the SNN graph at one resolution with one Louvain variant. The R
# half of the clustering sweep; 04b_run_clusteringleiden.py is the other half.
# One invocation per point of the sweep, three feature-set sizes by three PC
# counts by two algorithms by three resolutions:
#
#   Rscript 04a_run_clusteringlouvain.R v2000 25 0.8 1
#
# algo is passed straight through to FindClusters: 1 is the original Louvain, 2
# is Louvain with multilevel refinement. Both are in the sweep, and the value is
# recorded in the output name as .algo1 / .algo2.
#
# Requires Seurat and Matrix.
#
# Reads:
#   results/integration/clustering/SCAHN.snn_graph.<vargenes>.PC<nPCs>.rds
#     written by 03_run_SNN.R
#
# Writes:
#   results/integration/clustering/
#     SCAHN.clustering.<vargenes>.PC<nPCs>.algo<algo>.res<reso>.rds
#     read by 05_run_cluster_aggregate.R
# ==============================================================================

library(Seurat)
library(Matrix)

options(future.globals.maxSize = 20000 * 1024^2)

dir.results    = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/results/"
dir.clustering = paste0(dir.results, "integration/clustering/")


# -- Parameters ----------------------------------------------------------------

args = commandArgs(trailingOnly = TRUE)
vargenes = args[1]         # "v1000" "v2000" "v3000"
nPCs = as.integer(args[2]) # 20 25 30
reso = as.numeric(args[3]) # 0.5 0.8 1.2
algo = as.integer(args[4]) # 1 2


# -- Clustering ----------------------------------------------------------------

snn_graph = readRDS(paste0(dir.clustering, "SCAHN.snn_graph.", vargenes, ".PC", nPCs, ".rds"))

set.seed(42)
clusters = FindClusters(object = snn_graph, resolution = reso, algorithm = algo, verbose = TRUE)

saveRDS(clusters, paste0(dir.clustering, "SCAHN.clustering.", vargenes,
                         ".PC", nPCs, ".algo", algo, ".res", reso, ".rds"))
