# ==============================================================================
# Shared-nearest-neighbour graph from the harmony-corrected embedding
#
# Builds the SNN graph that the two clustering scripts partition, and writes it
# twice: as RDS for the Seurat/Louvain path and as matrix market for the
# Python/Leiden path. One invocation per point of the sweep, matching
# 02a_run_harmony_SCT.R:
#
#   Rscript 03_run_SNN.R v2000 25
#
# Requires Seurat, data.table and Matrix.
#
# Reads:
#   results/integration/harmony/SCAHN.harmony.Z_corr.<vargenes>.PC<nPCs>.csv
#     written by 02a_run_harmony_SCT.R
#   data/SCAHN.simplemeta.csv   supplies the cell barcodes, which the corrected
#                               embedding does not carry
#
# Writes:
#   results/integration/clustering/SCAHN.snn_graph.<vargenes>.PC<nPCs>.rds
#     read by 04a_run_clusteringlouvain.R
#   results/integration/clustering/SCAHN.snn_graph.<vargenes>.PC<nPCs>.mtx
#     read by 04b_run_clusteringleiden.py
# ==============================================================================

library(Seurat)
library(data.table)
library(Matrix)

options(future.globals.maxSize = 20000 * 1024^2)

dir.data       = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/data/"
dir.results    = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/results/"
dir.harmony    = paste0(dir.results, "integration/harmony/")
dir.clustering = paste0(dir.results, "integration/clustering/")


# -- Parameters ----------------------------------------------------------------

args = commandArgs(trailingOnly = TRUE)
vargenes = args[1] # "v1000" "v2000" "v3000"
nPCs = args[2]     # 20 25 30


# -- Corrected embedding -------------------------------------------------------

dtf = fread(paste0(dir.data, "SCAHN.simplemeta.csv"))
harmony.Z_corr = fread(paste0(dir.harmony, "SCAHN.harmony.Z_corr.", vargenes, ".PC", nPCs, ".csv"))
dim(harmony.Z_corr)

# Row order matches simplemeta, so the barcodes can be attached positionally.
harmony.Z_corr = as.matrix(harmony.Z_corr)
rownames(harmony.Z_corr) = dtf$cell.pid


# -- Neighbour graph -----------------------------------------------------------

# FindNeighbors builds the k = 20 nearest-neighbour graph and derives the
# shared-nearest-neighbour graph from it. Only the SNN is kept: it is what both
# clustering scripts partition, and the KNN graph is not used again.
set.seed(42)
nn_graph = FindNeighbors(
  object = harmony.Z_corr, k.param = 20, compute.SNN = TRUE, verbose = FALSE,
  return.neighbor = FALSE)

snn_graph = nn_graph$snn

saveRDS(snn_graph, paste0(dir.clustering, "SCAHN.snn_graph.", vargenes, ".PC", nPCs, ".rds"))

# The same graph in matrix market form, for the Python side
writeMM(snn_graph, paste0(dir.clustering, "SCAHN.snn_graph.", vargenes, ".PC", nPCs, ".mtx"))
