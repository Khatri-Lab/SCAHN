# ==============================================================================
# Harmony batch correction of the SCT embedding
#
# Corrects one PCA embedding for sample of origin and writes the corrected cell
# factors. One invocation per point of the sweep, three feature-set sizes by
# three PC counts, so nine in total:
#
#   Rscript 02a_run_harmony_SCT.R v2000 25
#
# Requires data.table and harmony.
#
# Reads:
#   data/SCAHN/SCAHN.pca.<vargenes>.emb.csv   written by 01_run_integration.R
#   data/SCAHN.simplemeta.csv                 one row per cell, in the same
#                                             order as the embedding
#
# Writes:
#   results/integration/harmony/SCAHN.harmony.Z_corr.<vargenes>.PC<nPCs>.csv
#     read by 03_run_SNN.R and 04b_run_clusteringleiden.py
# ==============================================================================

library(data.table)
library(harmony)

options(future.globals.maxSize = 20000 * 1024^2)

dir.data        = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/data/"
dir.results     = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/results/"
dir.integration = paste0(dir.data, "SCAHN/")
dir.harmony     = paste0(dir.results, "integration/harmony/")


# -- Parameters ----------------------------------------------------------------

args = commandArgs(trailingOnly = TRUE)
vargenes = args[1] # "v1000" "v2000" "v3000"
nPCs = args[2]     # 20 25 30


# -- Embedding and cell metadata -----------------------------------------------

Z_pca_ref = fread(paste0(dir.integration, "SCAHN.pca.", vargenes, ".emb.csv"))
dtf = fread(paste0(dir.data, "SCAHN.simplemeta.csv"))

# The embedding is written with all PCs; keep only as many as this run sweeps.
Z_pca_ref = Z_pca_ref[, 1:nPCs]


# -- Harmony -------------------------------------------------------------------

set.seed(8)
# max_iter = 10 is harmony's default and what produced the on-disk embeddings.
ref_harmObj = harmony::HarmonyMatrix(
  data_mat = Z_pca_ref, meta_data = dtf[, c("cell.pid", "sampleid.pid", "pid")],
  vars_use = c("sampleid.pid"), max_iter = 10, return_object = TRUE, nclust = 100,
  do_pca = FALSE)

# Z_corr comes out factors-by-cells, so transpose back to cells-by-factors
harmony.Z_corr = ref_harmObj$Z_corr

fwrite(t(harmony.Z_corr), paste0(dir.harmony, "SCAHN.harmony.Z_corr.", vargenes, ".PC", nPCs,
              ".csv"))
