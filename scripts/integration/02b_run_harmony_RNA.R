# ==============================================================================
# Harmony batch correction of the RNA embedding
#
# The RNA counterpart of 02a_run_harmony_SCT.R
#
# Requires data.table and harmony.
#
# Reads:
#   data/SCAHN/SCAHN.beforeintegration.pca.v2000.emb.csv
#     written by 01_run_integration.R
#   data/SCAHN.simplemeta.csv   one row per cell, in the same order as the
#                               embedding
#
# Writes:
#   results/integration/harmony/SCAHN.harmony_RNA.Z_corr.v2000.PC25.csv
# ==============================================================================

library(data.table)
library(harmony)

options(future.globals.maxSize = 20000 * 1024^2)

dir.data        = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/data/"
dir.results     = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/results/"
dir.integration = paste0(dir.data, "SCAHN/")
dir.harmony     = paste0(dir.results, "integration/harmony/")


# -- Parameters ----------------------------------------------------------------

vargenes = "v2000"
nPCs = 25


# -- Embedding and cell metadata -----------------------------------------------

Z_pca_ref = fread(paste0(dir.integration, "SCAHN.beforeintegration.pca.v2000.emb.csv"))
dtf = fread(paste0(dir.data, "SCAHN.simplemeta.csv"))

# The embedding is written with all PCs; keep only as many as this run uses.
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

fwrite(t(harmony.Z_corr), paste0(dir.harmony, "SCAHN.harmony_RNA.Z_corr.", vargenes, ".PC", nPCs,
              ".csv"))
