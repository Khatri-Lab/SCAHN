# ==============================================================================
# Symphony reference from the harmony-corrected SCT embedding
#
# Builds the symphony reference that query datasets are mapped onto, together
# with the UMAP model that places those queries in the published atlas
# coordinates. Last stage of the integration pipeline: it labels the reference
# with celltype, so it runs once the annotation has settled.
#
# This script and 08a_anno.R each read the other's output: the celltype column
# comes from 08a, and the UMAP coordinates 08a places its cells with come from
# here. They were run alternately, so neither reproduces the published files on
# its own. The a/b numbering records that: they are one step, not two.
#
# PCA and harmony are re-run here rather than read back from 01 and 02a, because
# symphony needs two things the CSV exports do not carry: the PCA feature
# loadings and the harmony object itself. The two verification loops check that
# the re-run lands on the embeddings already on disk.
#
# Expensive: reads the integrated Seurat object, re-runs PCA and harmony over
# all 1,053,128 cells, then fits a UMAP model. Run once, on its own.
#
# Requires Seurat, symphony, data.table, tibble, Matrix, harmony and uwot.
#
# Reads:
#   data/SCAHN/SCAHN.srt.integrated.rds         written by 01_run_integration.R
#   data/SCAHN/SCAHN.SCT.var.features.list.rds  written by 01_run_integration.R
#   data/SCAHN/SCAHN.pca.v2000.emb.csv          written by 01_run_integration.R,
#                                               for the sign-flip check only
#   results/integration/harmony/SCAHN.harmony.Z_corr.v2000.PC25.csv
#     written by 02a_run_harmony_SCT.R
#   results/annotation/SCAHN.cellmeta.csv       written by 08a_anno.R; one row
#                                               per cell, in the same order as
#                                               the Seurat object, and the
#                                               source of celltype
#
# Writes, all into results/symphony/:
#   SCAHN.harmony_SCT.v2000.PC25.uwot.model
#   SCAHN.harmony_SCT.v2000.PC25.uwot.emb.csv
#     read by 08a_anno.R, and by SCAHN.anno.tmp.Rmd, which still points at
#     the data/ copy
#   SCAHN.symphony.reference.rds
# ==============================================================================

library(Seurat)
library(symphony)
library(data.table)
library(tibble)
library(Matrix)
library(harmony)
library(uwot)

options(future.globals.maxSize = 20000 * 1024^2)

dir.data        = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/data/"
dir.results     = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/results/"
dir.integration = paste0(dir.data, "SCAHN/")
dir.annotation  = paste0(dir.results, "annotation/")
dir.harmony     = paste0(dir.results, "integration/harmony/")
dir.symphony    = paste0(dir.results, "symphony/")


# -- Parameters ----------------------------------------------------------------

vargenes = "v2000"
nPCs = 25

path.uwot = paste0(dir.symphony, "SCAHN.harmony_SCT.", vargenes, ".PC", nPCs, ".uwot.model")


# -- The integrated object and cell metadata -----------------------------------

srt.sct = readRDS(paste0(dir.integration, "SCAHN.srt.integrated.rds"))
dtf = fread(paste0(dir.annotation, "SCAHN.cellmeta.csv"))

# Everything below joins metadata to the embedding by position, so the two have
# to carry the same cells in the same order.
all.equal(dtf$cell.pid, as.character(srt.sct$cell.pid))


# -- PCA -----------------------------------------------------------------------

var.features.list = readRDS(paste0(dir.integration, "SCAHN.SCT.var.features.list.rds"))
var_genes = var.features.list[[vargenes]]

VariableFeatures(srt.sct) = var_genes
DefaultAssay(srt.sct) = "SCT"
srt.sct = RunPCA(srt.sct, features = var_genes, seed.use = 42)

### does the re-run match the embedding 01_run_integration.R wrote?

# The CSV has no cell names, so the comparison is on values only.
Z_pca_ref2 = as.matrix(fread(paste0(dir.integration, "SCAHN.pca.", vargenes, ".emb.csv")))

all.equal(Z_pca_ref2, srt.sct@reductions$pca@cell.embeddings, check.attributes = FALSE)
for (i in 1:50) {
  a = Z_pca_ref2[, i]
  b = unname(srt.sct@reductions$pca@cell.embeddings[, i])
  same = all.equal(a, b)
  flipped = all.equal(a, -b)
  cat(sprintf("PC%d: same=%s  flipped=%s\n", i, same, flipped))
}

# A mean relative difference of exactly 2 is the signature of a sign flip: PCA
# eigenvectors are defined only up to sign, so multiplying one by -1 is an
# equally valid solution. Solvers such as IRLBA do not fix that sign
# deterministically, and it can change with library version, BLAS backend or
# machine even at the same seed. Distances are sign-invariant, so harmony,
# clustering, UMAP and kNN downstream are all unaffected.


# -- Harmony -------------------------------------------------------------------

Z_pca_ref = srt.sct@reductions$pca@cell.embeddings[, 1:nPCs]
loadings = srt.sct@reductions$pca@feature.loadings[, 1:nPCs]

set.seed(8)
# max_iter = 10 is harmony's default and what produced the on-disk embeddings.
ref_harmObj = harmony::HarmonyMatrix(
  data_mat = Z_pca_ref, meta_data = dtf[, c("cell.pid", "sampleid.pid")],
  vars_use = c("sampleid.pid"), max_iter = 10, nclust = 100, return_object = TRUE,
  do_pca = FALSE)

### does the re-run match what 02a_run_harmony_SCT.R wrote?

# Z_corr comes out factors-by-cells, so transpose back to cells-by-factors
harmony.Z_corr = as.data.table(t(ref_harmObj$Z_corr))
harmony.Z_corr2 = fread(paste0(dir.harmony, "SCAHN.harmony.Z_corr.", vargenes, ".PC", nPCs, ".csv"))

all.equal(harmony.Z_corr2, harmony.Z_corr)
for (i in 1:nPCs) {
  a = harmony.Z_corr2[[i]]
  b = harmony.Z_corr[[i]]
  same = all.equal(a, b)
  flipped = all.equal(a, -b)
  cat(sprintf("PC%d: same=%s  flipped=%s\n", i, same, flipped))
}


# -- UMAP model ----------------------------------------------------------------

# Fitted on the corrected embedding as written to disk rather than on the re-run
# above, so that the model reproduces the published coordinates.
umap_model = uwot::umap(
  harmony.Z_corr2, n_neighbors = 30, n_components = 2, metric = "cosine",
  min_dist = 0.5, spread = 1, learning_rate = 1.0, set_op_mix_ratio = 1.0,
  local_connectivity = 1L, repulsion_strength = 1.0, negative_sample_rate = 5, init = "spectral",
  n_epochs = NULL, ret_model = TRUE, ret_nn = TRUE, n_sgd_threads = 1, verbose = TRUE, seed = 42)

uwot::save_uwot(umap_model, path.uwot)
fwrite(umap_model$embedding, paste0(dir.symphony, "SCAHN.harmony_SCT.", vargenes, ".PC", nPCs,
              ".uwot.emb.csv"))


# -- Gene means and standard deviations ----------------------------------------

# symphony scales query data with these, so they have to come from the same
# assay and the same gene set the reference is built on.
sct.data = srt.sct@assays[["SCT"]]@data[var_genes, ]

vargenes_means_sds = tibble(
  symbol = var_genes, mean = Matrix::rowMeans(sct.data))
vargenes_means_sds$stddev = symphony::rowSDs(sct.data, vargenes_means_sds$mean)


# -- Symphony reference --------------------------------------------------------

reference = buildReferenceFromHarmonyObj(
  ref_harmObj, dtf[, c("cell.pid", "sampleid.pid", "celltype")], vargenes_means_sds, loadings,
  verbose = TRUE, do_umap = FALSE, save_uwot_path = NULL)

# do_umap = FALSE above, so the model fitted here is attached by hand. To
# rebuild the reference without refitting, load the saved model instead:
#   umap_model = uwot::load_uwot(path.uwot)
reference$umap = umap_model
reference$save_uwot_path = path.uwot

saveRDS(reference, paste0(dir.symphony, "SCAHN.symphony.reference.rds"))
