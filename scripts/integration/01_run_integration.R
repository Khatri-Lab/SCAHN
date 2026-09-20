# ==============================================================================
# Integration: variable features, SCTransform, merge, and PCA
#
# First stage of the atlas pipeline. Turns the per-study Seurat objects into one
# merged object plus the PCA embeddings that 02a_run_harmony_SCT.R and
# 02b_run_harmony_RNA.R then batch-correct. Four steps:
#   1. RNA variable features, at each size of the sweep
#   2. SCTransform every study, then SCT variable features at the same sizes
#   3. merge into a single object, and export the raw counts as MTX
#   4. PCA on each feature set, once from SCT and once from RNA
#
# Expensive: reads the 1.4 GB per-study list and SCTransforms every study. Run
# once, on its own.
#
# Requires Seurat, data.table and Matrix. SCTransform also needs glmGamPoi,
# which it loads by namespace rather than attaching.
#
# Reads:
#   data/SCAHN/SCAHN.seurat.list.rds
#
# Writes:
#   data/SCAHN/SCAHN.RNA.var.features.list.rds
#   data/SCAHN/SCAHN.SCT.var.features.list.rds
#   data/SCAHN/SCAHN.pca.{v1000,v2000,v3000}.emb.csv
#     read by 02a_run_harmony_SCT.R
#   data/SCAHN/SCAHN.beforeintegration.pca.{v1000,v2000,v3000}.emb.csv
#     read by 02b_run_harmony_RNA.R, which by design takes only v2000: the RNA
#     arm is a single-configuration comparison, not a sweep
#   data/SCAHN/SCAHN.srt.integrated.rds
#     read by scores_and_pseudobulk.R and findmarkers.R
#   data/SCAHN/SCAHN.RNAcounts.mtx
#   data/SCAHN/SCAHN.RNAcounts.features.tsv
#   data/SCAHN/SCAHN.RNAcounts.barcodes.tsv
# ==============================================================================

library(Seurat)
library(data.table)
library(Matrix)

options(future.globals.maxSize = 20000 * 1024^2)

dir.data        = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/data/"
dir.integration = paste0(dir.data, "SCAHN/")


# -- Variable features, RNA ----------------------------------------------------

srtlist = readRDS(paste0(dir.integration, "SCAHN.seurat.list.rds"))

srtlist = lapply(X = srtlist, FUN = function(x) {
  DefaultAssay(x) = "RNA"
  x = NormalizeData(x)
  x = FindVariableFeatures(x, selection.method = "vst", nfeatures = 4000)
})

nfeatures.rna = c(3000, 2000, 1000)

var.features.list = setNames(
  lapply(nfeatures.rna, function(n) {
    SelectIntegrationFeatures(object.list = srtlist, nfeatures = n)
  }), paste0("v", nfeatures.rna))

saveRDS(var.features.list, paste0(dir.integration, "SCAHN.RNA.var.features.list.rds"))


# -- SCTransform, and variable features on SCT ---------------------------------

srtlist = lapply(
  X = srtlist, FUN = SCTransform, method = "glmGamPoi", seed.use = 42,
  assay = "RNA", return.only.var.genes = FALSE)

nfeatures.sct =  c(3000, 2000, 1000)

var.features.list = setNames(
  lapply(nfeatures.sct, function(n) {
    SelectIntegrationFeatures(object.list = srtlist, nfeatures = n)
  }), paste0("v", nfeatures.sct))

saveRDS(var.features.list, paste0(dir.integration, "SCAHN.SCT.var.features.list.rds"))


# -- Merge ---------------------------------------------------------------------

srt.sct = merge(x = srtlist[[1]], y = srtlist[2:length(srtlist)], merge.data = TRUE)

DefaultAssay(srt.sct) = "RNA"
srt.sct = NormalizeData(srt.sct)
saveRDS(srt.sct, paste0(dir.integration, "SCAHN.srt.integrated.rds"))


# -- Raw counts as MTX ---------------------------------------------------------

# The matrix market triplet plus the two name files, for tools that read counts
counts = GetAssayData(srt.sct, assay = "RNA", slot = "counts")
writeMM(counts, paste0(dir.integration, "SCAHN.RNAcounts.mtx"))
write.table(rownames(counts), paste0(dir.integration, "SCAHN.RNAcounts.features.tsv"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(colnames(counts), paste0(dir.integration, "SCAHN.RNAcounts.barcodes.tsv"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


# -- PCA over the feature-set sweep --------------------------------------------

# One PCA per feature set, from SCT and again from RNA. Each RunPCA overwrites
# the pca reduction, so each embedding is written out before the next iteration.
# Must be a subset of the names in the two var.features.list files above.
vargenes.pca = c("v1000", "v2000", "v3000")

### from SCT
var.features.list = readRDS(paste0(dir.integration, "SCAHN.SCT.var.features.list.rds"))
DefaultAssay(srt.sct) = "SCT"
for (v in vargenes.pca) {
  var_genes = var.features.list[[v]]
  VariableFeatures(srt.sct) = var_genes
  srt.sct = RunPCA(srt.sct, features = var_genes, seed.use = 42)
  fwrite(srt.sct@reductions$pca@cell.embeddings,
         paste0(dir.integration, "SCAHN.pca.", v, ".emb.csv"))
}

### from RNA, before any batch correction
var.features.list = readRDS(paste0(dir.integration, "SCAHN.RNA.var.features.list.rds"))
DefaultAssay(srt.sct) = "RNA"
srt.sct = ScaleData(srt.sct, features = unique(unlist(var.features.list[vargenes.pca])))

for (v in vargenes.pca) {
  var_genes = var.features.list[[v]]
  VariableFeatures(srt.sct) = var_genes
  srt.sct = RunPCA(srt.sct, features = var_genes, seed.use = 42)
  fwrite(srt.sct@reductions$pca@cell.embeddings,
         paste0(dir.integration, "SCAHN.beforeintegration.pca.", v, ".emb.csv"))
}
