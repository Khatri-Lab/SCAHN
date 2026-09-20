# ==============================================================================
# NeuMap mapped onto the human atlas by symphony
#
# Last stage of the mouse arm. Takes the human-symbol NeuMap object from
# NeuMap_orthologs.R, projects it into the atlas embedding with symphony, then
# transfers celltype labels from the reference by k nearest neighbours in that
# embedding: a plain majority vote, a distance-weighted vote, and the weighted
# vote's margin as a per-cell confidence.
#
#
# Requires Seurat, symphony and BiocNeighbors. uwot is called by namespace.
#
# Reads:
#   results/symphony/SCAHN.symphony.reference.rds
#   results/symphony/SCAHN.harmony_SCT.v2000.PC25.uwot.model
#     both written by 08b_symphony_reference.R
#   results/cross_species/NeuMap.135055.srt.rds
#     written by NeuMap_orthologs.R
#
# Writes, both into results/cross_species/:
#   NeuMap.signaturegenes.mat.rds
#   NeuMap.symphony.rds
#     both read by figures/03_cross_species.R
# ==============================================================================

library(Seurat)
library(symphony)
library(BiocNeighbors)

options(future.globals.maxSize = 20000 * 1024^2)

dir.results       = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/results/"
dir.symphony      = paste0(dir.results, "symphony/")
dir.cross.results = paste0(dir.results, "cross_species/")
SCAHN_FIGURES = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/scripts/figures"
source(file.path(SCAHN_FIGURES, "00_genesets.R"))

# -- Parameters ----------------------------------------------------------------

# Neighbours used for the label transfer. The vote matrix is reshaped to this
# width, so the two have to stay in step.
k.neighbours = 10

path.uwot = paste0(dir.symphony, "SCAHN.harmony_SCT.v2000.PC25.uwot.model")


# -- The atlas reference -------------------------------------------------------

reference = readRDS(paste0(dir.symphony, "SCAHN.symphony.reference.rds"))

umap_model = uwot::load_uwot(path.uwot)
reference$umap = umap_model
reference$save_uwot_path = path.uwot


# -- The query -----------------------------------------------------------------

srt_query = readRDS(paste0(dir.cross.results, "NeuMap.135055.srt.rds"))
srt_query = NormalizeData(srt_query)
srt_query = SCTransform(srt_query)


# -- Signature-gene matrix -----------------------------------------------------

mat = srt_query@assays$RNA@data[
  intersect(signaturegenes, rownames(srt_query@assays$RNA@data)), ]

saveRDS(mat, paste0(dir.cross.results, "NeuMap.signaturegenes.mat.rds"))


# -- Mapping -------------------------------------------------------------------

query = mapQuery(
  srt_query@assays$SCT@data, srt_query@meta.data, reference, vars = "sampleid",
  do_normalize = FALSE, do_umap = TRUE)


# -- Label transfer ------------------------------------------------------------

# BiocNeighbors wants cells as rows, and both embeddings come out factors-by-
# cells, so transpose each.
ref_emb = t(reference$Z_corr)
query_emb = t(query$Z)

knn_result = queryKNN(ref_emb, query_emb, k = k.neighbours, BNPARAM = AnnoyParam())

# One row per query cell, one column per neighbour, holding the neighbour's
# reference celltype
knn_labels = matrix(reference$meta_data$celltype[knn_result$index], ncol = k.neighbours)

### plain majority vote
query$cell_type_pred = apply(knn_labels, 1, function(x) {
  names(sort(table(x), decreasing = TRUE))[1]
})

### the same vote weighted by inverse distance, and its margin as a confidence
weights = 1 / (knn_result$distance + 1e-6)

query$cell_type_pred_weighted = sapply(seq_len(nrow(knn_labels)), function(i) {
  weighted_votes = tapply(weights[i, ], knn_labels[i, ], sum)
  names(which.max(weighted_votes))
})

query$cell_type_conf = sapply(seq_len(nrow(knn_labels)), function(i) {
  weighted_votes = tapply(weights[i, ], knn_labels[i, ], sum)
  max(weighted_votes) / sum(weighted_votes)
})

# The two votes should agree for all but the borderline cells
sort(table(query$cell_type_pred))
table(query$cell_type_pred_weighted)

saveRDS(query, paste0(dir.cross.results, "NeuMap.symphony.rds"))
