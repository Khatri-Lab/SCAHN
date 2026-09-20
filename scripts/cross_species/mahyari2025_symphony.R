# ==============================================================================
# mahyari2025 (NHP) mapped onto the human atlas by symphony
#
# The non-human-primate arm of the cross-species comparison, and the counterpart
# of NeuMap_orthologs.R plus NeuMap_symphony.R for the mouse. Four steps:
#   1. reprocess the published myeloid object of GSE277821 from scratch:
#      normalise, PCA, neighbours, Louvain and UMAP
#   2. export the per-cell metadata and UMAP
#   3. subset to the bone-marrow neutrophils
#   4. project them into the atlas embedding and transfer celltype labels by
#      k nearest neighbours
#
# The neutrophil subset is taken on ClusterNames_1.2, which is the authors' own
# published labelling carried in the GEO object, not the RNA_snn_res.0.8
# clustering computed above. That clustering is only exported for inspection, so
# nothing downstream of here depends on it.
#
# Expensive: a 1 GB Seurat object reprocessed end to end, then SCTransform and
# symphony on the subset. Run once, on its own.
#
# Requires Seurat, symphony, BiocNeighbors and data.table. uwot is called by
# namespace.
#
# Reads:
#   data/cross_species/GSE277821_RIRA.Myeloid.seurat.rds
#     the published myeloid object, as downloaded
#   results/symphony/SCAHN.symphony.reference.rds
#   results/symphony/SCAHN.harmony_SCT.v2000.PC25.uwot.model
#     both written by 08b_symphony_reference.R
#
# Writes:
#   results/cross_species/mahyari2025.metacell.csv
#     read back by the subset step below, and by figures/03_cross_species.R
#   results/cross_species/mahyari2025.signaturegenes.mat.rds
#   results/cross_species/mahyari2025.symphony.rds
#     both read by figures/03_cross_species.R
# ==============================================================================

library(Seurat)
library(symphony)
library(BiocNeighbors)
library(data.table)

options(future.globals.maxSize = 20000 * 1024^2)

dir.data          = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/data/"
dir.results       = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/results/"
dir.cross.data    = paste0(dir.data, "cross_species/")
dir.cross.results = paste0(dir.results, "cross_species/")
dir.symphony      = paste0(dir.results, "symphony/")
SCAHN_FIGURES = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/scripts/figures"
source(file.path(SCAHN_FIGURES, "00_genesets.R"))


# -- Parameters ----------------------------------------------------------------

nPCs = 20
reso = 0.8

### the neutrophil subset, on the authors' published ClusterNames_1.2 labels
clusters.neu = c(11, 13, 29, 23, 8, 21, 12, 17)
tissue.keep = "Bone marrow"

# Samples contributing fewer cells than this are dropped, so that per-sample
# summaries downstream are not built on a handful of cells.
min.cells.per.sample = 30

# Neighbours used for the label transfer. The vote matrix is reshaped to this
# width, so the two have to stay in step.
k.neighbours = 10

path.uwot = paste0(dir.symphony, "SCAHN.harmony_SCT.v2000.PC25.uwot.model")


# -- The atlas reference -------------------------------------------------------

reference = readRDS(paste0(dir.symphony, "SCAHN.symphony.reference.rds"))

# The uwot model is saved to its own path rather than inside the reference
# object, so attach it by hand
umap_model = uwot::load_uwot(path.uwot)
reference$umap = umap_model
reference$save_uwot_path = path.uwot


# -- Reprocessing the published object -----------------------------------------

srt.mahyari2025 = readRDS(paste0(dir.cross.data, "GSE277821_RIRA.Myeloid.seurat.rds"))

srt.mahyari2025 = NormalizeData(srt.mahyari2025)
srt.mahyari2025 = FindVariableFeatures(srt.mahyari2025, selection.method = "vst")
srt.mahyari2025 = ScaleData(object = srt.mahyari2025, features = VariableFeatures(srt.mahyari2025))
srt.mahyari2025 = RunPCA(srt.mahyari2025, features = VariableFeatures(object = srt.mahyari2025),
                         seed.use = 42)

srt.mahyari2025 = FindNeighbors(srt.mahyari2025, dims = 1:nPCs)
srt.mahyari2025 = FindClusters(srt.mahyari2025, resolution = reso)
srt.mahyari2025 = RunUMAP(srt.mahyari2025, dims = 1:nPCs, seed.use = 42, n.components = 2)


# -- Cell metadata -------------------------------------------------------------

dtf.mahyari2025 = as.data.frame(srt.mahyari2025@meta.data)
dtf.mahyari2025$cell = rownames(dtf.mahyari2025)

tmp = srt.mahyari2025@reductions$umap@cell.embeddings
colnames(tmp) = c("UMAP_1", "UMAP_2")
dtf.mahyari2025 = as.data.table(cbind(dtf.mahyari2025, tmp))

fwrite(dtf.mahyari2025, paste0(dir.cross.results, "mahyari2025.metacell.csv"))


# -- The bone-marrow neutrophils -----------------------------------------------

dtf.mahyari2025 = fread(paste0(dir.cross.results, "mahyari2025.metacell.csv"))
dtf.mahyari2025$sampleid = paste0(dtf.mahyari2025$SubjectId, ".", dtf.mahyari2025$Tissue)

dtf.mahyari2025.s = dtf.mahyari2025[ClusterNames_1.2 %in% clusters.neu & Tissue == tissue.keep]

samples.keep = names(which(table(dtf.mahyari2025.s$sampleid) >= min.cells.per.sample))
cells2include = dtf.mahyari2025.s[sampleid %in% samples.keep]$cell

srt.mahyari2025.s = subset(srt.mahyari2025, cells = cells2include)

srt.mahyari2025.s$tissue = gsub("Bone marrow", "bone marrow", srt.mahyari2025.s$Tissue)
srt.mahyari2025.s$sampleid = paste0(srt.mahyari2025.s$SubjectId, ".", srt.mahyari2025.s$Tissue)

srt.mahyari2025.s = NormalizeData(srt.mahyari2025.s)
srt.mahyari2025.s = SCTransform(srt.mahyari2025.s)


# -- Signature-gene matrix -----------------------------------------------------

mat = srt.mahyari2025.s@assays$RNA@data[
  intersect(signaturegenes, rownames(srt.mahyari2025.s@assays$RNA@data)), ]

saveRDS(mat, paste0(dir.cross.results, "mahyari2025.signaturegenes.mat.rds"))


# -- Mapping -------------------------------------------------------------------

query = mapQuery(
  srt.mahyari2025.s@assays$SCT@data, srt.mahyari2025.s@meta.data, reference, vars = "sampleid",
  do_normalize = FALSE, do_umap = TRUE)


# -- Label transfer ------------------------------------------------------------

# BiocNeighbors wants cells as rows, and both embeddings come out factors-by-
# cells, so transpose each.
ref_emb = t(reference$Z_corr)
query_emb = t(query$Z)

knn_result = queryKNN(ref_emb, query_emb, k = k.neighbours, BNPARAM = AnnoyParam())

# One row per query cell, one column per neighbour
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
sort(table(query$cell_type_pred_weighted))

saveRDS(query, paste0(dir.cross.results, "mahyari2025.symphony.rds"))
