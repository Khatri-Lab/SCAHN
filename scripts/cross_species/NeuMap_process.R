# ==============================================================================
# NeuMap: sample selection, QC, integration, harmony and clustering
#
# First stage of the mouse arm of the cross-species comparison, from the GEO
# downloads to the per-cell metadata NeuMap_orthologs.R takes on. Six steps:
#   1. read the fifteen per-sample Seurat objects of GSE266680
#   2. keep the nine that make up the NeuMap reference, and filter cells
#   3. normalise, then SCTransform each sample
#   4. integration features at seven sizes, then merge into one object
#   5. harmony on the SCT PCA, once per feature set
#   6. UMAP and Louvain at two resolutions, and export the cell metadata
#
# Six of the nine samples are taken whole. The other three come from knockout
# experiments, so only their wild-type arm is kept, identified by a _WT suffix
# on Sample_Name. The six remaining accessions are other experimental arms and
# are read but not used.
#
# Expensive: fifteen objects of 40 to 220 MB each, an SCTransform per sample, an
# 8.7 GB merged object, then harmony and clustering over 147,103 cells. Run
# once, on its own.
#
# Requires Seurat, data.table and harmony. SCTransform also needs glmGamPoi,
# which it loads by namespace rather than attaching.
#
# Reads:
#   data/cross_species/NeuMap/GSM8248*.rds
#     the fifteen per-sample objects as downloaded from GEO.
#
# Writes:
#   results/cross_species/NeuMap.srt.integrated.rds
#     read by NeuMap_orthologs.R
#   results/cross_species/NeuMap.SCT.var.features.list.rds
#   results/cross_species/NeuMap.harmony_SCT.{v2000,v3000}.embeddings.csv
#     written headerless, and read back by the clustering step below
#   results/cross_species/NeuMap.cellmeta.harmony_SCT.{v2000,v3000}.pc25.csv
#     the v2000 file is read by NeuMap_orthologs.R, which drops six clusters of
#     harmony_snn_res.0.8 before mapping
# ==============================================================================

library(Seurat)
library(data.table)
library(harmony)

options(future.globals.maxSize = 20000 * 1024^2)

dir.data          = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/data/"
dir.results       = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/results/"
dir.cross.data    = paste0(dir.data, "cross_species/")
dir.cross.results = paste0(dir.results, "cross_species/")
dir.neumap        = paste0(dir.cross.data, "NeuMap/")


# -- Parameters ----------------------------------------------------------------

### samples kept whole
samples.whole = c("GSM8248766_yo2_raw", "GSM8248765_yo1_raw",
                  "GSM8248762_pdac_raw", "GSM8248760_non_steady1_raw",
                  "GSM8248761_non_steady2_raw", "GSM8248759_lps_raw")

### knockout experiments, kept down to their wild-type arm
samples.wt = c("GSM8248752_cxcr4_raw", "GSM8248757_jb1_raw", "GSM8248758_jb2_raw")

### per-cell QC thresholds
min.counts   = 100
min.features = 100
max.mt       = 20

### integration feature-set sizes, and the subset that harmony is run on
nfeatures.sct = c(3000, 2500, 2000, 1500, 1000, 800, 500)
vargenes.harmony = c("v2000", "v3000")
nPCs.harmony = 25


# -- The per-sample objects ----------------------------------------------------

files = list.files(path = dir.neumap, pattern = "GSM8248", full.names = FALSE)

srtlist = list()
for (file in files) {
  pre = gsub("\\.rds$", "", file)
  srtlist[[pre]] = readRDS(paste0(dir.neumap, file))
}


# -- Sample selection ----------------------------------------------------------

srtlist.neumap = srtlist[samples.whole]

for (i in samples.wt) {
  obj = srtlist[[i]]
  cells_wt = colnames(obj)[grepl("_WT$", obj$Sample_Name)]
  srtlist.neumap[[i]] = subset(obj, cells = cells_wt)
}

length(unlist(sapply(srtlist.neumap, function(x) table(x$Sample_Name))))
sum(unlist(sapply(srtlist.neumap, function(x) table(x$Sample_Name))))

sapply(srtlist.neumap, function(x) summary(x$nCount_RNA))
sapply(srtlist.neumap, function(x) summary(x$nFeature_RNA))
sapply(srtlist.neumap, function(x) summary(x$mt_genes))


# -- Per-cell QC ---------------------------------------------------------------

srtlist.neumap.sub = list()
for (name in names(srtlist.neumap)) {
  srt = srtlist.neumap[[name]]
  metadata = as.data.frame(srt@meta.data)
  metadata$cell = rownames(metadata)
  metadata$sampleid = paste0(name, ".", metadata$Sample_Name)
  cells2keep = as.data.table(metadata)[nCount_RNA > min.counts & nFeature_RNA > min.features &
                                       mt_genes < max.mt]$cell
  srt@meta.data = metadata
  srt = subset(srt, cells = cells2keep)
  srtlist.neumap.sub[[name]] = srt
}


# -- Normalisation and variable features, RNA ----------------------------------

srtlist.neumap.sub = lapply(X = srtlist.neumap.sub, FUN = function(x) {
  DefaultAssay(x) = "RNA"
  x = NormalizeData(x)
  x = FindVariableFeatures(x, selection.method = "vst", nfeatures = 4000)
})


# -- SCTransform ---------------------------------------------------------------

srtlist.neumap.sub = lapply(
  X = srtlist.neumap.sub, FUN = SCTransform, method = "glmGamPoi", seed.use = 42,
  assay = "RNA", variable.features.n = 3000, return.only.var.genes = TRUE)


# -- Integration features, SCT -------------------------------------------------

var.features.list = setNames(
  lapply(nfeatures.sct, function(n) {
    SelectIntegrationFeatures(object.list = srtlist.neumap.sub, nfeatures = n)
  }), paste0("v", nfeatures.sct))

saveRDS(var.features.list, paste0(dir.cross.results, "NeuMap.SCT.var.features.list.rds"))


# -- Merge ---------------------------------------------------------------------

srt.sct = merge(x = srtlist.neumap.sub[[1]], y = srtlist.neumap.sub[2:length(srtlist.neumap.sub)],
                merge.data = TRUE)

saveRDS(srt.sct, paste0(dir.cross.results, "NeuMap.srt.integrated.rds"))


# -- Harmony -------------------------------------------------------------------

for (name in vargenes.harmony) {
  VariableFeatures(srt.sct) = var.features.list[[name]]
  DefaultAssay(srt.sct) = "SCT"
  srt.sct = RunPCA(srt.sct, features = var.features.list[[name]], seed.use = 42)

  srt.sct = RunHarmony(srt.sct, group.by.vars = "sampleid",
                       max_iter = 10, reduction.save = "harmony")

  fwrite(srt.sct@reductions$harmony@cell.embeddings,
         paste0(dir.cross.results, "NeuMap.harmony_SCT.", name, ".embeddings.csv"),
         col.names = FALSE)
}


# -- UMAP, clustering and cell metadata ----------------------------------------

for (name in vargenes.harmony) {
  for (nPC in nPCs.harmony) {
    cell.embeddings = fread(
      paste0(dir.cross.results, "NeuMap.harmony_SCT.", name, ".embeddings.csv"), header = FALSE)
    cell.embeddings = as.matrix(cell.embeddings)
    rownames(cell.embeddings) = colnames(srt.sct)

    harmony_dr = CreateDimReducObject(embeddings = cell.embeddings, key = "harmony_", assay = "SCT")
    srt.sct@reductions$harmony = harmony_dr

    srt.sct = RunUMAP(srt.sct, reduction = "harmony", reduction.name = "harmony.umap",
                      dims = 1:nPC, seed.use = 42, verbose = TRUE, n.components = 2)

    srt.sct = FindNeighbors(srt.sct, reduction = "harmony", dims = 1:nPC, verbose = FALSE,
                            graph.name = c("harmony_nn", "harmony_snn"))

    srt.sct = FindClusters(srt.sct, resolution = 0.8, graph.name = "harmony_snn")
    srt.sct = FindClusters(srt.sct, resolution = 1.2, graph.name = "harmony_snn")

    umapdf = as.data.table(
      srt.sct@reductions$harmony.umap@cell.embeddings[, 1:2])
    colnames(umapdf) = c("UMAP_1", "UMAP_2")

    umapdf$cell = srt.sct$cell
    umapdf$sampleid = srt.sct$sampleid
    umapdf$harmony_snn_res.0.8 = srt.sct$harmony_snn_res.0.8
    umapdf$harmony_snn_res.1.2 = srt.sct$harmony_snn_res.1.2

    fwrite(umapdf, paste0(dir.cross.results, "NeuMap.cellmeta.harmony_SCT.",
                          name, ".pc", nPC, ".csv"))
  }
}
