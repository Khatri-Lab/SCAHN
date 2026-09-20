# ==============================================================================
# NeuMap: mouse to human ortholog conversion
#
# Second stage of the mouse arm. Takes the merged NeuMap object and its
# clustering from NeuMap_process.R, drops the clusters that are not neutrophils,
# and re-expresses the counts on human gene symbols so that NeuMap_symphony.R
# can map the cells onto the human atlas.
#
# Counts are projected, not renamed. A mouse gene with several human orthologs
# has its counts split evenly between them, and several mouse genes mapping to
# one human gene are summed. convert_mouse_to_human below does both in a single
# sparse multiply by a weighted projection matrix.
#
# The mouse counts are kept alongside the human ones as a second assay, RNAmice,
# so nothing is lost by the conversion.
#
# Requires Seurat, data.table and Matrix.
#
# Reads:
#   data/cross_species/mousegenes_tohuman.rds
#     the biomaRt ortholog table, built once by the commented block below
#   data/cross_species/mousegenes_tohuman.manual.csv
#     the pairs biomaRt missed, appended to that table
#   results/cross_species/NeuMap.srt.integrated.rds
#   results/cross_species/NeuMap.cellmeta.harmony_SCT.v2000.pc25.csv
#     both written by NeuMap_process.R
#
# Writes:
#   results/cross_species/NeuMap.135055.srt.rds
#     read by NeuMap_symphony.R and by figures/03_cross_species.R.
# ==============================================================================

library(Seurat)
library(data.table)
library(Matrix)

options(future.globals.maxSize = 20000 * 1024^2)

dir.data          = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/data/"
dir.results       = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/results/"
dir.cross.data    = paste0(dir.data, "cross_species/")
dir.cross.results = paste0(dir.results, "cross_species/")
SCAHN_FIGURES = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/scripts/figures"
source(file.path(SCAHN_FIGURES, "00_genesets.R"))


# -- Parameters ----------------------------------------------------------------

vargenes = "v2000"
nPCs = 25

# Clusters of harmony_snn_res.0.8 that are not neutrophils, dropped before the
# conversion. Changing this changes the cell count in the output filename.
clusters.drop = c(10, 15, 20, 22, 23, 29)


# -- The ortholog table --------------------------------------------------------

# Built once with biomaRt against Ensembl 105 and cached to
# mousegenes_tohuman.rds, so this script does not need network access:
#
#   library(biomaRt)
#   mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl",
#                      version = 105, mirror = "uswest")
#   human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",
#                      version = 105, mirror = "uswest")
#   mouse_to_human = getLDS(
#     attributes = c("mgi_symbol"), filters = "mgi_symbol",
#     values = mouse_genes, mart = mouse,
#     attributesL = c("hgnc_symbol"), martL = human, uniqueRows = TRUE)

mouse_genes_human = readRDS(paste0(dir.cross.data, "mousegenes_tohuman.rds"))
mouse_genes_human = as.data.table(mouse_genes_human)
mouse_genes_human2 = fread(
  paste0(dir.cross.data, "mousegenes_tohuman.manual.csv"), header = FALSE)

colnames(mouse_genes_human) = c("mouse", "human")
colnames(mouse_genes_human2) = colnames(mouse_genes_human)

# The manual pairs are meant to be additions, so both of these should be empty
intersect(mouse_genes_human2$mouse, mouse_genes_human$mouse)
intersect(mouse_genes_human2$human, mouse_genes_human$human)

mouse_genes_human = rbind(mouse_genes_human, mouse_genes_human2)


# -- The conversion ------------------------------------------------------------

convert_mouse_to_human = function(expr_mat, mapping) {
  # expr_mat: genes (rows) x samples (cols), rownames = mouse gene symbols
  #           supports dense matrix, dgCMatrix, or any Matrix class
  # mapping:  data.table with columns "mouse" and "human"

  is_sparse = inherits(expr_mat, "sparseMatrix")

  mapping = copy(mapping)
  mapping = mapping[mouse != "" & human != ""]
  mapping = unique(mapping)

  ### how many genes map, and how many do so ambiguously
  mapping[, n_human := uniqueN(human), by = mouse]
  mapping[, n_mouse := uniqueN(mouse), by = human]

  cat(sprintf("Input mouse genes: %d\n", nrow(expr_mat)))
  cat(sprintf("Mapped: %d\n", sum(rownames(expr_mat) %in% mapping$mouse)))
  cat(sprintf("1:1 orthologs: %d\n", mapping[n_human == 1 & n_mouse == 1, uniqueN(mouse)]))
  cat(sprintf("1 mouse -> many human: %d mouse genes\n", mapping[n_human > 1, uniqueN(mouse)]))
  cat(sprintf("many mouse -> 1 human: %d human genes\n", mapping[n_mouse > 1, uniqueN(human)]))
  cat(sprintf("Unmapped: %d\n", sum(!rownames(expr_mat) %in% mapping$mouse)))

  ### one row per mouse-human pair, weighted so each mouse gene splits its
  ### counts evenly across the human genes it maps to
  m = unique(mapping[, .(mouse, human, n_human)])
  m[, weight := 1 / n_human]
  m = m[mouse %in% rownames(expr_mat)]

  ### the projection matrix, human genes by mouse genes. P %*% expr_mat does the
  ### divide-then-sum in one sparse multiply.
  human_genes = unique(m$human)
  mouse_genes = unique(m$mouse)

  human_idx = match(m$human, human_genes)
  mouse_idx = match(m$mouse, mouse_genes)

  P = sparseMatrix(
    i = human_idx, j = mouse_idx, x = m$weight, dims = c(length(human_genes), length(mouse_genes)),
    dimnames = list(human_genes, mouse_genes))

  out = P %*% expr_mat[mouse_genes, , drop = FALSE]

  if (!is_sparse) {
    out = as.matrix(out)
  }

  cat(sprintf("Output: %d human genes x %d samples\n", nrow(out), ncol(out)))

  return(out)
}


# -- NeuMap, less the non-neutrophil clusters ----------------------------------

srt.neumap = readRDS(paste0(dir.cross.results, "NeuMap.srt.integrated.rds"))
dtf.neumap = fread(paste0(dir.cross.results, "NeuMap.cellmeta.harmony_SCT.",
                          vargenes, ".pc", nPCs, ".csv"))

dtf.neumap$cell = rownames(srt.neumap@meta.data)

cells2keep = dtf.neumap[!harmony_snn_res.0.8 %in% clusters.drop]$cell
srt.neumap.s = subset(srt.neumap, cells = cells2keep)


# -- Human-symbol object -------------------------------------------------------

mat.human = convert_mouse_to_human(srt.neumap.s@assays$RNA@counts, mouse_genes_human)

# Which signature genes have no mouse ortholog in the table
setdiff(signaturegenes, rownames(mat.human))

srt_neumap = CreateSeuratObject(counts = mat.human, meta.data = srt.neumap.s@meta.data)

# Keep the mouse counts too, so the conversion loses nothing
srt_neumap[["RNAmice"]] = CreateAssayObject(
  counts = srt.neumap.s@assays$RNA@counts)

DefaultAssay(srt_neumap) = "RNA"
srt_neumap = NormalizeData(srt_neumap)

saveRDS(srt_neumap, paste0(dir.cross.results, "NeuMap.135055.srt.rds"))
