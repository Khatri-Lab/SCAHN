# ==============================================================================
# Marker genes per neutrophil subset and per subset category
#
# Seurat FindMarkers on the RNA assay of the integrated object, one identity
# against all remaining cells.
#
# Expensive: loads the full integrated object and fits one test per identity.
# Run once, on its own.
#
# Reads:
#   data/SCAHN/SCAHN.srt.integrated.rds
#   results/annotation/SCAHN.cellmeta.csv
#
# Writes:
#   results/annotation/SCAHN.subset.markers.csv
#     read by 04_infection_cancer.R, 05_sexdifference.R
#   results/annotation/SCAHN.subset_categ.markers.csv
#     read by scores_and_pseudobulk.R
# ==============================================================================

library(data.table)
library(Seurat)
library(future)

dir.data    = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/data/"
dir.results = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/results/"

# Each worker gets its own copy of the integrated object, so the 500 MB default
# globals limit is nowhere near enough.
options(future.globals.maxSize = 20000 * 1024^2)
plan("multicore", workers = 8)


# -- Integrated object and matching cell metadata ------------------------------

srt.sct = readRDS(paste0(dir.data, "SCAHN/SCAHN.srt.integrated.rds"))
dtf = fread(paste0(dir.results, "annotation/SCAHN.cellmeta.csv"))

# The identities come from dtf and the expression from srt.sct, so the two must
# be in the same cell order.
all.equal(dtf$cell.pid, as.character(srt.sct$cell.pid))

dtf[, celltype2 := ifelse(celltype %in% c("AZU1", "LTF"), "immature", "other")]
dtf[, celltype2 := ifelse(celltype %in% c("IFN1", "IFN2", "IFN3"), "IFN", celltype2)]

DefaultAssay(srt.sct) = "RNA"


# -- Markers per subset --------------------------------------------------------

srt.sct$cluster = dtf$celltype
Idents(srt.sct) = "cluster"

markers = NULL
for (i in unique(srt.sct$cluster)) {
  marker = FindMarkers(object = srt.sct, ident.1 = i, min.pct = 0.1, logfc.threshold = 0.25)
  marker$cluster = i
  marker$gene = rownames(marker)
  markers = rbind(markers, marker)
}
markers = as.data.table(markers)
fwrite(markers, paste0(dir.results, "annotation/SCAHN.subset.markers.csv"))


# -- Markers per subset category -----------------------------------------------

srt.sct$cluster = dtf$celltype2
Idents(srt.sct) = "cluster"

markers = NULL
for (i in c("IFN", "immature")) {
  marker = FindMarkers(object = srt.sct, ident.1 = i, min.pct = 0.1, logfc.threshold = 0.25)
  marker$cluster = i
  marker$gene = rownames(marker)
  markers = rbind(markers, marker)
}
markers = as.data.table(markers)
fwrite(markers, paste0(dir.results, "annotation/SCAHN.subset_categ.markers.csv"))
