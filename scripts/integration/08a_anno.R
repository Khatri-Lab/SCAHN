# ==============================================================================
# Subset annotation: marker dotplots, cluster UMAPs, labels and proportions
#
# Turns the consensus clusterings into the published neutrophil subset labels.
# Subsets a per-cell matrix to the selected genes, draws the dotplots and
# cluster UMAPs the resolution sweep was read from, cuts the 24 subset labels
# at res0.4, refines them cell by cell on marker expression, and writes the
# per-cell metadata and per-sample subset proportions the rest of the tree
# runs on.
#
# This script and 08b_symphony_reference.R each read the other's output: 08b
# labels the reference with the celltype column written here, and the UMAP
# coordinates used here come from the embedding 08b writes. They were run
# alternately -- an earlier embedding placed the cells while the labels were
# being cut, then 08b was re-run against the final labels -- so neither
# reproduces the published files on its own. The a/b numbering records that:
# they are one step, not two.
#
# Requires 00_setup.R (packages, paths, textsize) and Seurat
#
# Reads:
#   data/SCAHN/SCAHN.srt.integrated.rds  written by 01_run_integration.R
#   data/selectedgenes.rds
#   data/SCAHN.simplemeta.csv            written by 02a_run_harmony_SCT.R
#   data/samplemeta.csv
#   data/cellnumber.rds                  cells per sample over all cell types,
#                                        the denominator for neu.pct
#   results/symphony/SCAHN.harmony_SCT.v2000.PC25.uwot.emb.csv
#     written by 08b_symphony_reference.R
#   results/integration/clustering/consensus_results.<reso>.csv
#     written by 06b_run_consensuscluster.py
#
# Writes, all into results/annotation/:
#   SCAHN.selectedgenes.mat.expr.integrated.rds
#     read by figures/02_atlas.R and figures/03_cross_species.R
#   SCAHN.selectedgenes.mat.count.integrated.rds  read by this script only
#   SCAHN.cellmeta.csv
#     read by figures/01_load_data.R, 08b_symphony_reference.R,
#     findmarkers.R, scores_and_pseudobulk.R and slingshot.R
#   SCAHN.subsetproportion.rds
#     read by figures/01_load_data.R and the three meta_sc_*.R
#   SCAHN.concensus.<reso>.ht_markers.pdf
#   SCAHN.UMAP.pdf
# ==============================================================================

library(Seurat)

SCAHN_SCRIPTS = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/scripts/figures"
source(file.path(SCAHN_SCRIPTS, "00_setup.R"))     # packages, paths, textsize

dir.anno       = paste0(dir.results, "annotation/")
dir.clustering = paste0(dir.results, "integration/clustering/")

resolutions = c(0.1, 0.2, 0.4, 0.6, 0.8)

marker_genes = c("CST3", "LYZ", "CD14", "FCGR3A", "FCGR3B", "FCGR1A", "CD274", "S100A8", "CEACAM8",
                 "DEFA4", "DEFA3", "MARCO", "FABP4", "MCEMP1", "APOE", "APOC1", "NUPR1", "MRC1",
                 "C1QA", "SMIM25", "CSF1R", "CSF3R", "NAMPT", "UGCG", "SERPINF1", "CD1C", "FCER1A",
                 "CLEC4C", "PPBP", "PF4", "CD34", "SPINK2", "HBA1", "CD3E", "IL7R", "LTB", "CD8A",
                 "CD4", "CCR7", "SELL", "PTPRC", "ITGB1", "FOXP3", "KLRF1", "NKG7", "PRF1", "GNLY",
                 "GZMB", "MKI67", "STMN1", "IGHG1", "CD38", "CD27", "IGHM", "MS4A1", "CD79A",
                 "EPCAM", "KRT19", "CDH1", "DCN", "THY1", "COL1A1", "PECAM1", "CLDN5", "RAMP2",
                 "MS4A2", "GATA2", "CEACAM6", "TNFSF10", "TNFSF13B", "TNFSF14", "LGALS9", "CXCL10",
                 "GRN", "TNF", "CCL5", "CCL20", "IL1B", "ICAM1", "SPP1", "VEGFA", "VEGFB", "ITGAL",
                 "IL3RA", "OLIG1", "MOG", "CLDN11", "KIT", "SIGLEC5", "SIGLEC8", "ITGAM", "FLT1",
                 "KDR", "LTF", "CXCR2", "ALPL", "ANXA3", "CXCR1", "CYP4F3", "KCNJ15", "MME", "MMP9",
                 "PI3", "SLPI", "TNFRSF10C", "AZU1", "BPI", "CAMP", "CRISP3", "ELANE", "LCN2",
                 "MMP8", "MPO", "CD74", "HLA-DRA", "CCR1", "CCR2", "CXCR4", "CD99", "CD44", "CD47",
                 "CXCL8", "CCL3", "CD101", "CD36", "GBP2", "MX2", "ISG15", "CD177", "ARG1",
                 "S100A9", "OLFM4", "CTSG", "CD63", "LAMP1", "RETN", "PRTN3", "EGR1", "FOS",
                 "NFKBIA", "MT-CO2", "RPS18", "PTMA", "EPX", "CXCL2", "IGKC", "IGLC2", "IGHA1",
                 "IGLC3", "JCHAIN")


# -- Selected-gene matrices ----------------------------------------------------

srt = readRDS(paste0(dir.data, "SCAHN/SCAHN.srt.integrated.rds"))

genes = readRDS(paste0(dir.data, "selectedgenes.rds"))
setdiff(genes, rownames(srt@assays$RNA@data))
genes = intersect(genes, rownames(srt@assays$RNA@data))

saveRDS(srt@assays$RNA@data[genes, ],
        paste0(dir.anno, "SCAHN.selectedgenes.mat.expr.integrated.rds"))
saveRDS(srt@assays$RNA@counts[genes, ],
        paste0(dir.anno, "SCAHN.selectedgenes.mat.count.integrated.rds"))


# -- Per-cell table for the dotplots and UMAPs ---------------------------------

mat = readRDS(paste0(dir.anno, "SCAHN.selectedgenes.mat.count.integrated.rds"))

dtfplot = fread(paste0(dir.data, "SCAHN.simplemeta.csv"))
umap = fread(paste0(dir.results,
                    "symphony/SCAHN.harmony_SCT.v2000.PC25.uwot.emb.csv"))
mat = mat[, dtfplot$cell.pid]

dtfplot$UMAP_1 = umap$V1
dtfplot$UMAP_2 = umap$V2
dtfplot = as.data.table(cbind(dtfplot, umap))

dtfplot = cbind(dtfplot, as.data.table(t(mat)))


# -- Marker dotplots across the resolution sweep -------------------------------

for (reso in resolutions) {
  dfcluster = fread(paste0(dir.clustering, "consensus_results.", reso, ".csv"))
  dtfplot$cluster = dfcluster$consensus_label

  # one row of dots per cluster, so the canvas grows with the clustering
  fig_height = 2 + length(unique(dtfplot$cluster)) * 0.15

  p = plot_dotplot(dtfplot, features = intersect(marker_genes, colnames(dtfplot)),
                   group_by = "cluster", max_scale = 2, dot_scale = 4, title = NULL)
  ggsave(paste0(dir.anno, "SCAHN.concensus.", reso, ".ht_markers.pdf"),
         p, width = 24, height = fig_height)
}


# -- Cluster UMAPs -------------------------------------------------------------

p.umap.clusters = list()

for (reso in resolutions) {
  dfcluster = fread(paste0(dir.clustering, "consensus_results.", reso, ".csv"))
  dtfplot$cluster = dfcluster$consensus_label

  df.label = as.data.table(dtfplot)[, lapply(.SD, median), by = cluster,
                                    .SDcols = c("UMAP_1", "UMAP_2")]
  colnames(df.label)[2:3] = c("x", "y")
  df.label$labeltext = df.label[, "cluster", with = FALSE]

  clusters = sort(unique(dtfplot$cluster))
  col_cluster = expresso_colors("discrete_50", n = length(clusters))
  names(col_cluster) = clusters
  dtfplot$facetvar = as.character(reso)

  p.umap.clusters[[as.character(reso)]] =
    plot_scatter(dtfplot, x = "UMAP_1", y = "UMAP_2",
                 color_by = "cluster", color_type = "discrete", colors = col_cluster,
                 facet_by = "facetvar", facet_nrow = 5,
                 point_size = 0.001, alpha = 0.2, raster_dpi = 300,
                 shuffle = TRUE, seed = 13, label_df = df.label, title = NULL) +
    ggplot2::theme(
      text              = ggplot2::element_text(size = textsize),
      legend.key.width  = ggplot2::unit(0.2, "line"),
      legend.key.height = ggplot2::unit(0.2, "line")
    )
}

pdf(file = paste0(dir.anno, "SCAHN.UMAP.pdf"), width = 12, height = 8)
grid.arrange(grobs = p.umap.clusters, nrow = 2, as.table = TRUE, newpage = FALSE)
dev.off()


# -- Subset annotation ---------------------------------------------------------

# res0.1 is not carried: it is too coarse to be worth a column.
for (reso in c(0.8, 0.6, 0.4, 0.2)) {
  dfcluster = fread(paste0(dir.clustering, "consensus_results.", reso, ".csv"))
  dtfplot[[paste0("res", reso)]] = dfcluster$consensus_label
}

### Initial assignment, res0.4 cluster to subset

clustersanno = list(
  `AP-1` = c("4", "31"), AZU1 = c("6", "30"), `CCL3/4` = c("26", "41", "43"),
  CD74 = "10", CXCL = "36", EGR1 = "19", G0S2 = c("14", "24"), HSP = "29",
  IFN1 = "23", IFN2 = "13", IFN3 = c("1", "37", "22"), IL1B = "32",
  IL1R2 = "27", IL1RN = "34", lowDepth = c("15", "20"),
  LTF = c("8", "38", "33"), MME = c("2", "3", "28", "40"),
  MMP9 = c("0", "11", "17", "25", "35"), `NF-κB` = "5", PTGS2 = "16",
  S100A4 = "7", SLPI = c("18", "39"), TXNIP = c("9", "21", "44"),
  VEGFA = c("12", "42"))

setdiff(unique(dtfplot$res0.4), sort(unlist(clustersanno))) # character(0)

annocluster = list()
for (i in names(clustersanno)) {
  for (j in as.numeric(unlist(clustersanno[i]))) annocluster[[as.character(j)]] = i
}
dtfplot$celltype = unlist(annocluster)[as.character(dtfplot$res0.4)]

### Refinement on marker expression

# Granule and proliferation panel. A cell in cluster 6 or 30 with none of it
# detected has too little depth to place, rather than belonging to a subset.
granule = c("AZU1", "MPO", "ELANE", "CTSG", "PRTN3", "DEFA3", "DEFA4", "CEACAM8",
            "OLFM4", "LTF", "LCN2", "CAMP", "MMP8", "BPI", "SLPI", "STMN1", "MKI67")
dtfplot[, granule.sum := rowSums(.SD), .SDcols = granule]

dtfplot[, celltype := ifelse(res0.4 %in% c(6, 30) & granule.sum == 0,
                             "lowDepth", celltype)]
dtfplot[, celltype := ifelse(res0.4 %in% c(6, 30) & celltype %in% c("lowDepth") &
                               CD74 + `HLA-DRA` > 0, "CD74", celltype)]
dtfplot[, celltype := ifelse(res0.4 %in% c(6, 30) & celltype %in% c("lowDepth") &
                               DNAJB1 + HSPA1B > 0, "HSP", celltype)]

dtfplot[, celltype := ifelse(res0.4 %in% c(15) & AZU1 + MPO + ELANE + CTSG + PRTN3 > 0,
                             "AZU1", celltype)]
dtfplot[, celltype := ifelse(res0.4 %in% c(15) & CD74 + `HLA-DRA` > 0, "CD74", celltype)]

# The same two tests again without the lowDepth guard, which is not redundant:
# a cell positive for both panels was sent to CD74 above and is moved to HSP
# here, because this pair runs CD74 first and HSP second on the full set.
dtfplot[, celltype := ifelse(res0.4 %in% c(6, 30) & granule.sum == 0 &
                               CD74 + `HLA-DRA` > 0, "CD74", celltype)]
dtfplot[, celltype := ifelse(res0.4 %in% c(6, 30) & granule.sum == 0 &
                               DNAJB1 + HSPA1B > 0, "HSP", celltype)]

# AZU1 / LTF and the MMP9 / MME / LTF boundaries, each resolved on the markers
# that define the two subsets rather than on the cluster the cell landed in.
dtfplot[, celltype := ifelse(celltype %in% c("LTF") & LTF + LCN2 + CAMP + MMP8 == 0 &
                               AZU1 + MPO + ELANE + CTSG + PRTN3 > 0, "AZU1", celltype)]
dtfplot[, celltype := ifelse(celltype %in% c("AZU1") & LTF + LCN2 + CAMP + MMP8 > 0 &
                               AZU1 + MPO + ELANE + CTSG + PRTN3 == 0, "LTF", celltype)]

dtfplot[, celltype := ifelse(celltype %in% c("MMP9") & MMP9 == 0 &
                               LTF + LCN2 + CAMP + MMP8 + CEACAM8 > 0, "LTF", celltype)]
dtfplot[, celltype := ifelse(celltype %in% c("LTF") & MMP9 + MME > 0 &
                               LTF + LCN2 + CAMP + MMP8 + CEACAM8 == 0, "MMP9", celltype)]
dtfplot[, celltype := ifelse(celltype %in% c("MMP9") & MMP9 == 0 & MME > 0,
                             "MME", celltype)]
dtfplot[, celltype := ifelse(celltype %in% c("MME") & MMP9 > 2, "MMP9", celltype)]

dtfplot[, celltype := ifelse(celltype %in% c("VEGFA") & TXNIP > 0 &
                               VEGFA + CD83 + PPIF + SQSTM1 + C15orf48 + BHLHE40 +
                               CCRL2 + TNFAIP3 + TNFAIP6 + NFKBIA + IER3 == 0,
                             "TXNIP", celltype)]

dtfplot[, granule.sum := NULL]

fwrite(dtfplot[, c("UMAP_1", "UMAP_2", "pid", "sampleid.pid", "cell.pid",
                   "nCount_RNA", "nFeature_RNA", "percent.ribo", "percent.mt",
                   "res0.8", "res0.6", "res0.4", "res0.2", "celltype")],
       paste0(dir.anno, "SCAHN.cellmeta.csv"))


# -- Subset proportions --------------------------------------------------------

samplemeta = fread(paste0(dir.data, "samplemeta.csv"))
celltotalnumber = readRDS(paste0(dir.data, "cellnumber.rds"))

dtf.prop = cell_proportions(dtfplot, sample_col = "sampleid.pid",
                            cluster_col = "celltype")
dtf.prop[, number.allcell := ifelse(sampleid.pid %in% names(celltotalnumber),
                                    celltotalnumber[sampleid.pid], NA)]
dtf.prop[is.na(number.allcell)][, .N] # 159

metacols = c("pid", "sampleid", "subjectid", "group", "group.detail", "tissue",
             "tissue.detail")
samplemeta.m = samplemeta[match(dtf.prop$sampleid.pid, samplemeta$sampleid.pid)]
dtf.prop = cbind(dtf.prop, samplemeta.m[, metacols, with = FALSE])

table(dtf.prop[is.na(number.allcell)]$pid)
table(dtf.prop[tissue == "peripheral blood"]$pid)

# neutrophils as a fraction of all cells in the sample, not of the atlas
dtf.prop$neu.pct = dtf.prop$n_cells / dtf.prop$number.allcell

dtf.prop$source = samplemeta.m$source
dtf.prop$technique = samplemeta.m$technique
dtf.prop[, technique := ifelse(grepl("Chromium", technique), "10x Chromium", technique)]

# display form of source: drop the RBC-lysis qualifier and break the line
dtf.prop$source2 = gsub("_RBC.*", "", dtf.prop$source)
dtf.prop$source2 = gsub("fresh", "fresh\n", dtf.prop$source2)
dtf.prop$source2 = gsub("frozen", "frozen\n", dtf.prop$source2)
dtf.prop$source2 = gsub("WB", "whole blood", dtf.prop$source2)
dtf.prop$source2 = factor(dtf.prop$source2,
                          levels = c("frozen\nPBMC", "fresh\nPBMC", "fresh\nwhole blood"))

dtf.prop[, lapply(.SD, median), by = "tissue", .SDcols = "neu.pct"]
dtf.prop[, lapply(.SD, median), by = "source2", .SDcols = "neu.pct"]

saveRDS(dtf.prop, paste0(dir.anno, "SCAHN.subsetproportion.rds"))
