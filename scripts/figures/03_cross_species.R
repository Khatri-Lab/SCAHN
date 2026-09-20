# ==============================================================================
# Cross-species comparison
#
# Neutrophil subsets and marker expression compared across human (SCAHN), mouse
# (NeuMap, grieshaber-bouyer2021, xie2020) and non-human primate (mahyari2025,
# staupe2022) datasets: percent-expressing comparisons between species, NeuMap
# hub-gene and published mouse subset scores projected onto SCAHN, and symphony
# projections of mouse and NHP cells onto the SCAHN reference.
#
# Requires 00_setup.R (packages, figure settings, palettes) and 01_load_data.R
# (dtf, signaturegenes, col_celltype, col_tissue, textsize / textcolor,
# dir.data / dir.results / dir.fig).
#
# Reads, on top of what 01_load_data.R loads:
#   data/cross_species/mousegenes_tohuman.rds
#   data/cross_species/mousegenes_tohuman.manual.csv
#   data/cross_species/GSE277821_RIRA.Myeloid.seurat.rds
#   data/cross_species/NeuMap.samplemeta.csv
#   data/cross_species/human.PNG, data/cross_species/mouse.PNG,
#     data/cross_species/monkey.PNG
#   results/annotation/SCAHN.selectedgenes.mat.expr.integrated.rds
#   results/cross_species/NeuMap.135055.srt.rds
#   results/cross_species/dtf.neumap.csv     written by this script first, see below
#   results/cross_species/NeuMap.signaturegenes.mat.rds
#   results/cross_species/NeuMap.symphony.rds
#   results/cross_species/SCAHN.neumap.scores.rds
#   results/cross_species/srt.GSE165276.rds
#   results/cross_species/GSE165276.metacell.csv
#   results/cross_species/srt.xie2020mice.rds
#   results/cross_species/xie2020mice.metacell.csv
#   results/cross_species/srt.staupe2022.rds
#   results/cross_species/staupe2022.metacell.csv
#   results/cross_species/mahyari2025.metacell.csv
#   results/cross_species/mahyari2025.signaturegenes.mat.rds
#   results/cross_species/mahyari2025.symphony.rds
#
# Writes to results/cross_species/:
#   dtf.expr.pct.crossspecies.csv            per donor, human + mouse
#   dtf.expr.pct.sample.crossspecies.csv     per sample, human + mouse
#   dtf.expr.pct.sample.crossspeciesNHP.csv  per sample, human + NHP + mouse
#   dtf.neumap.csv                           NeuMap cell metadata after mapping,
#                                            written in the symphony section and
#                                            re-read by the percent-expressing one
#
# Writes to figures/original/:
#   Figure3.pdf   Figure 3
#   Fig.e4.pdf    Extended figure 4 (published mouse subset scores, NeuMap
#                 marker dotplot)
# ==============================================================================

SCAHN_SCRIPTS = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/scripts/figures"
source(file.path(SCAHN_SCRIPTS, "00_setup.R"))
source(file.path(SCAHN_SCRIPTS, "01_load_data.R"))

# "lowDepth" -> "Low depth" as in 02_atlas.R, here on the palette keys and on the
# two symphony predicted-label columns built below -- between them, every place
# this script reads a subset name from. dtf is left alone: nothing here draws its
# celltype column, it only rides along on hubgenes.
names(col_celltype)[names(col_celltype) == "lowDepth"] = "Low depth"


# -- Mouse-to-human gene mapping -----------------------------------------------

# The ortholog table was built once with biomaRt (Ensembl 105) and cached; the
# manual file adds pairs biomaRt missed. The conversion itself is done upstream:
# convert_mouse_to_human() in NeuMap_orthologs.R for the NeuMap object, while the
# xie2020, GSE165276 and staupe2022 objects arrive already converted and have no
# producer in this tree. Either way each mouse object carries human symbols in
# its RNA assay and the original mouse symbols in RNAmice; note RNAmice is the
# default assay, so the code below indexes @assays$RNA explicitly. The mapping is
# kept here only to label mouse gene names on the NeuMap gene UMAPs.
#
# library(biomaRt)
# mouse <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl",
#                     version = 105, mirror = "uswest")
# human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",
#                     version = 105, mirror = "uswest")
# mouse_to_human <- getLDS(
#     attributes = c("mgi_symbol"),
#     filters = "mgi_symbol",
#     values = mouse_genes,  # your vector of mouse gene names
#     mart = mouse,
#     attributesL = c("hgnc_symbol"),
#     martL = human,
#     uniqueRows = TRUE
# )

mouse_genes_human = readRDS(paste0(dir.data, "cross_species/mousegenes_tohuman.rds"))
mouse_genes_human = as.data.table(mouse_genes_human)
mouse_genes_human2 = fread(
  paste0(dir.data, "cross_species/mousegenes_tohuman.manual.csv"), header = FALSE)
colnames(mouse_genes_human) = c("mouse", "human")
colnames(mouse_genes_human2) = colnames(mouse_genes_human)

intersect(mouse_genes_human2$mouse, mouse_genes_human$mouse)
intersect(mouse_genes_human2$human, mouse_genes_human$human)

mouse_genes_human = rbind(mouse_genes_human, mouse_genes_human2)


# -- Signature gene panel ------------------------------------------------------

library(data.table)
library(Matrix)
library(Seurat)

signaturegenes = c("STMN1", "MKI67", "FUT4", "PRTN3", "CTSG", "ELANE",
                   "AZU1", "MPO",  "DEFA4", "DEFA3", "BPI", "CEACAM6", "CEACAM8",
                   "OLFM4", "LTF", "LCN2", "MMP8", "CAMP", "CRISP3", "TCN1", "HP", "FCN1", "PADI4",
                   "PGLYRP1", "ANXA3", "ARG1", "CD177", "MMP9",
                   "NQO2", "CYP4F3", "S100A8", "S100A9", "S100A12", "S100P", "TXN", "TSPO",
                   "PLAC8", "GPI", "OLR1",  "ORM1", "IL1R2", "S100A4",
                   "EGR1", "FOS", "FOSB", "JUN", "PTGS2", "MT2A",  "HES4",
                   "LY6E", "ISG15", "IFI44", "OASL", "IFI6", "OAS3",
                   "RSAD2", "HERC5", "MX1", "IFIT1", "IFIT2",  "IFIT5",
                   "DDX58", "XAF1", "IRF7", "EPSTI1", "APOL6", "GBP5", "GBP4", "GBP2",
                   "GBP1", "PARP14", "TNFSF10", "IFITM3", "CD274", "CXCL1", "CXCL2", "CXCL8",
                   "PPIF", "SPP1", "CDKN1A", "CD83", "SQSTM1", "PLAU", "CXCR4",
                   "CSTB", "LGALS3", "HMOX1", "VEGFA", "CCL3", "CCL4",
                   "TNFAIP3", "NFKBIA", "IL1B", "IL1RN", "CD69",
                   "SLPI", "PI3", "SIGLEC10", "TXNIP","G0S2",
                   "HSPA1A", "HSPA1B", "DNAJB1", "CD74", "HLA-DRA", "ADM",
                   "AQP9", "BCL6", "CR1", "APOBEC3A", "PTEN", "SELL",
                   "MME", "TNFRSF10C", "ALPL", "CMTM2", "DGAT2", "KCNJ15", "LRG1", "MGAM",
                   "CSF3R", "NAMPT", "FCGR3B", "CXCR2")


# -- NeuMap (mouse) projected onto SCAHN by symphony ---------------------------

query = readRDS(paste0(dir.results, "cross_species/NeuMap.symphony.rds"))
samplemeta.neumap = fread(paste0(dir.data, "cross_species/NeuMap.samplemeta.csv"))
sort(samplemeta.neumap$Raw_Neutrophil_counts - samplemeta.neumap$Neutrophil_counts)
samplemeta.neumap$tissue = gsub("blood", "peripheral blood", samplemeta.neumap$tissue)
samplemeta.neumap$tissue = gsub("bm", "bone marrow", samplemeta.neumap$tissue)
samplemeta.neumap$group = tolower(samplemeta.neumap$Disease)
samplemeta.neumap$group = gsub("Healthy", "Healthy", samplemeta.neumap$group)

dtf.neumap = as.data.table(query$umap)
colnames(dtf.neumap) = c("UMAP_1", "UMAP_2")
dtf.neumap$celltype = query$cell_type_pred
dtf.neumap$celltype_conf = query$cell_type_conf
mean(dtf.neumap$celltype_conf > 0.5)
which(table(dtf.neumap$celltype) <= 10)
table(dtf.neumap$celltype)

dtf.neumap$sampleid = query$meta_data$sampleid
dtf.neumap$nCount_RNA = query$meta_data$nCount_RNA
dtf.neumap$nFeature_RNA = query$meta_data$nFeature_RNA
dim(dtf.neumap)
dtf.neumap = cbind(
  dtf.neumap, samplemeta.neumap[match(dtf.neumap$sampleid, samplemeta.neumap$sampleid)][
    , c("tissue", "sex", "age", "group")])

fwrite(dtf.neumap, paste0(dir.results, "cross_species/dtf.neumap.csv"))

### UMAPs coloured by subset, tissue and disease group

dtf.neumap[celltype == "lowDepth", celltype := "Low depth"]
dtf.neumap$celltype = factor(dtf.neumap$celltype, levels = names(col_celltype))
df.label = as.data.table(dtf.neumap)[, lapply(.SD, median), by = celltype,
                                     .SDcols = c("UMAP_1", "UMAP_2")]
colnames(df.label)[2:3] = c("x", "y")
df.label$labeltext = df.label[, "celltype", with= FALSE]
df.label[,y:=ifelse(celltype == "IL1B", y-1, y)]
df.label[,y:=ifelse(celltype == "S100A4", y+0.5, y)]
df.label[,y:=ifelse(celltype == "IFN1", y-2, y)]
df.label[,x:=ifelse(celltype == "IFN1", x+1, x)]
df.label[,x:=ifelse(celltype == "CCL3/4", x+1, x)]
dtf.neumap$facetvar = "Subsets"

p.umap.celltype =
  plot_scatter(dtf.neumap, x = "UMAP_1", y = "UMAP_2",
               color_by = "celltype", color_type = "discrete", colors = col_celltype,
               facet_by = "facetvar", facet_nrow = 1,
               point_size = 0.01, point_alpha = 0.3, raster_dpi = 300,
               shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL,
               legend_ncol = 1, legend_point_size = 1, label = TRUE,
               label_df = df.label, label_size = 1.5, repel = FALSE) +
  theme_expresso(legend_position = "none", facet_label_face = "plain",
                 panel_background = "white", plot_background = "grey97",
                 legend_text_size = 8, text_size = 8,
                 legend_key_height = 0.1, show_axis = FALSE, grid = "none")

col_tissue2 = c(`peripheral blood` = "#EDD1D8",  `bone marrow` = "#73DAFF",
                spleen = "#8CC269", liver = "#7C1823", lung = "#DC3023", pancreas = "#FF0097",
                gut = "#FFDFB2", peritoneum = "#0E5FDB", skin = "#B0A4E3", brain = "grey17",
                `breast` = "#352A87", heart = "#6A1A99", placenta = "#F8D626")

dtf.neumap$facetvar = "Tissue"
p.umap.tissue =
  plot_scatter(dtf.neumap, x = "UMAP_1", y = "UMAP_2", color_by = "tissue", color_type = "discrete",
               colors = col_tissue2, facet_by = "facetvar", facet_nrow = 1,
               point_size = 0.01, point_alpha = 0.3, raster_dpi = 300,
               shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL,
               legend_ncol = 2, legend_point_size = 1.5) +
  scale_color_manual(values = col_tissue2, labels = cap_first, name = NULL,
                     na.value = "grey42") +
  theme_expresso(legend_position = c(1.05,0.6), facet_label_face = "plain",
                 panel_background = "white", plot_background = "grey97",
                 legend_text_size = 8, text_size = 8,
                 legend_key_height = 0.05, legend_key_spacing_y = 0.05,
                 show_axis = FALSE, grid = "none")

col_group2 = c(healthy = "#8CC269", flu = "cyan", `acute inflammation` = "#EFC000",
               `lung cancer` = "#DC3023", `breast cancer` = "#352A87",
               `pancreatic cancer` = "#FF0097", pancreatitis = "grey67",
               stroke = "#0E5FDB", infarction = "#7C1823",
               peritonitis = "#C5BB5C", fibrosis = "#06A5C7", `biliary damage` = "#7b72c5")

dtf.neumap$facetvar = "Group"
p.umap.group =
  plot_scatter(dtf.neumap, x = "UMAP_1", y = "UMAP_2", color_by = "group", color_type = "discrete",
               colors = col_group2, facet_by = "facetvar", facet_nrow = 1,
               point_size = 0.01, point_alpha = 0.3, raster_dpi = 300,
               shuffle = TRUE, seed = 42, na_color = "grey42",
               title = "Cells from NeuMap (mouse) projected to SCAHN",
               legend_ncol = 2, legend_point_size = 1.5) +
  scale_color_manual(values = col_group2, labels = cap_first, name = NULL,
                     na.value = "grey42") +
  theme_expresso(legend_position = c(1.05,0.6), text_size = 8,
                 facet_label_face = "plain", panel_background = "white",
                 plot_background = "grey97", legend_text_size = 8,
                 legend_key_height = 0.05, legend_key_spacing_y = 0.05,
                 show_axis = FALSE, grid = "none")

### marker dotplot
mat = readRDS(paste0(dir.results, "cross_species/NeuMap.signaturegenes.mat.rds"))
setdiff(signaturegenes, rownames(mat))

dtf.neumap.expr = as.data.table(cbind(dtf.neumap, as.matrix(t(mat))))

dtf.neumap.expr = cbind(dtf.neumap, as.matrix(t(mat)))
dtf.neumap.expr$DEFA4 = NULL # it is all 0 (+1)
dtf.neumap.expr$celltype = factor(dtf.neumap.expr$celltype, levels = rev(names(col_celltype)))

p.markerdot.neumap = plot_dotplot(
  dtf.neumap.expr[!is.na(celltype)],
  features = setdiff(intersect(signaturegenes, rownames(mat)), "DEFA4"),
  group_by = "celltype", max_scale = 2, dot_scale = 3, dot_stroke = 0.05,
  title = NULL, col_fontsize = 8, row_fontsize = 8, cluster_features = FALSE,
  cluster_groups = FALSE, show_axis = FALSE, feature_side = "top",
  colorbar_title = "Scaled expression")

### selected genes on the NeuMap UMAP, labelled with the mouse symbol
pgenes.NeuMap.mice = list()
for(id in c("STMN1", "MPO", "PRTN3", "LTF", "MMP8", "MMP9", "EGR1", "CXCL2",
            "CCL4", "IL1B", "SLPI", "HSPA1B")){
  id_mice = mouse_genes_human[human == id]$mouse[1]
  dtf.neumap.expr$facetvar = id_mice
  pgenes.NeuMap.mice[[id]] =
    plot_scatter(dtf.neumap.expr, x = "UMAP_1", y = "UMAP_2", color_by = id,
                 facet_by = "facetvar", color_type = "continuous", colors = c("grey99",
                            scales::dichromat_pal("DarkRedtoBlue.12")(12)[7:12]),
                 quantile_lower = 0.001, quantile_upper = 0.999,
                 shuffle = TRUE, seed = 42, point_size = 0.01,
                 point_alpha = 0.3, raster_dpi = 300, legend_ncol = 1,
                 facet_nrow = 1, title = NULL) +
    theme_expresso(text_size = 8, legend_position = "none",
                   show_axis = FALSE, grid = "none", facet_label_face = "italic",
                   panel_background = "white", plot_background = "grey97",
                   plot_margin = ggplot2::margin(0.02, 0.02, 0.02, 0.02, "line"))
}

### subset composition by tissue
dtf.neumap.prop = cell_proportions(dtf.neumap, sample_col = "sampleid", cluster_col = "celltype")
dtf.neumap.prop$`NA` = NULL
dtf.neumap.prop = cbind(
  dtf.neumap.prop, samplemeta.neumap[match(dtf.neumap.prop$sampleid, samplemeta.neumap$sampleid)][
    , c("tissue", "sex", "age", "group")])
table(dtf.neumap.prop$group)
dtf.neumap.prop[,group2:=ifelse(grepl("cancer", group), "cancer", group)]
dtf.neumap.prop[
  , tissue2 := ifelse(grepl("breast|liver|lung|pancreas|gut|heart|peritoneum|brain|placenta|skin",
                            tissue), "solid tissue", tissue)]
dtf.neumap.prop$tissue2 = factor(
  dtf.neumap.prop$tissue2, levels = c("peripheral blood", "bone marrow", "spleen", "solid tissue"))
table(dtf.neumap.prop$group2, dtf.neumap.prop$tissue2)

dtf.neumap.prop.m = melt(
  dtf.neumap.prop, measure.vars = intersect(colnames(dtf.neumap.prop), names(col_celltype)))
dtf.neumap.prop.m = as.data.table(dtf.neumap.prop.m)

dtf.neumap.prop.m.agg = dtf.neumap.prop.m[
  , lapply(.SD, function(x){mean(x)}), by = c("tissue2",  "variable"), .SDcols = c("value")]
dtf.neumap.prop.m.agg$tissue2 = factor(
  dtf.neumap.prop.m.agg$tissue2,
  levels = rev(c("peripheral blood", "bone marrow", "spleen", "solid tissue")))
dtf.neumap.prop.m.agg$variable = factor(dtf.neumap.prop.m.agg$variable,
                                        levels = rev(names(col_celltype)))
dtf.neumap.prop.m.agg$facetvar = "Subset proportions"

p.neumap.distribution =
  ggplot(dtf.neumap.prop.m.agg, aes(x = tissue2, y = value, fill = variable)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  theme_expresso(legend_position = "none") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = seq(0, 1, 0.2), position = "left") +
  # cap_first on the drawn tick labels only: tissue2 keeps the lower-case strings
  # its factor levels are built from just above.
  scale_x_discrete(position = "bottom", labels = cap_first) +
  facet_wrap(.~facetvar) +
  labs(x = NULL, y = NULL , title = NULL) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(),
        plot.title = element_text(size = textsize, face = "plain", hjust = 0,
                                  margin = margin(t = 0.1, r = 0.1, b = 0.1,
                                                  l = 0.1, unit = "line"))) +
  scale_fill_manual(values = col_celltype, name = "") + coord_flip()  +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(size = 2.3)))


# -- Percentage of cells expressing each gene, by dataset ----------------------

### SCAHN (human)
mat = readRDS(paste0(dir.results, "annotation/SCAHN.selectedgenes.mat.expr.integrated.rds"))
dtf.SCAHN.expr = cbind(dtf, as.matrix(t(mat[signaturegenes, dtf$cell.pid])))
dtf.SCAHN.expr$tissue0 = dtf.SCAHN.expr$tissue
dtf.SCAHN.expr[, tissue := ifelse(grepl("freshWB|neutrophil", source), "peripheral blood", "other")]
dtf.SCAHN.expr[, tissue := ifelse(grepl("bone marrow|spleen|lung|liver", tissue0), tissue0, tissue)]
table(dtf.SCAHN.expr$tissue)
dtf.SCAHN.expr$sampleid = dtf.SCAHN.expr$sampleid.pid

dtf.SCAHN.expr.pct = dtf.SCAHN.expr[nCount_RNA >= 400][,
  c(list(cellnumber = .N), lapply(.SD, function(x){mean(x>0)})), by = c("tissue", "pid"),
  .SDcols = signaturegenes]
dtf.SCAHN.expr.pct2 = dtf.SCAHN.expr[nCount_RNA >= 400][,
  c(list(cellnumber = .N), lapply(.SD, function(x){mean(x>0)})),
  by = c("tissue", "pid", "sampleid"), .SDcols = signaturegenes]
dtf.SCAHN.expr.pct$species = "human"
dtf.SCAHN.expr.pct2$species = "human"

### NeuMap (mouse)
srt.neumap = readRDS(paste0(dir.results, "cross_species/NeuMap.135055.srt.rds"))
dtf.neumap = fread(paste0(dir.results, "cross_species/dtf.neumap.csv"))
dtf.neumap.expr = cbind(dtf.neumap, as.matrix(t(
  srt.neumap@assays$RNA@data[intersect(signaturegenes, rownames(srt.neumap@assays$RNA@data)), ])))
dtf.neumap.expr$pid = "NeuMap"

dtf.neumap.expr.pct = dtf.neumap.expr[nCount_RNA >= 400][,
  c(list(cellnumber = .N), lapply(.SD, function(x){mean(x>0)})), by = c("tissue", "pid"),
  .SDcols = intersect(signaturegenes, rownames(srt.neumap@assays$RNA@data))]
dtf.neumap.expr.pct2 = dtf.neumap.expr[nCount_RNA >= 400][,
  c(list(cellnumber = .N), lapply(.SD, function(x){mean(x>0)})),
  by = c("tissue", "pid", "sampleid"),
  .SDcols = intersect(signaturegenes, rownames(srt.neumap@assays$RNA@data))]
dtf.neumap.expr.pct$species = "mouse"
dtf.neumap.expr.pct2$species = "mouse"

### grieshaber-bouyer2021 / GSE165276 (mouse)
srt.GSE165276 = readRDS(paste0(dir.results, "cross_species/srt.GSE165276.rds"))
meta.GSE165276 = fread(paste0(dir.results, "cross_species/GSE165276.metacell.csv"))
cells2include = meta.GSE165276[
  !RNAmice_snn_res.0.8 %in% c(3, 15, 14, 12, 18, 16, 17, 8)]$cell # 12962
srt.GSE165276 = subset(srt.GSE165276, cells = cells2include)
dtf.GSE165276 = as.data.table(srt.GSE165276@meta.data)
dtf.GSE165276.expr = cbind(dtf.GSE165276, as.matrix(t(
  srt.GSE165276@assays$RNA@data[
    intersect(signaturegenes, rownames(srt.GSE165276@assays$RNA@data)), ])))
dtf.GSE165276.expr$pid = "grieshaber-ouyer2021"

dtf.GSE165276.expr.pct = dtf.GSE165276.expr[nCount_RNA >= 400][,
  c(list(cellnumber = .N), lapply(.SD, function(x){mean(x>0)})), by = c("tissue", "pid"),
  .SDcols = intersect(signaturegenes, rownames(srt.GSE165276@assays$RNA@data))]
dtf.GSE165276.expr.pct2 = dtf.GSE165276.expr[nCount_RNA >= 400][,
  c(list(cellnumber = .N), lapply(.SD, function(x){mean(x>0)})),
  by = c("tissue", "pid", "sampleid"),
  .SDcols = intersect(signaturegenes, rownames(srt.GSE165276@assays$RNA@data))]
dtf.GSE165276.expr.pct$species = "mouse"
dtf.GSE165276.expr.pct2$species = "mouse"

### xie2020 (mouse)
srt.xie2020 = readRDS(paste0(dir.results, "cross_species/srt.xie2020mice.rds"))
meta.xie2020 = fread(paste0(dir.results, "cross_species/xie2020mice.metacell.csv"))
cells2include = meta.xie2020[
  RNAmice_snn_res.0.8 %in% c(17, 18, 19, 12, 11, 22, 28, 20, 4, 6, 15, 0, 2, 5, 7)]$cell # 31505
srt.xie2020 = subset(srt.xie2020, cells = cells2include)
dtf.xie2020 = as.data.table(srt.xie2020@meta.data)
dtf.xie2020.expr = cbind(dtf.xie2020, as.matrix(t(
  srt.xie2020@assays$RNA@data[
    intersect(signaturegenes, rownames(srt.xie2020@assays$RNA@data)), ])))
dtf.xie2020.expr$pid = "xie2020"

dtf.xie2020.expr.pct = dtf.xie2020.expr[nCount_RNA >= 400][,
  c(list(cellnumber = .N), lapply(.SD, function(x){mean(x>0)})), by = c("tissue", "pid"),
  .SDcols = intersect(signaturegenes, rownames(srt.xie2020@assays$RNA@data))]
dtf.xie2020.expr.pct2 = dtf.xie2020.expr[nCount_RNA >= 400][,
  c(list(cellnumber = .N), lapply(.SD, function(x){mean(x>0)})),
  by = c("tissue", "pid", "sampleid"),
  .SDcols = intersect(signaturegenes, rownames(srt.xie2020@assays$RNA@data))]
dtf.xie2020.expr.pct$species = "mouse"
dtf.xie2020.expr.pct2$species = "mouse"

### mahyari2025 (NHP)
srt.mahyari2025 = readRDS(paste0(dir.data, "cross_species/GSE277821_RIRA.Myeloid.seurat.rds"))

dtf.mahyari2025 = fread(paste0(dir.results, "cross_species/mahyari2025.metacell.csv"))
dtf.mahyari2025$sampleid = paste0(dtf.mahyari2025$SubjectId, "." , dtf.mahyari2025$Tissue)
dtf.mahyari2025.s = dtf.mahyari2025[
  ClusterNames_1.2 %in% c(11, 13, 29, 23, 8, 21, 12, 17) & Tissue == "Bone marrow"]
cells2include = dtf.mahyari2025.s[
  sampleid %in% names(which(table(dtf.mahyari2025.s$sampleid) >= 30))]$cell
length(cells2include)
srt.mahyari2025 = subset(srt.mahyari2025, cells = cells2include)
srt.mahyari2025$tissue = gsub("Bone marrow", "bone marrow", srt.mahyari2025$Tissue)
srt.mahyari2025 = NormalizeData(srt.mahyari2025)

dtf.mahyari2025 = as.data.table(srt.mahyari2025@meta.data)
dtf.mahyari2025$tissue = tolower(dtf.mahyari2025$Tissue)
dtf.mahyari2025$tissue = gsub("pbmc", "PBMC", dtf.mahyari2025$tissue)
dtf.mahyari2025$tissue = gsub("pln|mesln", "lymph node", dtf.mahyari2025$tissue)

dtf.mahyari2025.expr = cbind(dtf.mahyari2025, as.matrix(t(
  srt.mahyari2025@assays$RNA@data[
    intersect(signaturegenes, rownames(srt.mahyari2025@assays$RNA@data)), ])))
dtf.mahyari2025.expr$pid = "mahyari2025"
dtf.mahyari2025.expr$sampleid = paste0(dtf.mahyari2025.expr$SubjectId, dtf.mahyari2025.expr$tissue)

dtf.mahyari2025.expr.pct2 = dtf.mahyari2025.expr[nCount_RNA >= 400][,
  c(list(cellnumber = .N), lapply(.SD, function(x){mean(x>0)})),
  by = c("tissue", "pid", "sampleid"), .SDcols = intersect(signaturegenes,
                      rownames(srt.mahyari2025@assays$RNA@data))]
dtf.mahyari2025.expr.pct2$species = "NHP"

### staupe2022 (NHP)
srt.staupe2022 = readRDS(paste0(dir.results, "cross_species/srt.staupe2022.rds"))
dtf.staupe2022 = fread(paste0(dir.results, "cross_species/staupe2022.metacell.csv"))
cells2include = dtf.staupe2022[
  RNA_snn_res.0.8 %in% c(6, 9, 12, 13, 14, 40, 19)]$cell
srt.staupe2022 = subset(srt.staupe2022, cells = cells2include)

dtf.staupe2022 = as.data.table(srt.staupe2022@meta.data)

dtf.staupe2022.expr = cbind(dtf.staupe2022, as.matrix(t(
  srt.staupe2022@assays$RNA@data[
    intersect(signaturegenes, rownames(srt.staupe2022@assays$RNA@data)), ])))
dtf.staupe2022.expr$pid = "staupe2022"

dtf.staupe2022.expr.pct2 = dtf.staupe2022.expr[nCount_RNA >= 400][,
  c(list(cellnumber = .N), lapply(.SD, function(x){mean(x>0)})),
  by = c("tissue", "pid", "sampleid"), .SDcols = intersect(signaturegenes,
                      rownames(srt.staupe2022@assays$RNA@data))]
dtf.staupe2022.expr.pct2$species = "NHP"


# -- Combine the percent-expressing tables -------------------------------------

### per donor, human + mouse
colss = c("tissue", "pid", "species", "cellnumber", Reduce(intersect, list(
  signaturegenes, rownames(srt.xie2020@assays$RNA@data),
  rownames(srt.neumap@assays$RNA@data), rownames(srt.GSE165276@assays$RNA@data))))
dtf.expr.pct = rbind(dtf.SCAHN.expr.pct[,colss, with = FALSE],
                     dtf.neumap.expr.pct[,colss, with = FALSE],
                     dtf.GSE165276.expr.pct[,colss, with = FALSE],
                     dtf.xie2020.expr.pct[,colss, with = FALSE])
dtf.expr.pct = as.data.table(dtf.expr.pct)
fwrite(dtf.expr.pct, paste0(dir.results, "cross_species/dtf.expr.pct.crossspecies.csv"))

### per sample, human + mouse
colss = c("tissue", "pid", "species", "sampleid", "cellnumber", Reduce(intersect, list(
  signaturegenes, rownames(srt.xie2020@assays$RNA@data),
  rownames(srt.neumap@assays$RNA@data), rownames(srt.GSE165276@assays$RNA@data))))

dtf.expr.pct = rbind(dtf.SCAHN.expr.pct2[,colss, with = FALSE],
                     dtf.neumap.expr.pct2[,colss, with = FALSE],
                     dtf.GSE165276.expr.pct2[,colss, with = FALSE],
                     dtf.xie2020.expr.pct2[,colss, with = FALSE])
dtf.expr.pct = as.data.table(dtf.expr.pct)
fwrite(dtf.expr.pct, paste0(
  dir.results, "cross_species/dtf.expr.pct.sample.crossspecies.csv"))

### per sample, human + NHP + mouse
colss = c("tissue", "pid", "species", "sampleid", "cellnumber", Reduce(intersect, list(
  signaturegenes, rownames(srt.xie2020@assays$RNA@data),
  rownames(srt.neumap@assays$RNA@data), rownames(srt.GSE165276@assays$RNA@data),
  rownames(srt.mahyari2025@assays$RNA@data), rownames(srt.staupe2022@assays$RNA@data))))

dtf.expr.pct = rbind(dtf.SCAHN.expr.pct2[,colss, with = FALSE],
                     dtf.neumap.expr.pct2[,colss, with = FALSE],
                     dtf.GSE165276.expr.pct2[,colss, with = FALSE],
                     dtf.xie2020.expr.pct2[,colss, with = FALSE],
                     dtf.mahyari2025.expr.pct2[,colss, with = FALSE],
                     dtf.staupe2022.expr.pct2[,colss, with = FALSE])
dtf.expr.pct = as.data.table(dtf.expr.pct)
fwrite(dtf.expr.pct, paste0(dir.results, "cross_species/dtf.expr.pct.sample.crossspeciesNHP.csv"))


# -- Percent-expressing figure panels ------------------------------------------

### bone marrow: human vs NHP vs mouse
dtf.expr.pct = fread(paste0(dir.results, "cross_species/dtf.expr.pct.sample.crossspeciesNHP.csv"))
dtf.expr.pct.m = data.table::melt(
  dtf.expr.pct[cellnumber >= 100], id.vars = c("tissue", "pid", "species", "sampleid"),
  measure.vars = setdiff(colnames(dtf.expr.pct),
                         c("tissue", "pid", "species", "sampleid", "cellnumber")))
dtf.expr.pct.m$species = factor(dtf.expr.pct.m$species, levels = c("human", "NHP", "mouse"))
dtf.expr.pct.m.s = dtf.expr.pct.m[
  variable %in% c("MPO", "CTSG" ,"LTF", "LCN2", "MMP8", "MMP9", "CD177")]
dtf.expr.pct.m.s$variable = factor(
  dtf.expr.pct.m.s$variable, levels = c("MPO", "CTSG","LTF", "LCN2", "MMP8", "MMP9", "CD177"))

dtf.expr.pct.m[, .(median = median(value)), by = .(variable, species, tissue)]
dtf.expr.pct.m[variable %in% c("LTF","LCN2")][tissue == "bone marrow"]

p.species_compare.BM = ggplot(dtf.expr.pct.m.s[tissue %in% c("bone marrow")],
                              aes(x = species, y = value)) +
  ggbeeswarm::geom_quasirandom(aes(color = species), size = 0.2, width = 0.2,
                               show.legend = TRUE, alpha = 0.5, varwidth = FALSE, shape = 16) +
  geom_boxplot(color = textcolor, width = 0.2, size = 0.2,
               outlier.shape = NA, alpha = 1, fill = NA) +
  scale_y_continuous(breaks = seq(0,1,0.4), labels = scales::percent) +
  facet_wrap(.~ variable, nrow = 1) +
  labs(x = NULL, y = "Percentage of neutrophils\nexpressing each gene", title = "Bone marrow") +
  theme_expresso(text_size = 8, axis_text_size = 8, axis_title_size = 8,
                 legend_text_size = 8, legend_position = c(0.2,-2.05),
                 facet_label_face = "italic", panel_background = "white") +
  stat_compare_means(comparisons = list(c("human", "NHP"), c("NHP", "mouse"), c("human", "mouse")),
                     method = "wilcox.test", vjust = 1.2, size = 1.6,
                     label = "p.signif", color = "grey17", bracket.size = 0.1, tip.length = 0.02,
                     step.increase = c(0,0,0.06)) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(), panel.spacing = unit(0.05, "lines"),
        axis.title.y = element_text(hjust = 0.8, margin = margin(b = 0)),
        plot.title =  element_text(size = 8, color = textcolor,
                                   margin = margin(0.01, 0, 0.01,0, "line")),
        strip.text.x = element_text(size = 8, color = textcolor,
                                    margin = margin(0.05, 0, 0.05,0, "line"))) +
  scale_color_manual(values = c("human" = "#264EFF", "NHP" = "#2BAE85", "mouse" = "#F47B00")) +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 1.5)))

### peripheral blood, and lung or liver: human vs mouse
dtf.expr.pct = fread(paste0(dir.results, "cross_species/dtf.expr.pct.sample.crossspecies.csv"))
dtf.expr.pct.m = data.table::melt(
  dtf.expr.pct[cellnumber >= 100], id.vars = c("tissue", "pid", "species", "sampleid"),
  measure.vars = setdiff(colnames(dtf.expr.pct),
                         c("tissue", "pid", "species", "sampleid", "cellnumber")))
dtf.expr.pct.m$species = factor(dtf.expr.pct.m$species, levels = c("human", "NHP", "mouse"))

dtf.expr.pct.m.s = dtf.expr.pct.m[tissue %in% c("peripheral blood")][
  variable %in% c("S100A8","S100A4", "SLPI")]
dtf.expr.pct.m.s$variable = factor(dtf.expr.pct.m.s$variable, levels = c("S100A8","S100A4", "SLPI"))

p.species_compare.PB = ggplot(dtf.expr.pct.m.s[tissue %in% "peripheral blood"],
                              aes(x = species, y = value)) +
  ggbeeswarm::geom_quasirandom(aes(color = species), size = 0.05, width = 0.2,
                               show.legend = TRUE, alpha = 0.3, varwidth = FALSE, shape = 16) +
  geom_boxplot(color = textcolor, width = 0.2, size = 0.2,
               outlier.shape = NA, alpha = 1, fill = NA) +
  scale_y_continuous(breaks = seq(0,1,0.4), labels = scales::percent) +
  facet_wrap(.~ variable, nrow = 1) +
  labs(x = NULL, y = " \n ", title = "Peripheral blood") +
  theme_expresso(text_size = 8, axis_text_size = 8, axis_title_size = 8,
                 legend_position = "none", facet_label_face = "italic",
                 panel_background = "white") +
  stat_compare_means(comparisons = list(c("human", "mouse")),
                     method = "wilcox.test", vjust = 1.2, size = 1.6,
                     label = "p.signif", color = "grey17", bracket.size = 0.1, tip.length = 0.02) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(), panel.spacing = unit(0.05, "lines"),
        axis.title.y = element_text(hjust = 0.7, margin = margin(b = 0)),
        plot.title =  element_text(size = 8, color = textcolor,
                                   margin = margin(0.01, 0, 0.01,0, "line")),
        strip.text.x = element_text(size = 8, color = textcolor,
                                    margin = margin(0.05, 0, 0.05,0, "line"))) +
  scale_color_manual(values = c("human" = "#264EFF", "NHP" = "#2BAE85", "mouse" = "#F47B00"))

dtf.expr.pct.m.s = dtf.expr.pct.m[tissue %in% c( "liver", "lung")][
  variable %in% c("CXCL1", "CXCL2", "IL1B", "CCL4")]
dtf.expr.pct.m.s$variable = factor(
  dtf.expr.pct.m.s$variable, levels = c("CXCL1", "CXCL2", "IL1B", "CCL4"))

p.species_compare.liverorlung = ggplot(dtf.expr.pct.m.s, aes(x = species, y = value)) +
  ggbeeswarm::geom_quasirandom(aes(color = species), size = 0.05, width = 0.2,
                               show.legend = TRUE, alpha = 0.3, varwidth = FALSE, shape = 16) +
  geom_boxplot(color = textcolor, width = 0.2, size = 0.2,
               outlier.shape = NA, alpha = 1, fill = NA) +
  scale_y_continuous(breaks = seq(0,1,0.4), labels = scales::percent) +
  facet_wrap(.~ variable, nrow = 1) +
  labs(x = NULL, y = " \n ", title = "Lung or liver") +
  theme_expresso(text_size = 8, axis_text_size = 8, axis_title_size = 8,
                 legend_position = "none", facet_label_face = "italic",
                 panel_background = "white") +
  stat_compare_means(comparisons = list(c("human", "mouse")),
                     method = "wilcox.test", vjust = 1.2, size = 1.6,
                     label = "p.signif", color = "grey17", bracket.size = 0.1, tip.length = 0.02) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(), panel.spacing = unit(0.05, "lines"),
        axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), plot.title =  element_text(size = 8, color = textcolor,
                                   margin = margin(0.01, 0, 0.01,0, "line")),
        strip.text.x = element_text(size = 8, color = textcolor,
                                    margin = margin(0.05, 0, 0.05,0, "line"))) +
  scale_color_manual(values = c("human" = "#264EFF", "NHP" = "#B29DDA", "mouse" = "#F47B00"))


# -- NeuMap and published mouse subset scores in SCAHN -------------------------

hubgenes = readRDS(paste0(dir.results, "cross_species/SCAHN.neumap.scores.rds"))
colnames(hubgenes)

colnames(hubgenes) = c("cell.pid", "UMAP_1", "UMAP_2",
                       "immature\n(Cebpe+Mmp8+)", "IFN response\n(Ifit1+Cd274+)",
                       "IS-I\n(Cd14+Ptgs2+)", "preNeu\n(Ltf+mKi67+)",
                       "APC\n(H2+Cd74+)", "IS-II\n(Vegfa+Cd274+)",
                       "immuno-silent\n(Cd52+ circulating)",
                       "immuno-silent", "IS-I", "IFN response", "APC&IS-II", "immature", "preNeu",
                       "G0", "G1", "G2", "G3", "G4",
                       "G5a", "G5b", "G5c", "hG5a", "hG5b", "hG5c", "P1", "P2", "P3",
                       "P4", "neutrotime.pos", "neutrotime.neg", "T3", "T2", "T1")

hubgenes$celltype = dtf$celltype

### NeuMap hub-gene scores on the SCAHN UMAP

pscores.NeuMap.mice = list()
for(hub in c("preNeu\n(Ltf+mKi67+)", "immature\n(Cebpe+Mmp8+)",
             "IFN response\n(Ifit1+Cd274+)", "immuno-silent\n(Cd52+ circulating)",
             "IS-I\n(Cd14+Ptgs2+)", "IS-II\n(Vegfa+Cd274+)", "APC\n(H2+Cd74+)")){
  dtf1 = data.table(UMAP_1 = hubgenes$UMAP_1, UMAP_2 = hubgenes$UMAP_2,
                    id = hubgenes[[hub]] - 1, facetvar = gsub("\n.*", "", hub))
  pscores.NeuMap.mice[[hub]] =
    plot_scatter(dtf1, x = "UMAP_1", y = "UMAP_2", color_by = "id",
                 facet_by = "facetvar", color_type = "continuous", colors = c("grey99",
                            scales::dichromat_pal("DarkRedtoBlue.12")(12)[7:12]),
                 quantile_lower = 0.001, quantile_upper = 0.999,
                 shuffle = TRUE, seed = 42, point_size = 0.000001,
                 point_alpha = 0.2, raster_dpi = 400, legend_ncol = 1,
                 facet_nrow = 1, title = NULL) +
    theme_expresso(text_size = 9, legend_position = "none",
                   show_axis = FALSE, grid = "none", facet_label_face = "plain",
                   panel_background = "white", plot_background = "grey97",
                   plot_margin = ggplot2::margin(0.05, 0.05, 0.05, 0.05, "line"))
}

### published mouse subset scores on the SCAHN UMAP (Extended figure 4A-C)
pscores.mice = list()
for(hub in c("G0", "G1", "G2", "G3", "G4", "G5a", "G5b", "G5c",
"hG5a", "hG5b", "hG5c", "P1", "P2", "P3", "P4",  "T1", "T2", "T3")){
  dtf1 = data.table(UMAP_1 = hubgenes$UMAP_1, UMAP_2 = hubgenes$UMAP_2,
                    id = hubgenes[[hub]] - 1, facetvar = gsub("\n.*", "", hub))
  pscores.mice[[hub]] =
    plot_scatter(dtf1, x = "UMAP_1", y = "UMAP_2", color_by = "id",
                 facet_by = "facetvar", color_type = "continuous", colors = c("grey99",
                            scales::dichromat_pal("DarkRedtoBlue.12")(12)[7:12]),
                 quantile_lower = 0.001, quantile_upper = 0.999,
                 shuffle = TRUE, seed = 42, point_size = 0.000001,
                 point_alpha = 0.2, raster_dpi = 400, legend_ncol = 1,
                 facet_nrow = 1, title = NULL) +
    theme_expresso(text_size = 9, legend_position = "none",
                   show_axis = FALSE, grid = "none", facet_label_face = "plain",
                   panel_background = "white", plot_background = "grey97",
                   plot_margin = ggplot2::margin(0.01, 0.01, 0.01, 0.01, "line"))
}


# -- SCAHN-to-NeuMap subset correspondence -------------------------------------

### hand-curated mapping between SCAHN subsets and NeuMap subsets (Figure 3B)
mappings =  c(preNeu = "AZU1", immature = "LTF", immature = "MMP9",
              `immuno-silent\n or undefined` = "MMP9", `immuno-silent\n or undefined` = "S100A4",
              `immuno-silent\n or undefined` = "IL1R2", `IS-I` = "EGR1", `IS-I` = "AP-1",
              `IS-I` = "PTGS2", `IS-I` = "G0S2", `IFN response` = "IFN1", `IFN response` = "IFN2",
              `IFN response` = "IFN3",
              `IS-I` = "NF-κB", `IS-I` = "IL1B", `IS-I` = "IL1RN",
              `IS-II` = "CXCL", `IS-II` = "VEGFA", `IS-II` = "CCL3/4",
              `IS-II` = "NF-κB", `IS-II` = "IL1B", `IS-II` = "IL1RN",
              `immuno-silent\n or undefined`  = "SLPI", `IS-II` = "HSP", `APC` = "HSP",
              APC = "CD74",
              `immuno-silent\n or undefined` = "MME", `immuno-silent\n or undefined` = "TXNIP")

df = data.frame(SCAHN = mappings, NeuMap = names(mappings))

scahn_colors = c(AZU1 = "#DC3023", LTF = "#F47B00", MMP9 = "#F8D626",
                  S100A4 = "#FFDFB2", IL1R2 = "#acf300",
                  MME = "#80C684", TXNIP = "#D1BA58", EGR1 = "#1A5E1F",
                  `AP-1` = "#2BAE85", PTGS2 = "#AFBC65", G0S2 = "#E2E7BF", IFN1 = "#2E42B8",
                  IFN2 = "#0A73DC", IFN3 = "#264EFF", CXCL = "#E4C6D0", VEGFA = "#73DAFF",
                  `CCL3/4` = "#DCE318", `NF-κB` = "#5E34B1", IL1B = "#9E9AC8",
                  IL1RN = "#DCB0F2", SLPI = "#7CABB1", HSP = "#4D6D93", CD74 = "#C7E5C9")

neumap_colors = c(
  "preNeu" = "#8ABE75", "immature" = "#F0BB46",
  "IFN response" = "#d69aad", "IS-I" = "#40A1FF",
  "IS-II" = "#f19e73", "APC" = "#D1C4E9", "immuno-silent\n or undefined" = "grey99")

all_colors = c(scahn_colors, neumap_colors)

df$SCAHN  = factor(df$SCAHN,  levels = names(scahn_colors))
df$NeuMap = factor(df$NeuMap, levels = names(neumap_colors))
df$freq = 1 / ave(rep(1, nrow(df)), df$SCAHN, FUN = length)

p.alluvial = ggplot(df, aes(axis1 = SCAHN, axis2 = NeuMap, y = freq)) +
  geom_flow(aes(fill = NeuMap), width = 0.7) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.7, color = "white",
               linewidth = 0.01, size = 2) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 2.2, fontface = "plain", color = textcolor, lineheight = 0.8) +
  scale_x_discrete(limits = c("SCAHN", "NeuMap"), expand = c(0, 0), position = "top") +
  scale_fill_manual(values = all_colors, guide = "none") +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(color = "grey97", fill = "grey97"),
    panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
    axis.text.x.top = element_text(size = 8, face = "plain", color = textcolor,
                                   margin = margin(b = -8)),
    plot.margin = margin(0, 0, 0, 0, unit = "line"))


# -- mahyari2025 (NHP) projected onto SCAHN by symphony ------------------------

query = readRDS(paste0(dir.results, "cross_species/mahyari2025.symphony.rds"))

dtf.mahyari2025 = as.data.table(query$umap)
colnames(dtf.mahyari2025) = c("UMAP_1", "UMAP_2")
dtf.mahyari2025$celltype = query$cell_type_pred
dtf.mahyari2025$celltype_conf = query$cell_type_conf
table(query$cell_type_pred)
mean(dtf.mahyari2025$celltype_conf > 0.5)

sort(table(dtf.mahyari2025$celltype))
table(dtf.mahyari2025$tissue)

tt = query$meta_data[,c("SubjectId", "Tissue",  "nCount_RNA", "nFeature_RNA")]
colnames(tt)[1] = "subjectid"
colnames(tt)[2] = "tissue"

dtf.mahyari2025 = cbind(dtf.mahyari2025, tt)
dtf.mahyari2025$tissue = tolower(dtf.mahyari2025$tissue)
dtf.mahyari2025$tissue = gsub("pbmc", "peripheral blood", dtf.mahyari2025$tissue)
dtf.mahyari2025$tissue = gsub("^mesln|^pln", "lymph node", dtf.mahyari2025$tissue)
dtf.mahyari2025$sampleid = paste0(dtf.mahyari2025$subjectid, ".", dtf.mahyari2025$tissue)

table(query$cell_type_pred)
table(dtf.mahyari2025$celltype, dtf.mahyari2025$tissue )

### UMAP coloured by subset (Figure 3H)
# As on dtf.neumap above: the projected labels arrive keyed "lowDepth".
dtf.mahyari2025[celltype == "lowDepth", celltype := "Low depth"]
dtf.mahyari2025$celltype = factor(dtf.mahyari2025$celltype, levels = names(col_celltype))
df.label = as.data.table(dtf.mahyari2025)[, lapply(.SD, median), by = celltype,
                                          .SDcols = c("UMAP_1", "UMAP_2")]
colnames(df.label)[2:3] = c("x", "y")
df.label$labeltext = df.label[, "celltype", with= FALSE]
df.label[,y:=ifelse(celltype == "TXNIP",y+0.5,y)]
df.label[,x:=ifelse(celltype == "VEGFA",x-0.5,x)]
dtf.mahyari2025$facetvar = "Subsets"

p.umap.celltype.mahyari2025 =
  plot_scatter(dtf.mahyari2025, x = "UMAP_1", y = "UMAP_2",
               color_by = "celltype", color_type = "discrete", colors = col_celltype,
               facet_by = "facetvar", facet_nrow = 1,
               point_size = 0.1, point_alpha = 0.3, raster_dpi = 300,
               shuffle = TRUE, seed = 42, na_color = "grey42",
               title = "NHP bone marrow neutrophils mapped to SCAHN",
               legend_ncol = 1, legend_point_size = 1, label = TRUE,
               label_df = df.label, label_size = 1.5, repel = FALSE) +
  theme_expresso(legend_position = "none", panel_background = "white",
                 facet_label_face = "plain", legend_text_size = 8,
                 text_size = 8, legend_key_height = 0.1, show_axis = FALSE, grid = "none")

### selected genes on the mahyari2025 UMAP
mat = readRDS(paste0(dir.results, "cross_species/mahyari2025.signaturegenes.mat.rds"))
dtf.mahyari2025.expr = cbind(dtf.mahyari2025, as.matrix(t(mat)))

pgenes.mahyari2025 = list()
for(id in c("STMN1", "MPO", "PRTN3", "LTF", "MMP8", "MMP9")){
  dtf.mahyari2025.expr$facetvar = id
  pgenes.mahyari2025[[id]] =
    plot_scatter(dtf.mahyari2025.expr, x = "UMAP_1", y = "UMAP_2",
                 color_by = id, facet_by = "facetvar", color_type = "continuous",
                 colors = c("grey99", scales::dichromat_pal("DarkRedtoBlue.12")(12)[7:12]),
                 quantile_lower = 0.001, quantile_upper = 0.999,
                 shuffle = TRUE, seed = 42, point_size = 0.01,
                 point_alpha = 0.3, raster_dpi = 300, legend_ncol = 1,
                 facet_nrow = 1, title = NULL) +
    theme_expresso(text_size = 8, legend_position = "none",
                   show_axis = FALSE, grid = "none", facet_label_face = "italic",
                   panel_background = "white", plot_background = "grey97",
                   plot_margin = ggplot2::margin(0.02, 0.02, 0.02, 0.02, "line"))
}


# -- Figure 3 ------------------------------------------------------------------

library(png)
cairo_pdf(file = paste0(dir.fig, "Figure3.pdf"), width = 6.6, height = 6.2)
pushViewport(viewport(layout = grid.layout(nrow = 152, ncol = 146)))

### panel backgrounds
print(ggplot() + theme_void() + theme(plot.background = element_rect(color = "grey97",
                                            fill = "grey97")),
      vp = viewport(layout.pos.row = 2:32, layout.pos.col =  2:145))

print(ggplot() + theme_void() + theme(plot.background = element_rect(color = "grey97",
                                            fill = "grey97")),
      vp = viewport(layout.pos.row = 34:115, layout.pos.col =  2:34))

print(ggplot() + theme_void() + theme(plot.background = element_rect(color = "grey97",
                                            fill = "grey97")),
      vp = viewport(layout.pos.row = 34:115, layout.pos.col =  36:145))

print(ggplot() + theme_void() + theme(plot.background = element_rect(color = "grey97",
                                            fill = "grey97")),
      vp = viewport(layout.pos.row = 117:151, layout.pos.col =  2:64))

print(ggplot() + theme_void() + theme(plot.background = element_rect(color = "grey97",
                                            fill = "grey97")),
      vp = viewport(layout.pos.row = 117:151, layout.pos.col =  66:145))

### A: NeuMap hub-gene scores in SCAHN
pushViewport(viewport(layout.pos.row = 9:30, layout.pos.col = 6:144))
grid.arrange(grobs=pscores.NeuMap.mice, nrow = 1, as.table = TRUE, newpage = FALSE)
upViewport(1)

### C-E: NeuMap projected onto SCAHN
print(p.umap.group, vp = viewport(layout.pos.row = 36:65, layout.pos.col = 40:61))
print(p.umap.tissue, vp = viewport(layout.pos.row = 69:93, layout.pos.col =  40:61))

print(p.umap.celltype + theme(plot.margin = margin(0.01,0.05,0.01,0.12, "line")),
      vp = viewport(layout.pos.row = 40:74 , layout.pos.col = 110:144))
print(p.neumap.distribution, vp = viewport(layout.pos.row = 76:96, layout.pos.col = 95:144))

### B: SCAHN-to-NeuMap subset correspondence
print(p.alluvial, vp = viewport(layout.pos.row = 38:115, layout.pos.col = 3:33))

### F: marker genes on the NeuMap UMAP
pushViewport(viewport(layout.pos.row = 101:112, layout.pos.col = 37:144))
grid.arrange(grobs=pgenes.NeuMap.mice, nrow = 1, as.table = TRUE, newpage = FALSE)
upViewport(1)

### H: mahyari2025 (NHP) projected onto SCAHN
print(p.umap.celltype.mahyari2025, vp = viewport(layout.pos.row = 118:150, layout.pos.col = 6:31))

pushViewport(viewport(layout.pos.row = 123:150, layout.pos.col = 33:63))
grid.arrange(grobs=pgenes.mahyari2025, nrow = 2, as.table = TRUE, newpage = FALSE)
upViewport(1)

### I: percent-expressing across species
print(p.species_compare.BM, vp = viewport(layout.pos.row = 118:132, layout.pos.col = 66:144))
print(p.species_compare.PB, vp = viewport(layout.pos.row = 133:147, layout.pos.col = 66:105))
print(p.species_compare.liverorlung, vp = viewport(layout.pos.row = 133:147,
                                                   layout.pos.col = 106:144))

### species icons
img = readPNG(paste0(dir.data, "cross_species/human.PNG"))
img_grob = rasterGrob(img, interpolate = TRUE)
pushViewport(viewport(layout.pos.row = 2:6, layout.pos.col = 85:88))
grid.draw(img_grob)
popViewport()

img = readPNG(paste0(dir.data, "cross_species/human.PNG"))
img_grob = rasterGrob(img, interpolate = TRUE)
pushViewport(viewport(layout.pos.row = 36:40, layout.pos.col = 13:16))
grid.draw(img_grob)
popViewport()

img = readPNG(paste0(dir.data, "cross_species/mouse.PNG"))
img_grob = rasterGrob(img, interpolate = TRUE)
pushViewport(viewport(layout.pos.row = 37:40, layout.pos.col = 32:34))
grid.draw(img_grob)
popViewport()

img = readPNG(paste0(dir.data, "cross_species/mouse.PNG"))
img_grob = rasterGrob(img, interpolate = TRUE)
pushViewport(viewport(layout.pos.row = 35:38, layout.pos.col = 93:95))
grid.draw(img_grob)
popViewport()

img = readPNG(paste0(dir.data, "cross_species/monkey.PNG"))
img_grob = rasterGrob(img, interpolate = TRUE)
pushViewport(viewport(layout.pos.row = 117:120, layout.pos.col = 61:63))
grid.draw(img_grob)
popViewport()

### panel labels
grid.text(x = unit(0.02,"npc"), y=unit(0.97,"npc"), label = "a",
          gp=gpar(fontsize = 8,fontface = "bold"))
grid.text(x = unit(0.034,"npc"), y=unit(0.97,"npc"),
          label = "Expression of human orthologs of NeuMap (mouse) hub genes in SCAHN",
          gp=gpar(fontsize = 8, fontface = "plain"), hjust = 0)

grid.text(x = unit(0.02,"npc"), y=unit(0.77,"npc"), label = "b",
          gp=gpar(fontsize = 8,fontface = "bold"))
grid.text(x = unit(0.26,"npc"), y=unit(0.73,"npc"), label = "c",
          gp=gpar(fontsize = 8,fontface = "bold"))

grid.text(x = unit(0.26,"npc"), y=unit(0.545,"npc"), label = "d",
          gp=gpar(fontsize = 8,fontface = "bold"))
grid.text(x = unit(0.73,"npc"), y=unit(0.74,"npc"), label = "e",
          gp=gpar(fontsize = 8,fontface = "bold"))
grid.text(x = unit(0.26,"npc"), y=unit(0.35,"npc"), label = "f",
          gp=gpar(fontsize = 8,fontface = "bold"))
grid.text(x = unit(0.73,"npc"), y=unit(0.49,"npc"), label = "g",
          gp=gpar(fontsize = 8,fontface = "bold"))
grid.text(x = unit(0.02,"npc"), y=unit(0.225,"npc"), label = "h",
          gp=gpar(fontsize = 8,fontface = "bold"))
grid.text(x = unit(0.5,"npc"), y=unit(0.225,"npc"), label = "i",
          gp=gpar(fontsize = 8,fontface = "bold"))

dev.off()


# -- Extended figure 4 ---------------------------------------------------------

### bracket-and-label helper for grouping panels by source dataset
draw_bracket_label = function(bracket_rows, cols, label, gap = 2,
                                margin = c(0, 0, 0, 0),   # top, right, bottom, left in npc
                                lwd = 0.7, ticks = 0.5, h = 0.3, type = 4, curvature = 0.3,
                                bracket_col = "grey42",
                                textsize = 8, textface = "plain", textcol = "black") {
  text_row = min(bracket_rows) - gap
  all_rows = text_row:max(bracket_rows)
  n_all    = length(all_rows)
  n_brack  = length(bracket_rows)

  mt = margin[1]; mr = margin[2]; mb = margin[3]; ml = margin[4]

  # Outer viewport covers the full row/col region
  pushViewport(viewport(layout.pos.row = all_rows, layout.pos.col = cols))
  # Inner viewport applies margins
  pushViewport(viewport(
    x = unit(ml, "npc"), y = unit(mb, "npc"), width = unit(1 - ml - mr, "npc"),
    height = unit(1 - mt - mb, "npc"), just = c("left", "bottom")))

  bracket_y = (n_brack / 2) / n_all
  text_y    = 1 - 1 / (2 * n_all)

  grid.brackets(0, bracket_y, 1, bracket_y, lwd = lwd, ticks = ticks, h = h,
                type = type, curvature = curvature, col = bracket_col)
  grid.text(0.5, text_y, label = label,
            gp = gpar(fontsize = textsize, fontface = textface, col = textcol))

  upViewport(2)
}

cairo_pdf(file = paste0(dir.fig, "Fig.e4.pdf"), width = 10.5, height = 5.5)
pushViewport(viewport(layout=grid.layout(nrow = 155 , ncol = 450)))

### A: xie2020 G / hG subset scores
pushViewport(viewport(layout.pos.row = 11:32, layout.pos.col = 10:340))
grid.arrange(grobs=pscores.mice[c(paste0("G",c(0:4, "5a", "5b", "5c")), "hG5a", "hG5b", "hG5c")],
             nrow = 1, as.table = TRUE, newpage = FALSE)
upViewport(1)

### B: grieshaber-bouyer2021 P subset scores
pushViewport(viewport(layout.pos.row = 43:64, layout.pos.col = 10:130))
grid.arrange(grobs=pscores.mice[paste0("P",c(1:4))], nrow = 1, as.table = TRUE, newpage = FALSE)
upViewport(1)

### C: ng2024 T subset scores
pushViewport(viewport(layout.pos.row = 43:64, layout.pos.col = 160:250))
grid.arrange(grobs=pscores.mice[paste0("T",c(1:3))], nrow = 1, as.table = TRUE, newpage = FALSE)
upViewport(1)

draw_bracket_label(bracket_rows = 7:10, cols = 10:251, label = "Xie2020, mouse", gap = 2,
                   margin = c(0.01, 0.01, 0.01, 0.01), textsize = textsize, textcol = textcolor)
draw_bracket_label(bracket_rows = 7:10, cols = 251:340, label = "Xie2020, human", gap = 2,
                   margin = c(0.01, 0.01, 0.01, 0.01), textsize = textsize, textcol = textcolor)

draw_bracket_label(bracket_rows = 39:42, cols = 10:130,
                   label = "Grieshaber-bouyer2021, mouse", gap = 2,
                   margin = c(0.01, 0.01, 0.01, 0.01), textsize = textsize, textcol = textcolor)
draw_bracket_label(bracket_rows = 39:42, cols = 160:250, label = "Ng2024, mouse", gap = 2,
                   margin = c(0.01, 0.01, 0.01, 0.01), textsize = textsize, textcol = textcolor)

### D: NeuMap marker dotplot
print(p.markerdot.neumap, vp = viewport(layout.pos.row = 73:154, layout.pos.col = 2:450))

### panel labels
grid.text(x = unit(0.01,"npc"), y=unit(0.98,"npc"), label="a",
          gp=gpar(fontsize = textsize,fontface = "bold"))
grid.text(x = unit(0.01,"npc"), y=unit(0.75,"npc"), label="b",
          gp=gpar(fontsize = textsize,fontface = "bold"))
grid.text(x = unit(0.34,"npc"), y=unit(0.75,"npc"), label="c",
          gp=gpar(fontsize = textsize,fontface = "bold"))
grid.text(x = unit(0.01,"npc"), y=unit(0.53,"npc"), label="d",
          gp=gpar(fontsize = textsize,fontface = "bold"))
grid.text(x = unit(0.03,"npc"), y=unit(0.53,"npc"), label="NeuMap (mouse) mapped to SCAHN",
          gp=gpar(fontsize = textsize,fontface = "plain"), hjust = 0)
dev.off()
