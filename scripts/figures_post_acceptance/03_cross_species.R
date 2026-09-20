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
# Writes to figures/post_acceptance/:
#   Figure3.pdf   Figure 3
#   Fig.e4.pdf    Extended figure 4 (published mouse subset scores, NeuMap
#                 marker dotplot)
# ==============================================================================

SCAHN_SCRIPTS = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/scripts/figures_post_acceptance"
source(file.path(SCAHN_SCRIPTS, "00_setup.R"))
source(file.path(SCAHN_SCRIPTS, "01_load_data.R"))

# "lowDepth" -> "Low depth" as in 02_atlas.R, here on the palette keys and on the
# two symphony predicted-label columns built below -- between them, every place
# this script reads a subset name from. dtf is left alone: nothing here draws its
# celltype column, it only rides along on hubgenes.
names(col_celltype)[names(col_celltype) == "lowDepth"] = "Low depth"


# -- UMAP axis key -------------------------------------------------------------

p_umap_axis = function(xlab = "UMAP_1", ylab = "UMAP_2", title_size = textsize, line_size = 0.1,
                       arrow_size = 0.2, color = textcolor) {
  axis_line = element_line(color = color, linewidth = line_size, lineend = "round",
                           arrow = arrow(length = unit(arrow_size, "line"), type = "closed",
                                         angle = 20))
  ggplot(data.table(x = 0, y = 0), aes(x = x, y = y)) +
    geom_blank() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    labs(x = xlab, y = ylab) +
    theme_expresso_void(legend_position = "none") +
    theme(axis.line.x.bottom = axis_line, axis.line.y.left = axis_line,
          axis.title.x = element_text(size = title_size, color = color, hjust = 0,
                                      margin = margin(0.35, 0, 0, 0, "line")),
          axis.title.y = element_text(size = title_size, color = color, angle = 90, hjust = 0,
                                      margin = margin(0, 0.15, 0, 0, "line")),
          panel.background = element_blank(), plot.background = element_blank(),
          axis.ticks.length = unit(0, "pt"),
          plot.margin = margin(0, 0, 0, 0.1, "line"))
}

p.umap.axis = p_umap_axis()


# -- Mouse-to-human gene mapping -----------------------------------------------

# Built once with biomaRt (Ensembl 105) and cached; the manual file adds pairs
# biomaRt missed. Conversion happens upstream -- convert_mouse_to_human() in
# NeuMap_orthologs.R for NeuMap, while xie2020, GSE165276 and staupe2022 arrive
# converted and have no producer in this tree. Either way human symbols are in
# the RNA assay and mouse symbols in RNAmice, and RNAmice is the *default* assay,
# so everything below indexes @assays$RNA explicitly. The mapping is kept here
# only to label mouse gene names on the NeuMap gene UMAPs.
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

library(Matrix)
library(Seurat)

signaturegenes.atlas = signaturegenes
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
dtf.neumap$facetvar = "subsets"

p.umap.celltype =
  plot_scatter(dtf.neumap, x = "UMAP_1", y = "UMAP_2",
               color_by = "celltype", color_type = "discrete", colors = col_celltype,
               facet_nrow = 1,
               point_size = 0.01, point_alpha = 0.3, raster_dpi = 300,
               shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL,
               legend_ncol = 1, legend_point_size = 0.7, label = TRUE,
               label_df = df.label, label_size = pt2mm(textsize), repel = FALSE) +
  guides(color = guide_legend(ncol = 1, title = NULL, label.hjust = 0,
                              override.aes = list(size = 1.2, alpha = 0.7))) +
  theme_expresso(legend_position = c(1.01, 1), legend_justification = c(0, 1),
                 legend_title = FALSE, facet_label_face = "plain",
                 panel_background = "white", plot_background = "white",
                 legend_text_size = textsize, text_size = textsize,
                 legend_key_spacing_x = 0.01, legend_key_spacing_y = 0.01,
                 legend_key_height = 0.1, show_axis = FALSE, grid = "none") +
  theme(legend.title = element_text(margin = margin(b = 0.1, unit = "line")),
        legend.text  = element_text(margin = margin(-0.06, 0, -0.06, 0.1, "line")))


col_tissue2 = c(`peripheral blood` = "#EDD1D8",  `bone marrow` = "#73DAFF",
                spleen = "#8CC269", liver = "#7C1823", lung = "#DC3023", pancreas = "#FF0097",
                gut = "#FFDFB2", peritoneum = "#0E5FDB", skin = "#B0A4E3", brain = "grey17",
                `breast` = "#352A87", heart = "#6A1A99", placenta = "#F8D626")

dtf.neumap$facetvar = "Tissue"
p.umap.tissue =
  plot_scatter(dtf.neumap, x = "UMAP_1", y = "UMAP_2", color_by = "tissue",
               color_type = "discrete",
               colors = col_tissue2, facet_nrow = 1,
               point_size = 0.01, point_alpha = 0.3, raster_dpi = 300,
               shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL,
               legend_ncol = 1, legend_point_size = 1) +
  scale_color_manual(values = col_tissue2, labels = cap_first, name = NULL,
                     na.value = "grey42") +
  guides(color = guide_legend(ncol = 1, title = NULL, label.hjust = 0,
                              override.aes = list(size = 1, alpha = 0.7))) +
  theme_expresso(legend_position = c(1.01, 1), legend_justification = c(0, 1),
                 legend_title = FALSE, facet_label_face = "plain",
                 panel_background = "white", plot_background = "white",
                 legend_text_size = textsize, text_size = textsize,
                 legend_key_spacing_x = 0.01, legend_key_spacing_y = 0.01,
                 legend_key_height = 0.1, show_axis = FALSE, grid = "none") +
  theme(legend.title = element_text(margin = margin(b = 0.1, unit = "line")),
        legend.text  = element_text(margin = margin(-0.06, 0, -0.06, 0.1, "line")))

col_group2 = c(healthy = "#8CC269", flu = "cyan", `acute inflammation` = "#EFC000",
               `lung cancer` = "#DC3023", `breast cancer` = "#352A87",
               `pancreatic cancer` = "#FF0097", pancreatitis = "grey67",
               stroke = "#0E5FDB", infarction = "#7C1823",
               peritonitis = "#C5BB5C", fibrosis = "#06A5C7", `biliary damage` = "#7b72c5")

dtf.neumap$facetvar = "Group"
p.umap.group =
  plot_scatter(dtf.neumap, x = "UMAP_1", y = "UMAP_2", color_by = "group", color_type = "discrete",
               colors = col_group2, facet_nrow = 1,
               point_size = 0.01, point_alpha = 0.3, raster_dpi = 300,
               shuffle = TRUE, seed = 42, na_color = "grey42",
               title = NULL,
               legend_ncol = 1, legend_point_size = 1) +
  scale_color_manual(values = col_group2, labels = cap_first, name = NULL,
                     na.value = "grey42") +
  guides(color = guide_legend(ncol = 1, title = NULL, label.hjust = 0,
                              override.aes = list(size = 1, alpha = 0.7))) +
  theme_expresso(legend_position = c(1.01, 1), legend_justification = c(0, 1),
                 legend_title = FALSE, text_size = textsize,
                 facet_label_face = "plain", panel_background = "white",
                 plot_background = "white", legend_text_size = textsize,
                 legend_key_spacing_x = 0.01, legend_key_spacing_y = 0.01,
                 legend_key_height = 0.1, show_axis = FALSE, grid = "none") +
  theme(legend.title = element_text(margin = margin(b = 0.1, unit = "line")),
        legend.text  = element_text(margin = margin(-0.06, 0, -0.06, 0.1, "line")))

### marker dotplot
mat = readRDS(paste0(dir.results, "cross_species/NeuMap.signaturegenes.mat.rds"))
setdiff(signaturegenes, rownames(mat))

dtf.neumap.expr = as.data.table(cbind(dtf.neumap, as.matrix(t(mat))))
dtf.neumap.expr = cbind(dtf.neumap, as.matrix(t(mat)))
dtf.neumap.expr$DEFA4 = NULL # it is all 0 (+1)
dtf.neumap.expr$celltype = factor(dtf.neumap.expr$celltype, levels = rev(names(col_celltype)))

genes.markerdot.neumap = append(signaturegenes.atlas, "HSPA1A",
                                after = match("HSPA1B", signaturegenes.atlas) - 1L)

p.markerdot.neumap = plot_dotplot(
  dtf.neumap.expr[!is.na(celltype)],
  features = setdiff(intersect(genes.markerdot.neumap, rownames(mat)), "DEFA4"),
  group_by = "celltype", max_scale = 2, dot_scale = 2, dot_stroke = 0.05,
  title = NULL, col_fontsize = textsize, row_fontsize = textsize, cluster_features = FALSE,
  cluster_groups = FALSE, show_axis = FALSE, feature_side = "top",
  colorbar_title = "Scaled expression") +
  aes(x = .group, y = feature) +
  scale_x_discrete(position = "top", limits = rev) +
  scale_y_discrete(limits = rev) +
  theme(
    legend.position = c(1, 1), legend.justification = c(0, 1), legend.box.just = "left",
    axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5,
                                   face = "plain", size = textsize,
                                   margin = margin(b = 0.1, unit = "line")),
    axis.text.y = element_text(face = "italic", size = textsize,
                               margin = margin(r = 0.1, unit = "line")))

p.markerdot.neumap$layers[[1]]$aes_params$colour = NULL
p.markerdot.neumap = p.markerdot.neumap + aes(colour = expr) +
  scale_colour_gradientn(colours = expresso_colors("diverging_bwr2"),
                         limits = c(-2, 2), guide = "none")

### selected genes on the NeuMap UMAP, labelled with the mouse symbol

pgenes.NeuMap.mice = list()
for(id in c("STMN1", "MPO", "PRTN3", "LTF", "MMP8", "MMP9", "EGR1", "CXCL2",
            "CCL4", "IL1B", "SLPI", "CD74")){
  id_mice = mouse_genes_human[human == id]$mouse[1]
  dtf1 = data.table(UMAP_1 = dtf.neumap.expr$UMAP_1, UMAP_2 = dtf.neumap.expr$UMAP_2,
                    id = dtf.neumap.expr[[id]], facetvar = id_mice)
  pgenes.NeuMap.mice[[id]] =
    plot_scatter(dtf1, x = "UMAP_1", y = "UMAP_2", color_by = "id",
                 facet_by = "facetvar", color_type = "continuous", colors = c("grey99",
                            scales::dichromat_pal("DarkRedtoBlue.12")(12)[7:12]),
                 quantile_lower = 0.001, quantile_upper = 0.999,
                 shuffle = TRUE, seed = 42, point_size = 0.01,
                 point_alpha = 0.3, raster_dpi = 300, legend_ncol = 1,
                 facet_nrow = 1, title = NULL) +
    theme_expresso(text_size = textsize, legend_position = "none",
                   show_axis = FALSE, grid = "none", facet_label_face = "italic",
                   panel_background = "white", plot_background = "white",
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
  levels = c("peripheral blood", "bone marrow", "spleen", "solid tissue"))
dtf.neumap.prop.m.agg$variable = factor(dtf.neumap.prop.m.agg$variable,
                                        levels = names(col_celltype))
dtf.neumap.prop.m.agg$facetvar = "subset proportions"

p.neumap.distribution =
  ggplot(dtf.neumap.prop.m.agg, aes(x = tissue2, y = value, fill = variable)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  theme_expresso(legend_position = c(1.02,1), legend_text_size = textsize,
                 legend_key_spacing_y = 0.01, legend_key_height = 0.4) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = seq(0, 1, 0.2), position = "left") +
  scale_x_discrete(position = "bottom", labels = cap_first) +
  labs(x = NULL, y = "Subset proportion", title = NULL) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(),
        axis.line.x.bottom  = element_line(color = linecolor, linewidth = linesize),
        axis.line.y.left    = element_line(color = linecolor, linewidth = linesize),
        axis.ticks.x.bottom = element_line(color = linecolor, linewidth = linesize),
        axis.ticks.y.left   = element_line(color = linecolor, linewidth = linesize),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.justification = c(0,1),
        legend.text = element_text(margin = margin(-0.06, 0, -0.06, 0.1, "line")),
        plot.title = element_text(size = textsize, face = "plain", hjust = 0,
                                  margin = margin(t = 0.1, r = 0.1, b = 0.1,
                                                  l = 0.1, unit = "line"))) +
  scale_fill_manual(values = col_celltype, name = "") +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(size = 1.2)))


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
dtf.mahyari2025$sampleid = paste0(dtf.mahyari2025$SubjectId, ".", dtf.mahyari2025$Tissue)
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

p_exact = function(p) paste0(formatC(p, format = "e", digits = 1), "\n")

### bone marrow: human vs NHP vs mouse
dtf.expr.pct = fread(paste0(dir.results, "cross_species/dtf.expr.pct.sample.crossspeciesNHP.csv"))
dtf.expr.pct.m = data.table::melt(
  dtf.expr.pct[cellnumber >= 100], id.vars = c("tissue", "pid", "species", "sampleid"),
  measure.vars = setdiff(colnames(dtf.expr.pct),
                         c("tissue", "pid", "species", "sampleid", "cellnumber")))
dtf.expr.pct.m$species = factor(dtf.expr.pct.m$species, levels = c("human", "NHP", "mouse"))
dtf.expr.pct.m.s = dtf.expr.pct.m[
  variable %in% c("MPO", "CTSG", "PRTN3", "LTF", "LCN2", "MMP8", "MMP9", "CD177")]
dtf.expr.pct.m.s$variable = factor(dtf.expr.pct.m.s$variable,
  levels = c("MPO", "CTSG", "PRTN3", "LTF", "LCN2", "MMP8", "MMP9", "CD177"))

dtf.expr.pct.m[, .(median = median(value)), by = .(variable, species, tissue)]
dtf.expr.pct.m[variable %in% c("LTF","LCN2")][tissue == "bone marrow"]

genes.BM = c("MPO", "CTSG", "PRTN3", "LTF", "LCN2", "MMP8", "MMP9", "CD177")

p.species_compare.BM = list()
for (i in genes.BM) {
  p.species_compare.BM[[i]] = ggplot(
    dtf.expr.pct.m.s[tissue %in% "bone marrow"][variable %in% i],
    aes(x = species, y = value)) +
    ggbeeswarm::geom_quasirandom(aes(color = species), size = 0.2, width = 0.2,
                                 show.legend = TRUE, alpha = 0.5, varwidth = FALSE,
                                 shape = 16) +
    geom_boxplot(color = textcolor, width = boxwidth, size = 0.2,
                 outlier.shape = NA, alpha = 1, fill = NA) +
    labs(x = NULL, y =  paste0("*", i, "*<sup>+</sup>\u2009(%)"), title = NULL) +
    theme_expresso(text_size = textsize, axis_text_size = textsize,
                   axis_title_size = textsize, legend_position = "none",
                   grid = "none", panel_background = "white") +
    scale_y_continuous(breaks = c(0, 0.45, 0.9), labels = function(x) x * 100,
                       expand = expansion(mult = c(0.05, 0.2))) +
    ggsignif::geom_signif(
      comparisons = list(c("human", "NHP"), c("NHP", "mouse"), c("human", "mouse")),
      test = "wilcox.test",
      map_signif_level = p_exact,
      step_increase = c(0, 0.4, 0.43), extend_line = -0.08, tip_length = 0,
      vjust = 0.4, lineheight = 0.6, textsize = pt2mm(5), size = 0.1,
      color = "grey17") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = margin(0.1, 0.05, 0.1, 0.05, "line"),
          axis.title.y = element_markdown(size = textsize, color = textcolor,
                                          margin = margin(0, 0, 0, 0, "pt")),
          axis.text.y = element_text(margin = margin(0, 1, 0, 0, "pt"))) +
    scale_color_manual(values = c("human" = "#264EFF", "NHP" = "#2BAE85",
                                  "mouse" = "#F47B00")) +
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 1.5)))
}

### peripheral blood, and lung or liver: human vs mouse
dtf.expr.pct = fread(paste0(dir.results, "cross_species/dtf.expr.pct.sample.crossspecies.csv"))
dtf.expr.pct.m = data.table::melt(
  dtf.expr.pct[cellnumber >= 100], id.vars = c("tissue", "pid", "species", "sampleid"),
  measure.vars = setdiff(colnames(dtf.expr.pct),
                         c("tissue", "pid", "species", "sampleid", "cellnumber")))
dtf.expr.pct.m$species = factor(dtf.expr.pct.m$species, levels = c("human", "NHP", "mouse"))

dtf.expr.pct.m.s = dtf.expr.pct.m[tissue %in% c("peripheral blood")][
  variable %in% c("S100A4", "SLPI")]
dtf.expr.pct.m.s$variable = factor(dtf.expr.pct.m.s$variable, levels = c("S100A4", "SLPI"))

genes.PB = c( "S100A4", "SLPI")

p.species_compare.PB = list()
for (i in genes.PB) {
  p.species_compare.PB[[i]] = ggplot(
    dtf.expr.pct.m.s[tissue %in% "peripheral blood"][variable %in% i],
    aes(x = species, y = value)) +
    ggbeeswarm::geom_quasirandom(aes(color = species), size = 0.05, width = 0.2,
                                 show.legend = TRUE, alpha = 0.3, varwidth = FALSE,
                                 shape = 16) +
    geom_boxplot(color = textcolor, width = boxwidth, size = 0.2,
                 outlier.shape = NA, alpha = 1, fill = NA) +
    labs(x = NULL, y = paste0("*", i, "*<sup>+</sup>\u2009(%)"), title = NULL) +
    theme_expresso(text_size = textsize, axis_text_size = textsize,
                   axis_title_size = textsize, legend_position = "none",
                   grid = "none", panel_background = "white") +
    scale_y_continuous(breaks = c(0, 0.45, 0.9), labels = function(x) x * 100,
                       expand = expansion(mult = c(0.05, 0.18))) +
    ggsignif::geom_signif(
      comparisons = list(c("human", "mouse")), test = "wilcox.test",
      map_signif_level = p_exact,
      tip_length = 0, vjust = 0.4, lineheight = 0.6, textsize = pt2mm(5), size = 0.1,
      color = "grey17") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = margin(0.1, 0.05, 0.1, 0.05, "line"),
          axis.title.y = element_markdown(size = textsize, color = textcolor,
                                          margin = margin(0, 0, 0, 0, "pt")),
          axis.text.y = element_text(margin = margin(0, 1, 0, 0, "pt"))) +
    scale_color_manual(values = c("human" = "#264EFF", "NHP" = "#2BAE85",
                                  "mouse" = "#F47B00"))
}

dtf.expr.pct.m.s = dtf.expr.pct.m[tissue %in% c( "liver", "lung")][
  variable %in% c("CXCL1", "CXCL2", "IL1B", "CCL4")]
dtf.expr.pct.m.s$variable = factor(
  dtf.expr.pct.m.s$variable, levels = c("CXCL1", "CXCL2", "IL1B", "CCL4"))

genes.liverorlung = c("CXCL1", "CXCL2", "IL1B", "CCL4")

p.species_compare.liverorlung = list()
for (i in genes.liverorlung) {
  p.species_compare.liverorlung[[i]] = ggplot(
    dtf.expr.pct.m.s[variable %in% i], aes(x = species, y = value)) +
    ggbeeswarm::geom_quasirandom(aes(color = species), size = 0.05, width = 0.2,
                                 show.legend = TRUE, alpha = 0.3, varwidth = FALSE,
                                 shape = 16) +
    geom_boxplot(color = textcolor, width = boxwidth, size = 0.2,
                 outlier.shape = NA, alpha = 1, fill = NA) +
    labs(x = NULL, y = paste0("*", i, "*<sup>+</sup>\u2009(%)"), title = NULL) +
    theme_expresso(text_size = textsize, axis_text_size = textsize,
                   axis_title_size = textsize, legend_position = "none",
                   grid = "none", panel_background = "white") +
    scale_y_continuous(breaks = c(0, 0.45, 0.9), labels = function(x) x * 100,
                       expand = expansion(mult = c(0.05, 0.18))) +
    ggsignif::geom_signif(
      comparisons = list(c("human", "mouse")), test = "wilcox.test",
      map_signif_level = p_exact,
      tip_length = 0, vjust = 0.4, lineheight = 0.6, textsize = pt2mm(5), size = 0.1,
      color = "grey17") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = margin(0.1, 0.05, 0.1, 0.05, "line"),
          axis.title.y = element_markdown(size = textsize, color = textcolor,
                                          margin = margin(0, 0, 0, 0, "pt")),
          axis.text.y = element_text(margin = margin(0, 1, 0, 0, "pt"))) +
    scale_color_manual(values = c("human" = "#264EFF", "NHP" = "#B29DDA",
                                  "mouse" = "#F47B00"))
}


### species key, drawn as its own object because the panels above carry no legend
### (their x axis text is blanked, so species is otherwise unlabelled)
col_species = c(Human = "#264EFF", NHP = "#2BAE85", Mouse = "#F47B00")

df_leg.species = data.frame(
  label = factor(names(col_species), levels = names(col_species)), color = col_species)

p.speciesnames = ggplot(df_leg.species, aes(x = 1, y = label)) +
  geom_point(size = 0.5, aes(color = label)) +
  geom_text(aes(x = 1.1, label = label), hjust = 0, size = pt2mm(textsize),
            color = textcolor) +
  scale_color_manual(values = col_species) +
  scale_x_continuous(limits = c(0.95, 2)) +
  theme_void() +
  theme(
    legend.position = "none", plot.margin = margin(0, 0, 0, 0)) +
  scale_y_discrete(limits = rev(names(col_species)))


# -- NeuMap and published mouse subset scores in SCAHN -------------------------

hubgenes = readRDS(paste0(dir.results, "cross_species/SCAHN.neumap.scores.rds"))
colnames(hubgenes)

colnames(hubgenes) = c("cell.pid", "UMAP_1", "UMAP_2",
                       "Immature\n(Cebpe+Mmp8+)", "IFN response\n(Ifit1+Cd274+)",
                       "IS-I\n(Cd14+Ptgs2+)", "preNeu\n(Ltf+mKi67+)",
                       "APC\n(H2+Cd74+)", "IS-II\n(Vegfa+Cd274+)",
                       "immuno-silent\n(Cd52+ circulating)",
                       "immuno-silent", "IS-I", "IFN response", "APC&IS-II", "immature", "preNeu",
                       "G0", "G1", "G2", "G3", "G4",
                       "G5a", "G5b", "G5c", "hG5a", "hG5b", "hG5c", "P1", "P2", "P3",
                       "P4", "neutrotime.pos", "neutrotime.neg", "T3", "T2", "T1")

hubgenes$celltype = dtf$celltype

### NeuMap hub-gene scores on the SCAHN UMAP

italic_hubgenes = function(hub) {
  parts = strsplit(hub, "\n", fixed = TRUE)[[1]]
  if (parts[1] == "immuno-silent") return(paste0("Immuno-silent", "<br>"))
  genes = gsub("^\\(|\\)$", "", parts[2])
  paste0(parts[1], "<br>", gsub("([A-Za-z][A-Za-z0-9-]*)\\+", "*\\1*+", genes))
}

pscores.NeuMap.mice = list()
for(hub in c("preNeu\n(Ltf+mKi67+)", "Immature\n(Cebpe+Mmp8+)",
             "IFN response\n(Ifit1+Cd274+)", "immuno-silent\n(Cd52+ circulating)",
             "IS-I\n(Cd14+Ptgs2+)", "IS-II\n(Vegfa+Cd274+)", "APC\n(H2+Cd74+)")){
  dtf1 = data.table(UMAP_1 = hubgenes$UMAP_1, UMAP_2 = hubgenes$UMAP_2,
                    id = hubgenes[[hub]] - 1, facetvar = italic_hubgenes(hub))
  pscores.NeuMap.mice[[hub]] =
    plot_scatter(dtf1, x = "UMAP_1", y = "UMAP_2", color_by = "id",
                 facet_by = "facetvar", color_type = "continuous", colors = c("grey99",
                            scales::dichromat_pal("DarkRedtoBlue.12")(12)[7:12]),
                 quantile_lower = 0.001, quantile_upper = 0.999,
                 shuffle = TRUE, seed = 42, point_size = 0.000001,
                 point_alpha = 0.2, raster_dpi = 400, legend_ncol = 1,
                 facet_nrow = 1, title = NULL) +
    theme_expresso(text_size = textsize, legend_position = "none",
                   show_axis = FALSE, grid = "none", facet_label_face = "plain",
                   panel_background = "white", plot_background = "white",
                   plot_margin = ggplot2::margin(0.05, 0.05, 0.05, 0.05, "line")) +
    theme(strip.text = element_markdown(size = textsize, color = textcolor,
                                        lineheight = 0.8,
                                        margin = margin(0.05, 0, 0.05, 0, "line")))
}

### published mouse subset scores on the SCAHN UMAP
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
    theme_expresso(text_size = textsize, legend_position = "none",
                   show_axis = FALSE, grid = "none", facet_label_face = "plain",
                   panel_background = "white", plot_background = "white",
                   plot_margin = ggplot2::margin(0.01, 0.01, 0.01, 0.01, "line"))
}


# -- Cross-study subset correspondence table -----------------------------------
dt.correspond = data.table(
  "SCAHN" = c("AZU1", "LTF", "IFN1\u2013IFN3", "Cytokine/\ninflammation", "CD74"),
  "Xie2020" = c("G0\u2013G2", "G3", "G5c", "\u2013", "\u2013"),
  "Grieshaber-\nbouyer2021" = c("P1", "P2", "\u2013", "\u2013", "\u2013"),
  "Ng2024" = c("\u2013", "\u2013", "\u2013", "T3", "\u2013"),
  "NeuMap" = c("preNeu\n(Ltf+mKi67+)", "Immature\n(Cebpe+Mmp8+)",
               "IFN response\n(Ifit1+Cd274+)",
               "IS-I (Cd14+Ptgs2+),\nIS-II (Vegfa+Cd274+)", "APC\n(H2+Cd74+)"))

dt.correspond = tableGrob(
  dt.correspond, rows = NULL, theme = ttheme_default(
    core = list(
      fg_params = list(fontsize = textsize, fontface = "plain", lineheight = 0.8),
      bg_params = list(fill = c("white")), padding = unit(c(0.3, 0.4), "line")), colhead = list(
        fg_params = list(fontsize = textsize, fontface = "plain", lineheight = 0.8),
        bg_params = list(fill = c("white")), padding = unit(c(0.3, 0.4), "line"))))

### gene symbols in the NeuMap column set in italics
md_genes = function(x)
  gsub("\n", "<br>", gsub("([A-Za-z][A-Za-z0-9-]*)\\+", "*\\1*+", x), fixed = TRUE)

for (k in which(dt.correspond$layout$name == "core-fg" &
                dt.correspond$layout$l == ncol(dt.correspond))) {
  dt.correspond$grobs[[k]] = gridtext::richtext_grob(
    md_genes(dt.correspond$grobs[[k]]$label),
    gp = gpar(fontsize = textsize, col = textcolor, lineheight = 0.8),
    padding = unit(rep(0, 4), "pt"), box_gp = gpar(col = NA))
}

hrule = function(y) segmentsGrob(x0 = 0, x1 = 1, y0 = y, y1 = y,
                                 gp = gpar(lwd = 0.4, col = "grey42"))

for (i in seq_len(nrow(dt.correspond))) {
  dt.correspond = gtable_add_grob(
    dt.correspond, grobs = hrule(0), t = i, b = i, l = 1, r = ncol(dt.correspond),
    name = paste0("hrule.bottom.", i))
}
dt.correspond = gtable_add_grob(
  dt.correspond, grobs = hrule(1), t = 1, b = 1, l = 1, r = ncol(dt.correspond),
  name = "hrule.top")


# -- SCAHN-to-NeuMap subset correspondence -------------------------------------

### hand-curated mapping between SCAHN subsets and NeuMap subsets (Figure 3B)
mappings =  c(preNeu = "AZU1", Immature = "LTF", Immature = "MMP9",
              `Other` = "MMP9", `Other` = "S100A4",
              `Other` = "IL1R2", `IS-I` = "EGR1", `IS-I` = "AP-1",
              `IS-I` = "PTGS2", `IS-I` = "G0S2", `IFN response` = "IFN1", `IFN response` = "IFN2",
              `IFN response` = "IFN3",
              `IS-I` = "NF-κB", `IS-I` = "IL1B", `IS-I` = "IL1RN",
              `IS-II` = "CXCL", `IS-II` = "VEGFA", `IS-II` = "CCL3/4",
              `IS-II` = "NF-κB", `IS-II` = "IL1B", `IS-II` = "IL1RN",
              `Other`  = "SLPI", `IS-II` = "HSP", `APC` = "HSP",
              APC = "CD74",
              `Other` = "MME", `Other` = "TXNIP")

df = data.frame(SCAHN = mappings, NeuMap = names(mappings))

scahn_colors = c(AZU1 = "#DC3023", LTF = "#F47B00", MMP9 = "#F8D626",
                  S100A4 = "#FFDFB2", IL1R2 = "#acf300",
                  MME = "#80C684", TXNIP = "#D1BA58", EGR1 = "#1A5E1F",
                  `AP-1` = "#2BAE85", PTGS2 = "#AFBC65", G0S2 = "#E2E7BF", IFN1 = "#2E42B8",
                  IFN2 = "#0A73DC", IFN3 = "#264EFF", CXCL = "#E4C6D0", VEGFA = "#73DAFF",
                  `CCL3/4` = "#DCE318", `NF-κB` = "#5E34B1", IL1B = "#9E9AC8",
                  IL1RN = "#DCB0F2", SLPI = "#7CABB1", HSP = "#4D6D93", CD74 = "#C7E5C9")

neumap_colors = c(
  "preNeu" = "#8ABE75", "Immature" = "#F0BB46",
  "IFN response" = "#d69aad", "IS-I" = "#40A1FF",
  "IS-II" = "#f19e73", "APC" = "#D1C4E9", "Other" = "grey99")

all_colors = c(scahn_colors, neumap_colors)

df$SCAHN  = factor(df$SCAHN,  levels = names(scahn_colors))
df$NeuMap = factor(df$NeuMap, levels = names(neumap_colors))
df$freq = 1 / ave(rep(1, nrow(df)), df$SCAHN, FUN = length)

p.alluvial = ggplot(df, aes(axis1 = SCAHN, axis2 = NeuMap, y = freq)) +
  geom_flow(aes(fill = NeuMap), width = 0.7) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.7, color = "white",
               linewidth = 0.01, size = 2) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = pt2mm(textsize), fontface = "plain", color = textcolor, lineheight = 0.8) +
  scale_x_discrete(limits = c("SCAHN", "NeuMap"), expand = c(0, 0), position = "top") +
  scale_fill_manual(values = all_colors, guide = "none") +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(color = "white", fill = "white"),
    panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
    axis.text.x.top = element_text(size = textsize, face = "plain", color = textcolor,
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
table(dtf.mahyari2025$celltype, dtf.mahyari2025$tissue)

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
dtf.mahyari2025$facetvar = "subsets"

p.umap.celltype.mahyari2025 =
  plot_scatter(dtf.mahyari2025, x = "UMAP_1", y = "UMAP_2",
               color_by = "celltype", color_type = "discrete", colors = col_celltype,
               facet_nrow = 1,
               point_size = 0.1, point_alpha = 0.3, raster_dpi = 300,
               shuffle = TRUE, seed = 42, na_color = "grey42",
               title = NULL,
               legend_ncol = 1, legend_point_size = 1.2, label = TRUE,
               label_df = df.label, label_size = pt2mm(textsize), repel = FALSE) +
  guides(color = guide_legend(ncol = 1, title = NULL, label.hjust = 0,
                              override.aes = list(size = 1.2, alpha = 0.7))) +
  theme_expresso(legend_position = c(0.98, 1), legend_justification = c(0, 1),
                 legend_title = FALSE, panel_background = "white",
                 facet_label_face = "plain", legend_text_size = textsize,
                 text_size = textsize, legend_key_spacing_x = 0.01,
                 legend_key_spacing_y = 0.01, legend_key_height = 0.1,
                 show_axis = FALSE, grid = "none") +
  theme(legend.title = element_text(margin = margin(b = 0.1, unit = "line")),
        legend.text  = element_text(margin = margin(-0.06, 0, -0.06, 0.1, "line")))

### selected genes on the mahyari2025 UMAP (Figure 3H)
mat = readRDS(paste0(dir.results, "cross_species/mahyari2025.signaturegenes.mat.rds"))
dtf.mahyari2025.expr = cbind(dtf.mahyari2025, as.matrix(t(mat)))

pgenes.mahyari2025 = list()
for(id in c("STMN1", "MPO", "PRTN3", "LTF", "MMP8", "MMP9")){
  dtf1 = data.table(UMAP_1 = dtf.mahyari2025.expr$UMAP_1,
                    UMAP_2 = dtf.mahyari2025.expr$UMAP_2,
                    id = dtf.mahyari2025.expr[[id]], facetvar = id)
  pgenes.mahyari2025[[id]] =
    plot_scatter(dtf1, x = "UMAP_1", y = "UMAP_2",
                 color_by = "id", facet_by = "facetvar", color_type = "continuous",
                 colors = c("grey99", scales::dichromat_pal("DarkRedtoBlue.12")(12)[7:12]),
                 quantile_lower = 0.001, quantile_upper = 0.999,
                 shuffle = TRUE, seed = 42, point_size = 0.01,
                 point_alpha = 0.3, raster_dpi = 300, legend_ncol = 1,
                 facet_nrow = 1, title = NULL) +
    theme_expresso(text_size = textsize, legend_position = "none",
                   show_axis = FALSE, grid = "none", facet_label_face = "italic",
                   panel_background = "white", plot_background = "white",
                   plot_margin = ggplot2::margin(0.02, 0.02, 0.02, 0.02, "line"))
}


# -- Extended figure 4 ---------------------------------------------------------

draw_hline_label = function(line_rows, cols, label, gap = 2,
                            text_gap = 0.7,           # label to rule, in mm
                            margin = c(0, 0, 0, 0),   # top, right, bottom, left in npc
                            lwd = 0.7, line_col = "grey42",
                            textsize = 8, textface = "plain", textcol = textcolor) {
  text_row = min(line_rows) - gap
  all_rows = text_row:max(line_rows)
  n_all    = length(all_rows)
  n_line   = length(line_rows)

  mt = margin[1]; mr = margin[2]; mb = margin[3]; ml = margin[4]

  # Outer viewport covers the full row/col region
  pushViewport(viewport(layout.pos.row = all_rows, layout.pos.col = cols))
  # Inner viewport applies margins
  pushViewport(viewport(
    x = unit(ml, "npc"), y = unit(mb, "npc"), width = unit(1 - ml - mr, "npc"),
    height = unit(1 - mt - mb, "npc"), just = c("left", "bottom")))

  line_y = (n_line / 2) / n_all

  grid.lines(x = c(0, 1), y = c(line_y, line_y), gp = gpar(lwd = lwd, col = line_col))
  grid.text(unit(0.5, "npc"), unit(line_y, "npc") + unit(text_gap, "mm"), label = label,
            just = c("centre", "bottom"),
            gp = gpar(fontsize = textsize, fontface = textface, col = textcol))

  upViewport(2)
}

cairo_pdf(file = paste0(dir.fig, "Fig.e4.pdf"), width = 5.8, height = 8.8)
pushViewport(viewport(layout=grid.layout(nrow = 254, ncol = 380)))

### A: xie2020 G / hG subset scores
pushViewport(viewport(layout.pos.row = 11:55, layout.pos.col = 12:200))
grid.arrange(grobs=pscores.mice[c(paste0("G",c(0:4, "5a", "5b", "5c")))],
             nrow = 2, as.table = TRUE, newpage = FALSE)
upViewport(1)

pushViewport(viewport(layout.pos.row = 65:87, layout.pos.col = 12:153))
grid.arrange(grobs=pscores.mice[c( "hG5a", "hG5b", "hG5c")],
             nrow = 1, as.table = TRUE, newpage = FALSE)
upViewport(1)

# ### B: grieshaber-bouyer2021 P subset scores
pushViewport(viewport(layout.pos.row = 97:119, layout.pos.col = 12:200))
grid.arrange(grobs=pscores.mice[paste0("P",c(1:4))], nrow = 1, as.table = TRUE, newpage = FALSE)
upViewport(1)

### C: ng2024 T subset scores
pushViewport(viewport(layout.pos.row = 129:151, layout.pos.col = 12:153))
grid.arrange(grobs=pscores.mice[paste0("T",c(1:3))], nrow = 1, as.table = TRUE, newpage = FALSE)
upViewport(1)

draw_hline_label(line_rows = 7:10, cols = 22:190, label = "Xie2020, mouse", gap = 2,
                 margin = c(0.01, 0.01, 0.01, 0.01), textsize = textsize, textcol = textcolor)
draw_hline_label(line_rows = 61:64, cols = 22:143, label = "Xie2020, human", gap = 2,
                 margin = c(0.01, 0.01, 0.01, 0.01), textsize = textsize, textcol = textcolor)

draw_hline_label(line_rows = 93:96, cols = 22:190,
                 label = "Grieshaber-bouyer2021, mouse", gap = 2,
                 margin = c(0.01, 0.01, 0.01, 0.01), textsize = textsize, textcol = textcolor)
draw_hline_label(line_rows = 125:128, cols = 22:143, label = "Ng2024, mouse", gap = 2,
                 margin = c(0.01, 0.01, 0.01, 0.01), textsize = textsize, textcol = textcolor)

### E: cross-study subset correspondence table
pushViewport(viewport(layout.pos.row = 172:218, layout.pos.col = 12:200))
grid.draw(dt.correspond)
popViewport()

### F: NeuMap marker dotplot
print(p.markerdot.neumap, vp = viewport(layout.pos.row = 3:254, layout.pos.col = 202:350))

print(p.umap.axis + theme(plot.margin = margin(0, 0, 0, 0.1, "line")),
      vp = viewport(layout.pos.row = 48:60, layout.pos.col = 1:30))
print(p.umap.axis + theme(plot.margin = margin(0, 0, 0, 0.1, "line")),
      vp = viewport(layout.pos.row = 80:92, layout.pos.col = 1:30))
print(p.umap.axis + theme(plot.margin = margin(0, 0, 0, 0.1, "line")),
      vp = viewport(layout.pos.row = 112:124, layout.pos.col = 1:30))
print(p.umap.axis + theme(plot.margin = margin(0, 0, 0, 0.1, "line")),
      vp = viewport(layout.pos.row = 144:156, layout.pos.col = 1:30))

## panel labels
grid.text(x = unit(0.01,"npc"), y=unit(0.98,"npc"), label="a",
          gp=gpar(fontsize = panellabelsize,fontface = "bold", col = textcolor))

grid.text(x = unit(0.01,"npc"), y=unit(0.75,"npc"), label="b",
          gp=gpar(fontsize = panellabelsize,fontface = "bold", col = textcolor))

grid.text(x = unit(0.01,"npc"), y=unit(0.63,"npc"), label="c",
          gp=gpar(fontsize = panellabelsize,fontface = "bold", col = textcolor))

grid.text(x = unit(0.01,"npc"), y=unit(0.5,"npc"), label="d",
          gp=gpar(fontsize = panellabelsize,fontface = "bold", col = textcolor))


grid.text(x = unit(0.01,"npc"), y=unit(0.34,"npc"), label="e",
          gp=gpar(fontsize = panellabelsize,fontface = "bold", col = textcolor))

grid.text(x = unit(0.55,"npc"), y=unit(0.98,"npc"), label="f",
          gp=gpar(fontsize = panellabelsize,fontface = "bold", col = textcolor))

dev.off()


# -- Figure 3 ------------------------------------------------------------------
cairo_pdf(file = paste0(dir.fig, "Figure3.pdf"), width = 5.2, height = 7.3)
pushViewport(viewport(layout = grid.layout(nrow = 182, ncol = 141)))

### A: NeuMap hub-gene scores in SCAHN
pushViewport(viewport(layout.pos.row = 3:23, layout.pos.col = 6:140))
grid.arrange(grobs=pscores.NeuMap.mice, nrow = 1, as.table = TRUE, newpage = FALSE)
upViewport(1)

### C-E: NeuMap projected onto SCAHN
print(p.umap.group, vp = viewport(layout.pos.row = 31:57, layout.pos.col = 41:71))
print(p.umap.tissue, vp = viewport(layout.pos.row = 31:57, layout.pos.col =  93:119))

print(p.umap.celltype,
      vp = viewport(layout.pos.row = 59:101, layout.pos.col = 39:85))
print(p.neumap.distribution, vp = viewport(layout.pos.row = 60:108, layout.pos.col = 103:126))

### B: SCAHN-to-NeuMap subset correspondence
print(p.alluvial, vp = viewport(layout.pos.row = 31:108, layout.pos.col = 4:36))

# ### F: marker genes on the NeuMap UMAP
pushViewport(viewport(layout.pos.row = 110:121, layout.pos.col = 5:133))
grid.arrange(grobs=pgenes.NeuMap.mice, nrow = 1, as.table = TRUE, newpage = FALSE)
upViewport(1)

# ### H: mahyari2025 (NHP) projected onto SCAHN
print(p.umap.celltype.mahyari2025, vp = viewport(layout.pos.row = 125:167, layout.pos.col = 6:50))

pushViewport(viewport(layout.pos.row = 166:177, layout.pos.col = 5:60))
grid.arrange(grobs=pgenes.mahyari2025, nrow = 1, as.table = TRUE, newpage = FALSE)
upViewport(1)

# ### I: percent-expressing across species

pushViewport(viewport(layout.pos.row = 126:181, layout.pos.col = 62:100))
grid.arrange(grobs = p.species_compare.BM, nrow = 4, as.table = TRUE, newpage = FALSE)
upViewport(1)

pushViewport(viewport(layout.pos.row = 126:181, layout.pos.col = 104:120))
grid.arrange(grobs = p.species_compare.liverorlung, nrow = 4, as.table = TRUE, newpage = FALSE)
upViewport(1)

pushViewport(viewport(layout.pos.row = 154:181, layout.pos.col = 124:140))
grid.arrange(grobs = p.species_compare.PB, nrow = 2, as.table = TRUE, newpage = FALSE)
upViewport(1)

draw_hline_label(line_rows = 124:125, cols = 68:100, label = "Bone marrow", gap = 2,
                 margin = c(0.01, 0.01, 0.01, 0.01), textsize = textsize, textcol = textcolor)
draw_hline_label(line_rows = 124:125, cols = 110:120, label = "Liver/Lung", gap = 2,
                 margin = c(0.01, 0.01, 0.01, 0.01), textsize = textsize, textcol = textcolor)
draw_hline_label(line_rows = 152:153, cols = 130:140, label = "Peripheral\nblood", gap = 2,
                 text_gap = 0.4, margin = c(0.01, 0.01, 0.01, 0.01), textsize = textsize,
                 textcol = textcolor)

print(p.speciesnames, vp = viewport(layout.pos.row = 123:130, layout.pos.col = 124:140))


print(p.umap.axis + theme(plot.margin = margin(0, 0, 0, 0.1, "line")),
      vp = viewport(layout.pos.row = 171:182, layout.pos.col = 1:12))
print(p.umap.axis + theme(plot.margin = margin(0, 0, 0, 0.1, "line")),
      vp = viewport(layout.pos.row = 155:166, layout.pos.col = 1:12))
print(p.umap.axis + theme(plot.margin = margin(0, 0, 0, 0.1, "line")),
      vp = viewport(layout.pos.row = 114:125, layout.pos.col = 1:12))
print(p.umap.axis + theme(plot.margin = margin(0, 0, 0, 0.1, "line")),
      vp = viewport(layout.pos.row = 17:28, layout.pos.col = 1:12))


print(p.umap.axis + theme(plot.margin = margin(0, 0, 0, 0.1, "line")),
      vp = viewport(layout.pos.row = 91:102, layout.pos.col = 38:49))
print(p.umap.axis + theme(plot.margin = margin(0, 0, 0, 0.1, "line")),
      vp = viewport(layout.pos.row = 49:60, layout.pos.col = 38:49))
print(p.umap.axis + theme(plot.margin = margin(0, 0, 0, 0.1, "line")),
      vp = viewport(layout.pos.row = 49:60, layout.pos.col = 90:101))


### panel labels
grid.text(x = unit(0.01,"npc"), y=unit(0.99,"npc"), label = "a",
          gp=gpar(fontsize = panellabelsize,fontface = "bold", col = textcolor))

grid.text(x = unit(0.01,"npc"), y=unit(0.83,"npc"), label = "b",
          gp=gpar(fontsize = panellabelsize,fontface = "bold", col = textcolor))
grid.text(x = unit(0.28,"npc"), y=unit(0.83,"npc"), label = "c",
          gp=gpar(fontsize = panellabelsize,fontface = "bold", col = textcolor))
grid.text(x = unit(0.66,"npc"), y=unit(0.83,"npc"), label = "d",
          gp=gpar(fontsize = panellabelsize,fontface = "bold", col = textcolor))

grid.text(x = unit(0.28,"npc"), y=unit(0.66,"npc"), label = "e",
          gp=gpar(fontsize = panellabelsize,fontface = "bold", col = textcolor))
grid.text(x = unit(0.74,"npc"), y=unit(0.67,"npc"), label = "f",
          gp=gpar(fontsize = panellabelsize,fontface = "bold", col = textcolor))
grid.text(x = unit(0.01,"npc"), y=unit(0.41,"npc"), label = "g",
          gp=gpar(fontsize = panellabelsize,fontface = "bold", col = textcolor))

grid.text(x = unit(0.02,"npc"), y=unit(0.3,"npc"), label = "h",
          gp=gpar(fontsize = panellabelsize,fontface = "bold", col = textcolor))
grid.text(x = unit(0.45,"npc"), y=unit(0.33,"npc"), label = "i",
          gp=gpar(fontsize = panellabelsize,fontface = "bold", col = textcolor))

dev.off()
