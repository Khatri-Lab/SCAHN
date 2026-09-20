# ==============================================================================
# Atlas overview
#
# UMAPs, alluvial flow, subset distribution, QC, trajectory and marker dotplots
# for the SCAHN neutrophil atlas and for the peripheral-blood immune reference.
#
# Requires 00_setup.R (packages, figure settings, palettes) and 01_load_data.R
# (dtf, dtf.prop, signaturegenes, genes.*, dir.data / dir.results / dir.fig).
#
# Reads, on top of what 01_load_data.R loads:
#   results/annotation/SCAHN.selectedgenes.mat.expr.integrated.rds
#   results/annotation/SCAHN.cellscores.rds
#   results/annotation/SCAHN.slingshot.res.rds
#   data/peripheral_blood/blood.cellmeta.csv
#   data/peripheral_blood/blood.mat.selectedgenes.rds
#   data/peripheral_blood/blood.expr.sum.celltype.rds
#
# Writes to figures/original/:
#   fig1.umap_pid.pdf                Figure 1
#   fig1.umap_tissue_group_pct.pdf   Figure 1
#   Figure2.pdf                      Figure 2
#   Fig.e2.pdf                       Extended figure 2 (QC, trajectory)
#   Fig.s1.pdf                       Supplementary figure 1 (QC)
#   Fig.e3.pdf                       Extended figure 3 (peripheral blood)
# ==============================================================================

SCAHN_SCRIPTS = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/scripts/figures"
source(file.path(SCAHN_SCRIPTS, "00_setup.R"))
source(file.path(SCAHN_SCRIPTS, "01_load_data.R"))

# The subset named for sequencing depth is keyed "lowDepth" upstream and drawn
# "Low depth". Renamed on the keys themselves rather than at each label, so it
# reaches every panel at once: they all read the string off dtf$celltype,
# dtf.prop's columns or names(col_celltype), and cap_first() leaves it alone.
# Local to this script -- 00_setup.R and 01_load_data.R are shared, and
# 05_sexdifference.R still subsets on "lowDepth".
names(col_celltype)[names(col_celltype) == "lowDepth"] = "Low depth"
dtf[celltype == "lowDepth", celltype := "Low depth"]
setnames(dtf.prop, "lowDepth", "Low depth")


# -- SCAHN UMAPs and ggalluvial ------------------------------------------------

### tissue-by-group alluvial
dtf.summary = as.data.table(table(dtf$tissue, dtf$group))
dtf.summary = dtf.summary[N > 0]
colnames(dtf.summary)[1:2] = c("tissue", "group")
dtf.summary$tissue = factor(dtf.summary$tissue, levels = rev(names(col_tissue)))
dtf.summary$group  = factor(dtf.summary$group,  levels = rev(names(col_group)))

cap_tissue = function(x) sub("^Colonrectum$", "Colorectum", cap_first(x))

col_tissue_group = c(col_tissue, col_group)

dtf.summary_long = to_lodes_form(data.frame(dtf.summary),
                                 key = "category", value = "item", id = "index", axes = 2:1)
dtf.summary_long = as.data.table(dtf.summary_long)
dtf.summary_long[, colflow := ifelse(category == "group", item, NA)]

p.ggalluvial = ggplot(
  data = dtf.summary_long, aes(x = category, stratum = item, alluvium = index, y = N)) +
  geom_flow(aes(fill = item), width = 1/5) +
  geom_stratum(aes(fill = item), width = 1/5, size = .1) +
  scale_fill_manual(values = col_tissue_group) +
  theme_expresso_void(legend_position = "none") +
  coord_flip() +
  theme(plot.margin = unit(c(-0.8, -1.1, -0.8, -1.1), "line")) +
  labs(x = NULL, y = NULL)

### study
dtf$facetvar = "Study"
dtf$study = gsub("xue2022", "xue2023", dtf$pid)
dtf$study = gsub("myin2023", "myin2024", dtf$study)
dtf$study = gsub("reyfman2018", "reyfman2019", dtf$study)

p.umap.pid = plot_scatter(
  dtf, x = "UMAP_1", y = "UMAP_2", color_by = "study", color_type = "discrete", colors = col_study,
  facet_by = "facetvar", facet_nrow = 1, point_size = 0.01, point_alpha = 0.3, raster_dpi = 500,
  legend_ncol = 2, legend_point_size = 1.5,
  shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL) +
  theme_expresso(legend_position = c(1.03, 1.18), legend_justification = c(0, 1),
                 legend_key_width = 0.1, legend_key_height = 0.1,
                 legend_key_spacing_x = 0.01, legend_key_spacing_y = 0.03,
                 legend_text_size = 6, facet_label_face = "plain", show_axis = FALSE, grid = "none")

### tissue
dtf$facetvar = "Tissue"

p.umap.tissue = plot_scatter(
  dtf, x = "UMAP_1", y = "UMAP_2",
  color_by = "tissue", color_type = "discrete", colors = col_tissue,
  facet_by = "facetvar", facet_nrow = 1, point_size = .005, point_alpha = 0.3, raster_dpi = 500,
  legend_ncol = 4, legend_point_size = 2,
  shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL) +
  # Drawn labels only: adding a scale replaces the one plot_scatter() set, so its
  # values, name and na.value are repeated here.
  scale_color_manual(values = col_tissue, labels = cap_tissue, name = NULL,
                     na.value = "grey42") +
  theme_expresso(legend_position = c(1.01, 0.5), legend_justification = c(0, 0.5),
                 legend_key_spacing_x = 0.01, legend_key_spacing_y = 0.05,
                 facet_label_face = "plain", show_axis = FALSE, grid = "none")

### condition
dtf$facetvar = "Condition"

p.umap.group = plot_scatter(
  dtf, x = "UMAP_1", y = "UMAP_2", color_by = "group", color_type = "discrete", colors = col_group,
  facet_by = "facetvar", facet_nrow = 1, point_size = .005, point_alpha = 0.3, raster_dpi = 500,
  legend_ncol = 3, legend_point_size = 2,
  shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL) +
  # cap_severity(), not cap_first(): the qualifier trails the condition in the
  # data and leads it in the manuscript.
  scale_color_manual(values = col_group, labels = cap_severity, name = NULL,
                     na.value = "grey42") +
  theme_expresso(legend_position = c(1.01, 0.5), legend_justification = c(0, 0.5),
                 legend_key_spacing_x = 0.01, legend_key_spacing_y = 0.05,
                 facet_label_face = "plain", show_axis = FALSE, grid = "none")

### platform, whole blood only
dtf[, technique := ifelse(grepl("Chromium", technique), "10x Chromium", technique)]
dtf[, technique2 := ifelse(technique %in% c("10x Chromium", "BD Rhapsody"), technique, "other")]

p.umap.platform.wb = plot_scatter(
  dtf[sampleid.pid %in% dtf.prop[source2 %in% "fresh\nwhole blood"]$sampleid.pid],
  x = "UMAP_1", y = "UMAP_2", color_by = "celltype", color_type = "discrete", colors = col_celltype,
  facet_by = "technique2", facet_nrow = 3, point_size = 0.01, point_alpha = 0.3, raster_dpi = 300,
  legend_ncol = 1, legend_point_size = 2,
  shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL) +
  theme_expresso(legend_position = c(1, 0.5), legend_justification = c(0, 0.5),
                 legend_key_spacing_x = 0.01, legend_key_spacing_y = 0.01,
                 legend_text_size = 9, facet_label_face = "plain", show_axis = FALSE, grid = "none")

### neutrophil subsets
dtf$celltype = factor(dtf$celltype, levels = names(col_celltype))

df.label = as.data.table(dtf)[, lapply(.SD, median), by = celltype, .SDcols = c("UMAP_1", "UMAP_2")]
colnames(df.label)[2:3] = c("x", "y")
df.label$labeltext = df.label[, "celltype", with = FALSE]
df.label[, x := ifelse(celltype == "IL1B",  x + 1,   x)]
df.label[, x := ifelse(celltype == "IL1RN", x + 0.5, x)]
df.label[, y := ifelse(celltype == "IL1RN", y + 0.5, y)]
df.label[, y := ifelse(celltype == "IL1R2", y + 1,   y)]
df.label[, y := ifelse(celltype == "PTGS2", y - 1,   y)]

dtf$facetvar = NULL
dtf$facetvar = "Neutrophil subset"

p.umap.celltype = plot_scatter(
  dtf, x = "UMAP_1", y = "UMAP_2",
  color_by = "celltype", color_type = "discrete", colors = col_celltype,
  facet_by = "facetvar", facet_nrow = 1, point_size = 0.001, point_alpha = 0.3, raster_dpi = 300,
  shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL,
  label = TRUE, label_df = df.label, label_size = 2.7) +
  theme_expresso(legend_position = "none", facet_label_face = "plain",
                 show_axis = FALSE, grid = "none")

### mature vs immature
dtf[, celltype_matureornot := ifelse(celltype %in% c("AZU1", "LTF"), "immature", "mature")]

df.label2 = as.data.table(dtf)[, lapply(.SD, median), by = celltype_matureornot,
                               .SDcols = c("UMAP_1", "UMAP_2")]
colnames(df.label2)[2:3] = c("x", "y")
# Drawn labels only: celltype_matureornot stays lower case, being what the
# colors vector below is keyed by and what the x nudge tests against.
df.label2$labeltext = cap_first(df.label2$celltype_matureornot)
df.label2[, x := ifelse(celltype_matureornot == "immature", x + 0.5, x)]

dtf$facetvar = "Neutrophil subset"

p.umap.celltype.matureornot = plot_scatter(
  dtf, x = "UMAP_1", y = "UMAP_2", color_by = "celltype_matureornot", color_type = "discrete",
  colors = c("immature" = "#F47B00", "mature" = "grey87"), facet_by = "facetvar", facet_nrow = 1,
  point_size = 0.000001, point_alpha = 0.1, raster_dpi = 400,
  shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL,
  label = TRUE, label_df = df.label2, label_size = 2.7) +
  theme_expresso(legend_position = "none", show_axis = FALSE, grid = "none",
                 show_facet_label = FALSE)


# -- SCAHN UMAPs of genes and signature scores ---------------------------------

mat = readRDS(paste0(dir.results, "annotation/SCAHN.selectedgenes.mat.expr.integrated.rds"))
mat = as.matrix(mat)
mat = mat[, dtf$cell.pid]

dtfgenes = as.data.table(cbind(dtf, t(mat)))

dtf.score = readRDS(paste0(dir.results, "annotation/SCAHN.cellscores.rds"))
dtf.score = dtf.score[cell.pid %in% dtf$cell.pid]

dtfgenes$`granule\ngenes` = dtf.score$score.granules - 1
dtfgenes$immature         = dtf.score$score.immature.neu - 1
dtfgenes$degranulating    = dtf.score$score.degranulating.neu - 1
dtfgenes$antiprotease     = dtf.score$score.antiprotease.neu - 1
dtfgenes$mature           = dtf.score$score.mature.neu - 1
dtfgenes$total            = dtf.score$score.total.neu - 1
dtfgenes$`IFN-related`    = dtf.score$score.ifn - 1

### signature scores, one panel per signature

pscores.neu = list()
for (gene in c("immature", "degranulating", "antiprotease", "mature", "total", "IFN-related")) {
  dtf1 = data.table(UMAP_1 = dtfgenes$UMAP_1, UMAP_2 = dtfgenes$UMAP_2,
                    id = dtfgenes[[gene]], facetvar = cap_first(gene))

  pscores.neu[[gene]] = plot_scatter(
    dtf1, x = "UMAP_1", y = "UMAP_2", color_by = "id",
    facet_by = "facetvar", color_type = "continuous",
    colors = c("grey99", scales::dichromat_pal("DarkRedtoBlue.12")(12)[7:12]),
    quantile_lower = 0.001, quantile_upper = 0.999, shuffle = TRUE, seed = 42,
    point_size = 0.002, point_alpha = 0.2, raster_dpi = 400, legend_ncol = 1, facet_nrow = 1) +
    theme_expresso(legend_position = "none", show_axis = FALSE, grid = "none",
                   facet_label_face = "plain", panel_background = "white",
                   plot_background  = "grey97",
                   plot_margin = ggplot2::margin(0.05, 0.05, 0.05, 0.05, "line"))
}

### individual maturity markers
pgenes = list()
for (gene in c("FCGR3B", "CXCR2", "MME", "CEACAM8", "granule\ngenes")) {
  dtf1 = data.table(UMAP_1 = dtfgenes$UMAP_1, UMAP_2 = dtfgenes$UMAP_2,
                    id = dtfgenes[[gene]], facetvar = cap_first(gene))

  pgenes[[gene]] = plot_scatter(
    dtf1, x = "UMAP_1", y = "UMAP_2", color_by = "id",
    facet_by = "facetvar", color_type = "continuous",
    colors = c("grey99", scales::dichromat_pal("DarkRedtoBlue.12")(12)[7:12]),
    quantile_lower = 0.001, quantile_upper = 0.999, shuffle = TRUE, seed = 42,
    point_size = 0.0001, point_alpha = 0.1, raster_dpi = 400,
    legend_ncol = 1, colorbar_height = 0.9, colorbar_width = 0.05, facet_nrow = 1) +
    theme_expresso(legend_position = c(0.99, 0.48), legend_text_size = 6,
                   facet_label_face = "italic", show_axis = FALSE, grid = "none")
}

### subset-defining genes, small multiples
genes2display = c("STMN1", "AZU1", "MPO", "PRTN3", "CTSG", "ELANE", "DEFA4",
                  "LTF", "MMP8", "CAMP", "ARG1", "MMP9", "IL1R2",
                  "S100A8", "S100A4", "OLR1", "TXNIP", "EGR1", "FOS", "PTGS2",
                  "G0S2", "IFI6", "GBP5", "CD274",
                  "CXCL2", "CXCR4", "VEGFA", "CCL4", "NFKBIA", "IL1B", "IL1RN",
                  "SLPI", "PI3", "HSPA1B", "CD74")

pgenes2 = list()
for (gene in genes2display) {
  dtf1 = data.table(UMAP_1 = dtfgenes$UMAP_1, UMAP_2 = dtfgenes$UMAP_2,
                    id = dtfgenes[[gene]], facetvar = gene)

  pgenes2[[gene]] = plot_scatter(
    dtf1, x = "UMAP_1", y = "UMAP_2", color_by = "id",
    facet_by = "facetvar", color_type = "continuous",
    colors = c("grey99", scales::dichromat_pal("DarkRedtoBlue.12")(12)[7:12]),
    quantile_lower = 0.001, quantile_upper = 0.999, shuffle = TRUE, seed = 42,
    point_size = 0.00001, point_alpha = 0.1, raster_dpi = 400, legend_ncol = 1, facet_nrow = 1) +
    theme_expresso(legend_position = "none", facet_label_face = "italic",
                   show_axis = FALSE, grid = "none") +
    theme(
      plot.background  = element_rect(color = "transparent", fill = "transparent"),
      strip.background = element_rect(color = "transparent", fill = "transparent"),
      plot.margin = margin(0.001, 0.001, 0.001, 0.001, "line"),
      strip.text.x = element_text(size = 8, color = textcolor,
                                  margin = margin(0.02, 0.02, 0.05, 0.02, "line")))
}


# -- Neutrophil proportion in peripheral blood ---------------------------------

# wilhelm2024: not sure if it is fresh or frozen
# deng2021: it says PBMC in the paper, but very high proportion (>60%) of
#           neutrophils was observed in 2 out of 3 samples
dtf.prop.blood = dtf.prop[
  tissue %in% c("peripheral blood") & !is.na(number.allcell)][
  !tissue.detail %in% c("neutrophil", "granulocyte") &
    !pid %in% c("deng2021", "tabulasapiens2022", "han2020", "wilhelm2024")]

dtf.prop.blood[, lapply(.SD, function(x) { median(x) }), by = c("source2"), .SDcols = c("neu.pct")]
dtf.prop.blood[grepl("freshWB", source)][, lapply(.SD, function(x) { median(x) }),
                                         by = c("technique"), .SDcols = c("neu.pct")]

table(dtf.prop.blood$source2)
wilcox.test(dtf.prop.blood[source2 == "fresh\nPBMC"]$neu.pct,
            dtf.prop.blood[source2 == "frozen\nPBMC"]$neu.pct)

summary(dtf.prop.blood[source2 == "fresh\nPBMC"]$neu.pct)
summary(dtf.prop.blood[source2 == "frozen\nPBMC"]$neu.pct)
summary(dtf.prop.blood[source2 == "fresh\nwhole blood"]$neu.pct)

p.neu.pct.source2 = ggplot(dtf.prop.blood, aes(x = source2, y = neu.pct)) +
  ggbeeswarm::geom_quasirandom(size = 0.2, width = 0.2, show.legend = TRUE,
                               alpha = 0.7, varwidth = FALSE, color = "grey42", shape = 16) +
  geom_boxplot(width = 0.35, size = 0.3, outlier.shape = NA,
               alpha = 1, fill = NA, color = "grey17") +
  theme_expresso(grid = "x", axis_text_size = 10) +
  theme(legend.position = "none", plot.title = element_text(hjust = 1.5),
        plot.margin     = unit(c(.1, .1, .2, .1), "line"),
        axis.text.x     = element_text(angle = 0, hjust = 0.5)) +
  labs(x = NULL, y = NULL, title = "") +
  guides(color = guide_legend(ncol = 1, label.hjust = 0,
                              override.aes = list(size = 2, alpha = 1))) +
  scale_x_discrete(limits = rev, position = "bottom", labels = cap_first) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  geom_signif(annotation = c("ns"), y_position = .1, xmin = c(2), xmax = c(3),
              textsize = 3, tip_length = 0.05, size = 0.2) +
  coord_flip()

dtf.prop.blood$technique = as.character(dtf.prop.blood$technique)
dtf.prop.blood[, technique2 := ifelse(technique %in% c("10x Chromium", "BD Rhapsody"),
                                      technique, "other")]

p.neu.pct.technique2 = ggplot(dtf.prop.blood[grepl("freshWB", source)],
                              aes(x = technique2, y = neu.pct)) +
  ggbeeswarm::geom_quasirandom(size = 0.2, width = 0.2, show.legend = TRUE,
                               alpha = 0.7, varwidth = FALSE, shape = 16, color = "grey42") +
  geom_boxplot(width = 0.4, size = 0.3, outlier.shape = NA,
               alpha = 1, fill = NA, color = "grey17") +
  theme_expresso(axis_text_size = 10) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.margin = unit(c(.1, .1, .2, .23), "line")) +
  labs(x = NULL, y = NULL, title = NULL) +
  guides(color = guide_legend(ncol = 1, label.hjust = 0,
                              override.aes = list(size = 2, alpha = 1))) +
  scale_x_discrete(limits = rev, position = "bottom", labels = cap_first) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_flip()


# -- Neutrophil subset distribution, bar and boxplot ---------------------------

dtf.prop[, technique2 := ifelse(technique %in% c("10x Chromium", "BD Rhapsody"),
                                technique, "other")]
dtf.prop$immature = dtf.prop$AZU1 + dtf.prop$LTF
dtf.prop$`cytokine/inflammation subsets` =
  dtf.prop$CXCL + dtf.prop$`NF-κB` + dtf.prop$IL1B + dtf.prop$IL1RN +
  dtf.prop$VEGFA + dtf.prop$`CCL3/4`

### all tissues
dtf.prop.m = melt(
  dtf.prop, id.vars = c("group", "pid", "sampleid.pid", "tissue", "tissue.detail",
              "source", "technique2"), measure.vars = c(names(col_celltype), "immature",
                   "cytokine/inflammation subsets"))
dtf.prop.m = as.data.table(dtf.prop.m)
length(table(dtf.prop.m$sampleid.pid)) # 785

dtf.prop.m.aggmean = dtf.prop.m[, lapply(.SD, function(x) { mean(x) }),
                                by = c("tissue", "variable"), .SDcols = c("value")]
dtf.prop.m.aggmedian = dtf.prop.m[, lapply(.SD, function(x) { median(x) }),
                                  by = c("tissue", "variable"), .SDcols = c("value")]
dtf.prop.m.aggsd = dtf.prop.m[, lapply(.SD, function(x) { sd(x) }),
                              by = c("tissue", "variable"), .SDcols = c("value")]

### peripheral blood
dtf.prop.m.blood = dtf.prop.m[
  source %in% c("freshWB_RBClysis", "frozenPBMC", "freshWB_RBCdepletion",
                "freshPBMC", "neutrophil") & tissue %in% c("peripheral blood")]
dim(dtf.prop.m.blood)
dtf.prop.m.blood = as.data.table(dtf.prop.m.blood)
dtf.prop.m.blood[, tissue := ifelse(grepl("WB", source), "whole blood", source)]
dtf.prop.m.blood[, tissue := ifelse(grepl("neutrophil", source), "purified neutrophil", tissue)]
dtf.prop.m.blood[, tissue := ifelse(grepl("PBMC", source), "PBMC", tissue)]
table(dtf.prop.m.blood$tissue)

dtf.prop.m.blood.aggmean = dtf.prop.m.blood[, lapply(.SD, function(x) { mean(x) }),
                                            by = c("tissue", "variable"), .SDcols = c("value")]
dtf.prop.m.blood.aggmedian = dtf.prop.m.blood[, lapply(.SD, function(x) { median(x) }),
                                              by = c("tissue", "variable"), .SDcols = c("value")]
dtf.prop.m.blood.aggsd = dtf.prop.m.blood[, lapply(.SD, function(x) { sd(x) }),
                                          by = c("tissue", "variable"), .SDcols = c("value")]

summary(dtf.prop.m.blood[tissue == "whole blood"][variable == "AZU1"]$value +
          dtf.prop.m.blood[tissue == "whole blood"][variable == "LTF"]$value)
summary(dtf.prop.m.blood[tissue == "PBMC"][variable == "AZU1"]$value +
          dtf.prop.m.blood[tissue == "PBMC"][variable == "LTF"]$value)

### whole blood, by platform
dtf.prop.WB = dtf.prop[source %in% c("freshWB_RBClysis", "freshWB_RBCdepletion")]
table(dtf.prop.WB$group, dtf.prop.WB$technique2)

dtf.prop.BM = dtf.prop[tissue %in% c("bone marrow")]
table(dtf.prop.BM$group, dtf.prop.BM$technique2)

dtf.prop.m.blood_platform.aggmean = dtf.prop.m.blood[tissue %in% c("whole blood")][
  , lapply(.SD, function(x) { mean(x) }),
  by = c("pid", "technique2", "variable"), .SDcols = c("value")]
dtf.prop.m.blood_platform.aggmedian = dtf.prop.m.blood[tissue %in% c("whole blood")][
  , lapply(.SD, function(x) { median(x) }),
  by = c("pid", "technique2", "variable"), .SDcols = c("value")]
dtf.prop.m.blood_platform.aggsd = dtf.prop.m.blood[tissue %in% c("whole blood")][
  , lapply(.SD, function(x) { sd(x) }),
  by = c("pid", "technique2", "variable"), .SDcols = c("value")]

celltypeorder = c(rev(names(col_celltype)), "immature", "cytokine/inflammation subsets")
dtf.prop.m.blood_platform.aggmean$variable =
  factor(dtf.prop.m.blood_platform.aggmean$variable, levels = celltypeorder)
dtf.prop.m.blood_platform.aggmedian$variable =
  factor(dtf.prop.m.blood_platform.aggmedian$variable, levels = celltypeorder)
dtf.prop.m.blood_platform.aggsd$variable =
  factor(dtf.prop.m.blood_platform.aggsd$variable, levels = celltypeorder)

### all tissues and blood combined
dtf.prop.m.aggmean.combined =
  rbind(dtf.prop.m.aggmean, dtf.prop.m.blood.aggmean[, colnames(dtf.prop.m.aggmean), with = FALSE])
dtf.prop.m.aggmedian.combined =
  rbind(dtf.prop.m.aggmedian,
        dtf.prop.m.blood.aggmedian[, colnames(dtf.prop.m.aggmedian), with = FALSE])
dtf.prop.m.aggsd.combined =
  rbind(dtf.prop.m.aggsd, dtf.prop.m.blood.aggsd[, colnames(dtf.prop.m.aggsd), with = FALSE])

celltypeorder = c(rev(names(col_celltype)), "immature", "cytokine/inflammation subsets")
dtf.prop.m.aggmean.combined$variable =
  factor(dtf.prop.m.aggmean.combined$variable, levels = celltypeorder)
dtf.prop.m.aggmedian.combined$variable =
  factor(dtf.prop.m.aggmedian.combined$variable, levels = celltypeorder)
dtf.prop.m.aggsd.combined$variable =
  factor(dtf.prop.m.aggsd.combined$variable, levels = celltypeorder)

tissueorder = c("other", "BALF", "lung", "sputum", "other GI organs", "stomach", "pancreas",
                "brain", "colonrectum", "kidney", "liver", "spleen", "bone marrow", "cord blood",
                "PBMC", "purified neutrophil", "whole blood", "peripheral blood")

dtf.prop.m.aggmean.combined$tissue =
  factor(dtf.prop.m.aggmean.combined$tissue, levels = tissueorder)
dtf.prop.m.aggmedian.combined$tissue =
  factor(dtf.prop.m.aggmedian.combined$tissue, levels = tissueorder)
dtf.prop.m.aggsd.combined$tissue =
  factor(dtf.prop.m.aggsd.combined$tissue, levels = tissueorder)

### solid tissue vs blood
pct.tissue = dtf.prop.m[
  tissue %in% c("lung", "colonrectum", "other GI organs", "liver", "pancreas",
                "brain", "kidney", "stomach") |
    tissue.detail %in% c("bladder", "breast", "cervix", "heart", "ovarian", "ovary", "uterus")
][, c("tissue", "variable", "value")]
pct.tissue$tissue = "solid tissue"

pct.blood = dtf.prop.m.blood[grepl("WB|neutrophil", source)][
  , c("tissue", "variable", "value")]
pct.blood$tissue = "whole blood / purified neutrophil"

dtf.tmp = rbind(
  dtf.prop.m.blood[tissue == "PBMC"][, c("tissue", "variable", "value")], pct.blood,
  dtf.prop.m[tissue %in% c("cord blood", "bone marrow", "spleen")][
    , c("tissue", "variable", "value")], pct.tissue)

dtf.tmp$tissue = factor(dtf.tmp$tissue, levels = c("whole blood / purified neutrophil", "PBMC",
                                   "cord blood", "bone marrow", "spleen", "solid tissue"))

### bar, main figure
p.tissuedistribution = ggplot(
  dtf.prop.m.aggmean.combined[
    !tissue %in% c("other") & !variable %in% c("immature", "cytokine/inflammation subsets")
  ], aes(x = tissue, y = value, fill = variable)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  theme_expresso(axis_text_size = textsize, grid = "none",
                 legend_position = "none", legend_key_width = 0.2,
                 legend_key_height = 0.2, plot_margin = margin(.2, .2, .2, .2, "line")) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = seq(0, 1, 0.2), position = "left") +
  scale_x_discrete(position = "top", labels = cap_tissue) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(
    legend.text      = element_text(color = textcolor, face = "italic"),
    legend.spacing.y = unit(0.01, "line"), legend.spacing.x = unit(0.04, "line"),
    plot.title       = element_text(size = textsize, face = "plain", hjust = 0.5,
                                    margin = margin(t = 1, r = 0, b = 1, l = 0, unit = "pt"))) +
  scale_fill_manual(values = col_celltype, name = "") +
  coord_flip() +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE, size = 2.3))

### bar, by platform
dtf.prop.m.blood_platform.aggmean$technique_pid =
  paste0(dtf.prop.m.blood_platform.aggmean$technique2, ",", dtf.prop.m.blood_platform.aggmean$pid)

p.tissuedistribution.platform = ggplot(
  dtf.prop.m.blood_platform.aggmean[
    !variable %in% c("immature", "cytokine/inflammation subsets")],
  aes(x = technique_pid, y = value, fill = variable)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  theme_expresso(axis_text_size = textsize, grid = "none",
                 legend_position = "none", legend_key_width = 0.2,
                 legend_key_height = 0.2, plot_margin = margin(.2, .2, .2, .2, "line")) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = seq(0, 1, 0.2), position = "left") +
  scale_x_discrete(position = "top", limits = rev, labels = cap_first) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(
    legend.text      = element_text(color = textcolor, face = "italic"),
    legend.spacing.y = unit(0.01, "line"), legend.spacing.x = unit(0.04, "line"),
    plot.title       = element_text(size = textsize, face = "plain", hjust = 0.5,
                                    margin = margin(t = 1, r = 0, b = 1, l = 0, unit = "pt"))) +
  scale_fill_manual(values = col_celltype, name = "") +
  coord_flip() +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE, size = 2.3))

### heatmap, main figure
dtf.prop.m.agg.combined.mat = dcast(
  dtf.prop.m.aggmedian.combined[
    !tissue %in% c("other") & !variable %in% c("immature", "cytokine/inflammation subsets")
  ], tissue ~ variable, value.var = "value")
rowname = dtf.prop.m.agg.combined.mat$tissue
dtf.prop.m.agg.combined.mat$tissue = NULL
dtf.prop.m.agg.combined.mat = as.matrix(dtf.prop.m.agg.combined.mat)
rownames(dtf.prop.m.agg.combined.mat) = rowname
dtf.prop.m.agg.combined.mat = dtf.prop.m.agg.combined.mat[
  , intersect(names(col_celltype), colnames(dtf.prop.m.agg.combined.mat))]

dtf.prop.m.aggsd.combined.mat = dcast(
  dtf.prop.m.aggsd.combined[
    !tissue %in% c("other") & !variable %in% c("immature", "tissuesubset")],
  tissue ~ variable, value.var = "value")
rowname = dtf.prop.m.aggsd.combined.mat$tissue
dtf.prop.m.aggsd.combined.mat$tissue = NULL
dtf.prop.m.aggsd.combined.mat = as.matrix(dtf.prop.m.aggsd.combined.mat)
rownames(dtf.prop.m.aggsd.combined.mat) = rowname
dtf.prop.m.aggsd.combined.mat = dtf.prop.m.aggsd.combined.mat[
  , intersect(names(col_celltype), colnames(dtf.prop.m.aggsd.combined.mat))]

col_fun.tissueht = circlize::colorRamp2(
  seq(0, max(dtf.prop.m.agg.combined.mat), length = 8), c("grey99", pal_material("blue")(10)[4:10]))
dtf.prop.m.agg.combined.mat =
  dtf.prop.m.agg.combined.mat[rev(rownames(dtf.prop.m.agg.combined.mat)), ]
dtf.prop.m.aggsd.combined.mat =
  dtf.prop.m.aggsd.combined.mat[rev(rownames(dtf.prop.m.aggsd.combined.mat)), ]
mat.size = 1 - sqrt(dtf.prop.m.aggsd.combined.mat)

p.tissue.ht = Heatmap(
  dtf.prop.m.agg.combined.mat,
  rect_gp = gpar(type = "none"), col = col_fun.tissueht,
  na_col = "grey93", cluster_columns = FALSE, show_row_names = TRUE,
  cluster_rows    = FALSE, show_row_dend = TRUE, show_column_dend = TRUE,
  row_dend_side   = "left", column_dend_side = "top",
  clustering_method_rows      = "ward.D2", clustering_method_columns   = "ward.D2",
  clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
  row_names_side = "right", column_names_side = "top", column_names_rot = 90,
  row_labels = cap_tissue(rownames(dtf.prop.m.agg.combined.mat)),
  column_names_gp = gpar(fontsize = 10, fontface = "plain", col = col_celltype),
  row_names_gp = gpar(fontsize = 10, fontface = "plain"), column_title = NULL,
  column_title_side = "top", column_title_gp = gpar(fontsize = 8, fontface = "plain"),
  show_heatmap_legend = TRUE, cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(lwd = 0.2, col = "grey90", fill = NA, alpha = 0.2))
    grid.circle(x = x, y = y, r = 0.6 * mat.size[i, j] * min(unit.c(width, height)),
                gp = gpar(fill = col_fun.tissueht(dtf.prop.m.agg.combined.mat[i, j]),
                          col = NA, alpha = 0.9))
  }, heatmap_legend_param = list(
    title = "Median proportion", title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 9),
    legend_width = unit(4, "line"), legend_height = unit(.2, "line"),
    grid_width = unit(5, "line"), grid_height = unit(.1, "line"),
    by_row = FALSE, direction = "horizontal", title_position = "lefttop", at = c(0, 0.25, 0.5)))

### boxplot, solid tissue vs whole blood / purified neutrophil
p.tissuesubsets.pct = ggplot(
  dtf.tmp[variable %in% c("NF-κB", "IL1B", "IL1RN", "CXCL", "VEGFA", "CCL3/4")][
    tissue %in% c("solid tissue", "whole blood / purified neutrophil")],
  aes(x = tissue, y = value + 0.01)) +
  ggbeeswarm::geom_quasirandom(aes(color = tissue), size = .1, width = 0.2, show.legend = TRUE,
                               alpha = 0.2, varwidth = FALSE, shape = 16) +
  geom_boxplot(aes(color = tissue), width = 0.2, size = 0.2,
               outlier.shape = NA, alpha = 1, fill = NA) +
  labs(y = "Proportion in total neutrophils", x = NULL, title = NULL) +
  theme_expresso(axis_text_size = textsize, grid = "y",
                 legend_text_size = 8, legend_key_width = 0.3,
                 legend_key_height = 0.7, legend_justification = c(0.5, 0.5),
                 legend_position      = c(0.5, -0.08)) +
  scale_y_log10() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_text(hjust = 0.75), legend.key.size = unit(1, "lines"),
        plot.title      = element_text(size = textsize, face = "plain", hjust = 0,
                                       margin = margin(t = .1, r = 0, b = .1, l = 0,
                                                       unit = "pt"))) +
  scale_color_manual(values = c("solid tissue" = "grey17",
                                "whole blood / purified neutrophil" = "#EE6C00"),
                     labels = cap_first) +
  facet_wrap(. ~ variable, nrow = 3, scales = "free") +
  stat_compare_means(comparisons = list(c("whole blood / purified neutrophil", "solid tissue")),
                     method = "wilcox.test", vjust = 1.5, size = 2,
                     label = "p.format", color = "grey17", bracket.size = 0.1, tip.length = 0.01) +
  guides(color = guide_legend(nrow = 2, label.hjust = 0,
                              override.aes = list(size = 0.5, alpha = 1)))

p.tissuesubset.pct = ggplot(
  dtf.tmp[variable %in% c("cytokine/inflammation subsets")][
    tissue %in% c("solid tissue", "whole blood / purified neutrophil")],
  aes(x = tissue, y = value + 0.01)) +
  ggbeeswarm::geom_quasirandom(aes(color = tissue), size = .1, width = 0.2, show.legend = TRUE,
                               alpha = 0.2, varwidth = FALSE, shape = 16) +
  geom_boxplot(aes(color = tissue), width = 0.2, size = 0.2,
               outlier.shape = NA, alpha = 1, fill = NA) +
  labs(y = NULL, x = NULL, title = "Cytokine/inflammation\nsubsets") +
  theme_expresso(axis_text_size = textsize, grid = "y", legend_position = "none") +
  scale_y_log10() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.title   = element_text(size = textsize, face = "plain", hjust = 0,
                                    margin = margin(t = .1, r = 0, b = .1, l = 0, unit = "pt"))) +
  scale_color_manual(values = c("solid tissue" = "grey17",
                                "whole blood / purified neutrophil" = "#EE6C00")) +
  stat_compare_means(comparisons = list(c("whole blood / purified neutrophil", "solid tissue")),
                     method = "wilcox.test", vjust = 1.5, size = 2,
                     label = "p.format", color = "grey17", bracket.size = 0.1, tip.length = 0.01) +
  guides(color = guide_legend(nrow = 2, label.hjust = 0,
                              override.aes = list(size = 0.5, alpha = 1)))


# -- Neutrophil subset distribution, effect-size heatmap -----------------------

dtf.prop[, tissue2 := ifelse(grepl("WB|neutrophil", source) & tissue %in% c("peripheral blood"),
                             "whole blood / purified neutrophil", tissue)]
dtf.prop[, tissue2 := ifelse(source %in% c("freshPBMC", "frozenPBMC"), "PBMC", tissue2)]
dtf.prop[, tissue2 := ifelse(
  tissue %in% c("lung", "colonrectum", "other GI organs", "liver", "pancreas",
                "brain", "kidney", "stomach") |
    tissue.detail %in% c("bladder", "breast", "cervix", "heart", "ovarian", "ovary", "uterus"),
  "solid tissue", tissue2)]
dtf.prop = dtf.prop[tissue2 %in% c("bone marrow", "cord blood", "PBMC", "solid tissue",
                                   "whole blood / purified neutrophil", "spleen")]
dim(dtf.prop)
table(dtf.prop$tissue2)
colnames(dtf.prop)

subsetnames = names(col_celltype)

base_tissue   = "whole blood / purified neutrophil"
other_tissues = setdiff(unique(dtf.prop$tissue2), base_tissue)

prop_stats = rbindlist(lapply(subsetnames, function(ss) {
  rbindlist(lapply(other_tissues, function(tiss) {
    x = dtf.prop[tissue2 == tiss,        get(ss)]
    y = dtf.prop[tissue2 == base_tissue, get(ss)]
    x = x[is.finite(x)]
    y = y[is.finite(y)]
    if (length(x) < 2 || length(y) < 2) return(NULL)

    vals   = c(x, y)
    labels = c(rep(1L, length(x)), rep(0L, length(y)))
    mat    = matrix(vals, nrow = 1, dimnames = list(ss, NULL))
    hg     = expresso::gene_hedges_g(mat, labels)

    wt = wilcox.test(x, y, exact = FALSE)
    wt = t.test(x, y)

    data.table(
      subset = ss, tissue2 = tiss, n = length(x), n_base = length(y),
      median = median(x), median_base = median(y), hedges_g = hg$es[[ss]], p_value = wt$p.value)
  }))
}))

prop_stats[, FDR := p.adjust(p_value, method = "BH")]

### cell_fun: dot area scales with effect size, outline marks FDR <= cutoff
make_cellfun = function(val_mat, fdr_mat, cf, dot_scale = 0.35, fdr_cutoff = 0.1) {
  sz_mat = 1.2 - fdr_mat
  fd2 = fdr_mat
  fd2[fd2 > fdr_cutoff] = 1.2
  sz2_mat = 1.2 - fd2          # 0 when FDR > fdr_cutoff -> outline invisible

  function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height,
              gp = gpar(lwd = 0.5, col = "grey53", fill = "white", alpha = 0.2))
    v = val_mat[i, j]
    f = fdr_mat[i, j]
    if (is.na(v) || is.na(f)) return(invisible(NULL))
    r_base = dot_scale * min(unit.c(width, height))
    grid.circle(x, y, r = sz_mat[i, j] * r_base, gp = gpar(fill = cf(v), col = NA, alpha = 0.9))
    grid.circle(x, y, r = sz2_mat[i, j] * r_base, gp = gpar(fill = NA, col = "grey17", alpha = 1))
  }
}

hg_mat       = dcast(prop_stats, tissue2 ~ subset, value.var = "hedges_g")
fdr_mat_prop = dcast(prop_stats, tissue2 ~ subset, value.var = "FDR")

rn = hg_mat$tissue2
hg_mat = as.matrix(hg_mat[, -1])
rownames(hg_mat) = rn
fdr_mat_prop = as.matrix(fdr_mat_prop[, -1])
rownames(fdr_mat_prop) = rn

# keep column order consistent with subsetnames
cs = intersect(subsetnames, colnames(hg_mat))
hg_mat       = hg_mat[, cs, drop = FALSE]
fdr_mat_prop = fdr_mat_prop[, cs, drop = FALSE]

# fixed row order
row_order = c("bone marrow", "cord blood", "spleen", "PBMC", "solid tissue")
row_order = intersect(row_order, rownames(hg_mat))

col_order = c("AZU1", "LTF", "HSP", "CD74", "CXCL", "CCL3/4", "IL1B", "IL1RN",
              "MMP9", "NF-κB", "VEGFA", "SLPI", "G0S2", "IL1R2", "IFN3", "IFN1",
              "EGR1", "AP-1", "S100A4", "PTGS2", "IFN2", "Low depth", "TXNIP", "MME")

hg_mat       = hg_mat[row_order, col_order, drop = FALSE]
fdr_mat_prop = fdr_mat_prop[row_order, col_order, drop = FALSE]

hg_mat       = t(hg_mat)
fdr_mat_prop = t(fdr_mat_prop)

# color scale
clamp_hg = max(abs(hg_mat[is.finite(hg_mat)]), na.rm = TRUE)
clamp_hg = min(clamp_hg, 3)
cf_hg    = colorRamp2(seq(-clamp_hg, clamp_hg, length.out = 12),
                      dichromat::colorschemes$DarkRedtoBlue.12)

ht_prop = Heatmap(
  hg_mat, rect_gp = gpar(type = "none", fill = "white"),
  col = cf_hg, na_col = "grey", cluster_rows = FALSE, cluster_columns = FALSE,
  clustering_method_rows = "ward.D", clustering_method_columns = "ward.D",
  clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
  show_row_dend = FALSE, show_column_dend = TRUE, column_dend_side = "bottom",
  show_column_names = TRUE, column_names_side = "bottom", row_names_side = "left",
  column_names_rot = 30, column_names_gp = gpar(fontsize = textsize),
  column_labels = cap_first(colnames(hg_mat)),
  row_names_gp = gpar(fontsize = textsize), row_names_max_width = unit(40, "line"),
 column_title = expression("Proportion" ~ italic("vs.") ~ "WB"),
  column_title_side = "top", column_title_gp = gpar(fontsize = textsize, fontface = "plain"),
  show_heatmap_legend = TRUE, cell_fun = make_cellfun(hg_mat, fdr_mat_prop, cf_hg,
                          dot_scale = 0.8, fdr_cutoff = 0.05), heatmap_legend_param = list(
    title = "Effect size", title_gp = gpar(fontsize = textsize),
    labels_gp = gpar(fontsize = textsize),
    legend_widht = unit(0.8, "line"), legend_height = unit(0.2, "line"),
    grid_width = unit(0.8, "line"), grid_height = unit(.2, "line"),
    direction = "horizontal", title_position = "lefttop",
    at = round(seq(-clamp_hg, clamp_hg, length.out = 3), 1)))


# -- SCAHN signature genes dotplot ---------------------------------------------

mat = readRDS(paste0(dir.results, "annotation/SCAHN.selectedgenes.mat.expr.integrated.rds"))
dtfexpr = cbind(dtf, t(as.matrix(mat[signaturegenes, dtf$cell.pid])))

dtfexpr$celltype = factor(dtfexpr$celltype, levels = rev(names(col_celltype)))

# keyed by the color each panel is labelled in
gene_sets = list(
  "#a93434" = genes.immature.neu, "#F97D1C" = genes.degranulating.neu,
  "#7CABB1" = genes.antiprotease.neu, "#1A5E1F" = genes.mature.neu)

gene_color_map = unlist(lapply(
  names(gene_sets), function(col) setNames(rep(col, length(gene_sets[[col]])), gene_sets[[col]])))

label_map = setNames(signaturegenes, signaturegenes)
colored = intersect(names(gene_color_map), signaturegenes)
label_map[colored] = paste0('<span style="color:', gene_color_map[colored], '">',
                            colored, '</span>')

p.markerdot = plot_dotplot(
  dtfexpr, features = signaturegenes, group_by = "celltype",
  max_scale = 2, dot_scale = 3.2, dot_stroke = 0.05, title = NULL,
  col_fontsize = 8, row_fontsize = 8,
  cluster_features = FALSE, cluster_groups = FALSE, show_axis = FALSE,
  feature_side = "top", colorbar_title = "Scaled expression") +
  scale_x_discrete(labels = label_map, position = "top") +
  theme(
    legend.position = c(1, 0.4), axis.text.x.top = element_markdown(
      angle = 90, hjust = 0, vjust = 0.5, size = 8, margin = margin(b = 0.1, unit = "line")))

df.gene_sets = data.table(
  subset = c("Neutrophil-\nspecific\ngene signatures", "- Immature",
             "- Degranulating", "- Antiprotease", "- Mature"),
  color = c("#252525", names(gene_sets)), x = 1, y = (length(gene_sets) + 1):1)

p.gene_sets = ggplot(df.gene_sets, aes(x = x, y = y, label = subset, color = color)) +
  geom_text(angle = 0, size = 3, hjust = 0, lineheight = 0.8, vjust = 0) +
  scale_color_identity() +
  theme_void() +
  coord_cartesian(clip = "off")


# -- QC: nFeature, nCount, ribosome, MT ----------------------------------------

dtf[, technique2 := ifelse(grepl("10x Chromium", technique), "10x Chromium",
                           ifelse(grepl("BD Rhapsody", technique), "BD Rhapsody", "other"))]
table(dtf$technique2)

dtf.score = readRDS(paste0(dir.results, "annotation/SCAHN.cellscores.rds"))
dtf$`RPL/RPS expression` = dtf.score$score.ribosome
dtf$`Mitochondrial gene expression` = dtf$percent.mt

dtf.nFeaturenCount = dtf[
  source %in% c("freshPBMC", "frozenPBMC", "freshWB_RBClysis", "freshWB_RBCdepletion")][
  , lapply(.SD, function(x) { geom_mean(as.numeric(x)) }),
  by = c("source", "technique2", "sampleid.pid", "pid"),
  .SDcols = c("nCount_RNA", "nFeature_RNA", "RPL/RPS expression")]
dtf.nFeaturenCount$source = gsub("freshWB_RBClysis|freshWB_RBCdepletion",
                                 "freshWB", dtf.nFeaturenCount$source)
dtf.nFeaturenCount$source.tech = paste0(dtf.nFeaturenCount$source, ", ",
                                        dtf.nFeaturenCount$technique2)

dtf.nFeaturenCount$source.tech = factor(
  dtf.nFeaturenCount$source.tech, levels = c("freshPBMC, 10x Chromium", "freshPBMC, BD Rhapsody",
             "freshPBMC, other", "frozenPBMC, 10x Chromium", "frozenPBMC, BD Rhapsody",
             "freshWB, 10x Chromium", "freshWB, BD Rhapsody", "freshWB, other"))

# Drawn form of the source levels: "freshPBMC" is what the rows are subset,
# factored and grouped on above, "Fresh PBMC" is only what a label says. Splits
# the run-together pair and hands the rest to cap_first(), so "PBMC" and "WB"
# keep their case. Used at the two Fig.s1 sites that draw these strings: the
# p.nCount axis and the row names of the BD-vs-10x forest.
label_source = function(x) cap_first(sub("^(fresh|frozen)(PBMC|WB)", "\\1 \\2", x))

### meta-analysis: WB vs PBMC
table(dtf.nFeaturenCount$technique2)
dtf.nFeaturenCount[, class := ifelse(grepl("freshPBMC|frozenPBMC", source), 0,
                                     ifelse(grepl("freshWB", source), 1, NA))]
table(dtf.nFeaturenCount$class, dtf.nFeaturenCount$source)

features = c("nCount_RNA", "nFeature_RNA")
pids     = c("10x Chromium", "BD Rhapsody")

metaobj = lapply(pids, function(i) {
  sub = dtf.nFeaturenCount[technique2 == i]
  meta_dataset(
    expr = t(as.matrix(sub[, features, with = FALSE])), class = sub$class, label = i)
})

meta.nFeaturenCount.WBvsPBMC = meta_analysis(metaobj, outcome_type = "binary")

### meta-analysis: fresh vs frozen PBMC
table(dtf.nFeaturenCount$technique2)
dtf.nFeaturenCount[, class := ifelse(grepl("frozenPBMC", source), 0,
                                     ifelse(grepl("freshPBMC", source), 1, NA))]
table(dtf.nFeaturenCount$class, dtf.nFeaturenCount$source)

features = c("nCount_RNA", "nFeature_RNA")
pids     = c("10x Chromium", "BD Rhapsody")

metaobj = lapply(pids, function(i) {
  sub = dtf.nFeaturenCount[technique2 == i]
  meta_dataset(
    expr = t(as.matrix(sub[, features, with = FALSE])), class = sub$class, label = i)
})

meta.nFeaturenCount.freshvsfrozenPBMC = meta_analysis(metaobj, outcome_type = "binary")

### meta-analysis: BD Rhapsody vs 10x Chromium
table(dtf.nFeaturenCount$technique2)
dtf.nFeaturenCount[, class := ifelse(grepl("10x Chromium", technique2), 0,
                                     ifelse(grepl("BD Rhapsody", technique2), 1, NA))]
table(dtf.nFeaturenCount$class, dtf.nFeaturenCount$technique2)

features = c("nCount_RNA", "nFeature_RNA")
pids     = c("freshWB", "freshPBMC", "frozenPBMC")

metaobj = lapply(pids, function(i) {
  sub = dtf.nFeaturenCount[source == i]
  meta_dataset(
    expr = t(as.matrix(sub[, features, with = FALSE])), class = sub$class, label = i)
})

meta.nFeaturenCount.BDvs10x = meta_analysis(metaobj, outcome_type = "binary")

### per source and platform
p.nFeature = ggplot(dtf.nFeaturenCount, aes(x = source.tech, y = nFeature_RNA)) +
  ggbeeswarm::geom_quasirandom(aes(color = pid), size = 0.6, width = 0.2,
                               show.legend = TRUE, alpha = 0.7,
                               varwidth = FALSE, shape = 16, stroke = 0.1) +
  geom_boxplot(width = 0.2, size = 0.2, outlier.shape = NA,
               alpha = 1, fill = NA, color = "grey17") +
  labs(y = "nFeature_RNA", x = NULL, title = "") +
  coord_cartesian(clip = "off") +
  scale_y_log10() +
  scale_color_tableau("Tableau 20", labels = cap_first) +
  theme_expresso(legend_position = "none") +
  theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank(),
        plot.title = element_text(size = textsize, face = "plain", hjust = 0)) +
  guides(colour = guide_legend(override.aes = list(size = 2)))

p.nCount = ggplot(dtf.nFeaturenCount, aes(x = source.tech, y = nCount_RNA)) +
  ggbeeswarm::geom_quasirandom(aes(color = pid), size = 0.6, width = 0.2,
                               show.legend = TRUE, alpha = 0.7,
                               varwidth = FALSE, shape = 16, stroke = 0.1) +
  geom_boxplot(width = 0.2, size = 0.2, outlier.shape = NA,
               alpha = 1, fill = NA, color = "grey17") +
  labs(y = "nCount_RNA", x = NULL, title = "") +
  coord_cartesian(clip = "off") +
  scale_y_log10() +
  scale_x_discrete(labels = label_source) +
  scale_color_tableau("Tableau 20", labels = cap_first) +
  theme_expresso(legend_position = c(1, 1), legend_text_size = 8) +
  theme(axis.text.x = element_text(lineheight = 0.7, angle = -30, hjust = 0),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = textsize, face = "plain", hjust = 0)) +
  guides(colour = guide_legend(override.aes = list(size = 1.5)))

### per subset
dtf.nFeaturenCount2 = dtf[
  , lapply(.SD, function(x) geom_mean(as.numeric(x))), by = .(celltype, sampleid.pid, pid),
  .SDcols = c("nCount_RNA", "nFeature_RNA", "RPL/RPS expression", "Mitochondrial gene expression")]

medians = dtf.nFeaturenCount2[
  , lapply(.SD, median), by = "celltype",
  .SDcols = c("nCount_RNA", "nFeature_RNA", "RPL/RPS expression", "Mitochondrial gene expression")]

make_qr_plot = function(y_var) {
  ggplot(dtf.nFeaturenCount2,
         aes(x = factor(celltype, levels = medians[order(-get(y_var))]$celltype),
             y = .data[[y_var]])) +
    ggbeeswarm::geom_quasirandom(color = "grey77", size = 0.001, width = 0.2,
                                 alpha = 0.7, shape = 16) +
    geom_boxplot(width = 0.2, linewidth = 0.2, outlier.shape = NA, fill = NA, color = "grey17") +
    labs(y = y_var, x = NULL, title = NULL) +
    coord_cartesian(clip = "off") +
    theme_expresso(legend_position = "none") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.x = element_blank())
}

p.nFeature.subset = make_qr_plot("nFeature_RNA") + scale_y_log10()
p.nCount.subset   = make_qr_plot("nCount_RNA") + scale_y_log10()
p.ribosome.subset = make_qr_plot("RPL/RPS expression")
p.mt.subset       = make_qr_plot("Mitochondrial gene expression")


# -- Trajectory ----------------------------------------------------------------

sling.res = readRDS(paste0(dir.results, "annotation/SCAHN.slingshot.res.rds"))
curves_emb = sling.res$curves_emb
pseudotime_matrix = sling.res$pseudotime_matrix
dtf.s = dtf[cell.pid %in% rownames(pseudotime_matrix)]

lineages_df = lapply(names(curves_emb), function(lineage_name) {
  crv = curves_emb[[lineage_name]]
  data.frame(
    UMAP_1 = crv$s[, 1], UMAP_2 = crv$s[, 2], Lineage = lineage_name)
})

lineages_df = do.call(rbind, lineages_df)
lineages_df = as.data.table(lineages_df)

df.slingshot = data.table(
  UMAP_1 = dtf.s$UMAP_1, UMAP_2 = dtf.s$UMAP_2,
  celltype = dtf.s$celltype, sampleid.pid = dtf.s$sampleid.pid,
  facetvar = "Pseudotime", pseudotime = rowMeans2(pseudotime_matrix, na.rm = TRUE))

dtf.pseudotime = df.slingshot[, lapply(.SD, function(x) geom_mean(as.numeric(x))),
                              by = .(celltype, sampleid.pid), .SDcols = c("pseudotime")]
medians.pseudotime = dtf.pseudotime[, lapply(.SD, median), by = "celltype",
                                    .SDcols = c("pseudotime")]

p.pseudotime.subset = ggplot(
  dtf.pseudotime[pseudotime >= 30 & pseudotime <= 90], aes(x = factor(celltype,
                 levels = medians.pseudotime[order(pseudotime)]$celltype),
      y = .data[["pseudotime"]])) +
  ggbeeswarm::geom_quasirandom(color = "grey77", size = 0.01, width = 0.2,
                               alpha = 0.7, shape = 16, stroke = 0.1) +
  geom_boxplot(width = 0.2, linewidth = 0.2, outlier.shape = NA, fill = NA, color = "grey17") +
  labs(y = "Pseudotime", x = NULL, title = "Pseudotime in neutrophil subsets") +
  coord_cartesian(clip = "off") +
  theme_expresso(legend_position = "none") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.x = element_blank())

df.slingshot$celltype = factor(df.slingshot$celltype, levels = names(col_celltype))

df.label = as.data.table(df.slingshot)[, lapply(.SD, median), by = celltype,
                                       .SDcols = c("UMAP_1", "UMAP_2")]
colnames(df.label)[2:3] = c("x", "y")
df.label$labeltext = df.label[, "celltype", with = FALSE]
df.label[, x := ifelse(celltype == "IL1B",  x + 1,   x)]
df.label[, x := ifelse(celltype == "IL1RN", x - 0.5, x)]
df.label[, y := ifelse(celltype == "IL1RN", y + 0.5, y)]
df.label[, y := ifelse(celltype == "IL1R2", y + 1,   y)]
df.label[, y := ifelse(celltype == "PTGS2", y - 1,   y)]

p.trajectory = plot_scatter(
  df.slingshot, x = "UMAP_1", y = "UMAP_2", color_by = "pseudotime", facet_by = "facetvar",
  color_type = "continuous", colors = viridis_pal(option = "D")(10),
  quantile_lower = 0.001, quantile_upper = 0.999, shuffle = TRUE, seed = 42,
  point_size = 0.0001, point_alpha = 0.3, raster_dpi = 400, legend_ncol = 1, facet_nrow = 1,
  label = TRUE, label_df = df.label, label_size = 2.5) +
  geom_path(data = lineages_df, aes(group = Lineage), linewidth = 0.3) +
  theme_expresso(legend_position = c(1, 0.5), legend_key_width = 0.1,
                 legend_key_height = 0.3, legend_key_spacing_x = 0.01,
                 legend_key_spacing_y = 0.05, show_axis = FALSE, grid = "none") +
  scale_color_gradientn(colours = viridis_pal(option = "D")(10),
                        breaks = c(55, 100), labels = c("low", "high"))


# -- Immune cells in peripheral blood, UMAPs -----------------------------------

df.PB = fread(paste0(dir.data, "peripheral_blood/blood.cellmeta.csv"))
df.PB = df.PB[!grep("Neu-P|P-T|M-P|B-M|B-T|Mye-Lym|lowQual|doublet|Doublet", celltype)]
df.PB = df.PB[integrated_snn_res.1.2 != 27] # exclude Neu-Lym from the plot
df.PB$pid = gsub("[ab]$", "", df.PB$pid)

col_pid_PB = c("#F8D626", "#214DC8", "#B0A4E3", "#06A5C7", "#FFA631")
names(col_pid_PB) = sort(unique(df.PB$pid))

col_celltype_PB = c(
  `Immat Neutrophil` = "#0A73DC", `Mature Neutrophil` = "#2E42B8",
  `CD14 Monocyte` = "#8CC269", `CD16 Monocyte` = "#70F3FF", cDC = "#B0A4E3", pDC = "#CCA4E3",
  "Platelet" = "grey88", Erythrocyte = "#FD8CC1FF", HSPC = "#DC3023", "Eos/Baso" = "#EFC000",
  `CD4 T Naive` = "#7397AB", `CD8 T Naive` = "#D2F0F4",
  `CD4 T Eff/Mem` = "#33B7A0", `CD8 T Eff/Mem` = "#F86B1D",
  Treg = "#1A6840", `Prolif T/NK` = "#F9906F", NK = "#FFA631",
  B = "#F9FB0E", PB = "#FFDC91FF", "Neu-Lym" = "#003C67FF")


### study
df.PB$facetvar = "Study"
pointsizen = 0.001

p.umap.pid.PB = plot_scatter(
  df.PB, x = "UMAP_1", y = "UMAP_2", color_by = "pid", color_type = "discrete", colors = col_pid_PB,
  facet_by = "facetvar", facet_nrow = 1, point_size = 0.0001, point_alpha = 0.2, raster_dpi = 400,
  shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL,
  legend_ncol = 1, legend_point_size = 2) +
  # Drawn labels only, as on p.umap.tissue / p.umap.group above.
  scale_color_manual(values = col_pid_PB, labels = cap_first, name = NULL,
                     na.value = "grey42") +
  theme_expresso(legend_position = c(1.01, 0.5), legend_justification = c(0, 0.5),
                 legend_key_spacing_x = 0.01, legend_key_spacing_y = 0.01,
                 legend_key_height = 0.1, legend_text_size = 9,
                 facet_label_face = "plain", show_axis = FALSE, grid = "none")

### group
df.PB$facetvar = "Group"
pointsizen = 0.001

p.umap.group.PB = plot_scatter(
  df.PB, x = "UMAP_1", y = "UMAP_2",
  color_by = "group", color_type = "discrete", colors = col_group,
  facet_by = "facetvar", facet_nrow = 1, point_size = 0.0001, point_alpha = 0.2, raster_dpi = 400,
  shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL,
  legend_ncol = 1, legend_point_size = 2) +
  scale_color_manual(values = col_group, labels = cap_severity, name = NULL,
                     na.value = "grey42") +
  theme_expresso(legend_position = c(1.01, 0.5), legend_justification = c(0, 0.5),
                 legend_key_spacing_x = 0.01, legend_key_spacing_y = 0.01,
                 legend_key_height = 0.1, legend_text_size = 9,
                 facet_label_face = "plain", show_axis = FALSE, grid = "none")

### cell type
df.label = as.data.table(df.PB)[, lapply(.SD, median), by = celltype,
                                .SDcols = c("UMAP_1", "UMAP_2")]
colnames(df.label)[2:3] = c("x", "y")
df.label$labeltext = df.label[, "celltype", with = FALSE]
df.label = df.label[!is.na(celltype)]
df.label = as.data.table(df.label)
df.label[, y := ifelse(celltype == "CD4 T Naive",   y + 1,   y)]
df.label[, y := ifelse(celltype == "CD4 T Eff/Mem", y - 0.5, y)]
df.label[, y := ifelse(celltype == "Prolif T/NK",   y + 0.5, y)]
df.label[, y := ifelse(celltype == "CD16 Monocyte", y - 0.5, y)]
df.label[, x := ifelse(celltype == "CD16 Monocyte", x + 1,   x)]
df.label[, x := ifelse(celltype == "HSPC",          x - 0.5, x)]
df.label[, x := ifelse(celltype == "NK",            x - 1,   x)]

df.PB$facetvar = "Cell type"

p.umap.celltype.PB = plot_scatter(
  df.PB, x = "UMAP_1", y = "UMAP_2",
  color_by = "celltype", color_type = "discrete", colors = col_celltype_PB,
  facet_by = "facetvar", facet_nrow = 1, point_size = 0.001, point_alpha = 0.3, raster_dpi = 300,
  shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL,
  legend_ncol = 1, legend_point_size = 2, label = TRUE, label_df = df.label, label_size = 2.4) +
  theme_expresso(legend_position = c(1.01, 0.5), legend_justification = c(0, 0.5),
                 legend_key_spacing_x = 0.01, legend_key_spacing_y = 0.01,
                 legend_key_height = 0.1, legend_text_size = 9,
                 facet_label_face = "plain", show_axis = FALSE, grid = "none")

### neutrophils vs everything else
df.PB[, celltype.neu.main := ifelse(celltype == "Mature Neutrophil", "Mature Neu", celltype)]
df.PB[, celltype.neu.main := ifelse(celltype == "Immat Neutrophil", "Immat Neu", celltype.neu.main)]
df.PB[, celltype.neu.main := ifelse(celltype == "Neu-Lym", "Neu-Lym", celltype.neu.main)]
df.PB[, celltype.neu.main := ifelse(celltype == "Eos/Baso", "Eos/Baso", celltype.neu.main)]
df.PB[, celltype.neu.main := ifelse(grepl("Mono|cDC|pDC", celltype.neu.main),
                                    "Monocyte and DC", celltype.neu.main)]
df.PB[, celltype.neu.main := ifelse(grepl("CD4|CD8|Treg|NK", celltype.neu.main),
                                    "T and NK", celltype.neu.main)]
df.PB[, celltype.neu.main := ifelse(grepl("Erythrocyte|HSPC", celltype.neu.main),
                                    "other", celltype.neu.main)]
table(df.PB$celltype.neu.main)

col_celltype.neu.main = c(
  other = "grey80", `Monocyte and DC` = "grey80", `T and NK` = "grey80",
  `Mature Neu` = "#2E42B8", PB = "grey80", Platelet = "grey80",
  `Immat Neu` = "#0A73DC", B = "grey80", "Neu-Lym" = "grey80", "Eos/Baso" = "grey80")

df.label.PB = as.data.table(df.PB)[, lapply(.SD, median), by = celltype.neu.main,
                                   .SDcols = c("UMAP_1", "UMAP_2")]
colnames(df.label.PB)[2:3] = c("x", "y")
df.label.PB$labeltext = df.label.PB[, "celltype.neu.main", with = FALSE]
df.label.PB[, x := ifelse(celltype.neu.main == "Mature Neu", x + 1, x)]
df.label.PB[, x := ifelse(celltype.neu.main == "T and NK",   x - 1, x)]
df.label.PB[, y := ifelse(celltype.neu.main == "T and NK",   y - 1, y)]
df.label.PB = df.label.PB[!is.na(celltype.neu.main)][
  !celltype.neu.main %in% "other"]
df.label.PB = as.data.table(df.label.PB)

df.PB$facetvar = "Immune cells in\nperipheral blood"

p.umap.celltype.neu.main = plot_scatter(
  df.PB, x = "UMAP_1", y = "UMAP_2", color_by = "celltype.neu.main", color_type = "discrete",
  colors = col_celltype.neu.main, facet_by = "facetvar", facet_nrow = 1,
  point_size = 0.0001, point_alpha = 0.05, raster_dpi = 400,
  shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL,
  legend_ncol = 1, legend_point_size = 2, label = TRUE, label_df = df.label.PB, label_size = 2.3) +
  theme_expresso(legend_position = "none", facet_label_face = "plain",
                 show_axis = FALSE, grid = "none", panel_background = "white",
                 plot_background = "grey97", plot_margin = ggplot2::margin(0, 0, 0, 0, "line"))


# -- Immune cells in peripheral blood, gene expression -------------------------

mat = readRDS(paste0(dir.data, "peripheral_blood/blood.mat.selectedgenes.rds"))
setdiff(signaturegenes, rownames(mat))

mat = mat[, df.PB$cell.pid]

### gene dotplot
df.PB.genes = cbind(df.PB, as.data.frame(t(as.matrix(mat))))
df.PB.genes$celltype = factor(df.PB.genes$celltype, levels = rev(names(col_celltype_PB)))

p.markerdot.PB = plot_dotplot(
  df.PB.genes, features = signaturegenes, group_by = "celltype",
  max_scale = 2, dot_scale = 3.2, dot_stroke = 0.05, title = NULL,
  col_fontsize = 8, row_fontsize = 8,
  cluster_features = FALSE, cluster_groups = FALSE, show_axis = FALSE,
  feature_side = "top", colorbar_title = "Scaled expression")

### share of reads coming from neutrophils
dtfgenes.sum = readRDS(paste0(dir.data, "peripheral_blood/blood.expr.sum.celltype.rds"))
dtfgenes.sum.test = dtfgenes.sum[
  sampleid.pid %in% dtf.prop[neu.pct >= 0.1]$sampleid.pid][
  , lapply(.SD, function(x) { sum(x) }), by = c("celltype"), .SDcols = c("IL1A", "IL1B", "CAT")]

dtfgenes.sum[, celltype.v2 := ifelse(grepl("Neutrophil", celltype), "Neutrophil", "other")]
dtfgenes.sum.sum = dtfgenes.sum[
  sampleid.pid %in% dtf.prop[neu.pct >= 0.1]$sampleid.pid][
  , lapply(.SD, function(x) { sum(x) }), by = c("celltype.v2"), .SDcols = signaturegenes]

dtfgenes.sum.sum = melt(dtfgenes.sum.sum, id.vars = "celltype.v2",
                        variable.name = "gene", value.name = "reads")
dtfgenes.sum.sum = as.data.table(dtfgenes.sum.sum)
dtfgenes.sum.sum[, pct := reads / sum(reads) * 100, by = gene]
dtfgenes.sum.sum$celltype.v2 = factor(dtfgenes.sum.sum$celltype.v2,
                                      levels = c("other", "Neutrophil"))

p.neu.pct = ggplot(dtfgenes.sum.sum, aes(x = gene, y = pct, fill = celltype.v2)) +
  geom_col() +
  scale_fill_manual(values = c("Neutrophil" = "#2E42B8", "other" = "grey93")) +
  theme_expresso(show_axis = FALSE, grid = "none") +
  theme(legend.position = "none", axis.title.y = ggplot2::element_text(size = 8, angle = 0,
                                             vjust = 0.5, hjust = 1)) +
  labs(x = NULL, y = "Neu reads PCT")

### UMAPs of genes and signature scores
df.PB.genes = cbind(df.PB, as.data.frame(t(as.matrix(mat[
  c("S100A8", "S100A9", "PRTN3", "CTSG", "ELANE", "AZU1", "LTF",
    "ARG1", "ANXA3", "CD177", "MMP9", "PGLYRP1", "OLR1", "CD274",
    "SLPI", "PI3", "MME", "CYP4F3", "CXCL1", "IL1R2", "ORM1"), ]))))

mat = mat + 1

df.PB.genes$immature       = get_gene_scores(mat, genes.immature.neu, "") - 1
df.PB.genes$mature         = get_gene_scores(mat, genes.mature.neu, "") - 1
df.PB.genes$degranulating  = get_gene_scores(mat, genes.degranulating.neu, "") - 1
df.PB.genes$antiprotease   = get_gene_scores(mat, genes.antiprotease.neu, "") - 1
df.PB.genes$total          = get_gene_scores(mat, genes.total.neu, "") - 1
df.PB.genes$`IFN-related`  = get_gene_scores(mat, genes.ifn, "") - 1

pscores.neu.PB = list()
for (gene in c("immature", "degranulating", "antiprotease", "mature", "total", "IFN-related")) {
  # As the SCAHN score loop above: the strip is capitalised, the key is not.
  dtf1 = data.table(UMAP_1 = df.PB.genes$UMAP_1, UMAP_2 = df.PB.genes$UMAP_2,
                    id = df.PB.genes[[gene]], facetvar = cap_first(gene))

  pscores.neu.PB[[gene]] = plot_scatter(
    dtf1, x = "UMAP_1", y = "UMAP_2", color_by = "id",
    facet_by = "facetvar", color_type = "continuous",
    colors = c("grey99", scales::dichromat_pal("DarkRedtoBlue.12")(12)[7:12]),
    quantile_lower = 0.001, quantile_upper = 0.999, shuffle = TRUE, seed = 42,
    point_size = 0.002, point_alpha = 0.2, raster_dpi = 400, legend_ncol = 1, facet_nrow = 1) +
    theme_expresso(legend_position = "none", show_axis = FALSE, grid = "none",
                   facet_label_face = "plain", panel_background = "white",
                   plot_background  = "grey97",
                   plot_margin = ggplot2::margin(0.05, 0.05, 0.05, 0.05, "line"))
}

pgenes.blood = list()
for (gene in c("S100A8", "S100A9", "PRTN3", "CTSG", "ELANE", "AZU1", "LTF",
               "ARG1", "ANXA3", "CD177", "MMP9", "PGLYRP1", "OLR1", "CD274",
               "SLPI", "PI3", "MME", "CYP4F3", "CXCL1", "IL1R2", "ORM1")) {
  dtf1 = data.table(UMAP_1 = df.PB.genes$UMAP_1, UMAP_2 = df.PB.genes$UMAP_2,
                    id = df.PB.genes[[gene]], facetvar = gene)

  pgenes.blood[[gene]] = plot_scatter(
    dtf1, x = "UMAP_1", y = "UMAP_2", color_by = "id",
    facet_by = "facetvar", color_type = "continuous",
    colors = c("grey99", scales::dichromat_pal("DarkRedtoBlue.12")(12)[7:12]),
    quantile_lower = 0.001, quantile_upper = 0.999, shuffle = TRUE, seed = 42,
    point_size = 0.001, point_alpha = 0.2, raster_dpi = 400, legend_ncol = 1, facet_nrow = 1,
    colorbar_height = 1.2, colorbar_width = 0.05) +
    theme_expresso(legend_position = c(1, 0.48), legend_text_size = 6,
                   facet_label_face = "italic", show_axis = FALSE, grid = "none",
                   plot_margin = ggplot2::margin(0.05, 0.5, 0.05, 0.2, "line"))
}


# -- Figure 1 ------------------------------------------------------------------

cairo_pdf(file = paste0(dir.fig, "fig1.umap_pid.pdf"), width = 2.8, height = 2.4)
pushViewport(viewport(layout = grid.layout(nrow = 38, ncol = 41)))
print(p.umap.pid, vp = viewport(layout.pos.row = 2:22, layout.pos.col = 1:18))
dev.off()

cairo_pdf(file = paste0(dir.fig, "fig1.umap_tissue_group_pct.pdf"), width = 8.75, height = 3)
pushViewport(viewport(layout = grid.layout(nrow = 50, ncol = 140)))
print(p.ggalluvial, vp = viewport(layout.pos.row = 21:29, layout.pos.col = 1:85))
print(p.umap.tissue, vp = viewport(layout.pos.row = 1:20, layout.pos.col = 1:18))
print(p.umap.group, vp = viewport(layout.pos.row = 30:49, layout.pos.col = 1:18))
print(p.neu.pct.source2, vp = viewport(layout.pos.row = 4:29, layout.pos.col = 98:136))
print(p.neu.pct.technique2, vp = viewport(layout.pos.row = 30:47, layout.pos.col = 95:136))
dev.off()


# -- Figure 2, legends and bracket helpers -------------------------------------

df.twocircles = data.table(
  x = rep(1, 2), y = 1:2, color = c("grey88", "grey17"), labeltext = c("FDR > 0.05", "FDR ≤ 0.05"))

p.twocircles = ggplot(df.twocircles, aes(x = x, y = y, label = labeltext, color = color)) +
  geom_point(shape = 21, fill = "grey88", size = 2) +
  geom_text(angle = 0, size = 2.6, hjust = -0.2, color = textcolor) +
  theme(text = element_text(size = textsize, color = textcolor)) +
  scale_color_identity() +
  theme_void() +
  coord_cartesian(clip = "off")

circle_data = data.frame(
  x = c(2, 3.8, 6), y = c(2, 2, 2), xtext = c(2, 4, 8.5), radius = c(0.4, 0.65, 0.9),
  label = c("", "", "low"))

pthreecircle = ggplot(data = circle_data) +
  geom_circle(aes(x0 = x, y0 = y, r = radius), linewidth = 0.1, size = 0.1,
              fill = "grey83", color = "grey83") +
  geom_text(aes(x = xtext, y = y, label = label), color = textcolor) +
  theme_void() +
  labs(y = "sd  high", title = NULL) +
  scale_x_continuous(limits = c(1.6, 10)) +
  theme(text = element_text(size = 10, color = textcolor),
        plot.title = element_text(size = 10, hjust = 0.5, color = textcolor),
        axis.title.y = element_text(size = 10, color = textcolor))

df_leg = data.frame(
  label = factor(names(col_celltype), levels = names(col_celltype)), color = col_celltype)

p.subsetnames = ggplot(df_leg, aes(x = 1, y = label, color = label)) +
  geom_point(size = 1) +
  geom_text(aes(x = 1.05, label = label), hjust = 0, size = 3) +
  scale_color_manual(values = col_celltype) +
  scale_x_continuous(limits = c(0.95, 2)) +
  theme_void() +
  theme(
    legend.position = "none", plot.margin = margin(0, 0, 0, 0)) +
  scale_y_discrete(limits = rev(names(col_celltype)))

draw_bracket_label = function(bracket_rows, cols, label, gap = 2,
                              margin = c(0, 0, 0, 0),   # top, right, bottom, left in npc
                              lwd = 0.7, ticks = 0.5, h = 0.3, type = 4, curvature = 0.3,
                              bracket_col = "grey42", textsize = 8, textface = "plain",
                              textcol = "black") {
  text_row = min(bracket_rows) - gap
  all_rows = text_row:max(bracket_rows)
  n_all    = length(all_rows)
  n_brack  = length(bracket_rows)

  mt = margin[1]
  mr = margin[2]
  mb = margin[3]
  ml = margin[4]

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

draw_vbracket_label = function(rows, bracket_cols, text_cols, label, margin = c(0, 0, 0, 0),
                               lwd = 0.7, ticks = 0.5, h = 1.5, type = 4, curvature = 0.3,
                               bracket_col = "grey42", textsize = 8, textface = "plain",
                               textcol = "black", lineheight = 0.8) {
  mt = margin[1]
  mr = margin[2]
  mb = margin[3]
  ml = margin[4]
  text_y = mb + (1 - mt - mb) / 2  # vertical center adjusted for margins

  pushViewport(viewport(layout.pos.row = rows, layout.pos.col = bracket_cols))
  pushViewport(viewport(x = unit(ml, "npc"), y = unit(mb, "npc"), width = unit(1 - ml - mr, "npc"),
                        height = unit(1 - mt - mb, "npc"), just = c("left", "bottom")))
  grid.brackets(0.5, 0, 0.5, 1, lwd = lwd, ticks = ticks, h = unit(h, "mm"),
                type = type, curvature = curvature, col = bracket_col)
  upViewport(2)

  pushViewport(viewport(layout.pos.row = rows, layout.pos.col = text_cols))
  grid.text(1, text_y, label = label, just = "right",
            gp = gpar(fontsize = textsize, fontface = textface,
                      col = textcol, lineheight = lineheight))
  upViewport()
}

draw_nested_brackets = function(inner_rows, outer_rows, all_cols, inner_x = 0.25, outer_x = 0.75,
                                arm_x   = 0.0,   # how far left inner bracket arms extend
                                margin = c(0, 0, 0, 0), lwd = 0.7, col = "grey42") {
  all_rows = min(min(inner_rows), min(outer_rows)):max(max(inner_rows), max(outer_rows))
  n        = length(all_rows)
  top_row  = min(all_rows)

  ry = function(row, edge = "mid") {
    switch(edge, top = 1 - (row - top_row)       / n,
      bottom = 1 - (row - top_row + 1)   / n, 1 - (row - top_row + 0.5) / n)
  }

  inner_top    = ry(min(inner_rows), "top")
  inner_bottom = ry(max(inner_rows), "bottom")
  outer_top    = ry(min(outer_rows), "top")
  outer_bottom = ry(max(outer_rows), "bottom")

  mt = margin[1]
  mr = margin[2]
  mb = margin[3]
  ml = margin[4]
  gp = gpar(lwd = lwd, col = col)

  pushViewport(viewport(layout.pos.row = all_rows, layout.pos.col = all_cols))
  pushViewport(viewport(x = unit(ml, "npc"), y = unit(mb, "npc"), width = unit(1 - ml - mr, "npc"),
                        height = unit(1 - mt - mb, "npc"), just = c("left", "bottom")))

  # Inner ] bracket (arms point left)
  grid.lines(x = c(inner_x, inner_x), y = c(inner_bottom, inner_top),    gp = gp)
  grid.lines(x = c(arm_x,   inner_x), y = c(inner_top,    inner_top),    gp = gp)
  grid.lines(x = c(arm_x,   inner_x), y = c(inner_bottom, inner_bottom), gp = gp)

  # Outer ] bracket (arms point left to inner_x, connecting to inner bracket body)
  grid.lines(x = c(outer_x, outer_x), y = c(outer_bottom, outer_top),    gp = gp)
  grid.lines(x = c(inner_x, outer_x), y = c(outer_top,    outer_top),    gp = gp)
  grid.lines(x = c(inner_x, outer_x), y = c(outer_bottom, outer_bottom), gp = gp)

  upViewport(2)
}


# -- Figure 2 ------------------------------------------------------------------

cairo_pdf(file = paste0(dir.fig, "Figure2.pdf"), width = 15.4, height = 10.5)
pushViewport(viewport(layout = grid.layout(nrow = 146, ncol = 116)))

genes2display = c("STMN1", "AZU1", "MPO", "PRTN3", "CTSG", "ELANE", "DEFA4",
                  "LTF", "MMP8", "CAMP", "ARG1", "MMP9", "OLR1",
                  "S100A8", "S100A4", "IL1R2", "TXNIP", "EGR1", "FOS", "PTGS2",
                  "G0S2", "IFI6", "GBP5", "CD274",
                  "CXCL2", "CXCR4", "VEGFA", "CCL4", "NFKBIA", "IL1B", "IL1RN",
                  "SLPI", "PI3", "HSPA1B", "CD74")

pushViewport(viewport(layout.pos.row = 115:121, layout.pos.col = 2:115))
grid.arrange(grobs = pgenes2[genes2display], nrow = 1, as.table = TRUE, newpage = FALSE)
upViewport(1)

pushViewport(viewport(layout.pos.row = 132:145, layout.pos.col = 5:49))
grid.arrange(grobs = pscores.neu, nrow = 1, as.table = TRUE, newpage = FALSE)
upViewport(1)

pushViewport(viewport(layout.pos.row = 132:145, layout.pos.col = 66:110))
grid.arrange(grobs = pscores.neu.PB, nrow = 1, as.table = TRUE, newpage = FALSE)
upViewport(1)

print(p.umap.celltype.neu.main, vp = viewport(layout.pos.row = 124:145, layout.pos.col = 54:64))

draw_bracket_label(bracket_rows = 128:131, cols = 5:43,
                   label = "Neutrophil-specific gene signatures", gap = 2,
                   margin = c(0.1, 0.04, 0.1, 0.01), textsize = textsize, textcol = textcolor)

draw_bracket_label(bracket_rows = 128:131, cols = 66:104,
                   label = "Neutrophil-specific gene signatures", gap = 2,
                   margin = c(0.1, 0.04, 0.1, 0.01), textsize = textsize, textcol = textcolor)

print(pgenes$FCGR3B,  vp = viewport(layout.pos.row = 2:12,  layout.pos.col = 3:8))
print(pgenes$CXCR2,   vp = viewport(layout.pos.row = 13:23, layout.pos.col = 3:8))
print(pgenes$MME,     vp = viewport(layout.pos.row = 24:34, layout.pos.col = 3:8))
print(pgenes$CEACAM8, vp = viewport(layout.pos.row = 35:45, layout.pos.col = 3:8))
print(pgenes$`granule\ngenes`, vp = viewport(layout.pos.row = 46:59, layout.pos.col = 3:8))

print(p.umap.celltype, vp = viewport(layout.pos.row = 2:41, layout.pos.col = 13:36))
print(p.umap.celltype.matureornot, vp = viewport(layout.pos.row = 40:53, layout.pos.col = 12:20))

print(p.tissuedistribution + theme(axis.text.y = element_blank(),
              plot.margin = margin(0.6, 0.05, 0.4, 0.2, "line")),
      vp = viewport(layout.pos.row = 12:52, layout.pos.col = 42:50))

pushViewport(viewport(layout.pos.row = 3:52, layout.pos.col = 51:83))
draw(p.tissue.ht, merge_legend = TRUE, heatmap_legend_side = "bottom", newpage = FALSE)
upViewport(1)

pushViewport(viewport(layout.pos.row = 1:60, layout.pos.col = 84:99))
draw(ht_prop, merge_legend = TRUE, heatmap_legend_side = "bottom", newpage = FALSE)
upViewport(1)

print(p.subsetnames + theme(plot.margin = margin(0.1, 0.04, 0.01, 0.03, "line")),
      vp = viewport(layout.pos.row = 5:51, layout.pos.col = 36:40))

print(p.tissuesubset.pct + theme(plot.margin = margin(0.01, 0.01, 0.2, 0.35, "line")),
      vp = viewport(layout.pos.row = 2:16, layout.pos.col = 102:115))

print(p.tissuesubsets.pct + theme(plot.margin = margin(0.1, 0.01, 0.01, 0.03, "line")),
      vp = viewport(layout.pos.row = 17:53, layout.pos.col = 101:115))

print(p.twocircles + theme(plot.margin = margin(0, 0, 0, 0, "line")),
      vp = viewport(layout.pos.row = 54:55, layout.pos.col = 94:98))

print(p.markerdot + theme(plot.margin = margin(0.1, 0.6, 0.1, 0.1, "line")),
      vp = viewport(layout.pos.row = 62:113, layout.pos.col = 10:109))
print(p.gene_sets + theme(plot.margin = margin(0.1, 0.6, 0.1, 0.1, "line")),
      vp = viewport(layout.pos.row = 66:75, layout.pos.col = 108:111))

draw_bracket_label(bracket_rows = 62:64, cols = 16:35, label = "Granule genes", gap = 2,
                   margin = c(0.1, 0.04, 0.1, 0.01), textsize = textsize, textcol = textcolor)
draw_bracket_label(bracket_rows = 62:64, cols = 52:72, label = "Interferon-related genes", gap = 2,
                   margin = c(0.1, 0.02, 0.1, 0.025), textsize = textsize, textcol = textcolor)
draw_bracket_label(bracket_rows = 62:64, cols = 73:89,
                   label = "Cytokine/inflammation genes", gap = 2,
                   margin = c(0.1, 0.03, 0.1, 0.025), textsize = textsize, textcol = textcolor)

draw_vbracket_label(rows = 71:74, bracket_cols = 11, text_cols = 2:10,
                    margin = c(0.24, 0.01, 0.08, 0.01), label = "Immature",
                    textsize = textsize, textcol = textcolor)
draw_vbracket_label(rows = 72:76, bracket_cols = 6, text_cols = 1:5,
                    margin = c(0.05, 0.01, 0.1, 0.01), label = "Granule\nsubsets",
                    textsize = textsize, textcol = textcolor)
draw_vbracket_label(rows = 90:95, bracket_cols = 10, text_cols = 1:9,
                    margin = c(0.21, 0.01, 0.05, 0.01), label = "Interferon-\nrelated\nsubsets",
                    textsize = textsize, textcol = textcolor)
draw_vbracket_label(rows = 96:106, bracket_cols = 10, text_cols = 1:9,
                    margin = c(0.05, 0.01, 0.08, 0.01), label = "Cytokine/\ninflammation\nsubsets",
                    textsize = textsize, textcol = textcolor)

draw_nested_brackets(
  inner_rows = 16:20, outer_rows = 14:17, all_cols = 83:84, inner_x = 0.2, outer_x = 0.8)

grid.text(x = unit(0.01, "npc"), y = unit(0.985, "npc"), label = "a",
          gp = gpar(fontsize = textsize, fontface = "bold"))
grid.text(x = unit(0.1, "npc"),  y = unit(0.73, "npc"),  label = "b",
          gp = gpar(fontsize = textsize, fontface = "bold"))
grid.text(x = unit(0.1, "npc"),  y = unit(0.985, "npc"), label = "c",
          gp = gpar(fontsize = textsize, fontface = "bold"))
grid.text(x = unit(0.01, "npc"), y = unit(0.56, "npc"),  label = "d",
          gp = gpar(fontsize = textsize, fontface = "bold"))
grid.text(x = unit(0.01, "npc"), y = unit(0.23, "npc"),  label = "e",
          gp = gpar(fontsize = textsize, fontface = "bold"))

grid.text(x = unit(0.36, "npc"), y = unit(0.985, "npc"), label = "f",
          gp = gpar(fontsize = textsize, fontface = "bold"))
grid.text(x = unit(0.74, "npc"), y = unit(0.985, "npc"), label = "g",
          gp = gpar(fontsize = textsize, fontface = "bold"))
grid.text(x = unit(0.87, "npc"), y = unit(0.985, "npc"), label = "h",
          gp = gpar(fontsize = textsize, fontface = "bold"))
grid.text(x = unit(0.01, "npc"), y = unit(0.14, "npc"),  label = "i",
          gp = gpar(fontsize = textsize, fontface = "bold"))
grid.text(x = unit(0.45, "npc"), y = unit(0.14, "npc"),  label = "j",
          gp = gpar(fontsize = textsize, fontface = "bold"))

print(pthreecircle + theme(plot.margin = margin(0.1, 0.03, 0.1, 0, "line")),
      vp = viewport(layout.pos.row = 53:55, layout.pos.col = 62:71))

dev.off()


# -- Extended figure 2, QC and trajectory --------------------------------------

cairo_pdf(file = paste0(dir.fig, "Fig.e2.pdf"), width = 7.8, height = 7)
pushViewport(viewport(layout = grid.layout(nrow = 100, ncol = 79)))

print(p.nFeature.subset + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.4, "line")),
      vp = viewport(layout.pos.row = 3:30, layout.pos.col = 2:40))
print(p.nCount.subset + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "line")),
      vp = viewport(layout.pos.row = 3:30, layout.pos.col = 41:79))
print(p.ribosome.subset + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.3, "line")),
      vp = viewport(layout.pos.row = 32:59, layout.pos.col = 4:40))
print(p.mt.subset + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.5, "line")),
      vp = viewport(layout.pos.row = 32:59, layout.pos.col = 42:79))

print(p.trajectory, vp = viewport(layout.pos.row = 62:100, layout.pos.col = 8:33))
print(p.pseudotime.subset + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "line")),
      vp = viewport(layout.pos.row = 65:92, layout.pos.col = 43:79))

grid.text(x = unit(0.04, "npc"), y = unit(0.985, "npc"), label = "a",
          gp = gpar(fontsize = textsize, fontface = "bold"))
grid.text(x = unit(0.56, "npc"), y = unit(0.985, "npc"), label = "b",
          gp = gpar(fontsize = textsize, fontface = "bold"))
grid.text(x = unit(0.04, "npc"), y = unit(0.7, "npc"),  label = "c",
          gp = gpar(fontsize = textsize, fontface = "bold"))
grid.text(x = unit(0.56, "npc"), y = unit(0.7, "npc"),  label = "d",
          gp = gpar(fontsize = textsize, fontface = "bold"))
grid.text(x = unit(0.04, "npc"), y = unit(0.4, "npc"),   label = "e",
          gp = gpar(fontsize = textsize, fontface = "bold"))
grid.text(x = unit(0.56, "npc"), y = unit(0.38, "npc"),  label = "f",
          gp = gpar(fontsize = textsize, fontface = "bold"))

dev.off()


# -- Supplementary figure 1, QC ------------------------------------------------

cairo_pdf(file = paste0(dir.fig, "Fig.s1.pdf"), width = 6.5, height = 4)
pushViewport(viewport(layout = grid.layout(nrow = 44, ncol = 65)))

print(p.nFeature + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "line")),
      vp = viewport(layout.pos.row = 1:18, layout.pos.col = 3:30))
print(p.nCount + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.23, "line")),
      vp = viewport(layout.pos.row = 19:43, layout.pos.col = 2:30))

print(plot_forest(as_forest_data(meta.nFeaturenCount.BDvs10x,
                                 gene = "nCount_RNA", rowname_style = "name"),
                  boxsize = 0.5, plottitle = " "),
      vp = viewport(layout.pos.row = 4:15, layout.pos.col = 49:64))

print(plot_forest(as_forest_data(meta.nFeaturenCount.BDvs10x, gene = "nFeature_RNA"),
                  boxsize = 0.5, plottitle = "BD Rhapsody vs. 10X Chromium") +
        theme(axis.text.y = element_blank()),
      vp = viewport(layout.pos.row = 4:15, layout.pos.col = 40:48))

print(plot_forest(as_forest_data(meta.nFeaturenCount.WBvsPBMC,
                                 gene = "nCount_RNA", rowname_style = "name"),
                  boxsize = 0.5, plottitle = " ") + scale_x_continuous(breaks = seq(-3, 3, 1)),
      vp = viewport(layout.pos.row = 16:26, layout.pos.col = 49:65))

print(plot_forest(as_forest_data(meta.nFeaturenCount.WBvsPBMC, gene = "nFeature_RNA"),
                  boxsize = 0.5, plottitle = "WB vs. PBMC") + theme(axis.text.y = element_blank()),
      vp = viewport(layout.pos.row = 16:26, layout.pos.col = 40:48))

print(plot_forest(as_forest_data(meta.nFeaturenCount.freshvsfrozenPBMC,
                                 gene = "nCount_RNA", rowname_style = "name"),
                  boxsize = 0.5, plottitle = " "),
      vp = viewport(layout.pos.row = 27:37, layout.pos.col = 49:65))

print(plot_forest(as_forest_data(meta.nFeaturenCount.freshvsfrozenPBMC, gene = "nFeature_RNA"),
                  boxsize = 0.5, plottitle = "Fresh vs. frozen PBMC") +
        theme(axis.text.y = element_blank()),
      vp = viewport(layout.pos.row = 27:37, layout.pos.col = 40:48))


grid.text(x = unit(0.02, "npc"), y = unit(0.97, "npc"), label = "a",
          gp = gpar(fontsize = textsize, fontface = "bold"))
grid.text(x = unit(0.02, "npc"), y = unit(0.54, "npc"),  label = "b",
          gp = gpar(fontsize = textsize, fontface = "bold"))
grid.text(x = unit(0.6, "npc"),  y = unit(0.97, "npc"), label = "c",
          gp = gpar(fontsize = textsize, fontface = "bold"))

dev.off()


# -- Extended figure 3 ---------------------------------------------------------

cairo_pdf(file = paste0(dir.fig, "Fig.e3.pdf"), width = 14, height = 7)
pushViewport(viewport(layout = grid.layout(nrow = 195, ncol = 200)))

print(p.umap.pid.PB, vp = viewport(layout.pos.row = 4:42, layout.pos.col = 5:23))
print(p.umap.group.PB, vp = viewport(layout.pos.row = 45:83, layout.pos.col = 5:23))
print(p.umap.celltype.PB, vp = viewport(layout.pos.row = 4:83, layout.pos.col = 42:81))

pushViewport(viewport(layout.pos.row = 4:86, layout.pos.col = 102:198))
grid.arrange(grobs = pgenes.blood, nrow = 3, as.table = TRUE, newpage = FALSE)
upViewport(1)

print(p.markerdot.PB, vp = viewport(layout.pos.row = 85:183, layout.pos.col = 2:199))
print(p.neu.pct + theme(plot.margin = margin(0.04, 0.3, 0.04, 0.2, "line")),
      vp = viewport(layout.pos.row = 184:194, layout.pos.col = 3:193))

grid.text(x = unit(0.012, "npc"), y = unit(0.99, "npc"), label = "a",
          gp = gpar(fontsize = textsize, fontface = "bold"))
grid.text(x = unit(0.012, "npc"), y = unit(0.78, "npc"), label = "b",
          gp = gpar(fontsize = textsize, fontface = "bold"))
grid.text(x = unit(0.2, "npc"),   y = unit(0.99, "npc"), label = "c",
          gp = gpar(fontsize = textsize, fontface = "bold"))
grid.text(x = unit(0.5, "npc"),   y = unit(0.99, "npc"), label = "d",
          gp = gpar(fontsize = textsize, fontface = "bold"))

grid.text(x = unit(0.012, "npc"), y = unit(0.51, "npc"), label = "e",
          gp = gpar(fontsize = textsize, fontface = "bold"))
grid.text(x = unit(0.012, "npc"), y = unit(0.05, "npc"), label = "f",
          gp = gpar(fontsize = textsize, fontface = "bold"))

dev.off()
