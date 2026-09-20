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
# Writes to figures/post_acceptance/:
#   Figure1.pdf   Figure 1 (flowchart, UMAPs, alluvial, subset proportions)
#   Figure2.pdf   Figure 2
#   Fig.e2.pdf    Extended figure 2 (QC, trajectory)
#   Fig.s1.pdf    Supplementary figure 1 (QC)
#   Fig.e3.pdf    Extended figure 3 (peripheral blood)
# ==============================================================================

SCAHN_SCRIPTS = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/scripts/figures_post_acceptance"
source(file.path(SCAHN_SCRIPTS, "00_setup.R"))
source(file.path(SCAHN_SCRIPTS, "01_load_data.R"))
options(future.globals.maxSize= 20000 * 1024^2)

# The subset named for sequencing depth is keyed "lowDepth" upstream and drawn
# "Low depth". Renamed on the keys themselves rather than at each label, so it
# reaches every panel at once: they all read the string off dtf$celltype,
# dtf.prop's columns or names(col_celltype), and cap_first() leaves it alone.
# Local to this script -- 00_setup.R and 01_load_data.R are shared, and
# 05_sexdifference.R still subsets on "lowDepth".
names(col_celltype)[names(col_celltype) == "lowDepth"] = "Low depth"
dtf[celltype == "lowDepth", celltype := "Low depth"]
setnames(dtf.prop, "lowDepth", "Low depth")


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


# -- SCAHN UMAPs and ggalluvial ------------------------------------------------

### tissue-by-group alluvial
dtf.summary = as.data.table(table(dtf$tissue, dtf$group))
dtf.summary = dtf.summary[N > 0]
colnames(dtf.summary)[1:2] = c("tissue", "group")
dtf.summary$tissue = factor(dtf.summary$tissue, levels = rev(names(col_tissue)))
dtf.summary$group  = factor(dtf.summary$group,  levels = rev(names(col_group)))

col_tissue_group = c(col_tissue, col_group)

dtf.summary_long = to_lodes_form(data.frame(dtf.summary),
                                 key = "category", value = "item", id = "index", axes = 1:2)
dtf.summary_long = as.data.table(dtf.summary_long)
dtf.summary_long[, colflow := ifelse(category == "group", item, NA)]

nolabel = c("other", "other GI organs", "sputum", "brain", "kidney", "stomach",
            "other lung disease", "other cancer", "other disease")

twoline = c("Peripheral blood" = "Peripheral\nblood")

cap_tissue = function(x) sub("^Colonrectum$", "Colorectum", cap_first(x))

stratum_label = function(x, capitalise) {
  x = as.character(x)
  drawn = capitalise(x)
  drawn = ifelse(drawn %in% names(twoline), twoline[drawn], drawn)
  ifelse(x %in% nolabel, NA, drawn)
}

stratum_labels = function(cat, side, capitalise)
  geom_label_repel(data = dtf.summary_long[category == cat], stat = "stratum",
                   aes(label = after_stat(stratum_label(stratum, capitalise))),
                   size = pt2mm(textsize), color = textcolor, lineheight = 0.8,
                   fill = NA, linewidth = 0, label.padding = unit(0.5, "pt"),
                   direction = "y", hjust = if (side < 0) 1 else 0, nudge_x = side * 0.25,
                   xlim = if (side < 0) c(-Inf, NA) else c(NA, Inf),
                   min.segment.length = 0, segment.size = linesize, segment.color = "grey70",
                   box.padding = 0.02, seed = 42, na.rm = TRUE)

fmt_cells = function(x) ifelse(x == 0, "0", ifelse(x >= 1e6, paste0(x / 1e6, "M"),
                                                   paste0(x / 1e3, "k")))

p.ggalluvial = ggplot(
  data = dtf.summary_long, aes(x = category, stratum = item, alluvium = index, y = N)) +
  geom_flow(aes(fill = item), width = 1/5) +
  geom_stratum(aes(fill = item), width = 1/5, size = .1) +
  scale_fill_manual(values = col_tissue_group) +
  stratum_labels("tissue", -1, cap_tissue) +
  stratum_labels("group",  +1, cap_severity) +
  scale_x_discrete(labels = c(tissue = "Tissue", group = "Condition"),
                   expand = expansion(add = c(1.67, 2.45))) +
  scale_y_continuous(breaks = seq(0, 1e6, 2.5e5), labels = fmt_cells, expand = c(0, 0)) +
  labs(x = NULL, y = "Number of cells") +
  theme_expresso(legend_position = "none", show_axis = TRUE, grid = "none") +
  coord_cartesian(clip = "off")


### tissue
dtf$facetvar = "Tissue"

label_n = function(n, capitalise = cap_first) function(x) {
  k = n[x] / 1000
  paste0(capitalise(x), " (", ifelse(k < 10, sprintf("%.1fk", k), sprintf("%.0fk", k)), ")")
}

p.umap.tissue = plot_scatter(
  dtf, x = "UMAP_1", y = "UMAP_2", facet_by = "facetvar",
  color_by = "tissue", color_type = "discrete", colors = col_tissue,
  facet_nrow = 1, point_size = .005, point_alpha = 0.3, raster_dpi = 500,
  legend_ncol = 1, legend_point_size = 1,
  shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL) +
  scale_color_manual(values = col_tissue, labels = label_n(table(dtf$tissue), cap_tissue),
                     name = NULL, na.value = "grey42") +
  theme_expresso(legend_position = c(1.01, 1), legend_justification = c(0, 1),
                 legend_title = FALSE,
                 legend_key_spacing_x = 0.01, legend_key_spacing_y = 0.01,
                 legend_key_height = 0.1, facet_label_face = "plain", show_axis = FALSE,
                 grid = "none") +
  theme(legend.title = element_text(margin = margin(b = 0.1, unit = "line")),
        legend.text  = element_text(margin = margin(-0.04, 0, -0.04, 0.1, "line")))

### condition
dtf$facetvar = "Condition"

p.umap.group = plot_scatter(
  dtf, x = "UMAP_1", y = "UMAP_2", color_by = "group", facet_by = "facetvar",
  color_type = "discrete", colors = col_group,
  facet_nrow = 1, point_size = .005, point_alpha = 0.3, raster_dpi = 500,
  legend_ncol = 1, legend_point_size = 1,
  shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL) +
  scale_color_manual(values = col_group, labels = label_n(table(dtf$group), cap_severity),
                     name = NULL, na.value = "grey42") +
  theme_expresso(legend_position = c(1.01, 1), legend_justification = c(0, 1),
                 legend_title = TRUE,
                 legend_key_spacing_x = 0.01, legend_key_spacing_y = 0.01,
                 legend_key_height = 0.1, facet_label_face = "plain", show_axis = FALSE,
                 grid = "none") +
  theme(legend.title = element_text(margin = margin(b = 0.1, unit = "line")),
        legend.text  = element_text(margin = margin(-0.04, 0, -0.04, 0.1, "line")))

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
                 legend_text_size = textsize, facet_label_face = "plain", show_axis = FALSE,
                 grid = "none")

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

p.umap.celltype = plot_scatter(
  dtf, x = "UMAP_1", y = "UMAP_2",
  color_by = "celltype", color_type = "discrete", colors = col_celltype,
  facet_nrow = 1, point_size = 0.001, point_alpha = 0.3, raster_dpi = 300,
  shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL,
  label = TRUE, label_df = df.label, label_size = pt2mm(textsize)) +
  theme_expresso(legend_position = "none", facet_label_face = "plain",
                 show_axis = FALSE, grid = "none")

### mature vs immature
col_matureornot = c(Immature = "#F47B00", Mature = "grey87")

dtf[, celltype_matureornot := ifelse(celltype %in% c("AZU1", "LTF"), "Immature", "Mature")]

df.label2 = as.data.table(dtf)[, lapply(.SD, median), by = celltype_matureornot,
                               .SDcols = c("UMAP_1", "UMAP_2")]
colnames(df.label2)[2:3] = c("x", "y")
df.label2$labeltext = df.label2[, "celltype_matureornot", with = FALSE]
df.label2[, x := ifelse(celltype_matureornot == "Immature", x + 0.8, x)]

dtf$facetvar = "mature state"

p.umap.celltype.matureornot = plot_scatter(
  dtf, x = "UMAP_1", y = "UMAP_2", color_by = "celltype_matureornot", color_type = "discrete",
  colors = col_matureornot, facet_by = "facetvar", facet_nrow = 1,
  point_size = 0.000001, point_alpha = 0.1, raster_dpi = 400,
  shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL,
  label = TRUE, label_df = df.label2, label_size = pt2mm(textsize)) +
  theme_expresso(legend_position = "none", show_axis = FALSE, grid = "none",
                 show_facet_label = FALSE)


# -- SCAHN UMAPs of genes and signature scores ---------------------------------

mat = readRDS(paste0(dir.results, "annotation/SCAHN.selectedgenes.mat.expr.integrated.rds"))
mat = as.matrix(mat)
mat = mat[, dtf$cell.pid]

dtfgenes = as.data.table(cbind(dtf, t(mat)))

dtf.score = readRDS(paste0(dir.results, "annotation/SCAHN.cellscores.rds"))
dtf.score = dtf.score[cell.pid %in% dtf$cell.pid]

dtfgenes$`Granule\ngenes` = dtf.score$score.granules - 1
dtfgenes$immature         = dtf.score$score.immature.neu - 1
dtfgenes$degranulating    = dtf.score$score.degranulating.neu - 1
dtfgenes$antiprotease     = dtf.score$score.antiprotease.neu - 1
dtfgenes$mature           = dtf.score$score.mature.neu - 1
dtfgenes$total            = dtf.score$score.total.neu - 1
dtfgenes$`IFN-related`    = dtf.score$score.ifn - 1

### signature scores, one panel per signature

facetlabel.scores = c(
  "immature"      = "Immature\nNeu",
  "degranulating" = "Degranulating\nNeu",
  "antiprotease"  = "Antiprotease\nNeu",
  "mature"        = "Mature\nNeu",
  "total"         = "Total\nNeu",
  "IFN-related"   = "IFN\nrelated")

pscores.neu = list()
for (gene in c("immature", "degranulating", "antiprotease", "mature", "total", "IFN-related")) {
  dtf1 = data.table(UMAP_1 = dtfgenes$UMAP_1, UMAP_2 = dtfgenes$UMAP_2,
                    id = dtfgenes[[gene]], facetvar = facetlabel.scores[[gene]])

  pscores.neu[[gene]] = plot_scatter(
    dtf1, x = "UMAP_1", y = "UMAP_2", color_by = "id",
    facet_by = "facetvar", color_type = "continuous",
    colors = c("grey99", scales::dichromat_pal("DarkRedtoBlue.12")(12)[7:12]),
    quantile_lower = 0.001, quantile_upper = 0.999, shuffle = TRUE, seed = 42,
    point_size = 0.002, point_alpha = 0.2, raster_dpi = 400, legend_ncol = 1, facet_nrow = 1) +
    theme_expresso(legend_position = "none", show_axis = FALSE, grid = "none",
                   facet_label_face = "plain", panel_background = "white",
                   plot_background  = "white",
                   plot_margin = ggplot2::margin(0.05, 0.05, 0.05, 0.05, "line"))
}

### individual maturity markers
pgenes = list()
for (gene in c("FCGR3B", "CXCR2", "MME", "CEACAM8", "Granule\ngenes")) {
  dtf1 = data.table(UMAP_1 = dtfgenes$UMAP_1, UMAP_2 = dtfgenes$UMAP_2,
                    id = dtfgenes[[gene]], facetvar = gene)

  pgenes[[gene]] = plot_scatter(
    dtf1, x = "UMAP_1", y = "UMAP_2", color_by = "id",
    facet_by = "facetvar", color_type = "continuous",
    colors = c("grey99", scales::dichromat_pal("DarkRedtoBlue.12")(12)[7:12]),
    quantile_lower = 0.001, quantile_upper = 0.999, shuffle = TRUE, seed = 42,
    point_size = 0.0001, point_alpha = 0.1, raster_dpi = 400,
    legend_ncol = 1, colorbar_height = 0.9, colorbar_width = 0.05, facet_nrow = 1) +
    theme_expresso(legend_position = c(0.99, 0.48), legend_text_size = textsize,
                   facet_label_face = "italic", show_axis = FALSE, grid = "none")
}

### subset-defining genes, small multiples
genes2display = c("STMN1", "MKI67", "AZU1", "MPO", "PRTN3", "CTSG", "ELANE", "DEFA4",
                  "LTF", "MMP8", "CAMP", "ARG1", "MMP9", "IL1R2",
                  "CYP4F3", "S100A8", "S100A4", "OLR1", "TXNIP", "EGR1", "FOS", "PTGS2",
                  "G0S2", "IFI6", "GBP5", "CD274",
                  "CXCL2", "CXCL8", "CXCR4", "VEGFA", "CCL4", "NFKBIA", "IL1B", "IL1RN",
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
      strip.text.x = element_text(size = textsize, color = textcolor,
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
  geom_boxplot(width = boxwidth, size = 0.3, outlier.shape = NA,
               alpha = 1, fill = NA, color = "grey17") +
  theme_expresso(grid = "none", axis_text_size = textsize) +
  theme(legend.position = "none",
        plot.margin     = unit(c(.1, .1, .2, .1), "line"),
        axis.text.x     = element_text(angle = 30, hjust = 1)) +
  labs(x = NULL, y = "% of Neu", title = NULL) +
  guides(color = guide_legend(ncol = 1, label.hjust = 0,
                              override.aes = list(size = 2, alpha = 1))) +
  scale_x_discrete(labels = cap_first) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = function(x) x * 100) #+


dtf.prop.blood$technique = as.character(dtf.prop.blood$technique)
dtf.prop.blood[, technique2 := ifelse(technique %in% c("10x Chromium", "BD Rhapsody"),
                                      technique, "other")]

p.neu.pct.technique2 = ggplot(dtf.prop.blood[grepl("freshWB", source)],
                              aes(x = technique2, y = neu.pct)) +
  ggbeeswarm::geom_quasirandom(size = 0.2, width = 0.2, show.legend = TRUE,
                               alpha = 0.7, varwidth = FALSE, shape = 16, color = "grey42") +
  geom_boxplot(width = boxwidth, size = 0.3, outlier.shape = NA,
               alpha = 1, fill = NA, color = "grey17") +
  theme_expresso(axis_text_size = textsize) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1),
        plot.margin = unit(c(.1, .1, .2, .23), "line")) +
  labs(x = NULL, y = "% of Neu\nin whole blood", title = NULL) +
  guides(color = guide_legend(ncol = 1, label.hjust = 0,
                              override.aes = list(size = 2, alpha = 1))) +
  scale_x_discrete(labels = function(x) sub(" ", "\n", cap_first(x))) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = function(x) x * 100)


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
  scale_x_discrete(position = "top") +
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
  scale_x_discrete(position = "top", limits = rev) +
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
    !tissue %in% c("other", "peripheral blood") &
      !variable %in% c("immature", "cytokine/inflammation subsets")
  ], tissue ~ variable, value.var = "value")
rowname = dtf.prop.m.agg.combined.mat$tissue
dtf.prop.m.agg.combined.mat$tissue = NULL
dtf.prop.m.agg.combined.mat = as.matrix(dtf.prop.m.agg.combined.mat)
rownames(dtf.prop.m.agg.combined.mat) = rowname
dtf.prop.m.agg.combined.mat = dtf.prop.m.agg.combined.mat[
  , intersect(names(col_celltype), colnames(dtf.prop.m.agg.combined.mat))]

dtf.prop.m.aggsd.combined.mat = dcast(
  dtf.prop.m.aggsd.combined[
    !tissue %in% c("other", "peripheral blood") &
      !variable %in% c("immature", "tissuesubset")],
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

dtf.prop.m.agg.combined.mat   = t(dtf.prop.m.agg.combined.mat)
dtf.prop.m.aggsd.combined.mat = t(dtf.prop.m.aggsd.combined.mat)
mat.size = 1 - sqrt(dtf.prop.m.aggsd.combined.mat)

p.tissue.ht = Heatmap(
  dtf.prop.m.agg.combined.mat,
  rect_gp = gpar(type = "none"), col = col_fun.tissueht,
  na_col = "grey93", cluster_columns = FALSE, show_row_names = TRUE,
  cluster_rows    = FALSE, show_row_dend = TRUE, show_column_dend = TRUE,
  row_dend_side   = "left", column_dend_side = "top",
  clustering_method_rows      = "ward.D2", clustering_method_columns   = "ward.D2",
  clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
  row_names_side = "left", column_names_side = "bottom", column_names_rot = 45,
  # column_labels rather than renaming the matrix: cell_fun indexes the median and
  # size matrices by [i, j], and the tissue strings are the same ones the dcast
  # above keys on, so only what is drawn is capitalised and renamed.
  column_labels = cap_tissue(colnames(dtf.prop.m.agg.combined.mat)),
  column_names_gp = gpar(fontsize = textsize, fontface = "plain", col = textcolor),
  row_names_gp = gpar(fontsize = textsize, fontface = "plain", col = textcolor),
  row_title = "Neutrophil subset",
  row_title_gp = gpar(fontsize = textsize, fontface = "plain", col = textcolor),
  column_title = NULL,
  column_title_side = "top",
  column_title_gp = gpar(fontsize = textsize, fontface = "plain", col = textcolor),
  show_heatmap_legend = FALSE, cell_fun = function(j, i, x, y, width, height, fill) {
    grid.circle(x = x, y = y, r = 0.6 * mat.size[i, j] * min(unit.c(width, height)),
                gp = gpar(fill = col_fun.tissueht(dtf.prop.m.agg.combined.mat[i, j]),
                          col = NA, alpha = 0.9))
  })

lgd.tissue.ht = Legend(
  col_fun = col_fun.tissueht, at = c(0, 0.25, 0.5),
  title = "median proportion", title_gp = gpar(fontsize = textsize, col = textcolor),
  title_position = "leftcenter-rot", title_gap = unit(0.3, "line"),
  labels_gp = gpar(fontsize = textsize, col = textcolor),
  legend_width = unit(.2, "line"), legend_height = unit(3.2, "line"),
  grid_width = unit(.2, "line"), grid_height = unit(2, "line"),
  by_row = FALSE, direction = "vertical")

### boxplot, solid tissue vs whole blood / purified neutrophil

col_tissuetype = c("whole blood / purified neutrophil" = "#EE6C00", "solid tissue" = "grey17")

p.tissuesubsets.pct = list()
for (i in c("CXCL", "VEGFA", "CCL3/4", "NF-κB", "IL1B", "IL1RN")){
  p.tissuesubsets.pct[[i]] = ggplot(
  dtf.tmp[variable %in% i][
    tissue %in% c("solid tissue", "whole blood / purified neutrophil")],
  aes(x = tissue, y = 100*value + 1)) +
  ggbeeswarm::geom_quasirandom(aes(color = tissue), size = .05, width = 0.2, show.legend = TRUE,
                               alpha = 0.2, varwidth = FALSE, shape = 16) +
  geom_boxplot(aes(color = tissue), width = boxwidth, size = 0.2,
               outlier.shape = NA, alpha = 1, fill = NA) +
  labs(y = paste0(i, " (% of Neu)"), x = NULL, title = NULL) +
  theme_expresso(axis_text_size = textsize, grid = "none",
                 legend_position      = "none") +

  scale_y_log10(breaks = c(3, 30), expand = expansion(mult = c(0.05, 0.18))) +
  expand_limits(y = c(3, 30)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.key.size = unit(1, "lines"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.2, "line"),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
        plot.title = element_text(size = textsize, face = "plain", hjust = 0,
                                  margin = margin(t = 1.6, r = 0, b = .1, l = 0, unit = "pt"))) +
  scale_color_manual(values = col_tissuetype) +
  stat_compare_means(comparisons = list(c("whole blood / purified neutrophil", "solid tissue")),
                     method = "wilcox.test", vjust = 0, size = 2,
                     label = "p.format", color = "grey17", bracket.size = 0.1, tip.length = 0) +
  guides(color = guide_legend(nrow = 2, label.hjust = 0,
                              override.aes = list(size = 0.5, alpha = 1)))

}


p.tissuesubset.pct = ggplot(
  dtf.tmp[variable %in% c("cytokine/inflammation subsets")][
    tissue %in% c("solid tissue", "whole blood / purified neutrophil")],
  aes(x = tissue, y = 100*value + 1)) +
  ggbeeswarm::geom_quasirandom(aes(color = tissue), size = .05, width = 0.2, show.legend = TRUE,
                               alpha = 0.2, varwidth = FALSE, shape = 16) +
  geom_boxplot(aes(color = tissue), width = boxwidth, size = 0.2,
               outlier.shape = NA, alpha = 1, fill = NA) +
  labs(y = "Cytokine/\ninflammation\nsubsets (% of Neu)", x = NULL, title = NULL) +
  theme_expresso(axis_text_size = textsize, grid = "none", legend_position = "none") +
  scale_y_log10(expand = expansion(mult = c(0.05, 0.18))) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.2, "line"),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
        plot.title   = element_text(size = textsize, face = "plain", hjust = 0,
                                    margin = margin(t = .1, r = 0, b = .1, l = 0, unit = "pt"))) +
  scale_color_manual(values = col_tissuetype) +
  stat_compare_means(comparisons = list(c("whole blood / purified neutrophil", "solid tissue")),
                     method = "wilcox.test", vjust = 0, size = 2,
                     label = "p.format", color = "grey17", bracket.size = 0.1, tip.length = 0) +
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
# Kept as a separate object rather than narrowing dtf.prop in place. The tissues
# outside these six categories -- BALF and sputum keep tissue2 = tissue, so they are
# among them -- are needed by the panels built earlier from dtf.prop, and overwriting
# it here silently dropped their columns whenever an earlier block was re-run in a
# live session.
dtf.prop.tissue2 = dtf.prop[tissue2 %in% c("bone marrow", "cord blood", "PBMC", "solid tissue",
                                           "whole blood / purified neutrophil", "spleen")]
dim(dtf.prop.tissue2)
table(dtf.prop.tissue2$tissue2)
colnames(dtf.prop.tissue2)

subsetnames = names(col_celltype)

base_tissue   = "whole blood / purified neutrophil"
other_tissues = setdiff(unique(dtf.prop.tissue2$tissue2), base_tissue)

prop_stats = rbindlist(lapply(subsetnames, function(ss) {
  rbindlist(lapply(other_tissues, function(tiss) {
    x = dtf.prop.tissue2[tissue2 == tiss,        get(ss)]
    y = dtf.prop.tissue2[tissue2 == base_tissue, get(ss)]
    x = x[is.finite(x)]
    y = y[is.finite(y)]
    if (length(x) < 2 || length(y) < 2) return(NULL)

    # Hedges' g via expresso::gene_hedges_g (1-row matrix: tissue=1, base=0)
    vals   = c(x, y)
    labels = c(rep(1L, length(x)), rep(0L, length(y)))
    mat    = matrix(vals, nrow = 1, dimnames = list(ss, NULL))
    hg     = expresso::gene_hedges_g(mat, labels)

    wt = wilcox.test(x, y, exact = FALSE)
    #wt = t.test(x, y)

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
    v = val_mat[i, j]
    f = fdr_mat[i, j]
    if (is.na(v) || is.na(f)) return(invisible(NULL))
    r_base = dot_scale * min(unit.c(width, height))
    grid.circle(x, y, r = sz_mat[i, j] * r_base, gp = gpar(fill = cf(v), col = NA, alpha = 0.9))
    grid.circle(x, y, r = sz2_mat[i, j] * r_base,
                gp = gpar(fill = NA, col = "grey17", alpha = 1, lwd = 0.4))
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
  hg_mat, rect_gp = gpar(type = "none"),
  col = cf_hg, na_col = "grey", cluster_rows = FALSE, cluster_columns = FALSE,
  clustering_method_rows = "ward.D", clustering_method_columns = "ward.D",
  clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
  show_row_dend = FALSE, show_column_dend = TRUE, column_dend_side = "bottom",
  show_column_names = TRUE, column_names_side = "bottom", row_names_side = "left",
  column_names_rot = 45, column_names_gp = gpar(fontsize = textsize, col = textcolor),
  column_labels = cap_first(colnames(hg_mat)),
  row_names_gp = gpar(fontsize = textsize, col = textcolor), row_names_max_width = unit(40, "line"),
 #column_title = "proportion vs. WB",
  column_title = NULL, column_title_side = "top",
  column_title_gp = gpar(fontsize = textsize, fontface = "plain", col = textcolor),
  row_title = "Neutrophil subset",
  row_title_gp = gpar(fontsize = textsize, fontface = "plain", col = textcolor),
  show_heatmap_legend = FALSE, cell_fun = make_cellfun(hg_mat, fdr_mat_prop, cf_hg,
                          dot_scale = 1.3, fdr_cutoff = 0.05))

# Built here rather than through heatmap_legend_param for title_gap; see lgd.tissue.ht.
lgd.prop = Legend(
  col_fun = cf_hg, at = round(seq(-clamp_hg, clamp_hg, length.out = 3), 1),
  title = "Effect size",
  title_gp = gpar(fontsize = textsize, col = textcolor),
  title_position = "leftcenter-rot", title_gap = unit(0.3, "line"),
  labels_gp = gpar(fontsize = textsize, col = textcolor),
  legend_width = unit(.2, "line"), legend_height = unit(3.2, "line"),
  grid_width = unit(.2, "line"), grid_height = unit(2, "line"),
  direction = "vertical")


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
  col_fontsize = textsize, row_fontsize = textsize,
  cluster_features = FALSE, cluster_groups = FALSE, show_axis = FALSE,
  feature_side = "top", colorbar_title = "scaled expression") +
  #scale_x_discrete(labels = label_map, position = "top") +
  scale_x_discrete(position = "top") +
  theme(
    legend.position = c(1, 0.4), axis.text.x.top = element_markdown(
      angle = 90, hjust = 0, vjust = 0.5, size = textsize, margin = margin(b = 0.1, unit = "line")))

# The same dotplot with the axes swapped: genes on rows, subsets on columns.
# plot_dotplot() hardcodes x = feature, y = .group, so the mapping is re-pointed
# at the two columns it builds rather than flipped with coord_flip(), which would
# leave the constructor's rotated italic x-axis theme on the subset axis. Both
# scales take limits = rev because the constructor orders .group by the input
# factor levels (rev(col_celltype)) and feature by signaturegenes; reversing puts
# subsets left to right in col_celltype order and genes top to bottom in
# signaturegenes order, the reading order p.markerdot above has.
p.markerdot.t = plot_dotplot(
  dtfexpr, features = signaturegenes, group_by = "celltype",
  max_scale = 2, dot_scale = 2.3, dot_stroke = 0.05, title = NULL,
  col_fontsize = textsize, row_fontsize = textsize,
  cluster_features = FALSE, cluster_groups = FALSE, show_axis = FALSE,
  feature_side = "top", colorbar_title = "Scaled expression") +
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

# plot_dotplot() fixes the dot outline at colour = "grey67" in its single
# geom_point layer, so the outline goes by dropping that colour and mapping colour
# to expr through the same gradient as the fill -- each ring in the colour of its
# own dot. Not colour = NA (ggplot2 4.x reads it as missing and drops the row, dot
# and all), and not stroke = 0 or after_scale(fill) (both blank the "% expressed"
# keys, which come from this layer). The gradient is restated, so keep it in step
# with colors/max_scale in the plot_dotplot() call above.
p.markerdot.t$layers[[1]]$aes_params$colour = NULL
p.markerdot.t = p.markerdot.t + aes(colour = expr) +
  scale_colour_gradientn(colours = expresso_colors("diverging_bwr2"),
                         limits = c(-2, 2), guide = "none")

df.gene_sets = data.table(
  subset = c("neutrophil-\nspecific\ngene signatures", "- immature",
             "- degranulating", "- antiprotease", "- mature"),
  color = c("#252525", names(gene_sets)), x = 1, y = (length(gene_sets) + 1):1)

p.gene_sets = ggplot(df.gene_sets, aes(x = x, y = y, label = subset, color = color)) +
  geom_text(angle = 0, size = pt2mm(textsize), hjust = 0, lineheight = 0.8, vjust = 0) +
  scale_color_identity() +
  theme_void() +
  coord_cartesian(clip = "off")


# -- QC: nFeature, nCount, ribosome, MT ----------------------------------------

dtf[, technique2 := ifelse(grepl("10x Chromium", technique), "10x Chromium",
                           ifelse(grepl("BD Rhapsody", technique), "BD Rhapsody", "other"))]
table(dtf$technique2)

dtf.score = readRDS(paste0(dir.results, "annotation/SCAHN.cellscores.rds"))
dtf$`Ribosomal protein (RPL/RPS) gene expression` = dtf.score$score.ribosome
dtf$`Mitochondrial gene expression` = dtf$percent.mt

dtf.nFeaturenCount = dtf[
  source %in% c("freshPBMC", "frozenPBMC", "freshWB_RBClysis", "freshWB_RBCdepletion")][
  , lapply(.SD, function(x) { geom_mean(as.numeric(x)) }),
  by = c("source", "technique2", "sampleid.pid", "pid"),
  .SDcols = c("nCount_RNA", "nFeature_RNA", "Ribosomal protein (RPL/RPS) gene expression")]
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
    expr = t(as.matrix(sub[, features, with = FALSE])), class = sub$class,
    label = label_source(i))
})

meta.nFeaturenCount.BDvs10x = meta_analysis(metaobj, outcome_type = "binary")

### per source and platform
p.nFeature = ggplot(dtf.nFeaturenCount, aes(x = source.tech, y = nFeature_RNA)) +
  ggbeeswarm::geom_quasirandom(aes(color = pid), size = 0.6, width = 0.2,
                               show.legend = TRUE, alpha = 0.7,
                               varwidth = FALSE, shape = 16, stroke = 0.1) +
  geom_boxplot(width = boxwidth, size = 0.2, outlier.shape = NA,
               alpha = 1, fill = NA, color = "grey17") +
  labs(y = "nFeature_RNA", x = NULL, title = "") +
  coord_cartesian(clip = "off") +
  scale_y_log10() +
  scale_color_tableau("Tableau 20", labels = cap_first) +
  theme_expresso(legend_position = c(1.02,1)) +
  theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank(),
        legend.justification = c(0,1),
        plot.title = element_text(size = textsize, face = "plain", hjust = 0)) +
  guides(colour = guide_legend(override.aes = list(size = 1)))

p.nCount = ggplot(dtf.nFeaturenCount, aes(x = source.tech, y = nCount_RNA)) +
  ggbeeswarm::geom_quasirandom(aes(color = pid), size = 0.6, width = 0.2,
                               show.legend = TRUE, alpha = 0.7,
                               varwidth = FALSE, shape = 16, stroke = 0.1) +
  geom_boxplot(width = boxwidth, size = 0.2, outlier.shape = NA,
               alpha = 1, fill = NA, color = "grey17") +
  labs(y = "nCount_RNA", x = NULL, title = "") +
  coord_cartesian(clip = "off") +
  scale_y_log10() +
  scale_x_discrete(labels = label_source) +
  scale_color_tableau("Tableau 20", labels = cap_first) +
  theme_expresso(legend_position = "none", legend_text_size = textsize) +
  theme(axis.text.x = element_text(lineheight = 0.7, angle = -30, hjust = 0),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = textsize, face = "plain", hjust = 0)) +
  guides(colour = guide_legend(override.aes = list(size = 1.5)))

### per subset
dtf.nFeaturenCount2 = dtf[
  , lapply(.SD, function(x) geom_mean(as.numeric(x))), by = .(celltype, sampleid.pid, pid),
  .SDcols = c("nCount_RNA", "nFeature_RNA",
              "Ribosomal protein (RPL/RPS) gene expression", "Mitochondrial gene expression")]

medians = dtf.nFeaturenCount2[
  , lapply(.SD, median), by = "celltype",
  .SDcols = c("nCount_RNA", "nFeature_RNA",
              "Ribosomal protein (RPL/RPS) gene expression", "Mitochondrial gene expression")]

make_qr_plot = function(y_var) {
  ggplot(dtf.nFeaturenCount2,
         aes(x = factor(celltype, levels = medians[order(-get(y_var))]$celltype),
             y = .data[[y_var]])) +
    ggbeeswarm::geom_quasirandom(color = "grey77", size = 0.001, width = 0.2,
                                 alpha = 0.7, shape = 16) +
    geom_boxplot(width = boxwidth, linewidth = 0.2, outlier.shape = NA, fill = NA,
                 color = "grey17") +
    labs(y = y_var, x = NULL, title = NULL) +
    coord_cartesian(clip = "off") +
    theme_expresso(legend_position = "none") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.x = element_blank())
}

p.nFeature.subset = make_qr_plot("nFeature_RNA") + scale_y_log10()
p.nCount.subset   = make_qr_plot("nCount_RNA") + scale_y_log10()
p.ribosome.subset = make_qr_plot("Ribosomal protein (RPL/RPS) gene expression")
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
  facetvar = "pseudotime", pseudotime = rowMeans2(pseudotime_matrix, na.rm = TRUE))

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
  geom_boxplot(width = boxwidth, linewidth = 0.2, outlier.shape = NA, fill = NA, color = "grey17") +
  labs(y = "Pseudotime", x = NULL, title = NULL) +
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
  df.slingshot, x = "UMAP_1", y = "UMAP_2", color_by = "pseudotime",
  color_type = "continuous", colors = viridis_pal(option = "D")(10),
  quantile_lower = 0.001, quantile_upper = 0.999, shuffle = TRUE, seed = 42,
  point_size = 0.0001, point_alpha = 0.3, raster_dpi = 400, legend_ncol = 1, facet_nrow = 1,
  label = TRUE, label_df = df.label, label_size = pt2mm(textsize)) +
  geom_path(data = lineages_df, aes(group = Lineage), linewidth = 0.3) +
  theme_expresso(legend_position = c(0.9, 1), legend_justification = c(0, 1),
                 legend_title = TRUE, legend_key_width = 0.1,
                 legend_key_height = 0.3, legend_key_spacing_x = 0.01,
                 legend_key_spacing_y = 0.05, show_axis = FALSE, grid = "none") +
  scale_color_gradientn(name = "Pseudotime", colours = viridis_pal(option = "D")(10),
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
pointsizen = 0.001
p.umap.pid.PB = plot_scatter(
  df.PB, x = "UMAP_1", y = "UMAP_2", color_by = "pid", color_type = "discrete", colors = col_pid_PB,
  facet_nrow = 1, point_size = 0.0001, point_alpha = 0.2, raster_dpi = 400,
  shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL,
  legend_ncol = 1, legend_point_size = 2) +
  scale_color_manual(values = col_pid_PB, labels = cap_first, name = NULL,
                     na.value = "grey42") +
  guides(color = guide_legend(ncol = 1, title = "Study", label.hjust = 0,
                              override.aes = list(size = 2, alpha = 0.7))) +
  theme_expresso(legend_position = c(1.01, 1), legend_justification = c(0, 1),
                 legend_title = TRUE, legend_key_spacing_x = 0.01,
                 legend_key_spacing_y = 0.01, legend_key_height = 0.1,
                 legend_text_size = textsize, facet_label_face = "plain",
                 show_axis = FALSE, grid = "none") +
  theme(legend.title = element_text(margin = margin(b = 0.1, unit = "line")),
        legend.text  = element_text(margin = margin(-0.04, 0, -0.04, 0.1, "line")))

### group
pointsizen = 0.001
p.umap.group.PB = plot_scatter(
  df.PB, x = "UMAP_1", y = "UMAP_2",
  color_by = "group", color_type = "discrete", colors = col_group,
  facet_nrow = 1, point_size = 0.0001, point_alpha = 0.2, raster_dpi = 400,
  shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL,
  legend_ncol = 1, legend_point_size = 2) +
  scale_color_manual(values = col_group, labels = cap_first, name = NULL,
                     na.value = "grey42") +
  guides(color = guide_legend(ncol = 1, title = "Group", label.hjust = 0,
                              override.aes = list(size = 2, alpha = 0.7))) +
  theme_expresso(legend_position = c(1.01, 1), legend_justification = c(0, 1),
                 legend_title = TRUE, legend_key_spacing_x = 0.01,
                 legend_key_spacing_y = 0.01, legend_key_height = 0.1,
                 legend_text_size = textsize, facet_label_face = "plain",
                 show_axis = FALSE, grid = "none") +
  theme(legend.title = element_text(margin = margin(b = 0.1, unit = "line")),
        legend.text  = element_text(margin = margin(-0.04, 0, -0.04, 0.1, "line")))

### cell type
md_celltype_PB = function(x) {
  x = as.character(x)
  x[x == "Treg"]          = "T<sub>reg</sub>"
  x[x == "CD4 T Eff/Mem"] = "CD4<sup>+</sup> T<sub>EM</sub>"
  x[x == "CD8 T Eff/Mem"] = "CD8<sup>+</sup> T<sub>EM</sub>"
  x = sub("^CD([48]) T", "CD\\1<sup>+</sup> T", x)
  x = sub(" Naive", " naive", x, fixed = TRUE)
  x
}

md_legend_PB = function(x)
  paste0(md_celltype_PB(x),
         "<sup style='color:transparent'>+</sup><sub style='color:transparent'>M</sub>")

df.label = as.data.table(df.PB)[, lapply(.SD, median), by = celltype,
                                .SDcols = c("UMAP_1", "UMAP_2")]
colnames(df.label)[2:3] = c("x", "y")
df.label$labeltext = md_celltype_PB(df.label$celltype)
df.label = df.label[!is.na(celltype)]
df.label = as.data.table(df.label)
df.label[, y := ifelse(celltype == "CD4 T Naive",   y + 1,   y)]
df.label[, y := ifelse(celltype == "CD4 T Eff/Mem", y - 0.5, y)]
df.label[, y := ifelse(celltype == "Prolif T/NK",   y + 0.5, y)]
df.label[, y := ifelse(celltype == "CD16 Monocyte", y - 0.5, y)]
df.label[, x := ifelse(celltype == "CD16 Monocyte", x + 1,   x)]
df.label[, x := ifelse(celltype == "HSPC",          x - 0.5, x)]
df.label[, x := ifelse(celltype == "NK",            x - 1,   x)]

p.umap.celltype.PB = plot_scatter(
  df.PB, x = "UMAP_1", y = "UMAP_2",
  color_by = "celltype", color_type = "discrete", colors = col_celltype_PB,
  facet_nrow = 1, point_size = 0.001, point_alpha = 0.3, raster_dpi = 300,
  shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL,
  legend_ncol = 1, legend_point_size = 2, label = FALSE) +
  geom_richtext(data = df.label, aes(x = x, y = y, label = labeltext),
                size = pt2mm(textsize), color = textcolor, alpha = 0.8,
                fill = NA, label.color = NA, label.padding = unit(rep(0, 4), "pt"),
                inherit.aes = FALSE) +
  scale_color_manual(values = col_celltype_PB, name = NULL, na.value = "grey42",
                     labels = md_legend_PB) +
  guides(color = guide_legend(ncol = 1, title = "Cell type",
                              override.aes = list(size = 2, alpha = 0.7))) +
  theme_expresso(legend_position = c(1.01, 1), legend_justification = c(0, 1),
                 legend_title = TRUE, legend_key_spacing_x = 0.01,
                 legend_key_spacing_y = 0.01, legend_key_height = 0.1,
                 legend_text_size = textsize, facet_label_face = "plain",
                 show_axis = FALSE, grid = "none") +
  theme(legend.text  = element_markdown(hjust = 0,
                                        margin = margin(0.09, 0, -0.17, 0.1, "line")),
        legend.title = element_text(margin = margin(b = 0.1, unit = "line")))

### neutrophils vs everything else
df.PB[, celltype.neu.main := ifelse(celltype == "Mature Neutrophil", "Mature Neu", celltype)]
df.PB[, celltype.neu.main := ifelse(celltype == "Immat Neutrophil", "Immat Neu", celltype.neu.main)]
df.PB[, celltype.neu.main := ifelse(celltype == "Neu-Lym", "Neu-Lym", celltype.neu.main)]
df.PB[, celltype.neu.main := ifelse(celltype == "B", "B cell", celltype.neu.main)]
df.PB[, celltype.neu.main := ifelse(celltype == "Eos/Baso", "Eos/Baso", celltype.neu.main)]
df.PB[, celltype.neu.main := ifelse(grepl("Mono|cDC|pDC", celltype.neu.main),
                                    "Monocyte and DC", celltype.neu.main)]
df.PB[, celltype.neu.main := ifelse(grepl("CD4|CD8|Treg|NK", celltype.neu.main),
                                    "T/NK cell", celltype.neu.main)]
df.PB[, celltype.neu.main := ifelse(grepl("Erythrocyte|HSPC", celltype.neu.main),
                                    "other", celltype.neu.main)]

col_celltype.neu.main = c(
  `Immat Neu` = "#0A73DC", `Mature Neu` = "#2E42B8", "Eos/Baso" = "#EFC000",
  `Monocyte and DC` = "#8CC269", Platelet = "grey88",  `T/NK cell` = "#D2F0F4",
  PB = "#FFDC91FF",  `B cell` = "#F9FB0E",  other = "grey90")

df.label.PB = as.data.table(df.PB)[, lapply(.SD, median), by = celltype.neu.main,
                                   .SDcols = c("UMAP_1", "UMAP_2")]
colnames(df.label.PB)[2:3] = c("x", "y")
df.label.PB$labeltext = df.label.PB[, "celltype.neu.main", with = FALSE]
df.label.PB[, x := ifelse(celltype.neu.main == "B cell", x + 1, x)]
df.label.PB[, x := ifelse(celltype.neu.main == "Mature Neu", x + 1.5, x)]
df.label.PB[, y := ifelse(celltype.neu.main == "Mature Neu", y - 1.4, y)]
df.label.PB[, x := ifelse(celltype.neu.main == "T and NK",   x - 1, x)]
df.label.PB[, y := ifelse(celltype.neu.main == "T and NK",   y - 1, y)]
df.label.PB = df.label.PB[!is.na(celltype.neu.main)][
  !celltype.neu.main %in% "other"]
df.label.PB = as.data.table(df.label.PB)

p.umap.celltype.neu.main = plot_scatter(
  df.PB, x = "UMAP_1", y = "UMAP_2", color_by = "celltype.neu.main", color_type = "discrete",
  colors = col_celltype.neu.main,  facet_nrow = 1,
  point_size = 0.0001, point_alpha = 0.05, raster_dpi = 400,
  shuffle = TRUE, seed = 42, na_color = "grey42", title = NULL,
  legend_ncol = 1, legend_point_size = 2, label = TRUE, label_df = df.label.PB,
  label_size = pt2mm(textsize)) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 1.4, alpha = 1))) +
  theme_expresso(legend_position = c(1.02,0.98), legend_justification = c(0,1),
                 legend_text_size = textsize, legend_key_spacing_x = 0.01,
                 legend_key_spacing_y = 0.01, legend_key_height = 0.1,
                 facet_label_face = "plain",
                 show_axis = FALSE, grid = "none", panel_background = "white",
                 plot_background = "white", plot_margin = ggplot2::margin(0, 0, 0, 0, "line")) +
  theme(legend.text = element_text(margin = margin(0, 0, 0, 0.08, "line")))


p.umap.celltype.neu.main = suppressMessages(
  p.umap.celltype.neu.main +
    scale_x_continuous(breaks = scales::pretty_breaks(),
                       expand = expansion(mult = c(0.05, 0.15))) +
    scale_color_manual(values = col_celltype.neu.main, name = NULL, na.value = "grey42",
                       breaks = setdiff(names(col_celltype.neu.main), "other")))


# -- Immune cells in peripheral blood, gene expression -------------------------

mat = readRDS(paste0(dir.data, "peripheral_blood/blood.mat.selectedgenes.rds"))
setdiff(signaturegenes, rownames(mat))

mat = mat[, df.PB$cell.pid]

### gene dotplot
df.PB.genes = cbind(df.PB, as.data.frame(t(as.matrix(mat))))
df.PB.genes$celltype = factor(df.PB.genes$celltype, levels = rev(names(col_celltype_PB)))

p.markerdot.PB = plot_dotplot(
  df.PB.genes, features = signaturegenes, group_by = "celltype",
  max_scale = 2, dot_scale = 2, dot_stroke = 0.05, title = NULL,
  col_fontsize = textsize, row_fontsize = textsize,
  cluster_features = FALSE, cluster_groups = FALSE, show_axis = FALSE,
  feature_side = "bottom", colorbar_title = "Scaled expression") +
  aes(x = .group, y = feature) +
  scale_x_discrete(limits = rev, labels = md_celltype_PB) +
  scale_y_discrete(limits = rev) +
  theme(
    legend.position      = c(1.1,1),
    legend.justification = "top",
    legend.location      = "plot",
    axis.text.x = element_markdown(angle = 90,
                                   face = "plain", size = textsize,
                                   margin = margin(t = 0.1, unit = "line")),
    axis.text.y = element_text(face = "italic", size = textsize,
                               margin = margin(r = 0.1, unit = "line")))


p.markerdot.PB$layers[[1]]$aes_params$colour = NULL
p.markerdot.PB = p.markerdot.PB + aes(colour = expr) +
  scale_colour_gradientn(colours = expresso_colors("diverging_bwr2"),
                         limits = c(-2, 2), guide = "none")

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


p.neu.pct = ggplot(dtfgenes.sum.sum, aes(x = pct, y = gene, fill = celltype.v2,
                                         colour = celltype.v2)) +
  geom_col(linewidth = linesize / 2) +
  scale_fill_manual(values = c("Neutrophil" = "grey67", "other" = NA)) +
  scale_colour_manual(values = c("Neutrophil" = linecolor, "other" = NA), guide = "none") +
  scale_y_discrete(limits = rev) +
  theme_expresso(show_axis = FALSE, grid = "none") +
  theme(legend.position = "none",
        axis.line.x.bottom  = element_line(color = linecolor, linewidth = linesize),
        axis.ticks.x.bottom = element_line(color = linecolor, linewidth = linesize),
        axis.text.x         = element_text(size = textsize, color = textcolor),
        axis.title.x        = element_text(size = textsize, color = textcolor),
        axis.text.y         = element_text(face = "italic", size = textsize,
                                           color = textcolor,
                                           margin = margin(r = 0.1, unit = "line"))) +
  labs(y = NULL, x = "% (Neu reads)")

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
  dtf1 = data.table(UMAP_1 = df.PB.genes$UMAP_1, UMAP_2 = df.PB.genes$UMAP_2,
                    id = df.PB.genes[[gene]], facetvar = facetlabel.scores[[gene]])

  pscores.neu.PB[[gene]] = plot_scatter(
    dtf1, x = "UMAP_1", y = "UMAP_2", color_by = "id",
    facet_by = "facetvar", color_type = "continuous",
    colors = c("grey99", scales::dichromat_pal("DarkRedtoBlue.12")(12)[7:12]),
    quantile_lower = 0.001, quantile_upper = 0.999, shuffle = TRUE, seed = 42,
    point_size = 0.002, point_alpha = 0.2, raster_dpi = 400, legend_ncol = 1, facet_nrow = 1) +
    theme_expresso(legend_position = "none", show_axis = FALSE, grid = "none",
                   facet_label_face = "plain", panel_background = "white",
                   plot_background  = "white",
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
    theme_expresso(legend_position = c(1, 0.48), legend_text_size = textsize,
                   facet_label_face = "italic", show_axis = FALSE, grid = "none",
                   plot_margin = ggplot2::margin(0.05, 0.6, 0.05, 0.3, "line"))
}


# -- Figure 1a, study-design flow chart ----------------------------------------

library(grid)

lineheight = 1.45   # one text line, as a multiple of the font size
padx       = 0.10   # rule overhang past the widest string in a node, inches
padnode    = 0.05   # gap between a rule and the nearest line of text
gapnode    = 0.06   # gap between a rule and the arrow head or tail facing it
arrowpad   = 0.10   # arrow length beyond its label block, split over both ends
labelgap   = 0.05   # gap between the arrow and its label
rulelwd    = 0.8
arrowlwd   = 0.6

lh = textsize * lineheight / 72

# A node is a data frame of text lines with a per-line horizontal justification;
# an empty string is a blank line. row() recycles just over its text, so a node
# built from several row() calls can mix justifications -- which is what the
# atlas node needs, title flush with the left end of its rule and counts centred
# under it.
row  = function(text, just = "centre") data.frame(text = text, just = just)
node = function(...) do.call(rbind, list(...))

nodes = list(
  node(row(c("Public scRNA-seq", "datasets with mixed cell types"))),
  node(row(c("Neutrophils from", "multiple datasets"))),
  node(row(c("Single-Cell Atlas of Human Neutrophil", "(SCAHN)"), just = "center"),
       row(""),
       row(c("43 datasets", "40 studies", "9 countries", "7 platforms", "781 samples",
             "20+ tissues", "20+ diseases", "1M neutrophils"))))

# One arrow between consecutive nodes, so one fewer than there are nodes.
arrows = list(c("Collection", "Processing ", "Annotation"),
              c("Integration", "Clustering", "Annotation"))

# String widths need a device to measure against, and the device the panel is
# drawn on cannot be opened until the size the panel wants is known, so the
# measuring is done on a throwaway one first.
strwidth_in = function(x) {
  f = tempfile(fileext = ".pdf")
  pdf(f, width = 4, height = 4)
  pushViewport(viewport(gp = gpar(fontsize = textsize)))
  w = max(convertWidth(stringWidth(x), "in", valueOnly = TRUE))
  dev.off(); unlink(f)
  w
}

w.nodes  = sapply(nodes, function(n) strwidth_in(n$text) + 2 * padx)
h.nodes  = sapply(nodes, function(n) nrow(n) * lh + 2 * padnode)
h.arrows = sapply(arrows, function(a) length(a) * lh + 2 * arrowpad)

# Panel width is set by the widest node, unless an arrow label hanging off the
# centre line reaches further right than that node's right end does.
w.label = max(sapply(arrows, strwidth_in))
w.wf = max(max(w.nodes), 2 * (labelgap + w.label))
h.wf = sum(h.nodes) + sum(h.arrows) + 2 * length(arrows) * gapnode

flowchart = function(just = c("left", "top")) {
  free = c(convertWidth(unit(1, "npc"), "in", valueOnly = TRUE),
           convertHeight(unit(1, "npc"), "in", valueOnly = TRUE))
  if (free[1] < w.wf || free[2] < h.wf)
    warning(sprintf("flowchart needs %.2f x %.2f in, cell is %.2f x %.2f in",
                    w.wf, h.wf, free[1], free[2]))
  pushViewport(viewport(
    x = unit(switch(just[1], left = 0, centre = 0.5, right = 1), "npc"),
    y = unit(switch(just[2], bottom = 0, centre = 0.5, top = 1), "npc"), just = just,
    width = unit(w.wf, "in"), height = unit(h.wf, "in"),
    xscale = c(0, w.wf), yscale = c(0, h.wf)))
  cx = w.wf / 2
  ycur = h.wf
  for (i in seq_along(nodes)) {
    n = nodes[[i]]
    hw = w.nodes[i] / 2
    rule = function(y) grid.lines(unit(c(cx - hw, cx + hw), "native"), unit(c(y, y), "native"),
                                  gp = gpar(col = linecolor, lwd = rulelwd, lineend = "butt"))
    rule(ycur)
    y = ycur - padnode
    for (j in seq_len(nrow(n))) {
      left = n$just[j] == "left"
      grid.text(n$text[j], x = unit(if (left) cx - hw + padx else cx, "native"),
                y = unit(y - lh / 2, "native"),
                just = c(if (left) "left" else "centre", "centre"),
                gp = gpar(fontsize = textsize, col = textcolor))
      y = y - lh
    }
    ycur = y - padnode
    rule(ycur)
    if (i > length(arrows)) break
    ytail = ycur - gapnode
    yhead = ytail - h.arrows[i]
    grid.lines(unit(c(cx, cx), "native"), unit(c(ytail, yhead), "native"),
               gp = gpar(col = linecolor, fill = linecolor, lwd = arrowlwd),
               arrow = arrow(angle = 20, length = unit(0.06, "in"), type = "closed"))
    grid.text(paste(arrows[[i]], collapse = "\n"), x = unit(cx + labelgap, "native"),
              y = unit((ytail + yhead) / 2, "native"), just = c("left", "centre"),
              gp = gpar(fontsize = textsize, col = textcolor, lineheight = lineheight))
    ycur = yhead - gapnode
  }
  popViewport()
}


# -- Figure 1 ------------------------------------------------------------------

cairo_pdf(file = paste0(dir.fig, "Figure1.pdf"), width = 7.06, height = 4)
pushViewport(viewport(layout = grid.layout(nrow = 86, ncol = 146)))

pushViewport(viewport(layout.pos.row = 7:78, layout.pos.col = 2:36))
flowchart()
upViewport(1)

print(p.umap.tissue, vp = viewport(layout.pos.row = 8:35, layout.pos.col = 40:66))
print(p.umap.group, vp = viewport(layout.pos.row = 42:69, layout.pos.col = 40:66))

print(p.ggalluvial, vp = viewport(layout.pos.row = 3:58, layout.pos.col = 90:140))

print(p.neu.pct.source2, vp = viewport(layout.pos.row = 60:85, layout.pos.col = 92:118))
print(p.neu.pct.technique2, vp = viewport(layout.pos.row = 60:85, layout.pos.col = 119:146))

print(p.umap.axis + theme(plot.margin = margin(0, 0, 0, 0.2, "line")),
      vp = viewport(layout.pos.row = 28:38, layout.pos.col = 36:46))

print(p.umap.axis + theme(plot.margin = margin(0, 0, 0, 0.2, "line")),
      vp = viewport(layout.pos.row = 62:72, layout.pos.col = 36:46))

grid.text(x = unit(0.01, "npc"), y = unit(0.98, "npc"), label = "a",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.27, "npc"), y = unit(0.98, "npc"), label = "b",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.65, "npc"), y = unit(0.98, "npc"), label = "c",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))

grid.text(x = unit(0.65, "npc"), y = unit(0.31, "npc"), label = "d",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))

dev.off()


# -- Figure 2, legends and bracket helpers -------------------------------------

df.twocircles = data.table(
  x = rep(1, 2), y = 1:2, color = c("grey88", "grey17"), labeltext = c("FDR > 0.05", "FDR ≤ 0.05"))

p.twocircles = ggplot(df.twocircles, aes(x = x, y = y, label = labeltext, color = color)) +
  geom_point(shape = 21, fill = "grey88", size = 1.5) +
  geom_text(angle = 0, size = pt2mm(textsize), hjust = -0.12, color = textcolor) +
  theme(text = element_text(size = textsize, color = textcolor)) +
  scale_color_identity() +
  theme_void() +
  coord_cartesian(clip = "off")

p.twocircles.t = ggplot(df.twocircles, aes(x = y, y = x, label = labeltext, color = color)) +
  geom_point(shape = 21, fill = "grey88", size = 2, stroke = 0.3) +
  geom_text(angle = 90, size = pt2mm(textsize), hjust = -0.2, color = textcolor) +
  theme(text = element_text(size = textsize, color = textcolor)) +
  scale_color_identity() +
  scale_x_reverse() +
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
  theme(text = element_text(size = textsize, color = textcolor),
        plot.title = element_text(size = textsize, hjust = 0.5, color = textcolor),
        axis.title.y = element_text(size = textsize, color = textcolor))

circle_data.t = data.frame(
  x = c(2, 2, 2), y = c(9.6, 7.8, 5.6), ytext = c(12.2, 7.6, 2.9),
  radius = c(0.4, 0.65, 0.9), label = c("high", "", "low"))

pthreecircle.t = ggplot(data = circle_data.t) +
  geom_circle(aes(x0 = x, y0 = y, r = radius), linewidth = 0.1, size = 0.1,
              fill = "grey83", color = "grey83") +
  geom_text(aes(x = x, y = ytext, label = label), angle = 90, size = pt2mm(textsize),
            color = textcolor) +
  theme_void() +
  labs(y = "sd") +
  scale_y_continuous(limits = c(1.6, 13.3)) +
  scale_x_continuous(limits = c(1, 3)) +
  theme(text = element_text(size = textsize, color = textcolor),
        axis.title.y = element_text(size = textsize, color = textcolor, angle = 90))

df_leg = data.frame(
  label = factor(names(col_celltype), levels = names(col_celltype)), color = col_celltype)

p.subsetnames = ggplot(df_leg, aes(x = 1, y = label)) +
  geom_point(size = 0.8, aes(color = label)) +
  geom_text(aes(x = 1.1, label = label), hjust = 0, size = pt2mm(textsize), color = textcolor) +
  scale_color_manual(values = col_celltype) +
  scale_x_continuous(limits = c(0.95, 2)) +
  theme_void() +
  theme(
    legend.position = "none", plot.margin = margin(0, 0, 0, 0)) +
  scale_y_discrete(limits = rev(names(col_celltype)))

df_leg.matureornot = data.frame(
  label = factor(names(col_matureornot), levels = names(col_matureornot)),
  color = col_matureornot)

p.subsetnames.matureornot = ggplot(df_leg.matureornot, aes(x = 1, y = label)) +
  geom_point(size = 0.8, aes(color = label)) +
  geom_text(aes(x = 1.1, label = label), hjust = 0, size = pt2mm(textsize), color = textcolor) +
  scale_color_manual(values = col_matureornot) +
  scale_x_continuous(limits = c(0.95, 2)) +
  theme_void() +
  theme(
    legend.position = "none", plot.margin = margin(0, 0, 0, 0)) +
  scale_y_discrete(limits = rev(names(col_matureornot)))

col_tissuetype_tmp = c("Whole blood /\npurified neutrophil" = "#EE6C00", "Solid tissue" = "grey17")
df_leg.tissuetype = data.frame(
  label = factor(names(col_tissuetype_tmp), levels = names(col_tissuetype_tmp)),
  label_cap = cap_first(names(col_tissuetype_tmp)),
  color = col_tissuetype_tmp)

p.subsetnames.tissuetype = ggplot(df_leg.tissuetype, aes(x = 1, y = label)) +
  geom_point(size = 0.5, aes(color = label)) +
  geom_text(aes(x = 1.03, label = label_cap), hjust = 0, size = pt2mm(textsize),
            color = textcolor, lineheight = 0.8) +
  scale_color_manual(values = col_tissuetype_tmp) +
  scale_x_continuous(limits = c(0.95, 2)) +
  theme_void() +
  theme(
    legend.position = "none", plot.margin = margin(0, 0, 0, 0)) +
  scale_y_discrete(limits = rev(names(col_tissuetype_tmp)))

draw_bracket_label = function(bracket_rows, cols, label, gap = 2,
                              margin = c(0, 0, 0, 0),   # top, right, bottom, left in npc
                              lwd = 0.7, ticks = 0.5, h = 0.3, type = 4, curvature = 0.3,
                              bracket_col = "grey42", textsize = 8, textface = "plain",
                              textcol = textcolor) {
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
                               textcol = textcolor, lineheight = 0.8) {
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

draw_hline_label = function(line_rows, cols, label, gap = 2,
                            text_gap = 0.7,           # label to rule, in mm
                            lineheight = 1.2,
                            margin = c(0, 0, 0, 0),   # top, right, bottom, left in npc
                            lwd = 0.7, line_col = "grey42", textsize = 8,
                            textface = "plain", textcol = textcolor) {
  text_row = min(line_rows) - gap
  all_rows = text_row:max(line_rows)
  n_all    = length(all_rows)
  n_line   = length(line_rows)

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

  line_y = (n_line / 2) / n_all

  grid.lines(x = c(0, 1), y = c(line_y, line_y), gp = gpar(lwd = lwd, col = line_col))
  grid.text(unit(0.5, "npc"), unit(line_y, "npc") + unit(text_gap, "mm"), label = label,
            just = c("centre", "bottom"),
            gp = gpar(fontsize = textsize, fontface = textface, col = textcol,
                      lineheight = lineheight))

  upViewport(2)
}

draw_vline_label = function(rows, line_cols, label,
                            text_dir = c("horizontal", "vertical"),  # label orientation
                            text_gap = 0.7,           # label to rule, in mm
                            margin = c(0, 0, 0, 0),   # top, right, bottom, left in npc
                            lwd = 0.7, line_col = "grey42", textsize = 8, textface = "plain",
                            textcol = textcolor, lineheight = 0.8) {
  text_dir = match.arg(text_dir)
  mt = margin[1]
  mr = margin[2]
  mb = margin[3]
  ml = margin[4]

  pushViewport(viewport(layout.pos.row = rows, layout.pos.col = line_cols))
  pushViewport(viewport(x = unit(ml, "npc"), y = unit(mb, "npc"), width = unit(1 - ml - mr, "npc"),
                        height = unit(1 - mt - mb, "npc"), just = c("left", "bottom")))
  grid.lines(x = c(0.5, 0.5), y = c(0, 1), gp = gpar(lwd = lwd, col = line_col))
  grid.text(unit(0.5, "npc") - unit(text_gap, "mm"), unit(0.5, "npc"), label = label,
            rot = if (text_dir == "vertical") 90 else 0,
            just = if (text_dir == "vertical") c("centre", "bottom") else c("right", "centre"),
            gp = gpar(fontsize = textsize, fontface = textface,
                      col = textcol, lineheight = lineheight))
  upViewport(2)
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

cairo_pdf(file = paste0(dir.fig, "Figure2.pdf"), width = 7.6, height = 9.8)
pushViewport(viewport(layout = grid.layout(nrow = 208, ncol = 86)))

draw_vline_label(rows = 25:63, line_cols = 46, lwd = 0.4,
                 margin = c(0.03, 0.01, 0.02, 0.01), label = "Granule genes",
                 textsize = textsize, textcol = textcolor, text_dir = "vertical")

genes2display = c("STMN1", "MKI67", "AZU1", "MPO", "PRTN3", "CTSG", "ELANE", "DEFA4",
                  "LTF", "MMP8", "CAMP", "ARG1", "MMP9", "IL1R2",
                  "CYP4F3", "S100A8", "S100A4", "OLR1", "TXNIP", "EGR1", "FOS", "PTGS2",
                  "G0S2", "IFI6", "GBP5", "CD274",
                  "CXCL2", "CXCL8", "CXCR4", "VEGFA", "CCL4", "NFKBIA", "IL1B", "IL1RN",
                  "SLPI", "PI3", "HSPA1B", "CD74")

pushViewport(viewport(layout.pos.row = 12:202, layout.pos.col = 77:85))
grid.arrange(grobs = pgenes2[genes2display], ncol = 2, as.table = TRUE, newpage = FALSE)
upViewport(1)

print(p.umap.axis + theme(plot.margin = margin(0, 0, 0, 0.2, "line")),
      vp = viewport(layout.pos.row = 198:208, layout.pos.col = 74:79))

pushViewport(viewport(layout.pos.row = 154:168, layout.pos.col = 3:44))
grid.arrange(grobs = pscores.neu, nrow = 1, as.table = TRUE, newpage = FALSE)
upViewport(1)

pushViewport(viewport(layout.pos.row = 172:186, layout.pos.col = 3:44))
grid.arrange(grobs = pscores.neu.PB, nrow = 1, as.table = TRUE, newpage = FALSE)
upViewport(1)

print(p.umap.axis + theme(plot.margin = margin(0, 0, 0, 0.2, "line")),
      vp = viewport(layout.pos.row = 162:172, layout.pos.col = 1:6))
print(p.umap.axis + theme(plot.margin = margin(0, 0, 0, 0.2, "line")),
      vp = viewport(layout.pos.row = 180:190, layout.pos.col = 1:6))

print(p.umap.celltype.neu.main, vp = viewport(layout.pos.row = 187:205, layout.pos.col = 18:30))
print(p.umap.axis + theme(plot.margin = margin(0, 0, 0, 0.2, "line")),
      vp = viewport(layout.pos.row = 198:208, layout.pos.col = 16:21))

draw_hline_label(line_rows = 14:14, cols = 52:53, lwd = 0.4,
                 margin = c(0.5, 0.2, 0.1, 0.01), label = "Immature\nsubsets",
                 textsize = textsize, textcol = textcolor, text_gap = 0.3, lineheight = 0.7)
draw_hline_label(line_rows = 10:10, cols = 52:54, lwd = 0.4,
                 margin = c(0.5, 0.2, 0.1, 0.01), label = "Granule\nsubsets",
                 textsize = textsize, textcol = textcolor, text_gap = 0.3, lineheight = 0.7)
draw_hline_label(line_rows = 14:14, cols = 61:64, lwd = 0.4,
                 margin = c(0.1, 0.3, 0.1, 0.14), label = "IFN\nsubsets",
                 textsize = textsize, textcol = textcolor, text_gap = 0.3, lineheight = 0.7)
draw_hline_label(line_rows = 10:10, cols = 63:70, lwd = 0.4,
                 margin = c(0.5, 0.24, 0.1, 0.12), label = "Cytokine/\ninflammation\nsubsets",
                 textsize = textsize, textcol = textcolor, text_gap = 0.3, lineheight = 0.7)

print(pgenes$FCGR3B,  vp = viewport(layout.pos.row = 2:12,  layout.pos.col = 3:7))
print(pgenes$CXCR2,   vp = viewport(layout.pos.row = 13:23, layout.pos.col = 3:7))
print(pgenes$MME,     vp = viewport(layout.pos.row = 24:34, layout.pos.col = 3:7))
print(pgenes$CEACAM8, vp = viewport(layout.pos.row = 35:45, layout.pos.col = 3:7))
print(pgenes$`Granule\ngenes`, vp = viewport(layout.pos.row = 46:59, layout.pos.col = 3:8))
print(p.umap.axis + theme(plot.margin = margin(0, 0, 0, 0.2, "line")),
      vp = viewport(layout.pos.row = 53:63, layout.pos.col = 1:6))

print(p.umap.celltype, vp = viewport(layout.pos.row = 2:41, layout.pos.col = 13:36))
print(p.umap.celltype.matureornot, vp = viewport(layout.pos.row = 46:59, layout.pos.col = 14:21))
print(p.umap.axis, vp = viewport(layout.pos.row = 53:63, layout.pos.col = 12:17))
print(p.umap.axis, vp = viewport(layout.pos.row = 33:43, layout.pos.col = 12:17))

# ht_opt spacing, tightened for the two dot heatmaps only and restored below so the
# rest of the arm keeps the defaults set in 00_setup.R:
#   HEATMAP_LEGEND_PADDING  heatmap to right-side legend (default 2 mm)
#   TITLE_PADDING[1]        body side of the rotated row title (00_setup.R: 0.2 line)
#   DIMNAME_PADDING         both sides of the row names, and of the column names
#                           along the bottom (default 1 mm)
# Row names sit between the row title and the body, so TITLE_PADDING[1] and the
# two DIMNAME_PADDINGs are the whole adjustable part of that distance.
ht_opt_saved = list(HEATMAP_LEGEND_PADDING = ht_opt$HEATMAP_LEGEND_PADDING,
                    TITLE_PADDING          = ht_opt$TITLE_PADDING,
                    DIMNAME_PADDING        = ht_opt$DIMNAME_PADDING)
ht_opt$HEATMAP_LEGEND_PADDING = unit(0.1, "mm")
ht_opt$TITLE_PADDING          = unit(c(0.05, 0.2), "line")
ht_opt$DIMNAME_PADDING        = unit(0.1, "mm")

pushViewport(viewport(layout.pos.row = 64:117, layout.pos.col = 30:45))
draw(ht_prop, merge_legend = TRUE, heatmap_legend_side = "right",
     align_heatmap_legend = "heatmap_top", newpage = FALSE,
     background = "transparent", heatmap_legend_list = list(lgd.prop))
upViewport(1)

pushViewport(viewport(layout.pos.row = 64:120, layout.pos.col = 1:31))
draw(p.tissue.ht, merge_legend = TRUE, heatmap_legend_side = "right",
     align_heatmap_legend = "heatmap_top", newpage = FALSE,
     background = "transparent", heatmap_legend_list = list(lgd.tissue.ht))
upViewport(1)

print(p.twocircles + theme(plot.margin = margin(0, 0.1, 0, 0.1, "line")),
      vp = viewport(layout.pos.row = 83:84, layout.pos.col = 41:43))

print(pthreecircle.t + theme(plot.margin = margin(0.2, 0.4, 0.8, 0, "line")),
      vp = viewport(layout.pos.row = 82:96, layout.pos.col = 27:29))

ht_opt$HEATMAP_LEGEND_PADDING = ht_opt_saved$HEATMAP_LEGEND_PADDING
ht_opt$TITLE_PADDING          = ht_opt_saved$TITLE_PADDING
ht_opt$DIMNAME_PADDING        = ht_opt_saved$DIMNAME_PADDING

print(p.subsetnames + theme(plot.margin = margin(0.1, 0.01, 0.01, 0.03, "line")),
      vp = viewport(layout.pos.row = 4:50, layout.pos.col = 37:41))

print(p.subsetnames.matureornot + theme(plot.margin = margin(0.1, 0.01, 0.01, 0.03, "line")),
      vp = viewport(layout.pos.row = 47:50, layout.pos.col = 22:27))


print(p.subsetnames.tissuetype + theme(plot.margin = margin(0.1, 0.1, 0, 0.1, "line")),
      vp = viewport(layout.pos.row = 121:128, layout.pos.col = 36:45))

print(p.tissuesubset.pct,
     vp = viewport(layout.pos.row = 135:150, layout.pos.col = 33:45))

pushViewport(viewport(layout.pos.row = 121:151, layout.pos.col = 2:32))
grid.arrange(grobs = p.tissuesubsets.pct, nrow = 2, as.table = TRUE, newpage = FALSE)
upViewport(1)


print(p.markerdot.t + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "line")),
      vp = viewport(layout.pos.row = 12:207, layout.pos.col = 46:72))

grid.text(x = unit(0.01, "npc"), y = unit(0.99, "npc"), label = "a",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.15, "npc"),  y = unit(0.78, "npc"),  label = "b",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.15, "npc"),  y = unit(0.99, "npc"), label = "c",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))

grid.text(x = unit(0.55, "npc"), y = unit(0.98, "npc"), label = "d",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.88, "npc"), y = unit(0.98, "npc"), label = "e",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))

grid.text(x = unit(0.01, "npc"), y = unit(0.69, "npc"),  label = "f",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))

grid.text(x = unit(0.37, "npc"), y = unit(0.69, "npc"), label = "g",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.01, "npc"), y = unit(0.42, "npc"),  label = "h",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))

grid.text(x = unit(0.01, "npc"), y = unit(0.26, "npc"),  label = "i",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.01, "npc"), y = unit(0.17, "npc"),  label = "j",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))

dev.off()


# -- Extended figure 2, QC and trajectory --------------------------------------

cairo_pdf(file = paste0(dir.fig, "Fig.e2.pdf"), width = 6, height = 7)
pushViewport(viewport(layout = grid.layout(nrow = 100, ncol = 80)))

print(p.nFeature.subset + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.4, "line")),
      vp = viewport(layout.pos.row = 3:30, layout.pos.col = 2:39))
print(p.nCount.subset + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "line")),
      vp = viewport(layout.pos.row = 3:30, layout.pos.col = 40:79))
print(p.ribosome.subset + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.3, "line")),
      vp = viewport(layout.pos.row = 32:59, layout.pos.col = 4:39))
print(p.mt.subset + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.5, "line")),
      vp = viewport(layout.pos.row = 32:59, layout.pos.col = 41:79))

print(p.trajectory, vp = viewport(layout.pos.row = 62:100, layout.pos.col = 3:36))
print(p.pseudotime.subset + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.8, "line")),
      vp = viewport(layout.pos.row = 65:92, layout.pos.col = 40:79))
print(p.umap.axis, vp = viewport(layout.pos.row = 93:100, layout.pos.col = 2:8))

grid.text(x = unit(0.03, "npc"), y = unit(0.985, "npc"), label = "a",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.53, "npc"), y = unit(0.985, "npc"), label = "b",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.03, "npc"), y = unit(0.7, "npc"),  label = "c",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.53, "npc"), y = unit(0.7, "npc"),  label = "d",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.03, "npc"), y = unit(0.38, "npc"),   label = "e",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.53, "npc"), y = unit(0.38, "npc"),  label = "f",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))

dev.off()


# -- Supplementary figure 1, QC ------------------------------------------------

cairo_pdf(file = paste0(dir.fig, "Fig.s1.pdf"), width = 4, height = 6)
pushViewport(viewport(layout = grid.layout(nrow = 95, ncol = 70)))

print(p.nFeature + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "line")),
      vp = viewport(layout.pos.row = 1:18, layout.pos.col = 3:55))
print(p.nCount + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.23, "line")),
      vp = viewport(layout.pos.row = 19:43, layout.pos.col = 2:55))

print(plot_forest(as_forest_data(meta.nFeaturenCount.BDvs10x,
                                 gene = "nCount_RNA", rowname_style = "name"),
                  boxsize = 0.5, plottitle = "BD Rhapsody vs. 10X Chromium",
                  title_position = "center"),
      vp = viewport(layout.pos.row = 45:62, layout.pos.col = 3:30))

print(plot_forest(as_forest_data(meta.nFeaturenCount.BDvs10x, gene = "nFeature_RNA",
                                 rowname_style = "name"),
                  boxsize = 0.5, plottitle = "BD Rhapsody vs. 10X Chromium",
                  title_position = "center"),
      vp = viewport(layout.pos.row = 45:62, layout.pos.col = 36:63))

print(plot_forest(as_forest_data(meta.nFeaturenCount.WBvsPBMC,
                                 gene = "nCount_RNA", rowname_style = "name"),
                  boxsize = 0.5, plottitle = "WB vs. PBMC", title_position = "center") +
        scale_x_continuous(breaks = seq(-3, 3, 1)),
      vp = viewport(layout.pos.row = 64:78, layout.pos.col = 2:30))

print(plot_forest(as_forest_data(meta.nFeaturenCount.WBvsPBMC, gene = "nFeature_RNA",
                                 rowname_style = "name"),
                  boxsize = 0.5, plottitle = "WB vs. PBMC", title_position = "center"),
      vp = viewport(layout.pos.row = 64:78, layout.pos.col = 35:63))

print(plot_forest(as_forest_data(meta.nFeaturenCount.freshvsfrozenPBMC,
                                 gene = "nCount_RNA", rowname_style = "name"),
                  boxsize = 0.5, plottitle = "Fresh vs. frozen PBMC", title_position = "center"),
      vp = viewport(layout.pos.row = 80:94, layout.pos.col = 2:30))

print(plot_forest(as_forest_data(meta.nFeaturenCount.freshvsfrozenPBMC, gene = "nFeature_RNA",
                                 rowname_style = "name"),
                  boxsize = 0.5, plottitle = "Fresh vs. frozen PBMC", title_position = "center"),
      vp = viewport(layout.pos.row = 80:94, layout.pos.col = 35:63))

grid.text(x = unit(0.02, "npc"), y = unit(0.97, "npc"), label = "a",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.02, "npc"), y = unit(0.8, "npc"),  label = "b",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.02, "npc"),  y = unit(0.54, "npc"), label = "c",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))


dev.off()


# -- Extended figure 3 ---------------------------------------------------------

cairo_pdf(file = paste0(dir.fig, "Fig.e3.pdf"), width = 7.08, height = 9.2)
pushViewport(viewport(layout = grid.layout(nrow = 360, ncol = 186)))

print(p.umap.pid.PB, vp = viewport(layout.pos.row = 4:38, layout.pos.col = 5:25))
print(p.umap.group.PB, vp = viewport(layout.pos.row = 4:38, layout.pos.col = 44:64))

print(p.umap.celltype.PB, vp = viewport(layout.pos.row = 42:124, layout.pos.col = 5:64))

print(p.umap.axis, vp = viewport(layout.pos.row = 28:46, layout.pos.col = 1:12))
print(p.umap.axis, vp = viewport(layout.pos.row = 28:46, layout.pos.col = 40:51))
print(p.umap.axis, vp = viewport(layout.pos.row = 110:128, layout.pos.col = 1:12))
print(p.umap.axis, vp = viewport(layout.pos.row = 337:355, layout.pos.col = 1:12))

pushViewport(viewport(layout.pos.row = 150:348, layout.pos.col = 5:76))
grid.arrange(grobs = pgenes.blood, nrow = 7, as.table = TRUE, newpage = FALSE)
upViewport(1)

print(p.markerdot.PB, vp = viewport(layout.pos.row = 3:360, layout.pos.col = 83:150))
print(p.neu.pct + theme(plot.margin = margin(0.04, 0.3, 0.04, 0.2, "line")),
      vp = viewport(layout.pos.row = 3:345, layout.pos.col = 156:186))

grid.text(x = unit(0.012, "npc"), y = unit(0.99, "npc"), label = "a",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.23, "npc"), y = unit(0.99, "npc"), label = "b",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.012, "npc"),   y = unit(0.86, "npc"), label = "c",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.012, "npc"),   y = unit(0.58, "npc"), label = "d",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))

grid.text(x = unit(0.46, "npc"), y = unit(0.99, "npc"), label = "e",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.86, "npc"), y = unit(0.99, "npc"), label = "f",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))

dev.off()
