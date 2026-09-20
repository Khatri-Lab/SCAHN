# ==============================================================================
# Neutrophils in infectious disease and cancer
#
# Subset proportion and gene expression meta-analyses, correlation with the T/NK
# and mononuclear myeloid compartments, and outcome associations in whole blood
# (SUBSPACE + PUBLIC, SAVE-MORE) and in solid tumours (TCGA). The last section
# builds extended figure 7 from the per-cohort GLM and Cox matrices.

# Requires 00_setup.R (packages, figure settings, palettes) and 01_load_data.R
# (dtf, samplemeta, celltotalnumber, signaturegenes, gene_sets,
# dir.data / dir.results / dir.fig).
#
# Reads, on top of what 01_load_data.R loads:
#   results/annotation/SCAHN.subsetproportion.rds
#   results/annotation/SCAHN.subset.markers.csv
#   results/annotation/gene.exprpct.rds
#   results/meta_analysis/SCAHN.infect.exprmeta.selectedgenes.subsets.rds
#   results/meta_analysis/SCAHN.infect.exprmeta.selectedgenes.allneu.rds
#   results/meta_analysis/SCAHN.cancer.exprmeta.allgenes.allneu.rds
#   results/meta_analysis/dtfgenes.allneu.persample.geommean.rds
#   data/TNK.cellmeta.rds
#   data/myeloid.cellmeta.csv
#   data/myeloid.cellmeta_umap.csv
#   data/samplemeta.cancer.csv
#   results/clinical/SUBSPACE/SAVEMORE_processed.rds
#   results/clinical/SUBSPACE/cohort_glm_mv_matrices_PUBLIC_SUBSPACE.rds
#   results/clinical/SUBSPACE/cohort_glm_mv_matrices_PUBLIC_SUBSPACE.severity.rds
#   results/clinical/SUBSPACE/cohort_meta_re_mv_PUBLIC_SUBSPACE.rds
#   results/clinical/SUBSPACE/cohort_meta_re_mv_PUBLIC_SUBSPACE.severity.rds
#   results/clinical/TCGA/neut_rank.rds
#   results/clinical/TCGA/pan_cancer_cox_liu2018_mv_matrices.rds
#   results/clinical/TCGA/pan_cancer_meta_re_mv_neut.rds
#
# Writes to figures/post_acceptance/:
#   Figure4.pdf   Figure 4
#   Fig.e5.pdf    Extended figure 5 (infection forest plots, granule / IFN genes)
#   Fig.e6.pdf    Extended figure 6 (T/NK and myeloid UMAPs, all-compartment corr)
#   Fig.e7.pdf    Extended figure 7 (cohort and TCGA outcome heatmaps, forests)
# ==============================================================================

SCAHN_SCRIPTS = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/scripts/figures_post_acceptance"
source(file.path(SCAHN_SCRIPTS, "00_setup.R"))
source(file.path(SCAHN_SCRIPTS, "01_load_data.R"))

names(col_celltype)[names(col_celltype) == "lowDepth"] = "Low depth"

# Datasets are drawn as reference numbers, not study names. One table for the
# whole script, keyed on the lower-cased id so it catches a name whether or not
# cap_first() has been over it, and a function that falls back to cap_first() so
# an unlisted study still reads as it did.
ref_labels = c(
  combes2021       = "ref 14",
  schrepping2020   = "ref 5",
  zhang2023        = "ref 24",
  sinha2021        = "ref 15",
  wilk2021         = "ref 16",
  kaiser2024       = "ref 23",
  kwok2023         = "ref 6",
  `hu2022-oscc`    = "ref 27",
  hu2022           = "ref 27",
  `salcher2022-lc` = "LC, ref 26",
  `wu2024-gbc`     = "GBC, ref 7",
  `wu2024-hcc`     = "HCC, ref 7",
  `wu2024-other`   = "Other, ref 7",
  wu2024           = "ref 7",
  hu2023           = "ref 30",
  qian2020         = "ref 29",
  salcher2022      = "ref 26",
  wang2023         = "ref 28")

ref_label = function(x) {
  new = ref_labels[tolower(x)]
  ifelse(is.na(new), cap_first(x), new)
}


# -- Infection, neutrophil subset proportion -----------------------------------

col_infect = c(`COVID-19/flu,nonsevere` = "#7b72c5", `COVID-19,severe` = "#EFC000",
               `bacterial sepsis`       = "#F97D1C")

dtf.prop = readRDS(paste0(dir.results, "annotation/SCAHN.subsetproportion.rds"))
setnames(dtf.prop, "lowDepth", "Low depth")

### whole-blood samples from the infection studies
dtf.prop.infection =
  dtf.prop[pid %in% unique(dtf[grepl("COVID|sepsis|flu|conv", group)]$pid)][
    grepl("COVID|sepsis|flu|conv|healthy", group)][
    tissue %in% c("peripheral blood")][
    grepl("WB", source)]
table(dtf.prop.infection$group)

dtf.prop.infection[, condition := ifelse(
  group %in% c("COVID-19,nonsevere", "flu", "flu,pregnancy"), "COVID-19/flu,nonsevere",
  ifelse(group %in% c("healthy", "convalescent"), "healthy/convalescent", group))]

conditionorder = c("healthy/convalescent", "COVID-19/flu,nonsevere",
                   "COVID-19,severe", "bacterial sepsis")

dtf.prop.infection = dtf.prop.infection[neu.pct >= 0][n_cells >= 50]

dtf.prop.m = melt(dtf.prop.infection, id.vars = c("condition", "pid", "sampleid.pid"),
                  measure.vars = names(col_celltype))
dtf.prop.m = as.data.table(dtf.prop.m)

dtf.prop.m.aggmedian = dtf.prop.m[, lapply(.SD, function(x) {median(x)}),
                                  by = c("condition", "variable"), .SDcols = c("value")]
dtf.prop.m.aggsd = dtf.prop.m[, lapply(.SD, function(x) {sd(x)}), by = c("condition", "variable"),
                              .SDcols = c("value")]
dtf.prop.m.agg = dtf.prop.m[, lapply(.SD, function(x) {mean(x)}), by = c("condition", "variable"),
                            .SDcols = c("value")]
dtf.prop.m.agg[, lapply(.SD, function(x) {sum(x)}), by = c("condition"), .SDcols = c("value")]

dtf.prop.m.agg$variable = factor(dtf.prop.m.agg$variable, levels = names(col_celltype))
dtf.prop.m.aggmedian$variable = factor(dtf.prop.m.aggmedian$variable,
                                       levels = rev(names(col_celltype)))
dtf.prop.m.aggsd$variable = factor(dtf.prop.m.aggsd$variable, levels = rev(names(col_celltype)))

dtf.prop.m.agg$facetvar = "neutrophil subset in infectious diseases (whole blood)"

p.infectdistribution =
  ggplot(dtf.prop.m.agg, aes(x = condition, y = value, fill = variable)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  theme_expresso(axis_text_size = textsize, legend_text_size = textsize,
                 grid = "none", plot_background = "white",
                 legend_position = c(1.1, 1), legend_key_spacing_y = 0.01,
                 legend_key_height = 0.4,
                 plot_margin      = margin(5.5, 5.5, 5.5, 5.5, "pt")) +  # ggplot2 default
  theme(panel.grid.major = element_blank(), panel.border = element_blank(),
        axis.line.x.bottom  = element_line(color = linecolor, linewidth = linesize),
        axis.line.y.left    = element_line(color = linecolor, linewidth = linesize),
        axis.ticks.x.bottom = element_line(color = linecolor, linewidth = linesize),
        axis.ticks.y.left   = element_line(color = linecolor, linewidth = linesize),
        axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.y     = element_text(margin = margin(r = 0.5, l = 0, unit = "pt")),
        axis.text.y      = element_text(margin = margin(r = 0.5, l = 0, unit = "pt")),
        legend.justification = c(0, 1),
        legend.text      = element_text(color = textcolor,
                                        margin = margin(-0.06, 0, -0.06, 0.1, "line")),
        plot.title       = element_text(size = textsize, face = "plain", hjust = 0,
                                        margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1,
                                                        unit = "line"))) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = seq(0, 1, 0.2), position = "left") +
  scale_x_discrete(position = "bottom", limits = conditionorder, labels = cap_severity,
                   expand = expansion(add = 0.5)) +
  labs(x = NULL, y = "Subset proportion", title = NULL) +
  scale_fill_manual(values = col_celltype, name = "") +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(size = 1.2)))

### proportion meta-analysis, one contrast at a time
features = names(col_celltype)

### COVID-19,severe vs. healthy/convalescent
dtf.prop.infection[, class := ifelse(condition %in% "COVID-19,severe", 1, NA)]
dtf.prop.infection[, class := ifelse(condition %in% "healthy/convalescent", 0, class)]
pids = c("combes2021", "schrepping2020", "sinha2021", "wilk2021")

metaobj = lapply(pids, function(i) {
  sub = dtf.prop.infection[pid == i]
  meta_dataset(
    expr = t(as.matrix(sub[, features, with = FALSE])), class = sub$class, label = i)
})
meta.COVID19.severe = meta.prop.COVID19severe =
  meta_analysis(metaobj, outcome_type = "binary")

### bacterial sepsis vs. healthy/convalescent
dtf.prop.infection[, class := ifelse(condition %in% "bacterial sepsis", 1, NA)]
dtf.prop.infection[, class := ifelse(condition %in% "healthy/convalescent", 0, class)]
pids = c("combes2021", "kaiser2024", "kwok2023", "sinha2021")

metaobj = lapply(pids, function(i) {
  sub = dtf.prop.infection[pid == i]
  meta_dataset(
    expr = t(as.matrix(sub[, features, with = FALSE])), class = sub$class, label = i)
})
meta.sepsis = meta_analysis(metaobj, outcome_type = "binary")

### COVID-19/flu,nonsevere vs. healthy/convalescent
dtf.prop.infection[, class := ifelse(condition %in% "COVID-19/flu,nonsevere", 1, NA)]
dtf.prop.infection[, class := ifelse(condition %in% "healthy/convalescent", 0, class)]
pids = c("combes2021", "schrepping2020", "zhang2023")

metaobj = lapply(pids, function(i) {
  sub = dtf.prop.infection[pid == i]
  meta_dataset(
    expr = t(as.matrix(sub[, features, with = FALSE])), class = sub$class, label = i)
})
meta.nonsevere = meta_analysis(metaobj, outcome_type = "binary")

meta.COVID19.severe$summary$group = "COVID-19,severe"
meta.sepsis$summary$group         = "bacterial sepsis"
meta.nonsevere$summary$group      = "COVID-19/flu,nonsevere"

datadf = as.data.table(rbind(meta.COVID19.severe$summary, meta.sepsis$summary,
                             meta.nonsevere$summary))
datadf$celltype = factor(datadf$gene, levels = names(col_celltype))
datadf[, sig := ifelse(FDR <= 0.1, "FDR ≤ 0.1", "FDR > 0.1")]

nth = function(x, i) {
  x[(i - 1) %% length(x) + 1]
}
ci.value = -qnorm((1 - 0.95) / 2)

datadf$lower = datadf$pooled_ES - ci.value * datadf$SE
datadf$upper = datadf$pooled_ES + ci.value * datadf$SE

p.es.summary.infect =
  ggplot(data = datadf, aes(x = celltype, y = pooled_ES,
                            color = group, fill = group, alpha = sig)) +
  geom_hline(yintercept = 0, size = 0.2, alpha = 1, color = "grey67") +
  geom_tile(width = 0.7, height = 0.1, color = NA, position = position_dodge(0.7),
            show.legend = FALSE) +
  geom_errorbar(data = datadf, aes(x = celltype, ymin = lower, ymax = upper, group = group),
                position = position_dodge(0.7), width = .7, size = 0.2) +
  labs(x = " ", title = NULL,
       y = "Effect size\nvs. healthy/convalescent") +
  scale_color_manual(values = col_infect, breaks = names(col_infect), guide = "none") +
  scale_fill_manual(values = col_infect, breaks = names(col_infect), guide = "none") +
  scale_y_continuous(breaks = seq(-10, 2, 1)) +
  scale_alpha_manual(values = c("FDR > 0.1" = 0.2, "FDR ≤ 0.1" = 1)) +
  theme_expresso(axis_text_size = textsize, legend_text_size = textsize, grid = "none",
                 plot_background = "white", legend_position = c(1.02, 1),
                 legend_justification = c(0, 1), legend_key_width = 0.6,
                 legend_key_height = 0.1, plot_margin = margin(0, 0, 0, 0, "line")) +
  theme(plot.title       = element_text(margin = margin(0, 0, 0, 0, "line")),
        axis.title.x     = element_text(margin = margin(t = 0.05, unit = "line")),
        legend.box = "vertical", legend.spacing.y = unit(0.02, "line"),
        legend.text = element_text(margin = margin(l = 0.4, unit = "pt")),
        legend.key.size  = unit(1, "lines")) +
  coord_flip() +
  scale_x_discrete(limits = rev, position = "bottom") +
  guides(color = "none", shape = "none", size = "none",
         alpha = guide_legend(label.position = "right", nrow = 2, reverse = TRUE,
                              override.aes = list(alpha = 1, linewidth = 0.5,
                                                  color = c("grey17", "grey78"))))


# -- Infection, gene expression meta-analysis ----------------------------------

infect.exprmeta.sub = readRDS(paste0(
  dir.results, "meta_analysis/SCAHN.infect.exprmeta.selectedgenes.subsets.rds"))

### granule subsets
infect.exprmeta.sub$`granule subsets`$`COVID-19,severe`$summary$group =
  "COVID-19,severe"
infect.exprmeta.sub$`granule subsets`$`bacterial sepsis`$summary$group =
  "bacterial sepsis"
infect.exprmeta.sub$`granule subsets`$`COVID-19/flu,nonsevere`$summary$group =
  "COVID-19/flu,nonsevere"

datadf = as.data.table(rbind(
  infect.exprmeta.sub$`granule subsets`$`COVID-19,severe`$summary,
  infect.exprmeta.sub$`granule subsets`$`bacterial sepsis`$summary,
  infect.exprmeta.sub$`granule subsets`$`COVID-19/flu,nonsevere`$summary))
datadf$celltype = factor(datadf$gene, levels = names(col_celltype))

colorder = c("COVID-19/flu,nonsevere", "COVID-19,severe", "bacterial sepsis")

mat.es.granule = dcast(datadf, gene ~ group, value.var = "pooled_ES")
rowname = mat.es.granule$gene
mat.es.granule$gene = NULL
mat.es.granule = as.matrix(mat.es.granule)
rownames(mat.es.granule) = rowname

mat.pval.granule = dcast(datadf, gene ~ group, value.var = "p_value")
mat.pval.granule$gene = NULL
mat.pval.granule = as.matrix(mat.pval.granule)
rownames(mat.pval.granule) = rowname

genes2include = c("AZU1", "MPO", "PRTN3", "CTSG", "ELANE",
                  "DEFA4", "CEACAM6", "DEFA3", "BPI", "LYZ", "RETN", "ANXA1",
                  "CEACAM8", "OLFM4", "TCN1", "LTF", "LCN2", "MMP8", "CAMP",
                  "CRISP3", "HP", "PGLYRP1", "ANXA3", "FCN1", "ARG1", "PADI4",
                  "CD177", "MMP9", "NQO2", "S100A8", "S100A9", "S100A12",
                  "S100P", "TXN", "TSPO", "PLAC8", "GPI", "OLR1")

mat.es.granule   = mat.es.granule[genes2include, colorder]
mat.pval.granule = mat.pval.granule[genes2include, colorder]

max(mat.es.granule)
min(mat.es.granule)
ESmax = 2.5
mat.es.granule[mat.es.granule > ESmax]  = ESmax
mat.es.granule[mat.es.granule < -ESmax] = -ESmax

mat.fdr.granule = apply(mat.pval.granule, 2, function(x) {p.adjust(x, method = "fdr")})
mat.size.granule = 1 - mat.fdr.granule
mat.fdr.granule[mat.fdr.granule > 0.1] = 1
mat.size2.granule = 1 - mat.fdr.granule

col_fun = circlize::colorRamp2(
  seq(-ESmax, ESmax, length = 13), c(dichromat_pal("DarkRedtoBlue.12")(12)[1:6], "white",
    dichromat_pal("DarkRedtoBlue.12")(12)[7:12]))

colanno = HeatmapAnnotation(
  group = c("COVID-19/flu,nonsevere", "COVID-19,severe", "bacterial sepsis"),
  col = list(group = col_infect), annotation_legend_param = list(
    group = list(title = "Group", title_position = "topleft",
                 title_gp = gpar(fontsize = textsize, col = textcolor),
                 labels_gp = gpar(fontsize = textsize, col = textcolor),
                 legend_gp = gpar(fontsize = textsize), direction = "horizontal", nrow = 1,
                 legend_width = unit(.3, "line"), grid_width = unit(.3, "line"),
                 legend_height = unit(.3, "line"), grid_height = unit(.3, "line"))),
  show_legend = FALSE, show_annotation_name = FALSE,
  annotation_name_gp = gpar(fontsize = textsize, col = textcolor),
  annotation_name_side = "left", simple_anno_size = unit(0.3, "line"))

pht.granulesub =
  Heatmap(as.matrix(mat.es.granule), rect_gp = gpar(type = "none"), col = col_fun, na_col = "grey",
          cluster_columns = FALSE, cluster_rows = FALSE,
          show_row_dend = FALSE, show_column_dend = FALSE, bottom_annotation = colanno,
          column_dend_side = "bottom", clustering_method_rows = "ward.D2",
          clustering_method_columns = "ward.D2", clustering_distance_rows = "pearson",
          clustering_distance_columns = "pearson", show_column_names = FALSE,
          column_names_side = "top", column_names_rot = 45,
          column_names_gp = gpar(fontsize = textsize, col = textcolor),
          row_names_gp = gpar(fontsize = textsize, fontface = "italic", col = textcolor),
          row_names_max_width = unit(30, "line"),
          column_title = "      Granule subsets", column_title_side = "top",
          column_title_gp = gpar(fontsize = textsize, fontface = "plain", col = textcolor),
          show_heatmap_legend = FALSE, cell_fun = function(j, i, x, y, width, height, fill) {
            grid.circle(x = x, y = y, r = 4 * mat.size.granule[i, j] * min(unit.c(width, height)),
                        gp = gpar(fill = col_fun(mat.es.granule[i, j]), col = NA, alpha = 0.9))
            grid.circle(x = x, y = y, r = 4 * mat.size2.granule[i, j] * min(unit.c(width, height)),
                        gp = gpar(fill = NA, col = "grey17", alpha = 1))
          }, heatmap_legend_param = list(
            title = "ES", title_gp = gpar(fontsize = textsize, col = textcolor),
            labels_gp = gpar(fontsize = textsize, col = textcolor),
            legend_height = unit(3, "line"), legend_width = unit(.2, "line"),
            grid_height = unit(3, "line"), grid_width = unit(.1, "line"),
            at = seq(-2, 2, length = 3), by_row = FALSE,
            direction = "vertical", title_position = "topleft"))

### IFN subsets
infect.exprmeta.sub$`IFN subsets`$`COVID-19,severe`$summary$group =
  "COVID-19,severe"
infect.exprmeta.sub$`IFN subsets`$`bacterial sepsis`$summary$group =
  "bacterial sepsis"
infect.exprmeta.sub$`IFN subsets`$`COVID-19/flu,nonsevere`$summary$group =
  "COVID-19/flu,nonsevere"

datadf = as.data.table(rbind(
  infect.exprmeta.sub$`IFN subsets`$`COVID-19,severe`$summary,
  infect.exprmeta.sub$`IFN subsets`$`bacterial sepsis`$summary,
  infect.exprmeta.sub$`IFN subsets`$`COVID-19/flu,nonsevere`$summary))
datadf$celltype = factor(datadf$gene, levels = names(col_celltype))

colorder = c("COVID-19/flu,nonsevere", "COVID-19,severe", "bacterial sepsis")

mat.es.IFN = dcast(datadf, gene ~ group, value.var = "pooled_ES")
rowname = mat.es.IFN$gene
mat.es.IFN$gene = NULL
mat.es.IFN = as.matrix(mat.es.IFN)
rownames(mat.es.IFN) = rowname

mat.pval.IFN = dcast(datadf, gene ~ group, value.var = "p_value")
mat.pval.IFN$gene = NULL
mat.pval.IFN = as.matrix(mat.pval.IFN)
rownames(mat.pval.IFN) = rowname

genes2include = c("MT2A", "MT1X", "HES4", "LY6E", "ISG15", "IFI44L", "IFI44",
                  "OASL", "IFI6", "OAS3", "OAS2", "OAS1", "RSAD2", "HERC5",
                  "SAMD9", "SAMD9L", "MX1", "IFIT1", "IFIT2", "IFIT3", "IFIT5",
                  "DDX58", "XAF1", "IRF7", "EPSTI1", "APOL6", "GBP5", "GBP4",
                  "GBP2", "GBP1", "PARP14", "IRF1", "STAT1", "IFIH1", "TNFSF10",
                  "TNFSF13B", "IFITM1", "IFITM3") # 38

mat.es.IFN   = mat.es.IFN[genes2include, colorder]
mat.pval.IFN = mat.pval.IFN[genes2include, colorder]

max(mat.es.IFN)
min(mat.es.IFN)
ESmax = 2.5
mat.es.IFN[mat.es.IFN > ESmax]  = ESmax
mat.es.IFN[mat.es.IFN < -ESmax] = -ESmax

mat.fdr.IFN = apply(mat.pval.IFN, 2, function(x) {p.adjust(x, method = "fdr")})
mat.size.IFN = 1 - mat.fdr.IFN
mat.fdr.IFN[mat.fdr.IFN > 0.1] = 1
mat.size2.IFN = 1 - mat.fdr.IFN

col_fun = circlize::colorRamp2(
  seq(-ESmax, ESmax, length = 13), c(dichromat_pal("DarkRedtoBlue.12")(12)[1:6], "white",
    dichromat_pal("DarkRedtoBlue.12")(12)[7:12]))

colanno = HeatmapAnnotation(
  group = c("COVID-19/flu,nonsevere", "COVID-19,severe", "bacterial sepsis"),
  col = list(group = col_infect), annotation_legend_param = list(
    group = list(title = "Group", title_position = "topleft",
                 title_gp = gpar(fontsize = textsize, col = textcolor),
                 labels_gp = gpar(fontsize = textsize, col = textcolor),
                 legend_gp = gpar(fontsize = textsize), direction = "horizontal", nrow = 1,
                 legend_width = unit(.3, "line"), grid_width = unit(.3, "line"),
                 legend_height = unit(.3, "line"), grid_height = unit(.3, "line"))),
  show_legend = FALSE, show_annotation_name = FALSE,
  annotation_name_gp = gpar(fontsize = textsize, col = textcolor),
  annotation_name_side = "left", simple_anno_size = unit(0.3, "line"))

pht.IFNsub =
  Heatmap(as.matrix(mat.es.IFN), rect_gp = gpar(type = "none"), col = col_fun, na_col = "grey",
          cluster_columns = FALSE, cluster_rows = FALSE,
          show_row_dend = FALSE, show_column_dend = FALSE, bottom_annotation = colanno,
          column_dend_side = "bottom", clustering_method_rows = "ward.D2",
          clustering_method_columns = "ward.D2", clustering_distance_rows = "pearson",
          clustering_distance_columns = "pearson", show_column_names = FALSE,
          column_names_side = "top", column_names_rot = 45,
          column_names_gp = gpar(fontsize = textsize, col = textcolor),
          row_names_gp = gpar(fontsize = textsize, fontface = "italic", col = textcolor),
          row_names_max_width = unit(30, "line"),
          column_title = "IFN subsets", column_title_side = "top",
          column_title_gp = gpar(fontsize = textsize, fontface = "plain", col = textcolor),
          show_heatmap_legend = TRUE, cell_fun = function(j, i, x, y, width, height, fill) {
            grid.circle(x = x, y = y, r = 4 * mat.size.IFN[i, j] * min(unit.c(width, height)),
                        gp = gpar(fill = col_fun(mat.es.IFN[i, j]), col = NA, alpha = 0.9))
            grid.circle(x = x, y = y, r = 4 * mat.size2.IFN[i, j] * min(unit.c(width, height)),
                        gp = gpar(fill = NA, col = "grey17", alpha = 1))
          }, heatmap_legend_param = list(
            title = "Effect size vs.\nhealthy/convalescent",
            title_gp = gpar(fontsize = textsize, col = textcolor),
            labels_gp = gpar(fontsize = textsize, col = textcolor),
            legend_width = unit(.3, "line"), legend_height = unit(3, "line"),
            grid_width = unit(.3, "line"), grid_height = unit(3, "line"),
            at = seq(-2, 2, length = 3), by_row = FALSE,
            direction = "vertical", title_position = "leftcenter-rot"))


# -- Infection, gene expression meta-analysis, all neu -------------------------

### all neutrophils, five genes shown in figure 4
infect.exprmeta = readRDS(paste0(
  dir.results, "meta_analysis/SCAHN.infect.exprmeta.selectedgenes.allneu.rds"))

for (nm in names(infect.exprmeta))
  infect.exprmeta[[nm]]$per_dataset$dataset =
    ref_label(infect.exprmeta[[nm]]$per_dataset$dataset)

infect.exprmeta.summary = rbindlist(mapply(function(x, nm) {
  dt = as.data.table(x$summary)
  dt[, group := nm]
}, infect.exprmeta, names(infect.exprmeta), SIMPLIFY = FALSE))
infect.exprmeta.summary[gene == "CD274"]

colorder = c("COVID-19/flu,nonsevere", "COVID-19,severe", "bacterial sepsis")

mat.es.allneu.infect = dcast(infect.exprmeta.summary, gene ~ group, value.var = "pooled_ES")
rowname = mat.es.allneu.infect$gene
mat.es.allneu.infect$gene = NULL
mat.es.allneu.infect = as.matrix(mat.es.allneu.infect)
rownames(mat.es.allneu.infect) = rowname

mat.fdr.allneu = dcast(infect.exprmeta.summary, gene ~ group, value.var = "FDR")
mat.fdr.allneu$gene = NULL
mat.fdr.allneu = as.matrix(mat.fdr.allneu)
rownames(mat.fdr.allneu) = rowname

genes2include = c("ARG1", "OLR1", "S100A8", "S100A9", "CD274")

mat.es.allneu.infect = mat.es.allneu.infect[genes2include, colorder]
mat.fdr.allneu       = mat.fdr.allneu[genes2include, colorder]

max(mat.es.allneu.infect)
min(mat.es.allneu.infect)
ESmax = 2.5
mat.es.allneu.infect[mat.es.allneu.infect > ESmax]  = ESmax
mat.es.allneu.infect[mat.es.allneu.infect < -ESmax] = -ESmax

mat.size.allneu = 1.2 - mat.fdr.allneu
mat.fdr.allneu[mat.fdr.allneu > 0.1] = 1.2
mat.size2.allneu.infect = 1.2 - mat.fdr.allneu

mat.size.allneu         = sqrt(mat.size.allneu)
mat.size2.allneu.infect = sqrt(mat.size2.allneu.infect)

mat.es.allneu.infect    = t(mat.es.allneu.infect)
mat.size.allneu         = t(mat.size.allneu)
mat.size2.allneu.infect = t(mat.size2.allneu.infect)

col_fun = circlize::colorRamp2(
  seq(-ESmax, ESmax, length = 13), c(dichromat_pal("DarkRedtoBlue.12")(12)[1:6], "white",
    dichromat_pal("DarkRedtoBlue.12")(12)[7:12]))

pht.allneu.infect =
  Heatmap(as.matrix(mat.es.allneu.infect), rect_gp = gpar(type = "none", fill = "white"),
          col = col_fun, na_col = "grey", cluster_columns = FALSE, cluster_rows = FALSE,
          show_row_dend = FALSE, show_column_dend = FALSE,
          column_dend_side = "top", clustering_method_rows = "ward.D2",
          clustering_method_columns = "ward.D2", clustering_distance_rows = "pearson",
          clustering_distance_columns = "pearson", show_column_names = TRUE,
          column_names_side = "top", column_names_rot = 90,
          column_names_gp = gpar(fontsize = textsize, fontface = "italic", col = textcolor),
          row_names_side = "left",
          row_labels = cap_severity(rownames(mat.es.allneu.infect)),
          row_names_gp = gpar(fontsize = textsize, col = textcolor),
          row_names_max_width = unit(30, "line"), column_title = NULL, column_title_side = "top",
          column_title_gp = gpar(fontsize = textsize, fontface = "plain", col = textcolor),
          show_heatmap_legend = TRUE, cell_fun = function(j, i, x, y, width, height, fill) {
            grid.rect(x = x, y = y, width = width, height = height,
                      gp = gpar(col = NA, fill = "white", alpha = 0.2))
            grid.circle(x = x, y = y, r = .6 * mat.size.allneu[i, j] * min(unit.c(width, height)),
                        gp = gpar(fill = col_fun(mat.es.allneu.infect[i, j]),
                                  col = NA, alpha = 0.9))
            grid.circle(x = x, y = y, r = .6 * mat.size2.allneu.infect[i, j] *
                            min(unit.c(width, height)),
                        gp = gpar(fill = NA, col = "grey17", alpha = 1))
          }, heatmap_legend_param = list(
            title = "Effect size\nvs. healthy/\nconvalescent",
            title_gp = gpar(fontsize = textsize, col = textcolor, lineheight = 0.8),
            labels_gp = gpar(fontsize = textsize, col = textcolor),
            title_gap = unit(0.1, "line"), legend_width = unit(2, "line"),
            grid_height = unit(0.2, "line"),
            at = seq(-2, 2, length = 3), by_row = FALSE,
            direction = "horizontal", title_position = "topleft"))

pht.allneu.infect


# -- T/NK and mononuclear myeloid compartments ---------------------------------

# The T/NK names are set as in 02 (md_celltype_PB): a superscript plus on the CD4
# or CD8 lineage, and a subscript state where the manuscript abbreviates one. The
# labels here are drawn by ComplexHeatmap rather than by ggtext, so they are passed
# through gt_render(), which hands them to gridtext -- the same renderer, reached
# by ComplexHeatmap's own textGrob() wrapper, and it rotates richtext at arbitrary
# angles, which these column names need. Drawn labels only: the plain names are
# what the correlation matrices are built, subset and ordered on.
md_celltype = function(x) {
  x = as.character(x)
  x[x == "Treg"]            = "T<sub>reg</sub>"
  x[x == "CD4 T Eff/Mem"]   = "CD4<sup>+</sup> T<sub>EM</sub>"
  x[x == "CD8 T Eff/Mem"]   = "CD8<sup>+</sup> T<sub>EM</sub>"
  x[x == "CD4 T exhausted"] = "CD4<sup>+</sup> T<sub>EX</sub>"
  x[x == "CD8 T exhausted"] = "CD8<sup>+</sup> T<sub>EX</sub>"
  x = sub("^CD([48]) T", "CD\\1<sup>+</sup> T", x)
  x = sub(" Naive", " naive", x, fixed = TRUE)
  x
}

dtf.tnk = readRDS(paste0(dir.data, "TNK.cellmeta.rds"))
dtf.tnk$sampleid.pid = gsub("[ab]$", "", dtf.tnk$sampleid.pid)
samplemeta.cancer = fread(paste0(dir.data, "samplemeta.cancer.csv"))
length(table(dtf.tnk$sampleid.pid)) # 552

dtf.myeloid = fread(paste0(dir.data, "myeloid.cellmeta.csv"))
umap.myeloid = fread(paste0(dir.data, "myeloid.cellmeta_umap.csv"),
                   select = c("rowname", "UMAP_1", "UMAP_2", "sampleid.pid"))
dtf.myeloid = merge(dtf.myeloid, umap.myeloid, by = "rowname")
dtf.myeloid[, celltype := celltype.rescued2]
dtf.myeloid$sampleid.pid = gsub("[ab]$", "", dtf.myeloid$sampleid.pid)
dtf.myeloid$tissue = samplemeta[match(dtf.myeloid$sampleid.pid, samplemeta$sampleid.pid)]$tissue
dtf.myeloid$group = samplemeta[match(dtf.myeloid$sampleid.pid, samplemeta$sampleid.pid)]$group
length(table(dtf.myeloid$sampleid.pid)) # 591

### neu: remap IFN clusters, then compute proportions
dtf.neu = dtf
if (is.factor(dtf.neu$celltype)) {
  dtf.neu$celltype = unfactor(dtf.neu$celltype)
}
dtf.neu[, cluster := ifelse(grepl("IFN", celltype), "IFN", celltype)]
dtf.neu.s.pct = cell_proportions(dtf.neu, sample_col = "sampleid.pid", cluster_col = "cluster")
data.table::setnames(dtf.neu.s.pct, "n_cells", "neu.number")
dtf.neu.s.pct[, number.allcell := celltotalnumber[sampleid.pid]]
dtf.neu.s.pct[, neu.pct := neu.number / number.allcell]

### tnk: exclude doublets
dtf.tnk.s.pct = cell_proportions(dtf.tnk[!celltype %in% "doublet"], sample_col = "sampleid.pid",
                                 cluster_col = "celltype")
data.table::setnames(dtf.tnk.s.pct, "n_cells", "tnk.number")
dtf.tnk.s.pct[, tnk.pct := tnk.number / celltotalnumber[sampleid.pid]]

### myeloid: exclude other
dtf.myeloid.s.pct = cell_proportions(dtf.myeloid[!celltype %in% "other"],
                                   sample_col = "sampleid.pid", cluster_col = "celltype")
data.table::setnames(dtf.myeloid.s.pct, "n_cells", "myeloid.number")
dtf.myeloid.s.pct[, myeloid.pct := myeloid.number / celltotalnumber[sampleid.pid]]

dtfgenes.mean = readRDS(paste0(
  dir.results, "meta_analysis/dtfgenes.allneu.persample.geommean.rds"))
dtfgenes.mean$sampleid.pid = gsub("[ab]$", "", dtfgenes.mean$sampleid.pid)

length(intersect(dtf.tnk.s.pct$sampleid.pid, dtf.myeloid.s.pct$sampleid.pid)) # 552
length(intersect(dtf.neu.s.pct$sampleid.pid, dtf.myeloid.s.pct$sampleid.pid)) # 502

dtf.prop.neu.tnk = Reduce(
  function(x, y) merge(x, y, by = c("sampleid.pid"), all = FALSE), list(
    dtf.neu.s.pct, dtf.tnk.s.pct, dtf.myeloid.s.pct,
    dtfgenes.mean[, c("pid", "sampleid.pid", "CD274", "ARG1", "OLR1", "S100A8", "S100A9")]))

table(samplemeta[sampleid.pid %in% setdiff(dtf.neu.s.pct$sampleid.pid,
                                           dtf.myeloid.s.pct$sampleid.pid)]$pid)
table(samplemeta[sampleid.pid %in% setdiff(dtf.neu.s.pct$sampleid.pid,
                                           dtf.tnk.s.pct$sampleid.pid)]$pid)

setdiff(dtf.myeloid.s.pct$sampleid.pid, dtf.neu.s.pct$sampleid.pid)

dtf.prop.neu.tnk$pid = samplemeta[match(dtf.prop.neu.tnk$sampleid.pid, samplemeta$sampleid.pid)]$pid
dtf.prop.neu.tnk$group = samplemeta[match(dtf.prop.neu.tnk$sampleid.pid,
                                          samplemeta$sampleid.pid)]$group
dtf.prop.neu.tnk$tissue = samplemeta[match(dtf.prop.neu.tnk$sampleid.pid,
                                           samplemeta$sampleid.pid)]$tissue
dtf.prop.neu.tnk$tissue.detail = samplemeta[match(dtf.prop.neu.tnk$sampleid.pid,
                                                  samplemeta$sampleid.pid)]$tissue.detail

dtf.prop.neu.tnk = as.data.table(dtf.prop.neu.tnk) # 466  48
dim(dtf.prop.neu.tnk)

dtf.prop.neu.tnk$pid = gsub("a$|b$|c$|d$", "", dtf.prop.neu.tnk$pid)
dtf.prop.neu.tnk$pid = gsub("myin2024b", "yin2024", dtf.prop.neu.tnk$pid)
dtf.prop.neu.tnk$pid = gsub("reyfman2018", "reyfman2019", dtf.prop.neu.tnk$pid)
dtf.prop.neu.tnk$pid = gsub("xue2022", "xue2023", dtf.prop.neu.tnk$pid)

### infection samples with >= 50 cells in each of the three compartments
df.ntm.infect = dtf.prop.neu.tnk[sampleid.pid %in% dtf.prop.infection$sampleid.pid]
table(df.ntm.infect$pid)
table(dtf.prop.neu.tnk$pid)
table(dtf.prop.infection$pid)

df.ntm.infect = df.ntm.infect[
  tnk.number >= 50 & neu.number >= 50 & myeloid.number >= 50]

sample2exclude = ""
df.ntm.infect$combined = (df.ntm.infect$`CCL3/4.N` + df.ntm.infect$VEGFA.N) /
                         (df.ntm.infect$neu.number)

mat.infect = as.matrix(df.ntm.infect[, c(
  "CD4 T exhausted", "CD8 T exhausted", "Treg", "CD8 T cytotoxic", "CD4 T Eff/Mem", "T GZMK",
  "CD4 T Naive", "CD8 T Naive", "NK", "AZU1", "LTF", "IL1R2", "MMP9", "IFN",
  "CD274", "ARG1", "OLR1", "S100A8", "S100A9",
  "AP-1", "EGR1", "PTGS2", "TXNIP", "CCL3/4", "VEGFA", "IL1RN", "SLPI", "NF-κB", "IL1B")])

### single correlation
corr.infect = rcorr(mat.infect, type = "spearman")

rows = c("AZU1", "LTF", "MMP9", "CD274", "ARG1", "OLR1", "S100A8", "S100A9")
cols = c("CD4 T exhausted", "CD8 T exhausted", "Treg", "NK")

mat.corr.infect = corr.infect$r[rows, cols]
mat.corp.infect = corr.infect$P[rows, cols]

mat.corp.infect = matrix(p.adjust(as.vector(mat.corp.infect), method = "fdr"),
                         nrow = nrow(mat.corr.infect), ncol = ncol(mat.corr.infect))

mat.corr.infect
mat.corp.infect

NTM.Infect = function(celltype1, celltype2) {
  df.ntm.infect$celltype1 = df.ntm.infect[, celltype1, with = FALSE]
  df.ntm.infect$celltype2 = df.ntm.infect[, celltype2, with = FALSE]
  p =
    ggplot(df.ntm.infect[!sampleid.pid %in% sample2exclude], aes(x = celltype1, y = celltype2)) +
    labs(y = paste0(celltype2, " subset"), x = paste0(celltype1)) +
    geom_point(aes(color = group), size = 0.6, alpha = 0.8) +
    theme_expresso(plot_background = "white", legend_position = "none",
                   legend_justification = c(0, 0), legend_key_width = 0.1,
                   legend_key_height    = 0.1,
                   plot_margin          = margin(0.02, 0.05, 0.02, 0.05, "line")) +
    theme(legend.spacing.y = unit(0.04, "line"), legend.spacing.x = unit(0.04, "line"),
          legend.key       = element_blank(),
          axis.title       = element_text(margin = margin(0, 0, 0, 0)),
          axis.text        = element_text(margin = margin(0, 0, 0, 0)),
          axis.text.y      = element_text(hjust = 1)) +
    scale_color_manual(values = col_group) +
    scale_y_continuous(breaks = seq(0, 1, 0.1)) +
    stat_cor(color = textcolor, method = "spearman", size = 2.4, alpha = 0.5, label.sep = "\n") +
    geom_smooth(method = 'lm', alpha = 0.3, size = 0, span = 0.5, linetype = 0,
                color = "grey79", fill = "grey79") +
    stat_smooth(method = 'lm', geom = "line", alpha = 0.3,
                linewidth = 0.3, size = 1, span = 0.5, color = "#B71B1BFF")
  return(p)
}

### meta-correlation, selected neutrophil features vs. T/NK subsets
metacor.res = meta_corr(df.ntm.infect, vars = c("AZU1", "LTF", "MMP9", "IFN",
                                 "CD274", "ARG1", "OLR1", "S100A8", "S100A9",
                                 "CD4 T exhausted", "Treg", "NK"), cor_method = "pearson",
                        batch = "pid")

res.per_batch = as.data.table(metacor.res$per_batch)
res.per_batch = res.per_batch[, c("var1", "var2", "n", "r", "ci_lo", "ci_hi", "batch")]
colnames(res.per_batch) = c("var1", "var2", "samplesize", "estimate", "upper", "lower", "rowname")
pids = 1:length(unique(res.per_batch$rowname))
names(pids) = rev(unique(res.per_batch$rowname))
res.per_batch$index = pids[res.per_batch$rowname]
res.per_batch$pval  = NA
res.per_batch$label = NA

metacorr.infect.summary = as.data.table(metacor.res$summary)
metacorr.infect.summary = metacorr.infect.summary[
  , c("var1", "var2", "n_total", "r_random", "ci_random_lo", "ci_random_hi", "p_random")]
colnames(metacorr.infect.summary) = c("var1", "var2", "samplesize", "estimate",
                                      "upper", "lower", "pval")
metacorr.infect.summary$rowname = "summary"
metacorr.infect.summary$index   = 0
metacorr.infect.summary$label = paste0(
  "r:", round(metacorr.infect.summary$estimate, 2), "\n",
  "pval: ", formatC(metacorr.infect.summary$pval, format = "e", digits = 1))
res.per_batch = res.per_batch[, colnames(metacorr.infect.summary), with = FALSE]
metacorr.infect.summary = as.data.table(rbind(metacorr.infect.summary, res.per_batch))
metacorr.infect.summary = metacorr.infect.summary[order(-index)]

### meta-correlation heatmap, figure 4
roworders = c("AZU1", "LTF", "MMP9", "CD274", "OLR1", "ARG1", "S100A8", "S100A9")
colorders = c("CD4 T exhausted", "Treg", "NK")

corr.df = metacorr.infect.summary[rowname == "summary"][
  var1 %in% roworders & var2 %in% colorders]
corr.df$fdr = p.adjust(corr.df$pval, method = "fdr")
corr.mat = dcast(corr.df, var1 ~ var2, value.var = "estimate")
corr.mat.rowname = corr.mat$var1
corr.mat$var1 = NULL
corr.mat = as.matrix(corr.mat)
rownames(corr.mat) = corr.mat.rowname
fdr.mat = dcast(corr.df, var1 ~ var2, value.var = "fdr")
fdr.mat$var1 = NULL
fdr.mat = as.matrix(fdr.mat)
rownames(fdr.mat) = corr.mat.rowname

corr.mat = corr.mat[roworders, colorders]
fdr.mat  = fdr.mat[roworders, colorders]

rowname = c("AZU1 subset", "LTF subset", "MMP9 subset", "CD274", "OLR1", "ARG1", "S100A8", "S100A9")

make_labels = function(x, italicize = character()) {
  expr_txt = vapply(x, function(s) {
    if (s %in% italicize) sprintf("italic('%s')", s) else sprintf("'%s'", s)
  }, character(1))
  parse(text = expr_txt)
}

rowname = make_labels(rowname, italicize = c("CD274", "OLR1", "ARG1", "S100A8", "S100A9"))

rownames(corr.mat) = rowname
rownames(fdr.mat)  = rowname

size.mat = 1.2 - fdr.mat
fdr.mat[fdr.mat > 0.1] = 1.2
size2.mat = 1.2 - fdr.mat
label.mat = corr.mat
label.mat = round(label.mat, 2)
label.mat[fdr.mat > 0.1] = ""

col_fun_corr_infect = circlize::colorRamp2(
  seq(-0.5, 0.5, length = 13), c(dichromat_pal("DarkRedtoBlue.12")(12)[1:6], "white",
    dichromat_pal("DarkRedtoBlue.12")(12)[7:12]))

ht_opt$TITLE_PADDING = unit(c(0.2, 0.2), "line")

pht.metacorr.infect =
  Heatmap(corr.mat, rect_gp = gpar(type = "none", fill = "white"),
          col = col_fun_corr_infect, na_col = "grey", cluster_columns = FALSE, cluster_rows = FALSE,
          show_row_dend = FALSE, show_column_dend = FALSE, column_dend_side = "bottom",
          clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
          clustering_distance_rows = "pearson", clustering_distance_columns = "pearson",
          show_column_names = TRUE, column_names_side = "top", column_names_rot = 45,
          column_labels = gt_render(md_celltype(colnames(corr.mat))),
          column_names_gp = gpar(fontsize = textsize, col = textcolor),
          row_labels = rowname,
          row_names_gp = gpar(fontsize = textsize, col = textcolor),
          column_title = "                Proportion in total T/NK cells",
          column_title_side = "top",
          column_title_gp = gpar(fontsize = textsize, fontface = "plain", hjust = 1,
                                 col = textcolor),
          row_title = "Proportion or expression in neutrophils",
          row_title_gp = gpar(fontsize = textsize, fontface = "plain",
                              col = textcolor), show_heatmap_legend = FALSE,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.rect(x = x, y = y, width = width, height = height,
                      gp = gpar(col = NA, fill = "white", alpha = 0.2))
            grid.circle(x = x, y = y, r = 0.8 * size.mat[i, j] * min(unit.c(width, height)),
                        gp = gpar(fill = col_fun_corr_infect(corr.mat[i, j]),
                                  col = NA, alpha = 0.9))
            grid.circle(x = x, y = y, r = 0.8 * size2.mat[i, j] * min(unit.c(width, height)),
                        gp = gpar(fill = NA, col = "grey17", alpha = 1, lwd = 0.5))
            grid.text(x = x, y = y, label.mat[i, j],
            gp = gpar(fontsize = 5, col = textcolor))
          })

lgd.metacorr.infect = Legend(
  col_fun = col_fun_corr_infect, at = seq(-0.5, 0.5, length = 3),
  title = "Correlation", title_gp = gpar(fontsize = textsize, col = textcolor),
  title_position = "topcenter", title_gap = unit(0.3, "line"),
  labels_gp = gpar(fontsize = textsize, col = textcolor),
  legend_width = unit(2, "line"), legend_height = unit(.2, "line"),
  grid_width = unit(2, "line"), grid_height = unit(.2, "line"),
  by_row = FALSE, direction = "horizontal")

### meta-correlation across all T/NK and mononuclear myeloid subsets
subsets_neu_tnk_myeloid = c("AP-1", "AZU1", "CCL3/4", "CD74", "CXCL", "EGR1",
                          "G0S2", "HSP", "IFN", "IL1B", "IL1R2", "IL1RN", "LTF",
                          "MME", "MMP9", "NF-κB", "PTGS2", "S100A4", "SLPI", "TXNIP", "VEGFA",
                          "CD4 T Eff/Mem", "CD4 T IFN", "CD4 T Naive",
                          "CD4 T exhausted", "CD8 T Naive", "CD8 T cytotoxic",
                          "CD8 T exhausted", "NK", "NKT", "Prolif T/NK",
                          "T GZMK", "T/NK CCL3/4", "Treg", "other T", "γδ T",
                          "AP-1 MoMac", "Alveolar Mac", "C1Q MoMac",
                          "CCL2 MoMac", "CCL3/4 MoMac", "CD14 Mono",
                          "CD16 Mono", "CXCL10 MoMac", "DC3", "FOLR2 Mac",
                          "HSP MoMac", "IFN MoMac", "IL1B MoMac", "Int Mono",
                          "MT MoMac", "TREM2 Mac", "VEGFA MoMac", "cDC", "mregDC", "pDC")

setdiff(colnames(df.ntm.infect), subsets_neu_tnk_myeloid)

subset_pct_insamples = df.ntm.infect[, lapply(.SD, function(x) {mean(x > 0)}),
                                     .SDcols = subsets_neu_tnk_myeloid]

# if a subset is absent in more than 60% of samples, exclude it from analysis
subsets_exclude = names(subset_pct_insamples)[unlist(subset_pct_insamples[1]) < 0.4]

metacor.res.TNK_myeloid = meta_corr(df.ntm.infect, vars = setdiff(subsets_neu_tnk_myeloid,
                                                 subsets_exclude), cor_method = "pearson",
                                  batch = "pid")

res.TNK_myeloid.per_batch = as.data.table(metacor.res.TNK_myeloid$per_batch)
res.TNK_myeloid.per_batch = res.TNK_myeloid.per_batch[
  , c("var1", "var2", "n", "r", "ci_lo", "ci_hi", "batch")]
colnames(res.TNK_myeloid.per_batch) = c("var1", "var2", "samplesize", "estimate",
                                      "upper", "lower", "rowname")
pids = 1:length(unique(res.TNK_myeloid.per_batch$rowname))
names(pids) = rev(unique(res.TNK_myeloid.per_batch$rowname))
res.TNK_myeloid.per_batch$index = pids[res.TNK_myeloid.per_batch$rowname]
res.TNK_myeloid.per_batch$pval  = NA
res.TNK_myeloid.per_batch$label = NA

metacorr.TNK_myeloid.infect.summary = as.data.table(metacor.res.TNK_myeloid$summary)
metacorr.TNK_myeloid.infect.summary = metacorr.TNK_myeloid.infect.summary[
  , c("var1", "var2", "n_total", "r_random", "ci_random_lo", "ci_random_hi", "p_random")]
colnames(metacorr.TNK_myeloid.infect.summary) = c("var1", "var2", "samplesize",
                                                "estimate", "upper", "lower", "pval")
metacorr.TNK_myeloid.infect.summary$rowname = "summary"
metacorr.TNK_myeloid.infect.summary$index   = 0
metacorr.TNK_myeloid.infect.summary$label = paste0(
  "r:", round(metacorr.TNK_myeloid.infect.summary$estimate, 2), "\n",
  "pval: ", formatC(metacorr.TNK_myeloid.infect.summary$pval, format = "e", digits = 1))
res.TNK_myeloid.per_batch = res.TNK_myeloid.per_batch[
  , colnames(metacorr.TNK_myeloid.infect.summary), with = FALSE]
metacorr.TNK_myeloid.infect.summary = as.data.table(
  rbind(metacorr.TNK_myeloid.infect.summary, res.TNK_myeloid.per_batch))
metacorr.TNK_myeloid.infect.summary =
  metacorr.TNK_myeloid.infect.summary[order(-index)]

### all-compartment heatmap, extended figure 6
roworders = c("AZU1", "LTF", "MMP9", "S100A4", "IL1R2", "MME", "TXNIP", "EGR1",
              "AP-1", "PTGS2", "G0S2", "IFN", "CXCL", "VEGFA",
              "CCL3/4", "NF-κB", "IL1B", "IL1RN", "SLPI", "HSP", "CD74")

colorders = c("CD4 T Eff/Mem", "CD4 T exhausted", "CD4 T IFN", "CD4 T Naive",
              "CD8 T cytotoxic", "CD8 T Naive", "NK", "NKT",
              "other T", "Prolif T/NK", "T GZMK", "T/NK CCL3/4", "Treg", "γδ T",
              "CD14 Mono", "CD16 Mono", "Int Mono", "IL1B MoMac", "CCL2 MoMac",
              "CCL3/4 MoMac", "CXCL10 MoMac", "IFN MoMac", "MT MoMac",
              "C1Q MoMac", "HSP MoMac", "AP-1 MoMac", "cDC", "DC3", "pDC", "mregDC")

# column compartment annotation: T/NK vs myeloid (Mono, DC, MoMac)
colorders_myeloid = c("CD14 Mono", "CD16 Mono", "Int Mono", "IL1B MoMac",
                    "CCL2 MoMac", "CCL3/4 MoMac", "CXCL10 MoMac", "IFN MoMac",
                    "MT MoMac", "C1Q MoMac", "HSP MoMac", "AP-1 MoMac", "cDC",
                    "DC3", "pDC", "mregDC")

compartment.col = factor(ifelse(colorders %in% colorders_myeloid, "myeloid", "T/NK"),
                         levels = c("T/NK", "myeloid"))

col_compartment_TNK_myeloid = c("T/NK" = "#4F8F5B", "myeloid" = "#8C71B5")

colanno.TNK_myeloid = HeatmapAnnotation(
  compartment = compartment.col, col = list(compartment = col_compartment_TNK_myeloid),
  annotation_legend_param = list(
    compartment = list(title = NULL, title_position = "topleft",
                       title_gp = gpar(fontsize = textsize, col = textcolor),
                       labels_gp = gpar(fontsize = textsize, col = textcolor),
                       direction = "vertical", ncol = 1, legend_width = unit(.3, "line"),
                       grid_width = unit(.3, "line"), legend_height = unit(.3, "line"),
                       grid_height = unit(.3, "line"))), show_legend = TRUE,
  show_annotation_name = FALSE, annotation_name_gp = gpar(fontsize = textsize, col = textcolor),
  annotation_name_side = "left", simple_anno_size = unit(0.3, "line"))

corr.df = metacorr.TNK_myeloid.infect.summary[rowname == "summary"][
  var1 %in% roworders & var2 %in% colorders]
corr.df$fdr = p.adjust(corr.df$pval, method = "fdr")
corr.mat = dcast(corr.df, var1 ~ var2, value.var = "estimate")
corr.mat.rowname = corr.mat$var1
corr.mat$var1 = NULL
corr.mat = as.matrix(corr.mat)
rownames(corr.mat) = corr.mat.rowname
fdr.mat = dcast(corr.df, var1 ~ var2, value.var = "fdr")
fdr.mat$var1 = NULL
fdr.mat = as.matrix(fdr.mat)
rownames(fdr.mat) = corr.mat.rowname

corr.mat = corr.mat[roworders, colorders]
fdr.mat  = fdr.mat[roworders, colorders]

size.mat = 1.2 - fdr.mat
fdr.mat[fdr.mat > 0.1] = 1.2
size2.mat = 1.2 - fdr.mat
label.mat = corr.mat
label.mat = round(label.mat, 2)
label.mat[fdr.mat > 0.1] = ""

col_fun_corr_infect2 = circlize::colorRamp2(
  seq(-0.8, 0.8, length = 13), c(dichromat_pal("DarkRedtoBlue.12")(12)[1:6], "white",
    dichromat_pal("DarkRedtoBlue.12")(12)[7:12]))

ht_opt$TITLE_PADDING = unit(c(0.8, 0.2), "line")

pht.metacorr.infect.TNK_myeloid =
  Heatmap(corr.mat, rect_gp = gpar(type = "none", fill = "white"),
          col = col_fun_corr_infect2, na_col = "grey", cluster_columns = FALSE,
          cluster_rows = FALSE,
          show_row_dend = FALSE, show_column_dend = FALSE, top_annotation = colanno.TNK_myeloid,
          column_dend_side = "bottom",
          clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
          clustering_distance_rows = "pearson", clustering_distance_columns = "pearson",
          show_column_names = TRUE, column_names_side = "top", column_names_rot = 90,
          column_labels = gt_render(md_celltype(colnames(corr.mat))),
          column_names_gp = gpar(fontsize = textsize, lineheight = 0.7, col = textcolor),
          row_names_side = "left",
          row_names_gp = gpar(fontsize = textsize, col = textcolor),
          column_title = "Proportion in T/NK or mononuclear myeloid compartments in whole blood",
          column_title_side = "top",
          column_title_gp = gpar(fontsize = textsize, fontface = "plain", hjust = 1,
                                 col = textcolor),
          row_title = "Proportion in neutrophil compartment in whole blood",
          row_title_side = "left",
          row_title_gp = gpar(fontsize = textsize, fontface = "plain", col = textcolor),
          show_heatmap_legend = TRUE, cell_fun = function(j, i, x, y, width, height, fill) {
            grid.rect(x = x, y = y, width = width, height = height,
                      gp = gpar(col = NA, fill = "white", alpha = 0.2))
            grid.circle(x = x, y = y, r = 0.4 * size.mat[i, j] * min(unit.c(width, height)),
                        gp = gpar(fill = col_fun_corr_infect2(corr.mat[i, j]),
                                  col = NA, alpha = 0.9))
            grid.circle(x = x, y = y, r = 0.4 * size2.mat[i, j] * min(unit.c(width, height)),
                        gp = gpar(fill = NA, col = "grey17", alpha = 1, lwd = 0.5))
            grid.text(x = x, y = y, label.mat[i, j], gp = gpar(fontsize = 4,
                                col = ifelse(corr.mat[i, j] > 0.5, "white", textcolor)))
          }, heatmap_legend_param = list(
            title = "Correlation", title_gp = gpar(fontsize = textsize, col = textcolor),
            labels_gp = gpar(fontsize = textsize, col = textcolor),
            legend_width = unit(.2, "line"), legend_height = unit(3, "line"),
            grid_width = unit(.2, "line"), grid_height = unit(2, "line"),
            at = seq(-0.8, 0.8, length = 5), by_row = FALSE,
            direction = "vertical", title_position = "leftcenter-rot"))


# -- Infection, SUBSPACE + PUBLIC whole-blood mortality ------------------------

meta_mv = readRDS(paste0(dir.results, "clinical/SUBSPACE/cohort_meta_re_mv_PUBLIC_SUBSPACE.rds"))

# This summary spells the subset with an ASCII k where col_celltype, the atlas and
# the TCGA summary in the cancer section all use a Greek kappa. The panels below
# select rows with geneset %in% names(col_celltype), so the mismatch dropped the
# row silently -- and it is not a null one: OR 1.18 (1.00-1.39), FDR 0.083, so it
# draws at full opacity. Fixed on meta_mv itself because the score panel further
# down re-reads meta_mv$summary rather than the table built here.
meta_mv$summary$geneset = sub("NF-kB", "NF-κB", meta_mv$summary$geneset)

### subsets
dtf_mv = as.data.table(meta_mv$summary)
dtf_mv = dtf_mv[geneset %in% names(col_celltype)]
dtf_mv$FDR = p.adjust(dtf_mv$p_value, method = "BH")

# order gene sets by multivariate OR (descending)
gs_order = dtf_mv[order(dtf_mv$OR, decreasing = FALSE), geneset]
dtf_mv$geneset = factor(dtf_mv$geneset, levels = gs_order)

dtf_PUBLIC = dtf_mv

# FDR for line color; model for point fill
dtf_PUBLIC$sig = ifelse(!is.na(dtf_PUBLIC$FDR) & dtf_PUBLIC$FDR < 0.1, "FDR ≤ 0.1", "FDR > 0.1")
dtf_PUBLIC$sig = factor(dtf_PUBLIC$sig, levels = c("FDR ≤ 0.1", "FDR > 0.1"))

sig_colors  = c("FDR ≤ 0.1" = "grey17", "FDR > 0.1" = "grey78")
model_fills = c("multivariate" = "grey52", "univariate" = "white")

p.meta.public.subset = ggplot(dtf_PUBLIC, aes(x = log2(OR), y = geneset)) +
  geom_vline(xintercept = 0, size = 0.2, alpha = 1, color = "grey67") +
  geom_errorbarh(aes(xmin = log2(CI_lo), xmax = log2(CI_hi), color = sig), height = 0.6, size = 0.3,
                 position = position_dodgev(height = 0.5)) +
  geom_point(shape = 21, size = 0.8, color = "grey52", fill = "grey52",
             position = position_dodgev(height = 0.5)) +
  scale_color_manual(values = sig_colors, name = "FDR", guide = guide_legend(
                       override.aes = list(shape = NA, linetype = "solid", linewidth = 0.5),
                       label.position = "right", nrow = 2)) +
  scale_fill_manual(values = model_fills, name = "", guide = guide_legend(
                      reverse = TRUE, override.aes = list(shape = 21, color = "grey52",
                                          size = 2))) +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1)) +
  scale_y_discrete(labels = cap_first) +
  labs(
    title = NULL,
    x = "Log2 odds ratio\nD28-D30 mortality",
    y = " ") +
  theme_expresso(axis_title_size = textsize, legend_key_spacing_y = 0.05,
                 legend_key_height = 0.1, legend_key_width = 0.6, legend_text_size = textsize,
                 legend_position = c(1, 1), legend_justification = c(0, 1)) +
  theme(
    panel.grid.major.y = element_blank(), legend.spacing.y = unit(0.2, "line"),
    legend.text = element_text(margin = margin(l = 0.4, unit = "pt")),
    axis.title.x = element_text(margin = margin(t = 0.5, unit = "pt")),
    axis.title.y = element_text(margin = margin(r = 0)))

### signature scores
dtf_mv = as.data.table(meta_mv$summary)
dtf_mv$geneset = gsub(".neu", "", dtf_mv$geneset)
dtf_mv = dtf_mv[geneset %in% names(gene_sets)]
dtf_mv$FDR = p.adjust(dtf_mv$p_value, method = "BH")

# order gene sets by multivariate OR (descending)
gs_order = dtf_mv[order(dtf_mv$OR, decreasing = FALSE), geneset]
dtf_mv$geneset = factor(dtf_mv$geneset, levels = gs_order)

dtf_PUBLIC.score = dtf_mv

# FDR for line color; model for point fill
dtf_PUBLIC.score$sig = ifelse(
  !is.na(dtf_PUBLIC.score$FDR) & dtf_PUBLIC.score$FDR < 0.1, "FDR ≤ 0.1", "FDR > 0.1")
dtf_PUBLIC.score$sig = factor(dtf_PUBLIC.score$sig, levels = c("FDR ≤ 0.1", "FDR > 0.1"))

sig_colors  = c("FDR ≤ 0.1" = "grey17", "FDR > 0.1" = "grey78")
model_fills = c("multivariate" = "grey52", "univariate" = "white")

p.meta.public.score = ggplot(dtf_PUBLIC.score, aes(x = log2(OR), y = geneset)) +
  geom_vline(xintercept = 0, size = 0.2, alpha = 1, color = "grey67") +
  geom_errorbarh(aes(xmin = log2(CI_lo), xmax = log2(CI_hi), color = sig), height = 0.6, size = 0.3,
                 position = position_dodgev(height = 0.5)) +
  geom_point(shape = 21, size = 0.8, color = "grey52", fill = "grey52",
             position = position_dodgev(height = 0.5)) +
  scale_color_manual(values = sig_colors, name = "FDR", guide = guide_legend(
                       override.aes = list(shape = NA, linetype = "solid", linewidth = 0.5),
                       label.position = "right", nrow = 2)) +
  scale_fill_manual(values = model_fills, name = "", guide = guide_legend(
                      reverse = TRUE, override.aes = list(shape = 21, color = "grey52",
                                          size = 2))) +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1)) +
  scale_y_discrete(labels = cap_first) +
  labs(
    title = NULL,
    x = "Log2 odds ratio\nD28-D30 mortality",
    y = " \n \n") +
  theme_expresso(axis_title_size = textsize, legend_key_spacing_y = 0.05,
                 legend_key_height = 0.1, legend_key_width = 0.6, legend_text_size = textsize,
                 legend_position = c(1, 1), legend_justification = c(0, 1)) +
  theme(
    panel.grid.major.y = element_blank(), legend.spacing.y = unit(0.2, "line"),
    legend.text = element_text(margin = margin(l = 0.4, unit = "pt")),
    axis.title.x = element_text(margin = margin(t = 0.5, unit = "pt")),
    axis.title.y = element_text(margin = margin(r = 0)))


# -- Infection, SAVE-MORE anakinra x IL1B --------------------------------------

proc.sub = readRDS(paste0(dir.results, "clinical/SUBSPACE/SAVEMORE_processed.rds"))

hte_savemore = function(feature, data, subcohort = "savemore_1", tp = "baseline",
                        arm_keep = c("placebo", "anakinra"), ref_level = "placebo",
                        scale_score = FALSE, cohort_label = NULL) {

  # ── Subset the pre-computed score table once (shared across features) ──────
  d_all = as.data.table(data)
  if (!is.null(subcohort)) d_all = d_all[cohort %in% subcohort]
  if (!is.null(tp))        d_all = d_all[timepoint == tp]
  d_all = d_all[!duplicated(sampleid)]
  sc_names = unique(as.character(d_all$cohort))

  # ── Per-feature analysis (defined inline; closes over d_all above) ─────────
  run_one = function(f) {
    if (!f %in% names(d_all)) stop("feature '", f, "' not a column in `data`")
    score     = as.numeric(d_all[[f]])
    score_lab = f
    if (scale_score) score = as.numeric(scale(score))

    # 2. analysis frame
    d = data.table(sampleid = d_all$sampleid, score = score,
                   mort   = suppressWarnings(as.integer(d_all$d30_mort)),
                   arm    = trimws(as.character(d_all$treatment)),
                   age = suppressWarnings(as.numeric(d_all$age)), sex = factor(d_all$sex),
                   sofa = suppressWarnings(as.numeric(d_all$sofa)), cohort = d_all$cohort)
    d[arm %in% c("", "NA", "n/a", "N/A"), arm := NA_character_]
    d = d[arm %in% arm_keep]
    # Bar plot needs only mort/arm/score; glm additionally needs age/sex/sofa.
    d_bar = d[!is.na(mort) & !is.na(arm) & !is.na(score)]
    d     = d_bar[!is.na(age) & !is.na(sex) & !is.na(sofa)]
    d[,     arm := relevel(factor(arm), ref = ref_level)]
    d_bar[, arm := relevel(factor(arm), ref = ref_level)]
    arm_levels = levels(d$arm)
    if (length(arm_levels) != 2L) {
      tab = d[, .N, by = arm]
      stop("expected 2 arm levels, got ", length(arm_levels), ". Levels: ",
           paste(sprintf("'%s' (n=%d)", tab$arm, tab$N), collapse = ", "))
    }

    # 3. logistic regression — arm × score, adjusted for age + sex + sofa.
    # Pooled-cohort mode adds `cohort` as a fixed effect.
    pool = length(unique(d$cohort)) > 1L
    fml_full = if (pool) mort ~ arm * score + age + sex + sofa + cohort
               else      mort ~ arm * score + age + sex + sofa
    fml_red  = if (pool) mort ~ arm + score + age + sex + sofa + cohort
               else      mort ~ arm + score + age + sex + sofa
    fit  = glm(fml_full, data = d, family = binomial)
    fit0 = glm(fml_red,  data = d, family = binomial)
    cf    = summary(fit)$coefficients
    p_int = anova(fit0, fit, test = "LRT")[2, "Pr(>Chi)"]

    # 4. Panel E — median split + Fisher per stratum (uses bar-plot frame:
    # samples with missing covariates are still included here).
    d_bar[, score_bin := factor(ifelse(score > median(score), ">median", "≤median"),
                                levels = c("≤median", ">median"))]
    bar_df = d_bar[, .(mort_pct = mean(mort), n = .N, n_dead = sum(mort)),
                   keyby = .(score_bin, arm)]
    fish   = d_bar[, .(p = fisher.test(table(arm, mort))$p.value), by = score_bin]

    lab_gap = 0.18 * max(bar_df$mort_pct)
    bar_df[, lab_y := mort_pct + lab_gap * (diff(range(mort_pct)) < lab_gap) *
                                 (arm == arm_levels[1]), by = score_bin]

    rule_y = -0.14 * max(bar_df$mort_pct)

    pal = setNames(c("grey70", "firebrick3"), arm_levels)
    if (is.null(cohort_label))
      cohort_label = sprintf("%s (%s, n=%d)", if (length(sc_names) > 1L) "SAVE-MORE" else sc_names,
                              ifelse(is.null(tp), "all tp", tp), nrow(d_bar))

    pE = ggplot(bar_df, aes(score_bin, mort_pct, fill = arm)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.55) +
      geom_text(aes(y = lab_y, label = sprintf("%d/%d", n_dead, n)),
                position = position_dodge(width = 0.8), vjust = -0.3, size = pt2mm(textsize)) +
      geom_segment(data = fish, aes(x = as.numeric(score_bin) - 0.2,
                                    xend = as.numeric(score_bin) + 0.2,
                                    y = max(bar_df$mort_pct) * 1.3,
                                    yend = max(bar_df$mort_pct) * 1.3),
                   inherit.aes = FALSE, linewidth = 0.2) +
      geom_text(data = fish, aes(x = score_bin, y = max(bar_df$mort_pct) * 1.45,
                                 label = sprintf("%.2g", p)), inherit.aes = FALSE,
                                 size = pt2mm(textsize)) +

      geom_segment(data = fish, aes(x = as.numeric(score_bin) - 0.34,
                                    xend = as.numeric(score_bin) + 0.34,
                                    y = rule_y, yend = rule_y),
                   inherit.aes = FALSE, linewidth = 0.2) +
      scale_y_continuous(breaks = seq(0, 1, 0.05), labels = label_percent(suffix = ""),
                         expand = expansion(mult = c(0, 0.1))) +
      scale_fill_manual(values = pal, labels = cap_first) +
      coord_cartesian(ylim = c(0, NA), clip = "off") +
      scale_x_discrete(labels = c("≤median" = sprintf("*%s*<sup>lo</sup>", f),
                                  ">median"      = sprintf("*%s*<sup>hi</sup>", f))) +
      labs(x = NULL, y = "Mortality (%)", fill = NULL, title = NULL) +
      theme_expresso(text_size = textsize, axis_title_size = textsize,
      legend_text_size = textsize) +
      theme(legend.position = c(1.02, 1), legend.justification = c(0, 1),
            axis.ticks.x       = element_blank(),
            axis.text.x        = element_markdown(margin = margin(t = 6, b = 0)),
            axis.title.y       = element_text(margin = margin(r = 0)),
            panel.grid.major.x = element_blank()) +
      guides(fill = guide_legend(nrow = 2, override.aes = list(size = 2, alpha = 1)))

    # 5. Panel F — marginal effect (log-odds) by arm
    q    = quantile(d$score, c(0.025, 0.975), na.rm = TRUE)
    grid = CJ(score = seq(q[1], q[2], length.out = 100),
              arm   = factor(arm_levels, levels = arm_levels))
    grid[, `:=`(age = mean(d$age), sex = levels(d$sex)[1], sofa = mean(d$sofa))]
    if (pool) grid[, cohort := d$cohort[1]]
    pr = predict(fit, newdata = grid, type = "link", se.fit = TRUE)
    grid[, `:=`(fit = pr$fit, lo = pr$fit - 1.96 * pr$se.fit, hi = pr$fit + 1.96 * pr$se.fit)]

    pF = ggplot(grid, aes(score, fit, colour = arm, fill = arm, linetype = arm)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
      geom_line(linewidth = 0.8) +
      annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5,
               label = sprintf("p = %.2g", p_int), size = pt2mm(textsize)) +
      scale_colour_manual(values = pal) +
      scale_fill_manual(values = pal) +
      scale_y_continuous(breaks = seq(-10, 10, 2)) +
      scale_x_continuous(breaks = seq(0, 20, 2)) +
      scale_linetype_manual(values = setNames(c("dashed", "solid"), arm_levels)) +
      labs(x = score_lab, y = "Mortality (logit)", colour = NULL, fill = NULL,
           linetype = NULL, title = NULL) +
      theme_expresso(text_size = textsize, axis_title_size = textsize,
      legend_text_size = textsize) +
      theme(legend.position = "none", plot.title = element_text(face = "plain"),
            axis.title.x    = element_text(face = "italic",
                                           margin = margin(t = 0, b = 0)),
            axis.text.x     = element_text(margin = margin(t = 1, b = 0)),
            axis.ticks.length.x = unit(1, "pt"),
            axis.title.y    = element_text(margin = margin(r = 0)))

    list(fit = fit, p_interaction = p_int, summary = cf, n = nrow(d),
         n_events = sum(d$mort == 1L), bar_df = bar_df, fisher = fish,
         pred_grid = grid, plot_bar = pE, plot_marg = pF,
         feature = f, subcohort = sc_names, tp = tp, arm_levels = arm_levels)
  }

  if (length(feature) > 1L) {
    out = lapply(feature, run_one); names(out) = feature; out
  } else {
    run_one(feature)
  }
}

km_savemore = function(feature, data, time_col = "survdays", event_col = "d30_mort",
                       subcohort = NULL, tp = "baseline",
                       arm_keep = c("placebo", "anakinra"), ref_level = "placebo",
                       cohort_label = NULL, y_min = NULL) {
  if (!is.null(y_min) && (length(y_min) != 1L || !is.numeric(y_min) || y_min < 0 || y_min >= 1))
    stop("`y_min` must be a single numeric in [0, 1) or NULL")
  if (!requireNamespace("survival", quietly = TRUE))
    stop("install the 'survival' package")
  sv = asNamespace("survival")

  d_all = as.data.table(data)
  if (!is.null(subcohort)) d_all = d_all[cohort %in% subcohort]
  if (!is.null(tp))        d_all = d_all[timepoint == tp]
  d_all = d_all[!duplicated(sampleid)]
  sc_names = unique(as.character(d_all$cohort))

  # Fit treatment vs placebo within one score stratum; return survfit, p, ggplot.
  fit_one_stratum = function(d_sub, arm_levels, bin_label, f, pal, header, y_min = NULL) {
    sf = sv$survfit(sv$Surv(time, event) ~ arm, data = d_sub)
    p_logrank = sv$survdiff(sv$Surv(time, event) ~ arm, data = d_sub)$pvalue

    # Cox HR (treatment vs ref). Wrapped in tryCatch — singular fits (e.g. no
    # events in one arm) return NA so the plot still renders.
    hr_est = NA_real_; hr_lo = NA_real_; hr_hi = NA_real_
    if (length(unique(d_sub$arm)) == 2L && sum(d_sub$event) > 0L) {
      cx = tryCatch(sv$coxph(sv$Surv(time, event) ~ arm, data = d_sub),
                    error = function(e) NULL, warning = function(w) NULL)
      if (!is.null(cx) && length(coef(cx)) >= 1L) {
        s      = summary(cx)$conf.int
        hr_est = s[1, "exp(coef)"]
        hr_lo  = s[1, "lower .95"]
        hr_hi  = s[1, "upper .95"]
      }
    }

    ss = summary(sf, censored = TRUE)
    arm_clean = sub("^arm=", "", as.character(ss$strata))
    plot_df = data.table(time = ss$time, surv = ss$surv, lo = ss$lower, hi = ss$upper,
                         n_risk = ss$n.risk, n_event = ss$n.event, n_censor = ss$n.censor,
                         arm = factor(arm_clean, levels = arm_levels))
    init_df = data.table(time = 0, surv = 1, lo = 1, hi = 1,
                         n_risk = NA_integer_, n_event = NA_integer_, n_censor = NA_integer_,
                         arm = factor(arm_levels, levels = arm_levels))
    plot_df = rbind(init_df, plot_df, use.names = TRUE)
    setorder(plot_df, arm, time)
    censor_df = plot_df[!is.na(n_censor) & n_censor > 0]

    n_by_arm = d_sub[, .N, by = arm][match(arm_levels, arm), N]
    e_by_arm = d_sub[, sum(event), by = arm][match(arm_levels, arm), V1]
    n_by_arm[is.na(n_by_arm)] = 0L; e_by_arm[is.na(e_by_arm)] = 0L

    if (is.null(y_min)) {
      y_min = floor(min(plot_df$surv, na.rm = TRUE) / 0.05) * 0.05 - 0.02
      y_min = max(0, min(0.9, y_min))
    }

    p = ggplot(plot_df, aes(time, surv, colour = arm)) +
      geom_step(linewidth = 0.7) +
      geom_point(data = censor_df, shape = 3, size = 1.2, show.legend = FALSE) +
      scale_y_continuous(limits = c(y_min, 1), labels = label_percent(suffix = ""),
                         expand = expansion(mult = c(0.02, 0.04))) +
      scale_x_continuous(breaks = c(0, 10, 20), labels = c("D0", "D10", "D20"),
                         expand = expansion(mult = c(0.01, 0.04))) +
      scale_colour_manual(values = pal, labels = sprintf("%s (n=%d, d=%d)",
                                           arm_levels, n_by_arm, e_by_arm)) +
      annotate("text", x = -Inf, y = -Inf, hjust = -0.03, vjust = -0.3, label = if (is.na(hr_est))
                         sprintf("log-rank p = %.2g", p_logrank) else
                         sprintf("HR = %.2f\np = %.1g", hr_est, p_logrank),
                         size = pt2mm(textsize)) +
      labs(x = NULL, y = "Survival (%)", colour = NULL,
           title = sprintf("%s · %s · %s", header, f, bin_label)) +
      theme_expresso(text_size = textsize, axis_title_size = textsize,
      legend_text_size = textsize) +
      theme(legend.position = "none", axis.title.y = element_text(margin = margin(r = 0)),
            axis.title.x    = element_text(margin = margin(t = 0, b = 0)),
            axis.text.x     = element_text(margin = margin(t = 1, b = 0)),
            axis.ticks.length.x = unit(1, "pt"),
            plot.title      = element_markdown(face = "plain", hjust = 0.5,
                                               margin = margin(l = 0, b = 0)))

    list(survfit = sf, p_logrank = p_logrank, hr = hr_est, hr_lo = hr_lo, hr_hi = hr_hi,
         n = nrow(d_sub), n_events = sum(d_sub$event == 1L), plot_df = plot_df, plot = p)
  }

  run_one = function(f) {
    if (!f %in% names(d_all)) stop("feature '", f, "' not a column in `data`")

    d = data.table(sampleid = d_all$sampleid, score = as.numeric(d_all[[f]]),
                   time     = suppressWarnings(as.numeric(d_all[[time_col]])),
                   event    = suppressWarnings(as.integer(d_all[[event_col]])),
                   arm = trimws(as.character(d_all$treatment)), cohort = d_all$cohort)
    d[arm %in% c("", "NA", "n/a", "N/A"), arm := NA_character_]
    d = d[arm %in% arm_keep]
    d = d[!is.na(time) & !is.na(event) & !is.na(arm) & !is.na(score) & time >= 0]
    d[, arm := relevel(factor(arm), ref = ref_level)]
    arm_levels = levels(d$arm)
    if (length(arm_levels) != 2L)
      stop("expected 2 arm levels, got ", length(arm_levels))

    d[, score_bin := factor(ifelse(score > median(score), "> median", "≤ median"),
                            levels = c("≤ median", "> median"))]

    pal = setNames(c("grey60", "firebrick3"), arm_levels)
    header = if (is.null(cohort_label))
               sprintf("%s (%s)", if (length(sc_names) > 1L) "SAVE-MORE" else sc_names,
                       ifelse(is.null(tp), "all tp", tp))
             else cohort_label

    # y-axis floor. Caller-supplied y_min wins; otherwise data-driven, shared
    # across the two strata so they're visually comparable.
    y_min_shared = if (!is.null(y_min)) y_min else {
      surv_min = min(vapply(c("≤ median", "> median"), function(b) {
        sf_b = sv$survfit(sv$Surv(time, event) ~ arm, data = d[score_bin == b])
        min(summary(sf_b, censored = TRUE)$surv, 1, na.rm = TRUE)
      }, numeric(1)))
      max(0, min(0.9, floor(surv_min / 0.05) * 0.05 - 0.02))
    }

    low = fit_one_stratum(d[score_bin == "≤ median"], arm_levels, "score ≤ median", f, pal, header,
                           y_min = y_min_shared)
    high = fit_one_stratum(d[score_bin == "> median"], arm_levels, "score > median", f, pal, header,
                           y_min = y_min_shared)

    list(low = low, high = high, n = nrow(d), n_events = sum(d$event == 1L),
         feature = f, subcohort = sc_names, tp = tp, arm_levels = arm_levels)
  }

  if (length(feature) > 1L) {
    out = lapply(feature, run_one); names(out) = feature; out
  } else {
    run_one(feature)
  }
}

genesets = list(
  IL1B = "IL1B", IL1A = "IL1A")

sm_gs_score_tbl = rbindlist(lapply(c("savemore_1", "savemore_2"), function(nm) {
  el   = proc.sub[[nm]]
  expr = log2(as.matrix(el$tpm) + 1) + 1
  ph   = as.data.table(el$pheno)
  setnames(ph, "accession", "sampleid")
  ph[, cohort := nm]
  gs_scores = as.data.table(sapply(genesets, function(gs) {
    g_in = intersect(gs, rownames(expr))
    if (!length(g_in)) return(rep(NA_real_, ncol(expr)))
    apply(expr[g_in, , drop = FALSE], 2, geom_mean)
  }))
  gs_scores[, sampleid := colnames(expr)]
  merge(ph, gs_scores, by = "sampleid")
}), use.names = TRUE, fill = TRUE)

res_sm = hte_savemore(feature = c("IL1B", "IL1A"), data = sm_gs_score_tbl,
                      subcohort = NULL, tp = "baseline")

km_il1b = km_savemore("IL1B", data = sm_gs_score_tbl, subcohort = NULL,
                      tp = "baseline", y_min = 0.85)


# -- Cancer, neutrophil subset proportion --------------------------------------

dtf.prop = readRDS(paste0(dir.results, "annotation/SCAHN.subsetproportion.rds"))
setnames(dtf.prop, "lowDepth", "Low depth")
dtf.prop.cancertissue = dtf.prop[!grepl("peripheral blood", tissue)][
  !grepl("han2020", pid)]
dtf.prop.cancertissue = dtf.prop.cancertissue[grepl("cancer", group)]

samplemeta.cancer = fread(paste0(dir.data, "samplemeta.cancer.csv"))

dtf.prop.cancertissue$stage = samplemeta.cancer[
  match(dtf.prop.cancertissue$sampleid.pid, samplemeta.cancer$sampleid.pid)]$stage
dtf.prop.cancertissue$type = samplemeta.cancer[
  match(dtf.prop.cancertissue$sampleid.pid, samplemeta.cancer$sampleid.pid)]$type
dtf.prop.cancertissue$treatment = samplemeta.cancer[
  match(dtf.prop.cancertissue$sampleid.pid, samplemeta.cancer$sampleid.pid)]$treatment
dtf.prop.cancertissue$group = paste0(dtf.prop.cancertissue$type, ", ",
                                     dtf.prop.cancertissue$treatment)

dtf.prop.cancertissue = dtf.prop.cancertissue[n_cells >= 50]

dtf.prop.cancertissue[, condition := ifelse(group.detail == "PDAC", "PDAC", group)]
dtf.prop.cancertissue[, condition := ifelse(condition == "kidney cancer",
                                            "other cancer", condition)]
table(dtf.prop.cancertissue$condition)

dtf.prop.cancertissue$sex = samplemeta[
  match(dtf.prop.cancertissue$sampleid.pid, samplemeta$sampleid.pid)]$sex
dtf.prop.cancertissue$age = samplemeta[
  match(dtf.prop.cancertissue$sampleid.pid, samplemeta$sampleid.pid)]$age

dtf.prop.cancertissue = dtf.prop.cancertissue[!type %in% "metastasis"]
dtf.prop.cancertissue[, .N] # 155
table(dtf.prop.cancertissue$pid)
dtf.prop.cancertissue[pid %in% c("hu2022", "salcher2022", "wu2024")][
  treatment == "pre-treatment"][, .N] # 107

dtf.prop.cancertissue[pid %in% c("hu2022", "salcher2022", "wu2024", "hu2023",
                                 "wang2023", "zilionis2019")][
  treatment == "pre-treatment"][, .N] # 118

dtf.prop.cancertissue.m = melt(
  dtf.prop.cancertissue, id.vars = c("pid", "sampleid.pid", "condition", "type", "stage",
                   "treatment"), measure.vars = names(col_celltype))
dtf.prop.cancertissue.m = as.data.table(dtf.prop.cancertissue.m)
dtf.prop.cancertissue.m$condition = factor(
  dtf.prop.cancertissue.m$condition,
  levels = c("healthy", "lung cancer", "PDAC", "GI cancer", "other cancer"))

dtf.prop.cancertissue$pid0 = dtf.prop.cancertissue$pid
dtf.prop.cancertissue[, pid := ifelse(pid == "wu2024" & group.detail == "GBC", "wu2024-GBC", pid)]
dtf.prop.cancertissue[, pid := ifelse(pid == "wu2024" & group.detail == "HCC", "wu2024-HCC", pid)]
dtf.prop.cancertissue[, pid := ifelse(pid == "wu2024", "wu2024-other", pid)]
dtf.prop.cancertissue$pid = gsub("hu2022", "hu2022-OSCC", dtf.prop.cancertissue$pid)
dtf.prop.cancertissue$pid = gsub("salcher2022", "salcher2022-LC", dtf.prop.cancertissue$pid)

### bar
pids = c("hu2022", "salcher2022", "wu2024", "wang2023", "qian2020", "hu2023")

dtf.prop.m.agg = dtf.prop.cancertissue.m[pid %in% pids][
  , lapply(.SD, function(x) {mean(x)}),
  by = c("variable", "type", "treatment"), .SDcols = c("value")]

dtf.prop.m.agg$variable = factor(dtf.prop.m.agg$variable, levels = names(col_celltype))
dtf.prop.m.agg$group = paste0(dtf.prop.m.agg$type, ", ", dtf.prop.m.agg$treatment)
dtf.prop.m.agg$group = gsub("treatment", "tx", dtf.prop.m.agg$group)

grouporder.cancer = c("adjacent, pre-treatment", "cancer, pre-treatment",
                      "cancer, post-treatment")
dtf.prop.m.agg$group = factor(dtf.prop.m.agg$group,
                              levels = gsub("treatment", "tx", grouporder.cancer))

p.cancer.distribution =
  ggplot(dtf.prop.m.agg, aes(x = group, y = value, fill = variable)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  theme_expresso(axis_text_size = textsize, legend_text_size = textsize,
                 grid = "none", plot_background = "white",
                 legend_position = c(1.1, 1), legend_key_spacing_y = 0.01,
                 legend_key_height = 0.4, plot_margin = margin(.1, .1, .1, .1, "line")) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(),
        axis.line.x.bottom  = element_line(color = linecolor, linewidth = linesize),
        axis.line.y.left    = element_line(color = linecolor, linewidth = linesize),
        axis.ticks.x.bottom = element_line(color = linecolor, linewidth = linesize),
        axis.ticks.y.left   = element_line(color = linecolor, linewidth = linesize),
        axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.y     = element_text(margin = margin(r = 0.5, l = 0, unit = "pt")),
        axis.text.y      = element_text(margin = margin(r = 0.5, l = 0, unit = "pt")),
        legend.justification = c(0, 1),
        legend.text      = element_text(color = textcolor,
                                        margin = margin(-0.06, 0, -0.06, 0.1, "line")),
        plot.title       = element_text(size = textsize, face = "plain", hjust = 0,
                                        margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1,
                                                        unit = "line"))) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = seq(0, 1, 0.2), position = "left") +
  scale_x_discrete(position = "bottom", labels = cap_first,
                   expand = expansion(add = 0.5)) +
  labs(x = NULL, y = "Subset proportion", title = NULL) +
  scale_fill_manual(values = col_celltype, name = "") +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(size = 1.2)))

### proportion meta-analysis, cancer vs. adjacent
dtf.prop.cancertissue[, class := ifelse(type %in% "cancer", 1, NA)]
dtf.prop.cancertissue[, class := ifelse(type %in% "adjacent", 0, class)]
pids = c("hu2022-OSCC", "salcher2022-LC", "wu2024-GBC", "wu2024-HCC", "wu2024-other")
features = names(col_celltype)

metaobj = lapply(pids, function(i) {
  sub = dtf.prop.cancertissue[pid == i]
  meta_dataset(
    expr = t(as.matrix(sub[, features, with = FALSE])), class = sub$class, label = i)
})
cancer.propmeta = meta_analysis(metaobj, outcome_type = "binary")

cancer.propmeta$per_dataset$dataset = ref_label(cancer.propmeta$per_dataset$dataset)

datadf.cancer.adjacent.summary = as.data.table(cancer.propmeta$summary)
datadf.cancer.adjacent.summary[, sig := ifelse(FDR <= 0.1, "FDR ≤ 0.1", "FDR > 0.1")]
datadf.cancer.adjacent.summary$celltype = factor(
  datadf.cancer.adjacent.summary$gene, levels = names(col_celltype))

nth = function(x, i) {
  x[(i - 1) %% length(x) + 1]
}
ci.value = -qnorm((1 - 0.95) / 2)

datadf.cancer.adjacent.summary$lower =
  datadf.cancer.adjacent.summary$pooled_ES -
  ci.value * datadf.cancer.adjacent.summary$SE
datadf.cancer.adjacent.summary$upper =
  datadf.cancer.adjacent.summary$pooled_ES +
  ci.value * datadf.cancer.adjacent.summary$SE

p.es.summary.cancer.vs.adjacent =
  ggplot(data = datadf.cancer.adjacent.summary, aes(x = gene, y = pooled_ES, alpha = sig)) +
  geom_hline(yintercept = 0, size = 0.2, alpha = 1, color = "grey67") +
  geom_tile(width = 0.7, height = 0.1, position = position_dodge(0.7), color = NA,
            show.legend = FALSE) +
  geom_errorbar(data = datadf.cancer.adjacent.summary,
                aes(x = celltype, ymin = lower, ymax = upper, group = celltype),
                position = position_dodge(0.7), color = "grey7", width = .7, size = 0.2) +
  labs(x = " ", title = NULL,
       y = bquote("Effect size\npre-tx cancer vs. adjacent")) +
  scale_y_continuous(breaks = seq(-10, 2, 1)) +
  scale_alpha_manual(values = c("FDR > 0.1" = 0.2, "FDR ≤ 0.1" = 1)) +
  theme_expresso(axis_text_size = textsize, legend_text_size = textsize,
                 grid = "none", plot_background = "white",
                 legend_position = c(0.96, 1), legend_justification = c(0, 1),
                 legend_key_width = 0.6, legend_key_height = 0.1) +
  theme(plot.title       = element_text(size = textsize, face = "plain",
                                        margin = margin(0, 0, 0, 0, "line")),
        axis.title.x     = element_text(margin = margin(t = 0.05, unit = "line")),
        legend.box = "vertical", legend.spacing.y = unit(0.02, "line"),
        legend.text = element_text(margin = margin(l = 0.4, unit = "pt")),
        legend.key.size  = unit(1, "lines")) +
  coord_flip() +
  scale_x_discrete(limits = rev, position = "bottom") +
  guides(color = "none", shape = "none", size = "none",
         alpha = guide_legend(label.position = "right", nrow = 2, reverse = TRUE,
                              override.aes = list(alpha = 1, linewidth = 0.5,
                                                  color = c("grey17", "grey78"))))


# -- Cancer, gene expression meta-analysis -------------------------------------

cancer.exprmeta = readRDS(paste0(
  dir.results, "meta_analysis/SCAHN.cancer.exprmeta.allgenes.allneu.rds"))
cancer.exprmeta.summary = as.data.table(cancer.exprmeta$summary)
cancer.exprmeta.summary[is.na(pooled_ES)][, .N]
cancer.exprmeta.summary[, .N]
cancer.exprmeta.summary = cancer.exprmeta.summary[!is.na(pooled_ES)]

genes.pctexpr = readRDS(paste0(dir.results, "annotation/gene.exprpct.rds"))

markers = fread(paste0(dir.results, "annotation/SCAHN.subset.markers.csv"))
cancer.exprmeta.summary = cancer.exprmeta.summary[
  gene %in% c(signaturegenes, unique(markers[avg_log2FC >= 1]$gene),
              names(which(genes.pctexpr > 20)))]
cancer.exprmeta.summary[, .N]

cancer.exprmeta$per_dataset$dataset = gsub("hu2022", "hu2022-OSCC",
                                          cancer.exprmeta$per_dataset$dataset)
cancer.exprmeta$per_dataset$dataset = gsub("salcher2022", "salcher2022-LC",
                                          cancer.exprmeta$per_dataset$dataset)

cancer.exprmeta$per_dataset$dataset = ref_label(cancer.exprmeta$per_dataset$dataset)

cancer.exprmeta.summary$FDR = p.adjust(cancer.exprmeta.summary$p_value, method = "fdr")
cancer.exprmeta.summary[gene %in% c("CD274", "OLR1", "MMP9")]

cancer.exprmeta.summary[, label := ifelse(
  (FDR <= 0.05 & abs(pooled_ES) >= 1) | (FDR <= 0.12 & gene %in% signaturegenes), gene, NA)]
cancer.exprmeta.summary[, sig := ifelse(FDR <= 0.1, "FDR≤0.1", "FDR>0.1")]
cancer.exprmeta.summary$sig = factor(cancer.exprmeta.summary$sig, levels = c("FDR≤0.1", "FDR>0.1"))

cancer.exprmeta.summary[, label := ifelse(
  gene %in% c("LYZ", "TUBA1A", "CMTM2", "DHRS7", "CKAP4", "SH3KBP1", "CORO1A",
              "TGOLN2", "RTN3", "CSTA", "USP10", "FOS", "NQO2", "FCN1", "HES4",
              "RGS2", "EMB", "CAMP", "ICAM1"
              ), NA, label)]

cancer.exprmeta.summary[, label.must := ifelse(
  gene %in% c(genes.mature.neu, genes.degranulating.neu) & FDR <= 0.1, label, NA)]
cancer.exprmeta.summary[, label.rest := ifelse(is.na(label.must), label, NA)]

set.seed(42)
p.cancer.expr.meta =
  ggplot(cancer.exprmeta.summary, aes(x = pooled_ES, y = -log10(p_value))) +
  labs(title = NULL, x = "Effect size, pre-tx cancer vs. adjacent") +
  geom_point(size = .2, shape = 16, alpha = 0.8, aes(color = sig)) +
  geom_vline(xintercept = 0, size = 0.1, color = "grey47", alpha = 0.6) +
  geom_text_repel(aes(label = label.rest), size = pt2mm(5), alpha = 1,
                  force = 1, max.overlaps = 10, force_pull = 1, box.padding = 0,
                  fontface = "italic", segment.size = 0.03,
                  min.segment.length = 0.1, show.legend = FALSE) +
  geom_text_repel(aes(label = label.must), size = pt2mm(5), alpha = 1,
                  force = 1, max.overlaps = Inf, force_pull = 1, box.padding = 0,
                  fontface = "italic", segment.size = 0.03,
                  min.segment.length = 0.1, show.legend = FALSE) +
  theme_expresso(panel_background = "white", plot_background = "white",
                 axis_text_size = textsize, axis_title_size = textsize,
                 legend_text_size = textsize,
                 legend_position = c(0.8, 0.9), legend_justification = c(0, 0)) +
  theme(axis.title.x = element_text(margin = margin(t = 0.5, unit = "pt"))) +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  scale_color_manual(values = c("FDR≤0.1" = "indianred", "FDR>0.1" = "grey77")) +
  guides(color = guide_legend(override.aes = list(size = 0.8, alpha = 0.8),
                              label.position = "right", nrow = 2),
         shape = "none", alpha = "none")


# -- Cancer, T/NK and mononuclear myeloid correlation --------------------------

df.ntm.cancer = dtf.prop.neu.tnk[
  sampleid.pid %in% dtf.prop.cancertissue$sampleid.pid]

df.ntm.cancer = df.ntm.cancer[
  tnk.number >= 50 & neu.number >= 50 & myeloid.number >= 50]

df.ntm.cancer = df.ntm.cancer[
  pid %in% names(which(table(df.ntm.cancer$pid) >= 4))]

df.ntm.cancer$type = samplemeta.cancer[
  match(df.ntm.cancer$sampleid.pid, samplemeta.cancer$sampleid.pid)]$type
df.ntm.cancer$treatment = samplemeta.cancer[
  match(df.ntm.cancer$sampleid.pid, samplemeta.cancer$sampleid.pid)]$treatment

sample2exclude = ""
df.ntm.cancer$combined = (df.ntm.cancer$`CCL3/4.N` + df.ntm.cancer$VEGFA.N) /
                         (df.ntm.cancer$neu.number)
metacor.res = meta_corr(df.ntm.cancer, vars = c("SLPI", "MMP9", "VEGFA", "CCL3/4", "CD74",
                                 "IL1RN", "NF-κB", "IL1B",
                                 "CD274", "ARG1", "OLR1", "S100A8", "S100A9",
                                 "CD4 T exhausted", "CD8 T exhausted", "Treg", "NK"),
                        cor_method = "pearson", batch = "pid")

res.per_batch = as.data.table(metacor.res$per_batch)
res.per_batch = res.per_batch[, c("var1", "var2", "n", "r", "ci_lo", "ci_hi", "batch")]
colnames(res.per_batch) = c("var1", "var2", "samplesize", "estimate", "upper", "lower", "rowname")
pids = 1:length(unique(res.per_batch$rowname))
names(pids) = rev(unique(res.per_batch$rowname))
res.per_batch$index = pids[res.per_batch$rowname]
res.per_batch$pval  = NA
res.per_batch$label = NA
res.per_batch$rowname = ref_label(res.per_batch$rowname)

metacorr.cancer.summary = as.data.table(metacor.res$summary)
metacorr.cancer.summary = metacorr.cancer.summary[
  , c("var1", "var2", "n_total", "r_random", "ci_random_lo", "ci_random_hi", "p_random")]
colnames(metacorr.cancer.summary) = c("var1", "var2", "samplesize", "estimate",
                                      "upper", "lower", "pval")
metacorr.cancer.summary$rowname = "summary"
metacorr.cancer.summary$index   = 0
metacorr.cancer.summary$label = paste0(
  "r:", round(metacorr.cancer.summary$estimate, 2), "\n",
  "pval: ", formatC(metacorr.cancer.summary$pval, format = "e", digits = 1))
res.per_batch = res.per_batch[, colnames(metacorr.cancer.summary), with = FALSE]
metacorr.cancer.summary = as.data.table(rbind(metacorr.cancer.summary, res.per_batch))

metacorr.cancer.summary = metacorr.cancer.summary[order(-index)]

plot_forest(metacorr.cancer.summary[var1 == "IL1RN" & var2 == "Treg"][order(-index)],
            boxsize = 5, color = "#DD2F38", plottitle = NULL,
            facet = FALSE, xlab = "correlation")

### meta-correlation heatmap, figure 4
roworders = c("VEGFA", "CCL3/4", "NF-κB", "IL1RN", "SLPI", "CD74", "CD274", "OLR1")
colorders = c("CD4 T exhausted", "CD8 T exhausted", "Treg", "NK")

corr.df = metacorr.cancer.summary[rowname == "summary"][
  var1 %in% roworders & var2 %in% colorders]
corr.df$fdr = p.adjust(corr.df$pval, method = "fdr")
corr.mat = dcast(corr.df, var1 ~ var2, value.var = "estimate")
corr.mat.rowname = corr.mat$var1
corr.mat$var1 = NULL
corr.mat = as.matrix(corr.mat)
rownames(corr.mat) = corr.mat.rowname
fdr.mat = dcast(corr.df, var1 ~ var2, value.var = "fdr")
fdr.mat$var1 = NULL
fdr.mat = as.matrix(fdr.mat)
rownames(fdr.mat) = corr.mat.rowname

corr.mat = corr.mat[roworders, colorders]
fdr.mat  = fdr.mat[roworders, colorders]

rowname = roworders

make_labels = function(x, italicize = character()) {
  expr_txt = vapply(x, function(s) {
    if (s %in% italicize) sprintf("italic('%s')", s) else sprintf("'%s'", s)
  }, character(1))
  parse(text = expr_txt)
}

rowname = c("VEGFA subset", "CCL3/4 subset", "NF-κB subset", "IL1RN subset",
            "SLPI subset", "CD74 subset", "CD274", "OLR1")
rowname = make_labels(rowname, italicize = c("OLR1", "CD274"))

rownames(corr.mat) = rowname
rownames(fdr.mat)  = rowname

size.mat = 1 - fdr.mat
fdr.mat[fdr.mat > 0.1] = 1
size2.mat = 1 - fdr.mat
label.mat = corr.mat
label.mat = round(label.mat, 2)
label.mat[fdr.mat > 0.1] = ""

col_fun_corr_cancer = circlize::colorRamp2(
  seq(-1, 1, length = 13), c(dichromat_pal("DarkRedtoBlue.12")(12)[1:6], "white",
    dichromat_pal("DarkRedtoBlue.12")(12)[7:12]))

ht_opt$TITLE_PADDING = unit(c(0.2, 0.2), "line")

pht.metacorr.cancer =
  Heatmap(corr.mat, rect_gp = gpar(type = "none", fill = "white"),
          col = col_fun_corr_cancer, na_col = "grey", cluster_columns = FALSE, cluster_rows = FALSE,
          show_row_dend = FALSE, show_column_dend = FALSE, column_dend_side = "bottom",
          clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
          clustering_distance_rows = "pearson", clustering_distance_columns = "pearson",
          show_column_names = TRUE, column_names_side = "top", column_names_rot = 45,
          column_labels = gt_render(md_celltype(colnames(corr.mat))),
          column_names_gp = gpar(fontsize = textsize, col = textcolor),
          row_labels = rowname,
          row_names_gp = gpar(fontsize = textsize, col = textcolor),
          column_title = "                Proportion in total T/NK cells",
          column_title_side = "top",
          column_title_gp = gpar(fontsize = textsize, fontface = "plain", hjust = 1,
                                 col = textcolor),
          row_title = "Proportion or expression in neutrophils",
          row_title_gp = gpar(fontsize = textsize, fontface = "plain",
                              col = textcolor), show_heatmap_legend = FALSE,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.rect(x = x, y = y, width = width, height = height,
                      gp = gpar(col = NA, fill = "white", alpha = 0.2))
            grid.circle(x = x, y = y, r = 0.79 * size.mat[i, j] * min(unit.c(width, height)),
                        gp = gpar(fill = col_fun_corr_cancer(corr.mat[i, j]),
                                  col = NA, alpha = 0.9))
            grid.circle(x = x, y = y, r = 0.79 * size2.mat[i, j] * min(unit.c(width, height)),
                        gp = gpar(fill = NA, col = "grey17", alpha = 1, lwd = 0.5))
            grid.text(x = x, y = y, label.mat[i, j],
            gp = gpar(fontsize = 5, col = textcolor))
          })

lgd.metacorr.cancer = Legend(
  col_fun = col_fun_corr_cancer, at = seq(-1, 1, length = 3),
  title = "Correlation", title_gp = gpar(fontsize = textsize, col = textcolor),
  title_position = "topcenter", title_gap = unit(0.3, "line"),
  labels_gp = gpar(fontsize = textsize, col = textcolor),
  legend_width = unit(2, "line"), legend_height = unit(.2, "line"),
  grid_width = unit(2, "line"), grid_height = unit(.2, "line"),
  by_row = FALSE, direction = "horizontal")

### meta-correlation across all T/NK and mononuclear myeloid subsets
subsets_neu_tnk_myeloid = c("AP-1", "AZU1", "CCL3/4", "CD74", "CXCL", "EGR1",
                          "G0S2", "HSP", "IFN", "IL1B", "IL1R2", "IL1RN", "LTF",
                          "MME", "MMP9", "NF-κB", "PTGS2", "S100A4", "SLPI", "TXNIP", "VEGFA",
                          "CD4 T Eff/Mem", "CD4 T IFN", "CD4 T Naive",
                          "CD4 T exhausted", "CD8 T Naive", "CD8 T cytotoxic",
                          "CD8 T exhausted", "NK", "NKT", "Prolif T/NK",
                          "T GZMK", "T/NK CCL3/4", "Treg", "other T", "γδ T",
                          "AP-1 MoMac", "Alveolar Mac", "C1Q MoMac",
                          "CCL2 MoMac", "CCL3/4 MoMac", "CD14 Mono",
                          "CD16 Mono", "CXCL10 MoMac", "DC3", "FOLR2 Mac",
                          "HSP MoMac", "IFN MoMac", "IL1B MoMac", "Int Mono",
                          "MT MoMac", "TREM2 Mac", "VEGFA MoMac", "cDC", "mregDC", "pDC")

subset_pct_insamples = df.ntm.cancer[, lapply(.SD, function(x) {mean(x > 0)}),
                                     .SDcols = subsets_neu_tnk_myeloid]

# if a subset is absent in more than 60% of samples, exclude it from analysis
subsets_exclude = names(subset_pct_insamples)[unlist(subset_pct_insamples[1]) < 0.4]

metacor.res.TNK_myeloid = meta_corr(df.ntm.cancer, vars = setdiff(subsets_neu_tnk_myeloid,
                                                 subsets_exclude), cor_method = "pearson",
                                  batch = "pid")

res.TNK_myeloid.per_batch = as.data.table(metacor.res.TNK_myeloid$per_batch)
res.TNK_myeloid.per_batch = res.TNK_myeloid.per_batch[
  , c("var1", "var2", "n", "r", "ci_lo", "ci_hi", "batch")]
colnames(res.TNK_myeloid.per_batch) = c("var1", "var2", "samplesize", "estimate",
                                      "upper", "lower", "rowname")
pids = 1:length(unique(res.TNK_myeloid.per_batch$rowname))
names(pids) = rev(unique(res.TNK_myeloid.per_batch$rowname))
res.TNK_myeloid.per_batch$index = pids[res.TNK_myeloid.per_batch$rowname]
res.TNK_myeloid.per_batch$pval  = NA
res.TNK_myeloid.per_batch$label = NA

metacorr.TNK_myeloid.cancer.summary = as.data.table(metacor.res.TNK_myeloid$summary)
metacorr.TNK_myeloid.cancer.summary = metacorr.TNK_myeloid.cancer.summary[
  , c("var1", "var2", "n_total", "r_random", "ci_random_lo", "ci_random_hi", "p_random")]
colnames(metacorr.TNK_myeloid.cancer.summary) = c("var1", "var2", "samplesize",
                                                "estimate", "upper", "lower", "pval")
metacorr.TNK_myeloid.cancer.summary$rowname = "summary"
metacorr.TNK_myeloid.cancer.summary$index   = 0
metacorr.TNK_myeloid.cancer.summary$label = paste0(
  "r:", round(metacorr.TNK_myeloid.cancer.summary$estimate, 2), "\n",
  "pval: ", formatC(metacorr.TNK_myeloid.cancer.summary$pval, format = "e", digits = 1))
res.TNK_myeloid.per_batch = res.TNK_myeloid.per_batch[
  , colnames(metacorr.TNK_myeloid.cancer.summary), with = FALSE]
metacorr.TNK_myeloid.cancer.summary = as.data.table(
  rbind(metacorr.TNK_myeloid.cancer.summary, res.TNK_myeloid.per_batch))
metacorr.TNK_myeloid.cancer.summary =
  metacorr.TNK_myeloid.cancer.summary[order(-index)]

### all-compartment heatmap, extended figure 6
roworders = c("AZU1", "LTF", "MMP9", "S100A4", "IL1R2", "MME", "TXNIP", "EGR1",
              "AP-1", "PTGS2", "G0S2", "IFN", "CXCL", "VEGFA",
              "CCL3/4", "NF-κB", "IL1B", "IL1RN", "SLPI", "HSP", "CD74")

colorders = c("CD4 T Eff/Mem", "CD4 T exhausted", "CD4 T IFN", "CD4 T Naive",
              "CD8 T cytotoxic", "CD8 T exhausted", "CD8 T Naive", "NK", "NKT",
              "other T", "Prolif T/NK", "T GZMK", "T/NK CCL3/4", "Treg", "γδ T",
              "CD14 Mono", "CD16 Mono", "Int Mono", "IL1B MoMac", "CCL2 MoMac",
              "CCL3/4 MoMac", "VEGFA MoMac", "CXCL10 MoMac", "IFN MoMac",
              "MT MoMac", "C1Q MoMac", "HSP MoMac", "AP-1 MoMac", "cDC", "DC3", "pDC", "mregDC")

# column compartment annotation: T/NK vs myeloid (Mono, DC, MoMac)
colorders_myeloid = c("CD14 Mono", "CD16 Mono", "Int Mono", "IL1B MoMac",
                    "CCL2 MoMac", "CCL3/4 MoMac", "VEGFA MoMac",
                    "CXCL10 MoMac", "IFN MoMac", "MT MoMac", "C1Q MoMac",
                    "HSP MoMac", "AP-1 MoMac", "cDC", "DC3", "pDC", "mregDC")

compartment.col = factor(ifelse(colorders %in% colorders_myeloid, "myeloid", "T/NK"),
                         levels = c("T/NK", "myeloid"))

col_compartment_TNK_myeloid = c("T/NK" = "#4F8F5B", "myeloid" = "#8C71B5")

colanno.TNK_myeloid = HeatmapAnnotation(
  compartment = compartment.col, col = list(compartment = col_compartment_TNK_myeloid),
  annotation_legend_param = list(
    compartment = list(title = NULL, title_position = "topleft",
                       title_gp = gpar(fontsize = textsize, col = textcolor),
                       labels_gp = gpar(fontsize = textsize, col = textcolor),
                       direction = "vertical", ncol = 1, legend_width = unit(.3, "line"),
                       grid_width = unit(.3, "line"), legend_height = unit(.3, "line"),
                       grid_height = unit(.3, "line"))), show_legend = TRUE,
  show_annotation_name = FALSE, annotation_name_gp = gpar(fontsize = textsize, col = textcolor),
  annotation_name_side = "left", simple_anno_size = unit(0.3, "line"))

corr.df = metacorr.TNK_myeloid.cancer.summary[rowname == "summary"][
  var1 %in% roworders & var2 %in% colorders]
corr.df$fdr = p.adjust(corr.df$pval, method = "fdr")
corr.mat = dcast(corr.df, var1 ~ var2, value.var = "estimate")
corr.mat.rowname = corr.mat$var1
corr.mat$var1 = NULL
corr.mat = as.matrix(corr.mat)
rownames(corr.mat) = corr.mat.rowname
fdr.mat = dcast(corr.df, var1 ~ var2, value.var = "fdr")
fdr.mat$var1 = NULL
fdr.mat = as.matrix(fdr.mat)
rownames(fdr.mat) = corr.mat.rowname

corr.mat = corr.mat[roworders, colorders]
fdr.mat  = fdr.mat[roworders, colorders]

size.mat = 1.2 - fdr.mat
fdr.mat[fdr.mat > 0.1] = 1.2
size2.mat = 1.2 - fdr.mat
label.mat = corr.mat
label.mat = round(label.mat, 2)
label.mat[fdr.mat > 0.1] = ""

col_fun_corr_cancer2 = circlize::colorRamp2(
  seq(-0.8, 0.8, length = 13), c(dichromat_pal("DarkRedtoBlue.12")(12)[1:6], "white",
    dichromat_pal("DarkRedtoBlue.12")(12)[7:12]))

ht_opt$TITLE_PADDING = unit(c(0.8, 0.2), "line")

pht.metacorr.cancer.TNK_myeloid =
  Heatmap(corr.mat, rect_gp = gpar(type = "none", fill = "white"),
          col = col_fun_corr_cancer2, na_col = "grey", cluster_columns = FALSE,
          cluster_rows = FALSE,
          show_row_dend = FALSE, show_column_dend = FALSE, top_annotation = colanno.TNK_myeloid,
          column_dend_side = "bottom",
          clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
          clustering_distance_rows = "pearson", clustering_distance_columns = "pearson",
          show_column_names = TRUE, column_names_side = "top", column_names_rot = 90,
          column_labels = gt_render(md_celltype(colnames(corr.mat))),
          column_names_gp = gpar(fontsize = textsize, lineheight = 0.7, col = textcolor),
          row_names_side = "left",
          row_names_gp = gpar(fontsize = textsize, col = textcolor),
          column_title = "Proportion in T/NK or mononuclear myeloid compartments in tissue",
          column_title_side = "top",
          column_title_gp = gpar(fontsize = textsize, fontface = "plain", hjust = 1,
                                 col = textcolor),
          row_title = "Proportion in neutrophil compartment in tissue", row_title_side = "left",
          row_title_gp = gpar(fontsize = textsize, fontface = "plain",
                              col = textcolor), show_heatmap_legend = TRUE,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.rect(x = x, y = y, width = width, height = height,
                      gp = gpar(col = NA, fill = "white", alpha = 0.2))
            grid.circle(x = x, y = y, r = 0.4 * size.mat[i, j] * min(unit.c(width, height)),
                        gp = gpar(fill = col_fun_corr_cancer2(corr.mat[i, j]),
                                  col = NA, alpha = 0.9))
            grid.circle(x = x, y = y, r = 0.4 * size2.mat[i, j] * min(unit.c(width, height)),
                        gp = gpar(fill = NA, col = "grey17", alpha = 1, lwd = 0.5))
            grid.text(x = x, y = y, label.mat[i, j], gp = gpar(fontsize = 4,
                                col = ifelse(corr.mat[i, j] > 0.5, "white", textcolor)))
          }, heatmap_legend_param = list(
            title = "Correlation", title_gp = gpar(fontsize = textsize, col = textcolor),
            labels_gp = gpar(fontsize = textsize, col = textcolor),
            legend_width = unit(.2, "line"), legend_height = unit(3, "line"),
            grid_width = unit(.2, "line"), grid_height = unit(2, "line"),
            at = seq(-0.8, 0.8, length = 5), by_row = FALSE,
            direction = "vertical", title_position = "leftcenter-rot"))


# -- Cancer, TCGA pan-cancer survival ------------------------------------------

meta_mv = readRDS(paste0(dir.results, "clinical/TCGA/pan_cancer_meta_re_mv_neut.rds"))

### subsets
dtf_mv_tcga = as.data.table(meta_mv$summary)
dtf_mv_tcga = dtf_mv_tcga[gene_set %in% names(col_celltype)]

# order by multivariate HR (descending)
gs_order_tcga = dtf_mv_tcga[order(dtf_mv_tcga$pooled_hr, decreasing = FALSE), gene_set]
dtf_mv_tcga$gene_set = factor(dtf_mv_tcga$gene_set, levels = gs_order_tcga)

dtf_tcga.subset = dtf_mv_tcga

# FDR significance
dtf_tcga.subset$sig = ifelse(
  !is.na(dtf_tcga.subset$FDR) & dtf_tcga.subset$FDR < 0.1, "FDR ≤ 0.1", "FDR > 0.1")
dtf_tcga.subset$sig = factor(dtf_tcga.subset$sig, levels = c("FDR ≤ 0.1", "FDR > 0.1"))

p.meta.tcga.subset = ggplot(dtf_tcga.subset, aes(x = log2(pooled_hr), y = gene_set)) +
  geom_vline(xintercept = 0, size = 0.2, alpha = 1, color = "grey67") +
  geom_errorbarh(aes(xmin = log2(pooled_hr_lo), xmax = log2(pooled_hr_hi), color = sig),
                 height = 0.6, size = 0.4, position = position_dodgev(height = 0.5)) +
  geom_point(shape = 21, size = 0.8, color = "grey52", fill = "grey52",
             position = position_dodgev(height = 0.5)) +
  scale_color_manual(values = sig_colors, name = "FDR", guide = guide_legend(
                       override.aes = list(shape = NA, linetype = "solid", linewidth = 0.5),
                       label.position = "right", nrow = 2)) +
  scale_fill_manual(values = model_fills, name = "", guide = guide_legend(
                      reverse = TRUE, override.aes = list(shape = 21, color = "grey52",
                                          size = 2))) +
  scale_x_continuous(breaks = c(-0.2, 0, 0.2)) +
  labs(
    title = NULL, x = "Log2 hazard ratio\n(OS or PFI in TCGA)",
    y = " ") +
  theme_expresso(axis_title_size = textsize, legend_key_spacing_y = 0.05,
                 legend_key_height = 0.1, legend_key_width = 0.6,
                 legend_text_size = textsize,
                 legend_position = c(1, 1), legend_justification = c(0, 1)) +
  theme(
    panel.grid.major.y = element_blank(), legend.spacing.y = unit(0.2, "line"),
    legend.text = element_text(margin = margin(l = 0.4, unit = "pt")),
    axis.title.x = element_text(margin = margin(t = 0.5, unit = "pt")),
    axis.title.y = element_text(margin = margin(r = 0)))

### signature scores
dtf_mv_tcga = as.data.table(meta_mv$summary)
dtf_mv_tcga$gene_set = gsub(".neu", "", dtf_mv_tcga$gene_set)
dtf_mv_tcga = dtf_mv_tcga[gene_set %in% names(gene_sets)]

# order by multivariate HR (descending)
gs_order_tcga = dtf_mv_tcga[order(dtf_mv_tcga$pooled_hr, decreasing = FALSE), gene_set]
dtf_mv_tcga$gene_set = factor(dtf_mv_tcga$gene_set, levels = gs_order_tcga)

dtf_tcga.sccore = dtf_mv_tcga

# FDR significance
dtf_tcga.sccore$sig = ifelse(
  !is.na(dtf_tcga.sccore$FDR) & dtf_tcga.sccore$FDR < 0.1, "FDR ≤ 0.1", "FDR > 0.1")
dtf_tcga.sccore$sig = factor(dtf_tcga.sccore$sig, levels = c("FDR ≤ 0.1", "FDR > 0.1"))

p.meta.tcga.score = ggplot(dtf_tcga.sccore, aes(x = log2(pooled_hr), y = gene_set)) +
  geom_vline(xintercept = 0, size = 0.2, alpha = 1, color = "grey67") +
  geom_errorbarh(aes(xmin = log2(pooled_hr_lo), xmax = log2(pooled_hr_hi), color = sig),
                 height = 0.6, size = 0.4, position = position_dodgev(height = 0.5)) +
  geom_point(shape = 21, size = 0.8, color = "grey52", fill = "grey52",
             position = position_dodgev(height = 0.5)) +
  scale_color_manual(values = sig_colors, name = "FDR", guide = guide_legend(
                       override.aes = list(shape = NA, linetype = "solid", linewidth = 0.5),
                       label.position = "right", nrow = 2)) +
  scale_fill_manual(values = model_fills, name = "", guide = guide_legend(
                      reverse = TRUE, override.aes = list(shape = 21, color = "grey52",
                                          size = 2))) +
  scale_x_continuous(breaks = c(-0.2, 0, 0.2)) +
  scale_y_discrete(labels = cap_first) +
  labs(
    title = NULL, x = "Log2 hazard ratio\n(OS or PFI in TCGA)",
    y = " \n \n") +
  theme_expresso(axis_title_size = textsize, legend_key_spacing_y = 0.05,
                 legend_key_height = 0.1, legend_key_width = 0.6,
                 legend_text_size = textsize,
                 legend_position = c(1, 1), legend_justification = c(0, 1)) +
  theme(
    panel.grid.major.y = element_blank(), legend.spacing.y = unit(0.2, "line"),
    legend.text = element_text(margin = margin(l = 0.4, unit = "pt")),
    axis.title.x = element_text(margin = margin(t = 0.5, unit = "pt")),
    axis.title.y = element_text(margin = margin(r = 0)),
    plot.title = element_text(hjust = 0, margin = margin(l = -1)))


# -- Figure 4, sample tables and legend helpers --------------------------------

### sample numbers in infectious diseases
studyorder.infect = c("combes2021", "schrepping2020", "zhang2023", "sinha2021",
                      "wilk2021", "kaiser2024", "kwok2023")
sampletable.infect = as.data.table(table(dtf.prop.infection$pid, dtf.prop.infection$condition))
sampletable.infect = as.data.table(dcast(sampletable.infect, V1 ~ V2, value.var = "N"))
sampletable.infect = sampletable.infect[match(studyorder.infect, V1)]
sampletable.infect = sampletable.infect[, conditionorder, with = FALSE]
sampletable.infect = rbind(sampletable.infect,
                           as.list(as.integer(colSums(sampletable.infect))))
sampletable.infect[, study := ref_label(c(studyorder.infect, "total"))]

sampletable.infect = tableGrob(
  sampletable.infect, rows = NULL, cols = NULL, theme = ttheme_default(
    core = list(
      fg_params = list(fontsize = textsize, fontface = "plain", lineheight = 0.8),
      bg_params = list(fill = c("white")), padding = unit(c(0.5, 0.3), "line")), colhead = list(
      fg_params = list(fontsize = textsize, fontface = "plain", rot = 45,
                       hjust = 0.5, vjust = 1, lineheight = 0.8),
      bg_params = list(fill = c("white")), padding = unit(c(0.5, 0.3), "line")), rowhead = list(
      fg_params = list(fontsize = textsize, fontface = "plain"), padding = unit(c(0.5, 0.3),
      "line"))))

nc.infect = ncol(sampletable.infect)

title.infect = textGrob("Sample\nnumber",
                        gp = gpar(fontsize = textsize, col = textcolor,
                                  lineheight = 0.9))
titlewidth.infect = grobWidth(title.infect) + unit(0.3, "line")
sampletable.infect$widths = unit.c(
  rep((unit(1, "npc") - sampletable.infect$widths[nc.infect] - titlewidth.infect) *
        (1 / (nc.infect - 1)), nc.infect - 1),
  sampletable.infect$widths[nc.infect])

for (k in which(sampletable.infect$layout$name == "core-fg" &
                sampletable.infect$layout$l == nc.infect)) {
  sampletable.infect$grobs[[k]]$x     = unit(0.05, "npc")
  sampletable.infect$grobs[[k]]$hjust = 0
}

sampletable.infect = gtable_add_cols(sampletable.infect, widths = titlewidth.infect, pos = 0)
sampletable.infect = gtable_add_grob(
  sampletable.infect, title.infect, t = 1, b = nrow(sampletable.infect), l = 1, r = 1,
  name = "sampletitle")

for (j in 2:ncol(sampletable.infect)) {
  sampletable.infect = gtable_add_grob(
    sampletable.infect,
    segmentsGrob(x0 = 0, x1 = 0, y0 = 0, y1 = 1, gp = gpar(lwd = 0.3, col = "grey42")),
    t = 1, b = nrow(sampletable.infect), l = j, r = j, clip = "off",
    name = paste0("vrule.", j))
}

### sample numbers in cancer
studyorder.cancer = c("hu2022", "salcher2022", "wu2024", "wang2023",
                      "qian2020", "hu2023")
sampletable.cancer = as.data.table(table(dtf.prop.cancertissue$pid0, dtf.prop.cancertissue$group))
sampletable.cancer = as.data.table(dcast(sampletable.cancer, V1 ~ V2, value.var = "N"))
sampletable.cancer = sampletable.cancer[match(studyorder.cancer, V1)]
sampletable.cancer = sampletable.cancer[, grouporder.cancer, with = FALSE]
sampletable.cancer = rbind(sampletable.cancer,
                           as.list(as.integer(colSums(sampletable.cancer))))
sampletable.cancer[, study := ref_label(c(studyorder.cancer, "total"))]

sampletable.cancer = tableGrob(
  sampletable.cancer, rows = NULL, cols = NULL, theme = ttheme_default(
    core = list(
      fg_params = list(fontsize = textsize, fontface = "plain", lineheight = 0.8),
      bg_params = list(fill = c("white")), padding = unit(c(0.5, 0.3), "line")), colhead = list(
      fg_params = list(fontsize = textsize, fontface = "plain", rot = 45,
                       hjust = 0.5, vjust = 1, lineheight = 0.8),
      bg_params = list(fill = c("white")), padding = unit(c(0.5, 0.3), "line")), rowhead = list(
      fg_params = list(fontsize = textsize, fontface = "plain"), padding = unit(c(0.5, 0.3),
      "line"))))

nc.cancer = ncol(sampletable.cancer)
title.cancer = textGrob("Sample\nnumber",
                        gp = gpar(fontsize = textsize, col = textcolor,
                                  lineheight = 0.9))
titlewidth.cancer = grobWidth(title.cancer) + unit(0.3, "line")
sampletable.cancer$widths = unit.c(
  rep((unit(1, "npc") - sampletable.cancer$widths[nc.cancer] - titlewidth.cancer) *
        (1 / (nc.cancer - 1)), nc.cancer - 1),
  sampletable.cancer$widths[nc.cancer])

for (k in which(sampletable.cancer$layout$name == "core-fg" &
                sampletable.cancer$layout$l == nc.cancer)) {
  sampletable.cancer$grobs[[k]]$x     = unit(0.05, "npc")
  sampletable.cancer$grobs[[k]]$hjust = 0
}

sampletable.cancer = gtable_add_cols(sampletable.cancer, widths = titlewidth.cancer, pos = 0)
sampletable.cancer = gtable_add_grob(
  sampletable.cancer, title.cancer, t = 1, b = nrow(sampletable.cancer), l = 1, r = 1,
  name = "sampletitle")

for (j in 2:ncol(sampletable.cancer)) {
  sampletable.cancer = gtable_add_grob(
    sampletable.cancer,
    segmentsGrob(x0 = 0, x1 = 0, y0 = 0, y1 = 1, gp = gpar(lwd = 0.3, col = "grey42")),
    t = 1, b = nrow(sampletable.cancer), l = j, r = j, clip = "off",
    name = paste0("vrule.", j))
}

### subset names
df.subsetlist = data.table(
  subset = names(col_celltype), color = col_celltype, y = 1, x = 1:length(col_celltype))

p.subsetlist = ggplot(df.subsetlist, aes(x = x, y = y, label = subset, color = color)) +
  geom_text(angle = -90, size = pt2mm(textsize), hjust = 0) +
  scale_color_identity() +
  theme_void() +
  coord_cartesian(clip = "off")

### FDR labels
df.twocircles = data.table(
  x = rep(1, 2), y = 1:2, color = c("grey88", "grey17"), labeltext = c("FDR > 0.1", "FDR ≤ 0.1"))

p.twocircles = ggplot(df.twocircles, aes(x = x, y = y, label = labeltext, color = color)) +
  geom_point(shape = 21, fill = "grey88", size = 2) +
  geom_text(angle = 0, size = pt2mm(textsize), hjust = -0.135, color = textcolor) +
  theme(text = element_text(size = textsize, color = textcolor)) +
  scale_color_identity() +
  theme_void() +
  coord_cartesian(clip = "off")

p.twocircles.t = ggplot(df.twocircles, aes(x = y, y = x, label = labeltext, color = color)) +
  geom_point(shape = 21, fill = "grey88", size = 2, stroke = 0.3) +
  geom_text(angle = 90, size = pt2mm(textsize), hjust = -0.135, color = textcolor) +
  theme(text = element_text(size = textsize, color = textcolor)) +
  scale_color_identity() +
  scale_x_reverse() +
  theme_void() +
  coord_cartesian(clip = "off")

df.legend_fdr = data.table(
  y = c(1, 2), color = c("grey17", "grey88"), label = c("FDR ≤ 0.1", "FDR > 0.1"))

p.legend_fdr = ggplot(data = df.legend_fdr) +
  geom_text(aes(x = 0, y = y + 0.15, label = label),
            angle = 90, hjust = 0, size = pt2mm(textsize), color = textcolor) +
  ggplot2::geom_segment(aes(x = -0.3, xend = 0.3, y = y, yend = y, color = color), linewidth = 1) +
  scale_color_identity() +
  xlim(-0.5, 0.5) +
  ylim(0.5, 3.5) +
  theme_void() +
  coord_cartesian(clip = "off")

### shapes for cancer
df.shapes = data.table(
  y = rep(1, 2), x = 1:2, shape = c(16, 17), labeltext = c("adjacent", "cancer"))

p.shapes = ggplot(df.shapes, aes(x = x, y = y, label = labeltext, shape = shape)) +
  geom_point(color = "grey17", fill = "grey17", size = 2) +
  geom_text(angle = 0, size = pt2mm(textsize), hjust = -0.2, vjust = 0.5, color = textcolor) +
  theme(text = element_text(size = textsize, color = textcolor)) +
  scale_shape_identity() +
  theme_void() +
  coord_cartesian(clip = "off")


# -- Figure 4 ------------------------------------------------------------------

ht_opt$TITLE_PADDING = unit(c(0.2, 0.2), "line")
cairo_pdf(file = paste0(dir.fig, "Figure4.pdf"), width = 7.08, height = 9.45)
pushViewport(viewport(layout = grid.layout(nrow = 545, ncol = 470)))

print(p.infectdistribution + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "line")),
      vp = viewport(layout.pos.row = 5:120, layout.pos.col = 6:55))

print(p.cancer.distribution + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "line")),
      vp = viewport(layout.pos.row = 5:99, layout.pos.col = 262:302))

pushViewport(viewport(layout.pos.row = 121:178, layout.pos.col =  2:76))
grid.draw(sampletable.infect)
popViewport()

pushViewport(viewport(layout.pos.row = 121:170, layout.pos.col = 257:323))
grid.draw(sampletable.cancer)
popViewport()

### infection
print(p.es.summary.infect + theme(plot.margin = margin(0.1, 0.3, 0.1, 0.1, "line")),
      vp = viewport(layout.pos.row = 5:180, layout.pos.col = 93:200))

print(
  plot_forest(transform(as_forest_data(infect.exprmeta$`COVID-19,severe`, gene = "CD274",
                                       rowname_style = "name"),
                        gene = "Severe COVID-19"),
              boxsize = 0.5, plottitle = NULL, color = "#DD2F38",
              facet_font = "plain", y_position = "left", title_position = "center",
              summary_space = 1.4,
              xlab = "Effect size\nvs. healthy/convalescent") +
    theme(plot.margin = unit(c(0.42, 0.3, 0.3, 0.1), "line"),
          plot.background = element_rect(color = "white", fill = "white"),
          strip.text.x = element_markdown(size = textsize, lineheight = 0.9,
                                          margin = margin(0.08, 0, 0.08, 0, "line"))) +
    scale_x_continuous(breaks = seq(-1, 3, 1)),
  vp = viewport(layout.pos.row = 229:309, layout.pos.col = 8:101))

print(
  plot_forest(transform(as_forest_data(infect.exprmeta$`bacterial sepsis`, gene = "CD274",
                                       rowname_style = "name"),
                        gene = "Bacterial sepsis"),
              boxsize = 0.5, plottitle = NULL, color = "#DD2F38",
              facet_font = "plain", y_position = "left", title_position = "center",
              summary_space = 1.4,
              xlab = "Effect size\nvs. healthy/convalescent") +
    theme(plot.margin = unit(c(0.42, 0.3, 0.3, 0.1), "line"),
          plot.background = element_rect(color = "white", fill = "white"),
          strip.text.x = element_markdown(size = textsize, lineheight = 0.9,
                                          margin = margin(0.08, 0, 0.08, 0, "line"))) +
    scale_x_continuous(breaks = seq(-1, 3, 1)),
  vp = viewport(layout.pos.row = 229:309, layout.pos.col = 104:193))


pushViewport(viewport(layout.pos.row = 184:225, layout.pos.col = 8:160))
draw(pht.allneu.infect, merge_legend = FALSE, heatmap_legend_side = "right",
     align_heatmap_legend = "heatmap_top",
     newpage = FALSE, background = "white", padding = unit(c(0.05, 0.05, 0.05, 0.05), "line"))
popViewport()

heatmap_legend_padding.default = ht_opt$HEATMAP_LEGEND_PADDING
ht_opt$HEATMAP_LEGEND_PADDING = unit(0, "mm")

pushViewport(viewport(layout.pos.row = 310:420, layout.pos.col = 10:100))
draw(pht.metacorr.infect, newpage = FALSE, background = "white",
     padding = unit(c(0.05, 0.05, 0.05, 0.05), "line"))
popViewport()

pushViewport(viewport(layout.pos.row = 312:327, layout.pos.col = 93:120))
draw(lgd.metacorr.infect)
popViewport()

print(p.meta.public.subset + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "line")),
      vp = viewport(layout.pos.row = 315:459, layout.pos.col = 130:213))
print(p.meta.public.score + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "line")),
      vp = viewport(layout.pos.row = 463:524, layout.pos.col = 103:213))

print(res_sm$IL1B$plot_bar, vp = viewport(layout.pos.row = 430:485, layout.pos.col = 14:83))
print(res_sm$IL1B$plot_marg, vp = viewport(layout.pos.row = 493:540, layout.pos.col = 6:61))
print(km_il1b$low$plot + labs(title = "*IL1B*<sup>lo</sup>"),
      vp = viewport(layout.pos.row = 487:535, layout.pos.col = 63:117))

### cancer
print(p.es.summary.cancer.vs.adjacent + theme(plot.margin = margin(0.1, 0.3, 0.1, 0.1, "line")),
      vp = viewport(layout.pos.row = 5:150, layout.pos.col = 335:441))

pushViewport(viewport(layout.pos.row = 310:420, layout.pos.col = 246:350))
draw(pht.metacorr.cancer, newpage = FALSE, background = "white",
     padding = unit(c(0.05, 0.05, 0.05, 0.05), "line"))
popViewport()

pushViewport(viewport(layout.pos.row = 312:327, layout.pos.col = 331:358))
draw(lgd.metacorr.cancer)
popViewport()

ht_opt$HEATMAP_LEGEND_PADDING = heatmap_legend_padding.default

print(
  plot_forest(transform(as_forest_data(cancer.propmeta, gene = "CCL3/4",
                                       rowname_style = "name"),
                        gene = " "),
              boxsize = 0.5, plottitle = NULL, color = "#DD2F38",
              facet_font = "plain", title_position = "center", y_position = "left",
              summary_space = 1.4,
              xlab = "Effect size\npre-tx cancer vs. adjacent") +
    theme(plot.margin = unit(c(0.03, 0.1, 0.1, 0.1), "line"),
          plot.background = element_rect(color = "white", fill = "white"),
          strip.text.x = element_markdown(size = textsize, lineheight = 0.9,
                                          margin = margin(0.08, 0, 0.08, 0, "line"))) +
    scale_x_continuous(breaks = -3:3),
  vp = viewport(layout.pos.row = 151:225, layout.pos.col = 365:462))

print(
  plot_forest(transform(as_forest_data(cancer.exprmeta, gene = "CD274",
                                       rowname_style = "name"),
                        gene = " "),
              boxsize = 0.5, plottitle = NULL, color = "#DD2F38",
              facet_font = "plain", title_position = "center", y_position = "left",
              summary_space = 1.4,
              xlab = "Effect size\npre-tx cancer vs. adjacent") +
    theme(plot.margin = unit(c(0.03, 0.1, 0.1, 0.1), "line"),
          plot.background = element_rect(color = "white", fill = "white"),
          strip.text.x = element_markdown(size = textsize, lineheight = 0.9,
                                          margin = margin(0.08, 0, 0.08, 0, "line"))) +
    scale_x_continuous(breaks = -3:5),
  vp = viewport(layout.pos.row = 227:301, layout.pos.col = 365:462))


print(p.cancer.expr.meta + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "line")),
      vp = viewport(layout.pos.row = 187:302, layout.pos.col = 214:343))

print(p.meta.tcga.subset + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "line")),
      vp = viewport(layout.pos.row = 315:459, layout.pos.col = 358:438))
print(p.meta.tcga.score + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "line")),
      vp = viewport(layout.pos.row = 463:524, layout.pos.col = 332:438))

print(
  plot_forest(metacorr.cancer.summary[var1 == "VEGFA" & var2 == "CD8 T exhausted"][order(-index)],
              boxsize = 5, color = "#DD2F38",
              plottitle = NULL,
              facet = FALSE, y_position = "left", summary_space = 1.4,
              xlab = paste0("Correlation between<br>VEGFA subset &<br>",
                            "CD8<sup>+</sup> T<sub>EX</sub> proportion")) +
    theme(plot.margin = unit(c(0.03, 0.04, 0.04, 0.1), "line"),
          plot.background = element_rect(color = "white", fill = "white"),
          axis.title.x = element_markdown(size = textsize, color = textcolor,
                                          lineheight = 0.9, hjust = 0.5,
                                          margin = margin(t = 0.5, unit = "pt"))),
  vp = viewport(layout.pos.row = 458:532, layout.pos.col = 253:346))

### legends
pushViewport(viewport(layout.pos.row = 42:67, layout.pos.col = 205:244))
draw(Legend(title = NULL, title_gp = gpar(fontsize = textsize, col = textcolor),
            title_position = "topcenter", nrow = 3,
            labels = cap_severity(names(col_infect)),
            by_row = TRUE,
            labels_gp = gpar(fontsize = textsize, col = textcolor, fontface = "plain"),
            grid_height = unit(0.3, "line"), grid_width = unit(0.1, "line"), gap = unit(.1, "line"),
            legend_gp = gpar(alpha = 0.8, col = col_infect, fill = col_infect)))
upViewport(1)

print(p.twocircles + theme(plot.margin = margin(0, 0, 0, 0, "line")),
      vp = viewport(layout.pos.row = 195:200, layout.pos.col = 122:136))
print(p.twocircles + theme(plot.margin = margin(0, 0, 0, 0, "line")),
      vp = viewport(layout.pos.row = 330:336, layout.pos.col = 93:95))
print(p.twocircles + theme(plot.margin = margin(0, 0, 0, 0, "line")),
      vp = viewport(layout.pos.row = 330:336, layout.pos.col = 331:332))

### panel titles and labels
grid.text(x = unit(0.01, "npc"), y = unit(0.99, "npc"), label = "a",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.23, "npc"), y = unit(0.99, "npc"), label = "b",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.01, "npc"), y = unit(0.65, "npc"), label = "c",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.01, "npc"), y = unit(0.57, "npc"), label = "d",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.01, "npc"), y = unit(0.425, "npc"), label = "e",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.29, "npc"), y = unit(0.425, "npc"), label = "f",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.29, "npc"), y = unit(0.155, "npc"), label = "g",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.01, "npc"), y = unit(0.21, "npc"), label = "h",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))


grid.text(x = unit(0.55, "npc"), y = unit(0.99, "npc"), label = "i",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.73, "npc"), y = unit(0.99, "npc"), label = "j",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.78, "npc"), y = unit(0.72, "npc"), label = "k",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.49, "npc"), y = unit(0.67, "npc"), label = "l",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.78, "npc"), y = unit(0.58, "npc"), label = "m",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.54, "npc"), y = unit(0.425, "npc"), label = "n",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.54, "npc"), y = unit(0.17, "npc"), label = "o",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.78, "npc"), y = unit(0.425, "npc"), label = "p",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.78, "npc"), y = unit(0.155, "npc"), label = "q",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))

dev.off()


# -- Extended figure 5, infection forest plots ---------------------------------

cairo_pdf(file = paste0(dir.fig, "Fig.e5.pdf"), width = 7.08, height = 8)
pushViewport(viewport(layout = grid.layout(nrow = 110, ncol = 77)))

pushViewport(viewport(layout.pos.row = 6:80, layout.pos.col = 50:60))
draw(pht.granulesub, merge_legend = FALSE, heatmap_legend_side = "left",
     newpage = FALSE, padding = unit(c(0.05, 0.05, 0.05, 0.05), "line"))
upViewport(1)

pushViewport(viewport(layout.pos.row = 6:80, layout.pos.col = 61:76))
draw(pht.IFNsub, merge_legend = FALSE, heatmap_legend_side = "right",
     align_heatmap_legend = "heatmap_top", newpage = FALSE,
     padding = unit(c(0.05, 0.05, 0.05, 0.05), "line"))
upViewport(1)

print(p.twocircles + theme(plot.margin = margin(0, 0, 0, 0, "line")),
      vp = viewport(layout.pos.row = 20:22, layout.pos.col = 71:72))

pushViewport(viewport(layout.pos.row = 82:86, layout.pos.col = 46:71))
draw(Legend(title = NULL, title_gp = gpar(fontsize = textsize, col = textcolor),
            title_position = "topcenter", labels = cap_severity(names(col_infect)), nrow = 3,
            by_row = TRUE,
            labels_gp = gpar(fontsize = textsize, col = textcolor, fontface = "plain"),
            grid_height = unit(0.5, "line"), grid_width = unit(0.2, "line"), gap = unit(.2, "line"),
            legend_gp = gpar(alpha = 0.8, col = col_infect, fill = col_infect)))
upViewport(1)

cap_rownames = function(d) {
  d$rowname = ifelse(d$rowname == "summary", d$rowname, ref_label(d$rowname))
  d
}

### proportion, granule subsets
print(
  plot_forest(as_forest_data(meta.nonsevere, gene = "MMP9", rowname_style = "n1n0"),
              boxsize = 0.5, plottitle = " ", color = "grey67") +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "line"), axis.text.y = element_blank()),
  vp = viewport(layout.pos.row = 4:19, layout.pos.col = 35:44))

print(
  plot_forest(as_forest_data(meta.nonsevere, gene = "LTF", rowname_style = "n1n0"),
              boxsize = 0.5, plottitle = " ", color = "grey67") +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "line"), axis.text.y = element_blank()),
  vp = viewport(layout.pos.row = 4:19, layout.pos.col = 25:34))

print(
  plot_forest(cap_rownames(as_forest_data(meta.nonsevere, gene = "AZU1", rowname_style = "n1n0")),
              boxsize = 0.5, plottitle = "Non-severe COVID-19/flu vs. healthy/convalescent",
              color = "grey67") +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.3), "line")),
  vp = viewport(layout.pos.row = 4:19, layout.pos.col = 2:24))

print(
  plot_forest(as_forest_data(meta.COVID19.severe, gene = "MMP9", rowname_style = "n1n0"),
              boxsize = 0.5, plottitle = " ", color = "#DD2F38") +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "line"), axis.text.y = element_blank()),
  vp = viewport(layout.pos.row = 21:37, layout.pos.col = 35:44))

print(
  plot_forest(as_forest_data(meta.COVID19.severe, gene = "LTF", rowname_style = "n1n0"),
              boxsize = 0.5, plottitle = " ", color = "#DD2F38") +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "line"), axis.text.y = element_blank()),
  vp = viewport(layout.pos.row = 21:37, layout.pos.col = 25:34))

print(
  plot_forest(cap_rownames(as_forest_data(meta.COVID19.severe, gene = "AZU1",
                                          rowname_style = "n1n0")),
              boxsize = 0.5, plottitle = "Severe COVID-19 vs. healthy/convalescent",
              color = "#DD2F38") +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "line")),
  vp = viewport(layout.pos.row = 21:37, layout.pos.col = 2:24))

print(
  plot_forest(as_forest_data(meta.sepsis, gene = "MMP9", rowname_style = "n1n0"),
              boxsize = 0.5, plottitle = " ", color = "#DD2F38") +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "line"), axis.text.y = element_blank()),
  vp = viewport(layout.pos.row = 39:55, layout.pos.col = 35:44))

print(
  plot_forest(as_forest_data(meta.sepsis, gene = "LTF", rowname_style = "n1n0"),
              boxsize = 0.5, plottitle = " ", color = "#DD2F38") +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "line"), axis.text.y = element_blank()),
  vp = viewport(layout.pos.row = 39:55, layout.pos.col = 25:34))

print(
  plot_forest(cap_rownames(as_forest_data(meta.sepsis, gene = "AZU1", rowname_style = "n1n0")),
              boxsize = 0.5, plottitle = "Bacterial sepsis vs. healthy/convalescent",
              color = "#DD2F38") +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.85), "line")),
  vp = viewport(layout.pos.row = 39:55, layout.pos.col = 2:24))

### proportion, IFN subsets
print(
  plot_forest(as_forest_data(meta.nonsevere, gene = "IFN3", rowname_style = "n1n0"),
              boxsize = 0.5, plottitle = " ", color = "#DD2F38") +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "line"), axis.text.y = element_blank()),
  vp = viewport(layout.pos.row = 58:73, layout.pos.col = 35:44))

print(
  plot_forest(as_forest_data(meta.nonsevere, gene = "IFN2", rowname_style = "n1n0"),
              boxsize = 0.5, plottitle = " ", color = "grey67") +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "line"), axis.text.y = element_blank()),
  vp = viewport(layout.pos.row = 58:73, layout.pos.col = 25:34))

print(
  plot_forest(cap_rownames(as_forest_data(meta.nonsevere, gene = "IFN1", rowname_style = "n1n0")),
              boxsize = 0.5, plottitle = "Non-severe COVID-19/flu vs. healthy/convalescent",
              color = "#DD2F38") +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.3), "line")),
  vp = viewport(layout.pos.row = 58:73, layout.pos.col = 2:24))

print(
  plot_forest(as_forest_data(meta.COVID19.severe, gene = "IFN3", rowname_style = "n1n0"),
              boxsize = 0.5, plottitle = " ", color = "#DD2F38") +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "line"), axis.text.y = element_blank()) +
    scale_x_continuous(breaks = seq(-2, 2, 1)),
  vp = viewport(layout.pos.row = 75:91, layout.pos.col = 35:44))

print(
  plot_forest(as_forest_data(meta.COVID19.severe, gene = "IFN2", rowname_style = "n1n0"),
              boxsize = 0.5, plottitle = " ", color = "grey67") +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "line"), axis.text.y = element_blank()),
  vp = viewport(layout.pos.row = 75:91, layout.pos.col = 25:34))

print(
  plot_forest(cap_rownames(as_forest_data(meta.COVID19.severe, gene = "IFN1",
                                          rowname_style = "n1n0")),
              boxsize = 0.5, plottitle = "Severe COVID-19 vs. healthy/convalescent",
              color = "#DD2F38") +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "line")) +
    scale_x_continuous(breaks = seq(-2, 2, 1)),
  vp = viewport(layout.pos.row = 75:91, layout.pos.col = 2:24))

print(
  plot_forest(as_forest_data(meta.sepsis, gene = "IFN3", rowname_style = "n1n0"),
              boxsize = 0.5, plottitle = " ", color = "grey67") +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "line"), axis.text.y = element_blank()) +
    scale_x_continuous(breaks = seq(-2, 2, 1)),
  vp = viewport(layout.pos.row = 93:109, layout.pos.col = 35:44))

print(
  plot_forest(as_forest_data(meta.sepsis, gene = "IFN2", rowname_style = "n1n0"),
              boxsize = 0.5, plottitle = " ", color = "grey67") +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "line"), axis.text.y = element_blank()),
  vp = viewport(layout.pos.row = 93:109, layout.pos.col = 25:34))

print(
  plot_forest(cap_rownames(as_forest_data(meta.sepsis, gene = "IFN1", rowname_style = "n1n0")),
              boxsize = 0.5, plottitle = "Bacterial sepsis vs. healthy/convalescent",
              color = "grey67") +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.85), "line")),
  vp = viewport(layout.pos.row = 93:109, layout.pos.col = 2:24))

grid.text(x = unit(0.02, "npc"), y = unit(0.98, "npc"), label = "a",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.02, "npc"), y = unit(0.48, "npc"), label = "b",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.62, "npc"), y = unit(0.98, "npc"), label = "c",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))

dev.off()


# -- Extended figure 6, T/NK and myeloid UMAPs ---------------------------------

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


dtf.tnk = readRDS(paste0(dir.data, "TNK.cellmeta.rds"))
dtf.tnk$sampleid.pid = gsub("[ab]$", "", dtf.tnk$sampleid.pid)
setdiff(unique(dtf.tnk$sampleid.pid), samplemeta$sampleid.pid)

dtf.tnk$tissue = samplemeta[match(dtf.tnk$sampleid.pid, samplemeta$sampleid.pid)]$tissue
dtf.tnk$group = samplemeta[match(dtf.tnk$sampleid.pid, samplemeta$sampleid.pid)]$group

dtf.tnk$pid = gsub("a$|b$|c$|d$", "", dtf.tnk$pid)
dtf.tnk$pid = gsub("myin2024b", "yin2024", dtf.tnk$pid)
dtf.tnk$pid = gsub("reyfman2018", "reyfman2019", dtf.tnk$pid)
dtf.tnk$pid = gsub("xue2022", "xue2023", dtf.tnk$pid)
dtf.tnk = dtf.tnk[!grepl("doublet", celltype)]
dtf.tnk[, tissue2 := ifelse(tissue %in% c("peripheral blood", "lung"), tissue, "other")]

dtf.tnk = dtf.tnk[sampleid.pid %in% dtf.prop.neu.tnk$sampleid.pid]
dtf.tnk[, .N]
length(unique(dtf.tnk$sampleid.pid))

pointsizen       = 0.01
pointsizen.small = 0.0001

### study
dtf.tnk$facetvar = "Study"
p.umap.pid = plot_scatter(
  as.data.frame(dtf.tnk), x = "UMAP_1", y = "UMAP_2",
  color_by = "pid", color_type = "discrete", colors = col_study, na_color = "grey42",
  facet_by = "facetvar", facet_nrow = 5,
  point_size = pointsizen.small, point_alpha = 0.2, raster_dpi = 400,
  shuffle = TRUE, seed = 13, show_axis = FALSE, grid = "none", title = "") +
  theme_expresso(
    axis_text_size = textsize, legend_text_size = textsize, show_axis = FALSE, grid = "none",
    plot_background = "white", legend_position = "none", legend_justification = c(0, 0.5),
    legend_key_width = 0.2, legend_key_height = 0.3,
    legend_key_spacing_x = 0.05, legend_key_spacing_y = 0.05, facet_label_face = "plain") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol = 2, label.hjust = 0, override.aes = list(size = 2, alpha = 1)))

### condition
dtf.tnk$facetvar = "Condition"
dtf.tnk$group = factor(dtf.tnk$group, levels = names(col_group))

p.umap.group = plot_scatter(
  as.data.frame(dtf.tnk), x = "UMAP_1", y = "UMAP_2",
  color_by = "group", color_type = "discrete", colors = col_group, na_color = "grey42",
  facet_by = "facetvar", facet_nrow = 5,
  point_size = pointsizen.small, point_alpha = 0.2, raster_dpi = 400,
  shuffle = TRUE, seed = 13, show_axis = FALSE, grid = "none", title = "") +
  theme_expresso(
    axis_text_size = textsize, legend_text_size = textsize, show_axis = FALSE, grid = "none",
    plot_background = "white", legend_position = "none", legend_justification = c(0, 0.5),
    legend_key_width = 0.2, legend_key_height = 0.3,
    legend_key_spacing_x = 0.05, legend_key_spacing_y = 0.05, facet_label_face = "plain") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol = 2, label.hjust = 0, override.aes = list(size = 2, alpha = 1)))

### tissue
col_tissue2 = c(
  "peripheral blood" = "#73DAFF", "lung" = "#DC3023", "other" = "grey77")
length(unique(dtf.tnk$tissue))
dtf.tnk$facetvar = "Tissue"
dtf.tnk$tissue2 = factor(dtf.tnk$tissue2, levels = names(col_tissue2))

p.umap.tissue = plot_scatter(
  as.data.frame(dtf.tnk), x = "UMAP_1", y = "UMAP_2",
  color_by = "tissue2", color_type = "discrete", colors = col_tissue2, na_color = "grey42",
  facet_by = "facetvar", facet_nrow = 5,
  point_size = pointsizen.small, point_alpha = 0.2, raster_dpi = 400,
  shuffle = TRUE, seed = 13, show_axis = FALSE, grid = "none", title = "") +
  theme_expresso(
    axis_text_size = textsize, legend_text_size = textsize, show_axis = FALSE, grid = "none",
    plot_background = "white", legend_position = "none", legend_justification = c(0, 0.5),
    legend_key_width = 0.2, legend_key_height = 0.3,
    legend_key_spacing_x = 0.05, legend_key_spacing_y = 0.05, facet_label_face = "plain") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol = 1, label.hjust = 0, override.aes = list(size = 2, alpha = 1)))

### T/NK subsets
col_celltype_TNK = c(`CD4 T Eff/Mem` = "#FFA631", `CD4 T exhausted` = "#F9906F",
                     `CD4 T IFN` = "#214DC8", `CD4 T Naive` = "#55BB8A",
                     `CD8 T cytotoxic` = "#06A5C7", `CD8 T exhausted` = "#DC3023",
                     `CD8 T Naive` = "#8CC269", NK = "#B0A4E3",
                     NKT = "#EAFF56", `other T` = "grey87",
                     `Prolif T/NK` = "#F8D626", `T GZMK` = "#FFF143",
                     `T/NK CCL3/4` = "#acf300", Treg = "#E4C6D0", `γδ T` = "#F3D3E7")

df.label = as.data.table(dtf.tnk)[, lapply(.SD, median), by = celltype,
                                  .SDcols = c("UMAP_1", "UMAP_2")]
colnames(df.label)[2:3] = c("x", "y")
df.label$labeltext = df.label[, "celltype", with = FALSE]
df.label[, y := ifelse(celltype == "CD8 T exhausted", y + 1, y)]
df.label[, y := ifelse(celltype == "T/NK CCL3/4", y + 0.5, y)]
df.label[, y := ifelse(celltype == "CD8 T cytotoxic", y - 0.5, y)]
df.label[, x := ifelse(celltype == "CD8 T Naive", x + 2, x)]
df.label[, labeltext.md := md_celltype(as.character(celltype))]

dtf.tnk$facetvar = "T/NK subsets"
p.umap.celltype = plot_scatter(
  as.data.frame(dtf.tnk), x = "UMAP_1", y = "UMAP_2",
  color_by = "celltype", color_type = "discrete", colors = col_celltype_TNK, na_color = "grey42",
  facet_by = "facetvar", facet_nrow = 5,
  point_size = pointsizen, point_alpha = 0.2, raster_dpi = 300,
  shuffle = TRUE, seed = 13, show_axis = FALSE, grid = "none", title = "",
  label = FALSE) +
  geom_richtext(data = df.label, aes(x = x, y = y, label = labeltext.md),
                inherit.aes = FALSE, size = pt2mm(5), colour = textcolor,
                alpha = 0.8, fontface = "plain", fill = NA, label.colour = NA,
                label.padding = unit(rep(0, 4), "pt")) +
  theme_expresso(
    axis_text_size = textsize, legend_text_size = textsize, show_axis = FALSE, grid = "none",
    plot_background = "white", legend_position = "none", legend_justification = c(0, 0.5),
    legend_key_width = 0.2, legend_key_height = 0.7,
    legend_key_spacing_x = 0.05, legend_key_spacing_y = 0.05, facet_label_face = "plain") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol = 1, label.hjust = 0, override.aes = list(size = 2, alpha = 1)))

### mononuclear myeloid cells
dtf.myeloid = fread(paste0(dir.data, "myeloid.cellmeta.csv"))
umap.myeloid = fread(paste0(dir.data, "myeloid.cellmeta_umap.csv"),
                   select = c("rowname", "UMAP_1", "UMAP_2", "sampleid.pid"))
dtf.myeloid = merge(dtf.myeloid, umap.myeloid, by = "rowname")
dtf.myeloid[, celltype := celltype.rescued2]

dtf.myeloid$sampleid.pid = gsub("[ab]$", "", dtf.myeloid$sampleid.pid)
setdiff(unique(dtf.myeloid$sampleid.pid), samplemeta$sampleid.pid)

dtf.myeloid$tissue = samplemeta[match(dtf.myeloid$sampleid.pid, samplemeta$sampleid.pid)]$tissue
dtf.myeloid$group = samplemeta[match(dtf.myeloid$sampleid.pid, samplemeta$sampleid.pid)]$group

dtf.myeloid$pid = gsub("a$|b$|c$|d$", "", dtf.myeloid$pid)
dtf.myeloid$pid = gsub("myin2024b", "yin2024", dtf.myeloid$pid)
dtf.myeloid$pid = gsub("reyfman2018", "reyfman2019", dtf.myeloid$pid)
dtf.myeloid$pid = gsub("xue2022", "xue2023", dtf.myeloid$pid)
dtf.myeloid = dtf.myeloid[!grepl("doublet", celltype)]
dtf.myeloid[, tissue2 := ifelse(tissue %in% c("peripheral blood", "lung"), tissue, "other")]

dtf.myeloid = dtf.myeloid[sampleid.pid %in% dtf.prop.neu.tnk$sampleid.pid]
dtf.myeloid[, .N]
length(unique(dtf.myeloid$sampleid.pid))

pointsizen       = 0.01
pointsizen.small = 0.0001

### study
dtf.myeloid$facetvar = "Study"
p.umap.pid.myeloid = plot_scatter(
  as.data.frame(dtf.myeloid), x = "UMAP_1", y = "UMAP_2",
  color_by = "pid", color_type = "discrete", colors = col_study, na_color = "grey42",
  facet_by = "facetvar", facet_nrow = 5,
  point_size = pointsizen.small, point_alpha = 0.2, raster_dpi = 400,
  shuffle = TRUE, seed = 13, show_axis = FALSE, grid = "none", title = "") +
  theme_expresso(
    axis_text_size = textsize, legend_text_size = textsize, show_axis = FALSE, grid = "none",
    plot_background = "white", legend_position = "none", legend_justification = c(0, 0.5),
    legend_key_width = 0.2, legend_key_height = 0.3,
    legend_key_spacing_x = 0.05, legend_key_spacing_y = 0.05, facet_label_face = "plain") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol = 2, label.hjust = 0, override.aes = list(size = 2, alpha = 1)))

### condition
dtf.myeloid$facetvar = "Condition"
dtf.myeloid$group = factor(dtf.myeloid$group, levels = names(col_group))

p.umap.group.myeloid = plot_scatter(
  as.data.frame(dtf.myeloid), x = "UMAP_1", y = "UMAP_2",
  color_by = "group", color_type = "discrete", colors = col_group, na_color = "grey42",
  facet_by = "facetvar", facet_nrow = 5,
  point_size = pointsizen.small, point_alpha = 0.2, raster_dpi = 400,
  shuffle = TRUE, seed = 13, show_axis = FALSE, grid = "none", title = "") +
  theme_expresso(
    axis_text_size = textsize, legend_text_size = textsize, show_axis = FALSE, grid = "none",
    plot_background = "white", legend_position = "none", legend_justification = c(0, 0.5),
    legend_key_width = 0.2, legend_key_height = 0.3,
    legend_key_spacing_x = 0.05, legend_key_spacing_y = 0.05, facet_label_face = "plain") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol = 2, label.hjust = 0, override.aes = list(size = 2, alpha = 1)))

### tissue
col_tissue2 = c(
  "peripheral blood" = "#73DAFF", "lung" = "#DC3023", "other" = "grey77")
length(unique(dtf.myeloid$tissue))
dtf.myeloid$facetvar = "Tissue"
dtf.myeloid$tissue2 = factor(dtf.myeloid$tissue2, levels = names(col_tissue2))

p.umap.tissue.myeloid = plot_scatter(
  as.data.frame(dtf.myeloid), x = "UMAP_1", y = "UMAP_2",
  color_by = "tissue2", color_type = "discrete", colors = col_tissue2, na_color = "grey42",
  facet_by = "facetvar", facet_nrow = 5,
  point_size = pointsizen.small, point_alpha = 0.2, raster_dpi = 400,
  shuffle = TRUE, seed = 13, show_axis = FALSE, grid = "none", title = "") +
  theme_expresso(
    axis_text_size = textsize, legend_text_size = textsize, show_axis = FALSE, grid = "none",
    plot_background = "white", legend_position = "none", legend_justification = c(0, 0.5),
    legend_key_width = 0.2, legend_key_height = 0.3,
    legend_key_spacing_x = 0.05, legend_key_spacing_y = 0.05, facet_label_face = "plain") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol = 1, label.hjust = 0, override.aes = list(size = 2, alpha = 1)))

### mononuclear myeloid subsets
col_celltype_myeloid = c(
  `CD14 Mono` = "#8CC269", `CD16 Mono` = "#AFBC65",
  `Int Mono` = "#55BB8A", `IL1B MoMac` = "#F0BB46",
  `CCL2 MoMac` = "#FFA631", `CCL3/4 MoMac` = "#acf300",
  `CXCL10 MoMac` = "#DC3023", `IFN MoMac` = "#214DC8",
  `C1Q MoMac` = "#EAFF56", `FOLR2 Mac` = "#E4C6D0",
  `VEGFA MoMac` = "#73DAFF", `HSP MoMac` = "#4D6D93", `AP-1 MoMac` = "#2BAE85", `MT MoMac` = "cyan",
  `Alveolar Mac` = "#06A5C7", `TREM2 Mac` = "#B36D61", cDC = "#F3D3E7", DC3 = "#9D4EDD",
  pDC = "#C883C2", mregDC = "#F9906F", other = "grey87")

dtf.myeloid$celltype = factor(dtf.myeloid$celltype, levels = names(col_celltype_myeloid))

df.label = as.data.table(dtf.myeloid)[, lapply(.SD, median), by = celltype,
                                    .SDcols = c("UMAP_1", "UMAP_2")]
colnames(df.label)[2:3] = c("x", "y")
df.label$labeltext = df.label[, "celltype", with = FALSE]
df.label[, x := ifelse(celltype == "AP-1 MoMac", x + 1.2, x)]
df.label[, y := ifelse(celltype == "AP-1 MoMac", y - 0.2, y)]
df.label[, y := ifelse(celltype == "HSP MoMac", y + 0.2, y)]
df.label[, y := ifelse(celltype == "FOLR2 Mac", y - 1, y)]
df.label[, y := ifelse(celltype == "CCL2 MoMac", y - 0.5, y)]
df.label[, y := ifelse(celltype == "TREM2 Mac", y + 0.5, y)]
df.label[, y := ifelse(celltype == "MT MoMac", y - 0.5, y)]
df.label[, labeltext.md := md_celltype(as.character(celltype))]

dtf.myeloid$facetvar = "Mononuclear myeloid subsets"
p.umap.celltype.myeloid = plot_scatter(
  as.data.frame(dtf.myeloid), x = "UMAP_1", y = "UMAP_2",
  color_by = "celltype", color_type = "discrete", colors = col_celltype_myeloid,
  na_color = "grey42", facet_by = "facetvar", facet_nrow = 5,
  point_size = pointsizen, point_alpha = 0.2, raster_dpi = 300,
  shuffle = TRUE, seed = 13, show_axis = FALSE, grid = "none", title = "",
  label = FALSE) +
  geom_richtext(data = df.label, aes(x = x, y = y, label = labeltext.md),
                inherit.aes = FALSE, size = pt2mm(5), colour = textcolor,
                alpha = 0.8, fontface = "plain", fill = NA, label.colour = NA,
                label.padding = unit(rep(0, 4), "pt")) +
  theme_expresso(
    axis_text_size = textsize, legend_text_size = textsize, show_axis = FALSE, grid = "none",
    plot_background = "white", legend_position = "none", legend_justification = c(0, 0.5),
    legend_key_width = 0.2, legend_key_height = 0.7,
    legend_key_spacing_x = 0.05, legend_key_spacing_y = 0.05, facet_label_face = "plain") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol = 2, label.hjust = 0, override.aes = list(size = 2, alpha = 1)))


# -- Extended figure 6 ---------------------------------------------------------

### assemble
ht_opt$TITLE_PADDING = unit(c(0.8, 0.2), "line")
cairo_pdf(file = paste0(dir.fig, "Fig.e6.pdf"), width = 7.08, height = 9.4)
pushViewport(viewport(layout = grid.layout(nrow = 64, ncol = 70)))

print(p.umap.celltype, vp = viewport(layout.pos.row = 1:15, layout.pos.col = 3:24))
print(p.umap.pid, vp = viewport(layout.pos.row = 17:22, layout.pos.col = 3:9))
print(p.umap.tissue, vp = viewport(layout.pos.row = 17:22, layout.pos.col = 10:16))
print(p.umap.group, vp = viewport(layout.pos.row = 17:22, layout.pos.col = 17:23))

print(p.umap.celltype.myeloid, vp = viewport(layout.pos.row = 41:55, layout.pos.col = 3:24))
print(p.umap.pid.myeloid, vp = viewport(layout.pos.row = 57:62, layout.pos.col = 3:9))
print(p.umap.tissue.myeloid, vp = viewport(layout.pos.row = 57:62, layout.pos.col = 10:16))
print(p.umap.group.myeloid, vp = viewport(layout.pos.row = 57:62, layout.pos.col = 17:23))

pushViewport(viewport(layout.pos.row = 1:32, layout.pos.col = 26:70))
draw(pht.metacorr.infect.TNK_myeloid, merge_legend = TRUE,
     heatmap_legend_side = "right", align_heatmap_legend = "heatmap_top",
     newpage = FALSE, background = "white",
     padding = unit(c(0.1, 0.1, 0.1, 0.1), "line"))
popViewport()

pushViewport(viewport(layout.pos.row = 33:64, layout.pos.col = 26:70))
draw(pht.metacorr.cancer.TNK_myeloid, merge_legend = TRUE,
     heatmap_legend_side = "right", align_heatmap_legend = "heatmap_top",
     newpage = FALSE, background = "white",
     padding = unit(c(0.1, 0.1, 0.1, 0.1), "line"))
popViewport()

print(p.twocircles.t, vp = viewport(layout.pos.row = 17:21, layout.pos.col = 68:69))
print(p.twocircles.t, vp = viewport(layout.pos.row = 49:52, layout.pos.col = 68:69))

pushViewport(viewport(layout.pos.row = 23:33, layout.pos.col = 2:20))
draw(Legend(title = "Study", title_gp = gpar(fontsize = textsize, col = textcolor),
            title_position = "topleft", title_gap = unit(0.1, "line"),
            labels = cap_first(names(col_study[sort(unique(dtf.prop.neu.tnk$pid))])), nrow = 11,
            by_row = TRUE,
            type = "points", pch = 16, size = unit(1.5, "mm"), background = "transparent",
            labels_gp = gpar(fontsize = textsize, col = textcolor, fontface = "plain"),
            grid_height = unit(0.5, "line"), grid_width = unit(0.2, "line"), gap = unit(.2, "line"),
            legend_gp = gpar(
              alpha = 0.8, col = col_study[sort(unique(dtf.prop.neu.tnk$pid))],
              fill = col_study[sort(unique(dtf.prop.neu.tnk$pid))])),
     x = unit(0, "npc"), just = "left")
upViewport(1)

pushViewport(viewport(layout.pos.row = 33:34, layout.pos.col = 2:20))
draw(Legend(title = "Tissue", title_gp = gpar(fontsize = textsize, col = textcolor),
            title_position = "topleft", labels = cap_first(names(col_tissue2)), nrow = 1,
            by_row = TRUE,
            title_gap = unit(0.1, "line"), type = "points", pch = 16, size = unit(1.5, "mm"),
            background = "transparent",
            labels_gp = gpar(fontsize = textsize, col = textcolor, fontface = "plain"),
            grid_height = unit(0.5, "line"), grid_width = unit(0.2, "line"), gap = unit(.2, "line"),
            legend_gp = gpar(alpha = 0.8, col = col_tissue2, fill = col_tissue2)),
     x = unit(0, "npc"), just = "left")
upViewport(1)

pushViewport(viewport(layout.pos.row = 35:39, layout.pos.col = 2:20))
draw(Legend(title = "Condition", title_gp = gpar(fontsize = textsize, col = textcolor),
            title_position = "topleft",
            labels = cap_severity(names(col_group[sort(unique(dtf.prop.neu.tnk$group))])),
            nrow = 6, by_row = TRUE, title_gap = unit(0.1, "line"),
            type = "points", pch = 16, size = unit(1.5, "mm"), background = "transparent",
            labels_gp = gpar(fontsize = textsize, col = textcolor, fontface = "plain"),
            grid_height = unit(0.5, "line"), grid_width = unit(0.2, "line"), gap = unit(.2, "line"),
            legend_gp = gpar(
              alpha = 0.8, col = col_group[sort(unique(dtf.prop.neu.tnk$group))],
              fill = col_group[sort(unique(dtf.prop.neu.tnk$group))])),
     x = unit(0, "npc"), just = "left")
upViewport(1)

print(p.umap.axis + theme(plot.margin = margin(0.08, 0, 0, 0.05, "line")),
      vp = viewport(layout.pos.row = 13:16, layout.pos.col = 1:5))
print(p.umap.axis + theme(plot.margin = margin(0.08, 0, 0, 0.05, "line")),
      vp = viewport(layout.pos.row = 20:23, layout.pos.col = 1:5))

print(p.umap.axis + theme(plot.margin = margin(0.08, 0, 0, 0.05, "line")),
      vp = viewport(layout.pos.row = 53:56, layout.pos.col = 1:5))
print(p.umap.axis + theme(plot.margin = margin(0, 0, 0.08, 0.05, "line")),
      vp = viewport(layout.pos.row = 61:64, layout.pos.col = 1:5))

grid.text(x = unit(0.01, "npc"), y = unit(0.985, "npc"), label = "a",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.01, "npc"), y = unit(0.37, "npc"), label = "b",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.41, "npc"), y = unit(0.985, "npc"), label = "c",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.41, "npc"), y = unit(0.48, "npc"), label = "d",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))

dev.off()


# -- Extended figure 7 ---------------------------------------------------------

### univariate + multivariate forest
# Built from the per-cohort meta-analysis outputs, filtered to genesets in
# col_celltype and ordered by multivariate OR. meta_mv leads and meta_uni is
# optional because only the multivariate half is drawn.
make_or_forest = function(meta_mv, meta_uni = NULL, gs_keep = names(col_celltype),
                          title = "", xlab = "Log2 odds ratio",
                          x_breaks = seq(-4, 4, 0.5)) {
  d_mv = as.data.table(meta_mv$summary)[geneset %in% gs_keep][, model := "multivariate"]
  d_mv$FDR = p.adjust(d_mv$p_value, method = "BH")

  gs_order = d_mv[order(OR, decreasing = FALSE), geneset]
  d_mv$geneset = factor(d_mv$geneset, levels = gs_order)

  if (!is.null(meta_uni)) {
    d_uni = as.data.table(meta_uni$summary)[geneset %in% gs_keep][, model := "univariate"]
    d_uni$FDR = p.adjust(d_uni$p_value, method = "BH")
    d_uni$geneset = factor(d_uni$geneset, levels = gs_order)
  }

  d = d_mv
  d$sig = factor(ifelse(!is.na(d$FDR) & d$FDR < 0.1, "FDR ≤ 0.1", "FDR > 0.1"),
                 levels = c("FDR ≤ 0.1", "FDR > 0.1"))

  sig_colors  = c("FDR ≤ 0.1" = "grey17",     "FDR > 0.1" = "grey78")
  model_fills = c("multivariate" = "grey52", "univariate" = "white")

  ggplot(d, aes(x = log2(OR), y = geneset, group = model)) +
    geom_vline(xintercept = 0, size = 0.2, color = "grey67") +
    geom_errorbarh(aes(xmin = log2(CI_lo), xmax = log2(CI_hi), color = sig),
                   height = 0.6, size = 0.4,
                   position = position_dodgev(height = 0.5)) +
    geom_point(shape = 21, size = 0.8, color = "grey52", fill = "grey52",
               position = position_dodgev(height = 0.5)) +
    scale_color_manual(values = sig_colors, name = "",
                       guide = guide_legend(override.aes = list(
                         shape = NA, linetype = "solid", linewidth = 1))) +
    scale_fill_manual(values = model_fills, name = "",
                      guide = guide_legend(reverse = TRUE, override.aes = list(
                        shape = 21, color = "grey52", size = 2))) +
    scale_x_continuous(breaks = x_breaks) +
    # cap_first() on the drawn gene-set names only; the column they come from is
    # what the rows are filtered and ordered on.
    scale_y_discrete(labels = cap_first) +
    labs(title = title, x = xlab, y = NULL) +
    # grid = "none" blanks both major gridlines, which theme_expresso() otherwise
    # draws dashed grey67.
    #
    # "right" + "top" puts the key outside the plotting area with a layout column
    # of its own. Deliberately not the overlay 06_framingham_survival.R and
    # 07_other.R use (position c(1, 1)): that takes no layout space and so needs a
    # wide right margin to show in, and these two forests are printed with
    # plot.margin 0.1 line on the right, which clips such a key away completely.
    theme_expresso(axis_title_size = textsize, legend_text_size = textsize,
                   grid = "none", legend_position = "right",
                   legend_justification = "top") +
    theme(text         = element_text(color = textcolor),
          axis.text    = element_text(color = textcolor),
          axis.title   = element_text(color = textcolor),
          plot.title   = element_text(color = textcolor),
          legend.text  = element_text(color = textcolor),
          legend.title = element_text(color = textcolor))
}


### cell function factory
# Circle size scales with FDR (smaller FDR = larger circle).
# A grey outline ring appears only when FDR <= fdr_cutoff.
make_cellfun = function(val_mat, fdr_mat, cf, dot_scale = 0.35, fdr_cutoff = 0.1) {
  sz_mat  = 1.2 - fdr_mat
  fd2     = fdr_mat; fd2[fd2 > fdr_cutoff] = 1.2
  sz2_mat = 1.2 - fd2          # 0 when FDR > fdr_cutoff -> outline invisible

  function(j, i, x, y, width, height, fill) {
    # col = NA, so no cell border: the grey53 lwd 0.5 outline this carried drew a
    # grid over all six heatmaps, and 05_sexdifference.R's make_cellfun has never
    # had one. The faint white fill stays, being what lifts a cell off the panel.
    grid.rect(x, y, width, height, gp = gpar(col = NA, fill = "white", alpha = 0.2))
    v = val_mat[i, j]; f = fdr_mat[i, j]
    if (is.na(v) || is.na(f)) return(invisible(NULL))
    r_base = dot_scale * min(unit.c(width, height))
    grid.circle(x, y, r = sz_mat[i, j]  * r_base,
                gp = gpar(fill = cf(v), col = NA, alpha = 0.9))
    # The ring is the only stroke in the cell -- the filled circle under it is
    # col = NA -- and it carried grid's default lwd = 1, which is 0.75 pt and about
    # three times the weight of everything else in the arm. linesize is the arm's
    # line width, so this draws the ring at the same weight as an axis line; the
    # same expression sets the table rules in 06_framingham_survival.R.
    grid.circle(x, y, r = sz2_mat[i, j] * r_base,
                gp = gpar(fill = NA, col = "grey17", alpha = 1,
                          lwd = linesize * ggplot2::.pt))
  }
}


### cohort x gene-set effect-size heatmap (multivariate)
# Turns the p-value matrix into BH FDR for circle size. Defaults match the
# SUBSPACE OR use case (value_field = "OR", p_field = "p_lrt", legend_step = 1);
# pass value_field = "HR" / p_field = "p_value" / legend_step = 0.5 for TCGA HR
# matrices. row_keep subsets rows, e.g. cancer types with neutrophil presence.
make_or_heatmap = function(glm_mv, gs_keep = names(col_celltype),
                           row_keep    = NULL,
                           value_field = "OR",
                           p_field     = "p_lrt",
                           title = "Log2 OR", legend_title = "Log2 OR",
                           cf = cf_or, dot_scale = 0.9, clamp = clamp_or,
                           legend_step = 1) {
  if (!value_field %in% names(glm_mv))
    stop(sprintf("value_field '%s' not found; available: %s",
                 value_field, paste(names(glm_mv), collapse = ", ")))
  if (!p_field %in% names(glm_mv))
    stop(sprintf("p_field '%s' not found; available: %s",
                 p_field, paste(names(glm_mv), collapse = ", ")))

  rows = if (is.null(row_keep)) rownames(glm_mv[[value_field]])
  else intersect(row_keep, rownames(glm_mv[[value_field]]))
  cols = intersect(gs_keep, colnames(glm_mv[[value_field]]))

  v_mat   = log2(glm_mv[[value_field]][rows, cols, drop = FALSE])
  p_mat   = glm_mv[[p_field]][rows, cols, drop = FALSE]
  fdr_mat = matrix(p.adjust(as.vector(p_mat), method = "BH"),
                   nrow = nrow(p_mat), dimnames = dimnames(p_mat))
  v_clust = v_mat; v_clust[!is.finite(v_clust)] = 0

  Heatmap(
    v_clust,
    rect_gp           = gpar(type = "none", fill = "white"),
    col               = cf, na_col = "grey",
    clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
    clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
    cluster_rows = FALSE, cluster_columns = TRUE,
    show_row_dend = FALSE, show_column_dend = TRUE,
    column_dend_side  = "bottom",
    show_column_names = TRUE, column_names_side = "top",
    column_names_rot  = 45,
    column_names_gp   = gpar(fontsize = textsize, col = textcolor),
    column_labels     = cap_first(colnames(v_clust)),
    # Row names stay upper case rather than going through cap_first(): they are
    # cancer-type acronyms on the TCGA matrices and cohort ids on the others, and
    # cap_first() would turn "gse101702" into "Gse101702".
    row_labels        = toupper(rownames(v_clust)),
    row_names_gp      = gpar(fontsize = textsize, col = textcolor),
    row_names_max_width = unit(40, "line"),
    column_title      = title, column_title_side = "top",
    column_title_gp   = gpar(fontsize = textsize, fontface = "plain", col = textcolor),
    show_heatmap_legend = TRUE,
    cell_fun          = make_cellfun(v_mat, fdr_mat, cf, dot_scale = dot_scale),
    heatmap_legend_param = list(
      title = legend_title, title_gp = gpar(fontsize = textsize, col = textcolor),
      labels_gp = gpar(fontsize = textsize, col = textcolor),
      legend_height = unit(3, "line"), grid_width = unit(0.1, "line"),
      direction = "vertical", title_position = "leftcenter-rot",
      at = seq(-clamp, clamp, by = legend_step)
    )
  )
}

dir.clinical.tcga     = paste0(dir.results, "clinical/TCGA/")
dir.clinical.subspace = paste0(dir.results, "clinical/SUBSPACE/")

### matrices
glm_mv = readRDS(paste0(dir.clinical.subspace,
                        "cohort_glm_mv_matrices_PUBLIC_SUBSPACE.rds"))
tcga_mv_mat = readRDS(paste0(dir.clinical.tcga,
                             "pan_cancer_cox_liu2018_mv_matrices.rds"))
neut_rank = readRDS(paste0(dir.clinical.tcga, "neut_rank.rds"))

### colour scales, separate for OR and HR
pal12    = dichromat::colorschemes$DarkRedtoBlue.12
clamp_or = 2
clamp_hr = 1
cf_or    = colorRamp2(seq(-clamp_or, clamp_or, length.out = 12), pal12)
cf_hr    = colorRamp2(seq(-clamp_hr, clamp_hr, length.out = 12), pal12)

### FDR outline legend, two circles
df.twocircles = data.table(
  x = rep(1, 2), y = 1:2, color = c("grey17", "grey88"),
  labeltext = c("FDR ≤ 0.1", "FDR > 0.1"))

p.twocircles = ggplot(df.twocircles, aes(x = x, y = y, label = labeltext, color = color)) +
  geom_point(shape = 21, fill = "grey88", size = 2) +
  geom_text(angle = 0, size = pt2mm(textsize), hjust = -0.2, color = textcolor) +
  theme(text = element_text(size = textsize, color = textcolor)) +
  scale_color_identity() +
  theme_void() +
  coord_cartesian(clip = "off")

### all subsets, SUBSPACE + PUBLIC and TCGA
neut_ct = neut_rank$cancer_type[neut_rank$pct_nonzero > 20]

# The SUBSPACE and PUBLIC objects spell this subset "NF-kB" with an ASCII k, where
# col_celltype and the TCGA matrices use the Greek letter. intersect() and %in% on
# the raw names therefore dropped it from every panel of this figure.
fix_nfkb = function(x) gsub("NF-kB", "NF-κB", x, fixed = TRUE)

blank = "          "
pad_ptgs2 = function(x) gsub("PTGS2", paste0("PTGS2", blank), x, fixed = TRUE)

colnames(glm_mv$OR)           = pad_ptgs2(fix_nfkb(colnames(glm_mv$OR)))
colnames(glm_mv$p_lrt)        = pad_ptgs2(fix_nfkb(colnames(glm_mv$p_lrt)))
colnames(tcga_mv_mat$HR)      = pad_ptgs2(fix_nfkb(colnames(tcga_mv_mat$HR)))
colnames(tcga_mv_mat$p_value) = pad_ptgs2(fix_nfkb(colnames(tcga_mv_mat$p_value)))

# Gene sets shared between SUBSPACE and TCGA, so both heatmaps line up
gs_keep_disp = pad_ptgs2(names(col_celltype))
gs_shared    = Reduce(intersect, list(gs_keep_disp,
                                      colnames(glm_mv$OR),
                                      colnames(tcga_mv_mat$HR)))

ht_subspace = make_or_heatmap(
  glm_mv,
  gs_keep      = gs_shared,
  title        = NULL,
  legend_title = "Log2 odds ratio\ncritical illness\n28-30-day mortality",
  cf           = cf_or,
  dot_scale    = 0.8,
  clamp        = clamp_or
)

ht_tcga = make_or_heatmap(
  tcga_mv_mat,
  gs_keep      = gs_shared,
  row_keep     = neut_ct,
  value_field  = "HR",
  p_field      = "p_value",
  title        = NULL,
  legend_title = "Log2 hazard ratio\nTCGA OS or PFI",
  cf           = cf_hr,
  dot_scale    = 0.4,
  clamp        = clamp_hr,
  legend_step  = 0.5
)

### neutrophil-subset specific signatures
colnames(tcga_mv_mat$HR)      = gsub(".neu", "", colnames(tcga_mv_mat$HR))
colnames(tcga_mv_mat$p_value) = gsub(".neu", "", colnames(tcga_mv_mat$p_value))

gs_neu        = c("immature", "degranulating", "mature", "antiprotease", "total")
gs_neu_shared = Reduce(intersect, list(gs_neu,
                                       colnames(glm_mv$OR),
                                       colnames(tcga_mv_mat$HR)))

ht_subspace2 = make_or_heatmap(
  glm_mv,
  gs_keep      = gs_neu_shared,
  title        = NULL,
  legend_title = "Log2 odds ratio\ncritical illness\n28-30-day mortality",
  cf           = cf_or,
  dot_scale    = 0.65,
  clamp        = clamp_or
)

ht_tcga2 = make_or_heatmap(
  tcga_mv_mat,
  gs_keep      = gs_neu_shared,
  row_keep     = neut_ct,
  value_field  = "HR",
  p_field      = "p_value",
  title        = NULL,
  legend_title = "Log2 hazard ratio\nTCGA OS or PFI",
  cf           = cf_hr,
  dot_scale    = 1.5,
  clamp        = clamp_hr,
  legend_step  = 0.5
)


### severe vs. nonsevere matrices

meta_mv_sev = readRDS(paste0(dir.clinical.subspace,
                             "cohort_meta_re_mv_PUBLIC_SUBSPACE.severity.rds"))
glm_mv_sev  = readRDS(paste0(dir.clinical.subspace,
                             "cohort_glm_mv_matrices_PUBLIC_SUBSPACE.severity.rds"))

colnames(glm_mv_sev$OR)    = pad_ptgs2(fix_nfkb(colnames(glm_mv_sev$OR)))
colnames(glm_mv_sev$p_lrt) = pad_ptgs2(fix_nfkb(colnames(glm_mv_sev$p_lrt)))

meta_mv_sev$summary$geneset = fix_nfkb(as.character(meta_mv_sev$summary$geneset))

gs_shared_sev = intersect(gs_keep_disp, colnames(glm_mv_sev$OR))

### all subsets
p.meta.severity = make_or_forest(
  meta_mv_sev,
  title    = NULL,
  xlab     = "Log2 odds ratio (severe vs. nonsevere)",
  x_breaks = seq(-2, 2, 1)
)

clamp_sev = 2.5
cf_sev    = colorRamp2(seq(-clamp_sev, clamp_sev, length.out = 12), pal12)

ht_severity = make_or_heatmap(
  glm_mv_sev,
  gs_keep      = gs_shared_sev,
  title        = NULL,
  legend_title = "Log2 odds ratio\nsevere vs. non-severe",
  cf           = cf_sev,
  dot_scale    = 0.76,
  clamp        = clamp_sev,
  legend_step  = clamp_sev
)

### neutrophil-subset specific signatures
p.meta.severity2 = make_or_forest(
  meta_mv_sev,
  gs_keep = gs_neu,
  title   = NULL,
  xlab    = "Log2 odds ratio (severe vs. non-severe)"
)

ht_severity2 = make_or_heatmap(
  glm_mv_sev,
  gs_keep      = gs_neu,
  title        = NULL,
  legend_title = "Log2 odds ratio\nsevere vs. non-severe",
  dot_scale    = 0.57
)

### assemble
cairo_pdf(file = paste0(dir.fig, "Fig.e7.pdf"), width = 6, height = 9.44)
pushViewport(viewport(layout = grid.layout(nrow = 160, ncol = 1270)))

### all subsets
pushViewport(viewport(layout.pos.row = 1:34, layout.pos.col = 3:845))
draw(ht_subspace, merge_legend = TRUE, heatmap_legend_side = "right",
     align_heatmap_legend = "heatmap_top", newpage = FALSE)
upViewport(1)

pushViewport(viewport(layout.pos.row = 106:160, layout.pos.col = 3:774))
draw(ht_tcga, merge_legend = TRUE, heatmap_legend_side = "right",
     align_heatmap_legend = "heatmap_top", newpage = FALSE)
upViewport(1)

### neutrophil subsets
pushViewport(viewport(layout.pos.row = 1:34, layout.pos.col = 844:1249))
draw(ht_subspace2, merge_legend = TRUE, heatmap_legend_side = "right",
     align_heatmap_legend = "heatmap_top", newpage = FALSE)
upViewport(1)

pushViewport(viewport(layout.pos.row = 106:160, layout.pos.col = 844:1176))
draw(ht_tcga2, merge_legend = TRUE, heatmap_legend_side = "right",
     align_heatmap_legend = "heatmap_top", newpage = FALSE)
upViewport(1)

### severe vs. nonsevere, forest over heatmap

pushViewport(viewport(layout.pos.row = 71:105, layout.pos.col = 3:845))
draw(ht_severity, merge_legend = TRUE, heatmap_legend_side = "right",
     align_heatmap_legend = "heatmap_top", newpage = FALSE)
upViewport(1)

pushViewport(viewport(layout.pos.row = 71:105, layout.pos.col = 844:1270))
draw(ht_severity2, merge_legend = TRUE, heatmap_legend_side = "right",
     align_heatmap_legend = "heatmap_top", newpage = FALSE)
upViewport(1)

print(p.meta.severity + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "line")),
      vp = viewport(layout.pos.row = 36:72, layout.pos.col = 151:598))

print(p.meta.severity2 + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "line")),
      vp = viewport(layout.pos.row = 50:72, layout.pos.col = 652:1099))


### FDR keys, one per heatmap in the top two rows
print(p.twocircles + theme(plot.margin = margin(0, 0, 0, 0, "line")),
      vp = viewport(layout.pos.row = 25:26, layout.pos.col = 638:818))
print(p.twocircles + theme(plot.margin = margin(0, 0, 0, 0, "line")),
      vp = viewport(layout.pos.row = 25:26, layout.pos.col = 1090:1170))


print(p.twocircles + theme(plot.margin = margin(0, 0, 0, 0, "line")),
      vp = viewport(layout.pos.row = 94:95, layout.pos.col = 636:816))

print(p.twocircles + theme(plot.margin = margin(0, 0, 0, 0, "line")),
      vp = viewport(layout.pos.row = 94:95, layout.pos.col = 1110:1190))

print(p.twocircles + theme(plot.margin = margin(0, 0, 0, 0, "line")),
      vp = viewport(layout.pos.row = 130:131, layout.pos.col = 590:770))
print(p.twocircles + theme(plot.margin = margin(0, 0, 0, 0, "line")),
      vp = viewport(layout.pos.row = 130:131, layout.pos.col = 981:1151))

### panel labels
grid.text(x = unit(0.02, "npc"), y = unit(0.99, "npc"), label = "a",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.02, "npc"), y = unit(0.78, "npc"), label = "b",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.02, "npc"), y = unit(0.54, "npc"), label = "c",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))
grid.text(x = unit(0.02, "npc"), y = unit(0.33, "npc"), label = "d",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", col = textcolor))

dev.off()
