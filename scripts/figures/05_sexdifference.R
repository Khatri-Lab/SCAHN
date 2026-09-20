# ==============================================================================
# Sex differences in the neutrophil compartment
#
# Male vs. female contrasts across four layers: single-cell subset proportion and
# gene expression, purified neutrophil bulk RNA-seq (E-GEAD-397), a whole-blood
# and PBMC bulk transcriptome meta-analysis, and Framingham GLMs of the signature
# scores against age, sex, BMI and smoking. Every bulk contrast is age-adjusted
# and logistic (sex[M = 1] ~ gene + age), so effect sizes are log odds ratios
# with >0 meaning higher in males.
#
# Requires 00_setup.R (packages, figure settings, palettes) and 01_load_data.R
# (dtf, samplemeta, gene_sets, genes.* signature vectors,
# dir.data / dir.results / dir.fig). Also uses margins and car, loaded below.
#
# Reads, on top of what 01_load_data.R loads:
#   results/annotation/SCAHN.subsetproportion.rds
#   results/annotation/SCAHN.subset.markers.csv
#   results/meta_analysis/SCAHN.sex.exprmeta.allgenes.allneu.rds
#   results/Framingham/FH.score_pheno.csv
#   results/sex_difference/ota2021.Neu_LDG.sex.ageAdj.logit.DEG.csv
#   results/sex_difference/dtf.sexMeta.ageAdj.shared.rds
#
# Writes to figures/original/:
#   Fig.s2.pdf    Supplementary figure 2 (GLM coefficients and marginal effects)
#   Figure5.pdf   Figure 5
# ==============================================================================

SCAHN_SCRIPTS = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/scripts/figures"
source(file.path(SCAHN_SCRIPTS, "00_setup.R"))
source(file.path(SCAHN_SCRIPTS, "01_load_data.R"))

gene_sets.gsea = gene_sets
gene_sets.gsea$IFN = c(genes.ifn, "DDX58", "IRF7", "SAMD9L")


# -- Single cell, subset proportion --------------------------------------------

dtf.prop = readRDS(paste0(dir.results, "annotation/SCAHN.subsetproportion.rds"))

### healthy whole-blood samples, with wigerblad2022 restricted by hand
dtf.prop[pid == "wigerblad2022"]$sampleid
dtf.prop.sex.a = dtf.prop[group == "healthy"][
  tissue %in% c("peripheral blood")][
  grepl("WB", source) | grepl("gupta2020|montaldo2022", pid)][
  pid != "wigerblad2022"]
dtf.prop.sex.b = dtf.prop[group == "healthy"][grepl("wigerblad2022", pid)][
  sampleid %in% c("Gran_F1", "Gran_F2", "Gran_F3", "Gran_M1", "Gran_M2",
                  "Gran_M3", "Neut_M0", "WB_oF4", "WB_oF5", "WB_oM4")]
dtf.prop.sex = rbind(dtf.prop.sex.a, dtf.prop.sex.b)
dtf.prop.sex[, .N]
dtf.prop.sex$sex = samplemeta[match(dtf.prop.sex$sampleid.pid, samplemeta$sampleid.pid)]$sex
dtf.prop.sex$age = samplemeta[match(dtf.prop.sex$sampleid.pid, samplemeta$sampleid.pid)]$age
dtf.prop.sex = dtf.prop.sex[!sex %in% c("unknown")]

### proportion meta-analysis, male vs. female
dtf.prop.sex[, class := ifelse(sex %in% "M", 1, NA)]
dtf.prop.sex[, class := ifelse(sex %in% "F", 0, class)]

pids = c("gupta2020", "kaiser2024", "kwok2023", "wigerblad2022", "sinha2021")
features = names(col_celltype)
metaobj = lapply(pids, function(i) {
  sub = dtf.prop.sex[pid == i]
  meta_dataset(
    expr = t(as.matrix(sub[, features, with = FALSE])), class = sub$class, label = i)
})

meta.prop.sex = meta_analysis(metaobj, outcome_type = "binary", var_method = "DeltaMethod")

datadf.sex.summary = as.data.table(meta.prop.sex$summary)
datadf.sex.summary$FDR = p.adjust(datadf.sex.summary$p_value, method = "fdr")

datadf.sex.summary = datadf.sex.summary[gene != "lowDepth"][order(-pooled_ES)]
datadf.sex.summary$celltype = factor(datadf.sex.summary$gene,
                                     levels = datadf.sex.summary$gene)

ci.value = -qnorm((1 - 0.95) / 2)

datadf.sex.summary$lower = datadf.sex.summary$pooled_ES -
                           ci.value * datadf.sex.summary$SE
datadf.sex.summary$upper = datadf.sex.summary$pooled_ES +
                           ci.value * datadf.sex.summary$SE

p.es.summary.sex =
  ggplot(data = datadf.sex.summary[p_value<0.05], aes(x = celltype, y = pooled_ES)) +
  theme_expresso(
    legend_position = "none", legend_justification = c(0, 0),
    legend_key_width = 0.1, legend_key_height = 0.1,
    legend_key_spacing_y = 0.02, panel_background = backgroundcol, plot_background = "grey97") +
  theme(
    plot.title = element_text(face = 1, margin = margin(0, 0, 0, 0, "line")),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 8, face = "plain", color = textcolor)) +
  geom_hline(yintercept = 0, size = 0.3, alpha = 1, color = "grey67") +
  geom_tile(width = 0.7, height = 0.1, position = position_dodge(0.7), alpha = 1, fill = "grey42") +
  geom_errorbar(data = datadf.sex.summary[p_value<0.05],
                aes(x = celltype, ymin = lower, ymax = upper, group = celltype),
                position = position_dodge(0.7), color = "grey7", width = 0.7, size = 0.2) +
  labs(x = NULL, y = bquote("Effect size"), title = "Proportion") +
  scale_y_continuous(breaks = seq(-20, 20, 1)) +
  scale_x_discrete(limits = rev, position = "top") +
  coord_flip() +
  guides(color = "none", shape = "none", size = "none",
         alpha = guide_legend(label.theme = element_text(angle = 90, hjust = 0.5,
                                                             vjust = 0.5), label.position = "top",
                              direction      = "vertical"))


# -- Single cell, gene expression GSEA -----------------------------------------

meta.expr.sex = readRDS(paste0(dir.results, "meta_analysis/SCAHN.sex.exprmeta.allgenes.allneu.rds"))
generanks = meta.expr.sex$summary$pooled_ES

names(generanks) = meta.expr.sex$summary$gene
sum(is.na(meta.expr.sex$summary$pooled_ES))
generanks = generanks[is.finite(generanks)]

### subset marker gene sets, ranked by the pooled sex effect size
markers = fread(paste0(dir.results, "annotation/SCAHN.subset.markers.csv"))
markers.s = markers[avg_log2FC >= 1.6]
subset.genesets = list()
for (i in unique(markers.s$cluster)) {
  tmp = markers.s[cluster == i]
  if (nrow(tmp) > 0) {
    subset.genesets[[i]] = tmp$gene
  }
}

gsea.expr.sex = fgseaMultilevel(pathways = subset.genesets, stats = generanks, nPermSimple = 10000)
gsea.expr.sex = as.data.table(gsea.expr.sex)
gsea.expr.sex = gsea.expr.sex[order(NES)]
gsea.expr.sex$FDR = p.adjust(gsea.expr.sex$pval, method = "fdr")

mat.NES = t(as.matrix(gsea.expr.sex[, "NES"]))
colnames(mat.NES) = gsea.expr.sex$pathway
mat.size = t(as.matrix(gsea.expr.sex[, "FDR"]))
colnames(mat.size) = gsea.expr.sex$pathway
mat.size = 1 - mat.size

mat.size2 = mat.size
mat.size2[mat.size2 < 0.9] = NA

maxNES = max(abs(mat.NES))

col_fun.NES = circlize::colorRamp2(
  seq(-maxNES, maxNES, length = 13), c(dichromat_pal("DarkRedtoBlue.12")(12)[1:6], "white",
    dichromat_pal("DarkRedtoBlue.12")(12)[7:12]))

p.expr.summary.sex =
  Heatmap(mat.NES, rect_gp = gpar(type = "none"), col = col_fun.NES, na_col = "grey93",
          cluster_columns = FALSE, show_row_names = FALSE,
          cluster_rows = FALSE, show_row_dend = TRUE, show_column_dend = TRUE,
          row_dend_side = "left", column_dend_side = "top", clustering_method_rows = "ward.D2",
          clustering_method_columns = "ward.D2", clustering_distance_rows = "euclidean",
          clustering_distance_columns = "euclidean",
          column_names_side = "top", column_names_rot = 90,
          column_names_gp = gpar(fontsize = 8, fontface = "plain"),
          row_names_gp = gpar(fontsize = 8, fontface = "plain"),
          column_title_rot = 0, column_title = NULL, column_title_side = "top",
          column_title_gp = gpar(jfontsize = 9, fontface = "plain"), show_heatmap_legend = FALSE,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.rect(x = x, y = y, width = width, height = height,
                      gp = gpar(lwd = 0.2, col = "grey90", fill = NA, alpha = 0.2))
            grid.circle(x = x, y = y, r = 7 * mat.size[i, j] * min(unit.c(width, height)),
                        gp = gpar(fill = col_fun.NES(mat.NES[i, j]), col = NA, alpha = 0.9))
            grid.circle(x = x, y = y, r = 7 * mat.size2[i, j] * min(unit.c(width, height)),
                        gp = gpar(color = "grey7", fill = NA, alpha = 0.9))
          }, heatmap_legend_param = list(
            title = "NES", title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8),
            legend_width = unit(2, "line"), legend_height = unit(.2, "line"),
            grid_width = unit(2, "line"), grid_height = unit(.1, "line"),
            by_row = FALSE, direction = "horizontal", title_position = "lefttop", at = c(-2, 0, 2)))


# -- Purified neutrophils, E-GEAD-397 ------------------------------------------

res.ota = fread(paste0(dir.results, "sex_difference/ota2021.Neu_LDG.sex.ageAdj.logit.DEG.csv"))
res.ota = res.ota[!is.na(p_value)]
res.ota[, log2_or := log_or / log(2)]

res.gsea = list()

generanks = res.ota$log2_or
names(generanks) = res.ota$gene
set.seed(1)
res.gsea[["sex.ota.ageAdj"]] = fgseaMultilevel(pathways = gene_sets.gsea, stats = generanks,
                                               nPermSimple = 10000)

### volcano
res.ota[, inset := ifelse(
  gene %in% genes.immature.neu, "immature",
  ifelse(gene %in% genes.degranulating.neu, "degranulating",
         ifelse(gene %in% genes.antiprotease.neu, "antiprotease",
                ifelse(gene %in% genes.mature.neu, "mature",
                       ifelse(gene %in% genes.ifn, "IFN", "other")))))]

res.ota = rbind(
  res.ota[inset == "other"], res.ota[inset != "other"])

res.ota[, p_value := ifelse(p_value < 1e-04, 1e-04, p_value)]
res.ota[, label := ifelse(inset != "other", gene, NA)]
res.ota$inset = factor(res.ota$inset, levels = names(col_set))

p.otz.ageAdj = ggplot(res.ota[!is.na(p_value)], aes(x = log2_or, y = -log10(p_value))) +
  theme_expresso(
    legend_position = c(0.35, 0.65), legend_justification = c(0, 0),
    legend_key_width = 0.1, legend_key_height = 0.1,
    legend_key_spacing_x = 0.04, legend_key_spacing_y = 0.01,
    legend_text_size = 8, panel_background = backgroundcol,
    plot_background = "grey97", axis_title_size = 8) +
  theme(
    plot.title = element_text(face = "plain", margin = margin(0, 0, 0, 0, "line"))) +
  labs(x = "Log2 odds ratio", title = expression("Male" ~ italic("vs.") ~ "female (neutrophils)")) +
  geom_point(shape = 16, size = 0.5, alpha = 0.5, aes(color = inset)) +
  geom_text_repel(aes(label = label, color = inset), size = 1.3, alpha = 1,
                  force = 1.5, max.overlaps = 40, force_pull = 1.5, box.padding = 0.02,
                  fontface = "italic", segment.size = 0.03, min.segment.length = 0.3,
                  show.legend = FALSE) +
  scale_color_manual(values = col_set, labels = cap_first) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 6, 1)) +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1)),
         shape = "none", alpha = "none")

### rank bar
res.ota = res.ota[order(-log2_or)]
res.ota$rank = 1:nrow(res.ota)

res.ota[, inset := ifelse(
  gene %in% genes.immature.neu, "immature",
  ifelse(gene %in% genes.degranulating.neu, "degranulating",
         ifelse(gene %in% genes.antiprotease.neu, "antiprotease",
                ifelse(gene %in% genes.mature.neu, "mature",
                       ifelse(gene %in% genes.ifn, "IFN", "other")))))]

res.ota = rbind(
  res.ota[inset == "other"], res.ota[inset != "other"])

p.sex.bar.otz.ageAdj =
  ggplot(res.ota, aes(x = -rank, y = log2_or)) +
  geom_bar(aes(color = inset, fill = inset), stat = "identity", width = 1, linewidth = 0) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = col_set) +
  scale_fill_manual(values = col_set)

### NES table
ota.table = as.data.table(res.gsea[["sex.ota.ageAdj"]])[
  !pathway %in% c("total")][, c("pathway", "NES", "pval")]
ota.table$NES = round(ota.table$NES, 2)
ota.table$pval = formatC(ota.table$pval, digits = 1, format = "e")
ota.table = ota.table[match(names(col_set), pathway)][!is.na(pathway)]

nes_mat = matrix(ota.table$NES, ncol = 1, dimnames = list(ota.table$pathway, "NES"))
col_fun = colorRamp2(c(-3, 0, 3), c("#2A0BD9", "white", "#A60021"))

ht.NES.ota.ageAdj = Heatmap(
  nes_mat, col = col_fun, name = "NES", cluster_rows = FALSE,
  row_names_side = "left", row_labels = cap_first(rownames(nes_mat)), column_title = NULL,
  show_column_names = FALSE, cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", nes_mat[i, j]), x, y, gp = gpar(fontsize = 8))
  }, right_annotation = rowAnnotation(
    pval = anno_text(ota.table$pval, gp = gpar(fontsize = 8)),
    annotation_name_gp = gpar(fontsize = 8), annotation_name_side = "top",
    annotation_name_rot = 0, show_annotation_name = TRUE), row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8), show_heatmap_legend = FALSE)


# -- Bulk transcriptome, age adjusted ------------------------------------------

dtf.sexMeta.ageAdj.shared = readRDS(paste0(
  dir.results, "sex_difference/dtf.sexMeta.ageAdj.shared.rds"))

# Number of studies is now `k.adj.WB` / `k.adj.PBMC` (was numStudies.WB / numStudies.PBMC).
dtf.sexMeta.ageAdj.shared[k.adj.PBMC > 2 & k.adj.WB > 1][, .N]
dtf.sexMeta.ageAdj.shared = dtf.sexMeta.ageAdj.shared[
  k.adj.PBMC > 1 & k.adj.WB > 1][!chr %in% c("X", "Y") & !is.na(chr)]

dtf.sexMeta.ageAdj.shared[, effectSize.WB       := log2(pooled_or.adj.WB)]
dtf.sexMeta.ageAdj.shared[, effectSize.PBMC     := log2(pooled_or.adj.PBMC)]
dtf.sexMeta.ageAdj.shared[, effectSizePval.WB   := pooled_p.adj.WB]
dtf.sexMeta.ageAdj.shared[, effectSizePval.PBMC := pooled_p.adj.PBMC]

maxES = 2.16

dtf.sexMeta.ageAdj.shared[, effectSize.WB2 := ifelse(effectSize.WB > maxES, maxES, effectSize.WB)]
dtf.sexMeta.ageAdj.shared[, effectSize.WB2 := ifelse(effectSize.WB2 < -maxES,
                                                     -maxES, effectSize.WB2)]

dtf.sexMeta.ageAdj.shared[, effectSize.PBMC2 := ifelse(effectSize.PBMC > maxES,
                                                       maxES, effectSize.PBMC)]
dtf.sexMeta.ageAdj.shared[, effectSize.PBMC2 := ifelse(effectSize.PBMC2 < -maxES,
                                                       -maxES, effectSize.PBMC2)]

dtf.sexMeta.ageAdj.shared$effectSizeFDR.WB =
  p.adjust(dtf.sexMeta.ageAdj.shared$effectSizePval.WB, method = "fdr")
dtf.sexMeta.ageAdj.shared$effectSizeFDR.PBMC =
  p.adjust(dtf.sexMeta.ageAdj.shared$effectSizePval.PBMC, method = "fdr")

dtf.sexMeta.ageAdj.shared[, sig := ifelse(
  effectSizeFDR.WB < 0.1 & effectSizeFDR.PBMC < 0.1, "significant.both",
  ifelse(effectSizeFDR.WB < 0.1 & effectSizeFDR.PBMC > 0.1, "significant.WB",
         ifelse(effectSizeFDR.WB > 0.1 & effectSizeFDR.PBMC < 0.1,
                "significant.PBMC", "not significant")))]

dtf.sexMeta.ageAdj.shared$effectSize.diff =
  dtf.sexMeta.ageAdj.shared$effectSize.WB -
  dtf.sexMeta.ageAdj.shared$effectSize.PBMC

dtf.sexMeta.ageAdj.shared[, inset := ifelse(
  gene %in% genes.immature.neu, "immature",
  ifelse(gene %in% genes.degranulating.neu, "degranulating",
         ifelse(gene %in% genes.antiprotease.neu, "antiprotease",
                ifelse(gene %in% genes.mature.neu, "mature",
                       ifelse(gene %in% genes.ifn, "IFN", "other")))))]

dtf.sexMeta.ageAdj.shared[, inset3 := ifelse(!is.na(inset) & inset != "other", "yes", "no")]
dtf.sexMeta.ageAdj.shared[, genelabel := ifelse(inset3 == "yes", gene, NA)]

sigsize = c(
  "significant.both" = 0.3, "significant.WB" = 0.3,
  "significant.PBMC" = 0.3, "not significant" = 0.01)

sigshape = c(
  "significant.both" = 16, "significant.WB" = 15, "significant.PBMC" = 17, "not significant" = 1)

dtf.sexMeta.ageAdj.shared = rbind(
  dtf.sexMeta.ageAdj.shared[sig == "not significant" & is.na(genelabel)],
  dtf.sexMeta.ageAdj.shared[sig != "not significant" & is.na(genelabel)],
  dtf.sexMeta.ageAdj.shared[sig == "not significant" & !is.na(genelabel)],
  dtf.sexMeta.ageAdj.shared[sig != "not significant" & !is.na(genelabel)])

### scatter, WB vs. PBMC
v = prcomp(cbind(dtf.sexMeta.ageAdj.shared$effectSize.WB2,
                 dtf.sexMeta.ageAdj.shared$effectSize.PBMC2))$rotation
beta = v[2, 1] / v[1, 1]

dtf.sexMeta.ageAdj.shared$inset = factor(dtf.sexMeta.ageAdj.shared$inset, levels = names(col_set))

p.sex.PBMC.WB.ageAdj = ggplot(dtf.sexMeta.ageAdj.shared,
                              aes(x = effectSize.WB2, y = effectSize.PBMC2)) +
  theme_expresso(
    legend_position = c(0.1, 0.01), legend_justification = c(0, 0),
    legend_key_width = 0.1, legend_key_height = 0.1,
    legend_key_spacing_x = 0.01, legend_key_spacing_y = 0.01,
    legend_text_size = 8, panel_background = backgroundcol,
    plot_background = "grey97", axis_title_size = 8) +
  theme(
    plot.title = element_text(face = "plain", margin = margin(0, 0, 0, 0, "line")),
    legend.spacing.y = unit(0.04, "line"), legend.spacing.x = unit(-0.02, "line"),
    legend.box = "horizontal", legend.box.just = "bottom") +
  labs(x = "Log2 odds ratio, WB", y = "Log2 odds ratio, PBMC",
       title = expression("Male" ~ italic("vs.") ~ "female (WB & PBMC)")) +
  geom_hline(yintercept = 0, size = 0.2, color = linecolor, alpha = 0.5) +
  geom_vline(xintercept = 0, size = 0.2, color = linecolor, alpha = 0.5) +
  geom_abline(intercept = 0, slope = beta, color = linecolor, size = 0.2, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = linecolor, size = 0.2, alpha = 0.5) +
  geom_point(stroke = 0.2, aes(color = inset, shape = sig, size = sig, alpha = inset3)) +
  geom_text_repel(aes(label = genelabel, color = inset), size = 1.3, alpha = 1,
                  force = 1.5, max.overlaps = 40, force_pull = 1.5, box.padding = 0.02,
                  fontface = "italic", segment.size = 0.03, min.segment.length = 0.3,
                  show.legend = FALSE) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_color_manual(values = col_set, labels = cap_first) +
  scale_shape_manual(values = sigshape, labels = cap_first) +
  scale_alpha_manual(values = c("yes" = 0.9, "no" = 0.3)) +
  scale_size_manual(values = sigsize) +
  stat_cor(color = textcolor, method = "spearman", size = 2, alpha = 0.8) +
  guides(
    shape = guide_legend(order = 1, override.aes = list(size = 1.5)),
    color = guide_legend(order = 2), alpha = "none", size = "none")

### rank bars
dtf.sexMeta.ageAdj.shared = dtf.sexMeta.ageAdj.shared[order(-effectSize.WB)]
dtf.sexMeta.ageAdj.shared$rank.WB = 1:nrow(dtf.sexMeta.ageAdj.shared)
dtf.sexMeta.ageAdj.shared = dtf.sexMeta.ageAdj.shared[order(-effectSize.PBMC)]
dtf.sexMeta.ageAdj.shared$rank.PBMC = 1:nrow(dtf.sexMeta.ageAdj.shared)

dtf.sexMeta.ageAdj.shared[, inset := ifelse(
  gene %in% genes.immature.neu, "immature",
  ifelse(gene %in% genes.degranulating.neu, "degranulating",
         ifelse(gene %in% genes.antiprotease.neu, "antiprotease",
                ifelse(gene %in% genes.mature.neu, "mature",
                       ifelse(gene %in% genes.ifn, "IFN", "other")))))]

dtf.sexMeta.ageAdj.shared = rbind(
  dtf.sexMeta.ageAdj.shared[inset == "other"], dtf.sexMeta.ageAdj.shared[inset != "other"])

generanks = dtf.sexMeta.ageAdj.shared$effectSize.WB
names(generanks) = dtf.sexMeta.ageAdj.shared$gene
set.seed(1)
res.gsea[["sex.WB.ageAdj"]] = fgseaMultilevel(pathways = gene_sets.gsea, stats = generanks,
                                              nPermSimple = 10000)

generanks = dtf.sexMeta.ageAdj.shared$effectSize.PBMC
names(generanks) = dtf.sexMeta.ageAdj.shared$gene
set.seed(1)
res.gsea[["sex.PBMC.ageAdj"]] = fgseaMultilevel(pathways = gene_sets.gsea, stats = generanks,
                                                nPermSimple = 10000)

p.sex.WB.ageAdj = ggplot(dtf.sexMeta.ageAdj.shared, aes(x = -rank.WB, y = effectSize.WB2)) +
  geom_bar(aes(color = inset, fill = inset), stat = "identity", width = 1, linewidth = 0) +
  labs(x = NULL, y = "Log2(OR)\nWB", title = NULL) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = c(col_set)) +
  scale_fill_manual(values = c(col_set))

p.sex.PBMC.ageAdj = ggplot(dtf.sexMeta.ageAdj.shared, aes(x = -rank.PBMC, y = -effectSize.PBMC2)) +
  geom_bar(aes(color = inset, fill = inset), stat = "identity", width = 1, linewidth = 0) +
  labs(x = NULL, y = "Log2(OR)\nPBMC", title = NULL) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = c(col_set)) +
  scale_fill_manual(values = c(col_set)) +
  coord_flip()

### NES tables
WB.table = as.data.table(res.gsea[["sex.WB.ageAdj"]])[
  !pathway %in% c("total")][, c("pathway", "NES", "pval")]
WB.table$NES = round(WB.table$NES, 2)
WB.table$pval = formatC(WB.table$pval, digits = 1, format = "e")
WB.table = WB.table[match(names(col_set), pathway)][!is.na(pathway)]
PBMC.table = as.data.table(res.gsea[["sex.PBMC.ageAdj"]])[
  !pathway %in% c("total")][, c("pathway", "NES", "pval")]
PBMC.table$NES = round(PBMC.table$NES, 2)
PBMC.table$pval = formatC(PBMC.table$pval, digits = 1, format = "e")
PBMC.table = PBMC.table[match(names(col_set), pathway)][!is.na(pathway)]

nes_mat.WB = matrix(WB.table$NES, ncol = 1, dimnames = list(WB.table$pathway, "NES"))
nes_mat.PBMC = matrix(PBMC.table$NES, ncol = 1, dimnames = list(PBMC.table$pathway, "NES"))

ht.NES.WB.ageAdj = Heatmap(
  nes_mat.WB, col = col_fun, name = "NES", cluster_rows = FALSE,
  row_names_side = "left", column_title = NULL, show_column_names = FALSE, show_row_names = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", nes_mat.WB[i, j]), x, y, gp = gpar(fontsize = 8))
  }, right_annotation = rowAnnotation(
    pval = anno_text(WB.table$pval, gp = gpar(fontsize = 8)),
    annotation_name_gp = gpar(fontsize = 8), annotation_name_side = "top",
    annotation_name_rot = 0, show_annotation_name = TRUE), row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8), show_heatmap_legend = FALSE)

ht.NES.PBMC.ageAdj = Heatmap(
  nes_mat.PBMC, col = col_fun, name = "NES", cluster_rows = FALSE,
  row_names_side = "left", row_labels = cap_first(rownames(nes_mat.PBMC)),
  column_title = NULL,
  show_column_names = FALSE, cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", nes_mat.PBMC[i, j]), x, y, gp = gpar(fontsize = 8))
  }, right_annotation = rowAnnotation(
    pval = anno_text(PBMC.table$pval, gp = gpar(fontsize = 8)),
    annotation_name_gp = gpar(fontsize = 8), annotation_name_side = "top",
    annotation_name_rot = 0, show_annotation_name = TRUE), row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8), show_heatmap_legend = FALSE)

### hypergeometric test: immature genes among the top 1% by effect size
genes2test = genes.immature.neu
N = length(unique(dtf.sexMeta.ageAdj.shared$gene))  # Total number of unique genes in List B
K = nrow(dtf.sexMeta.ageAdj.shared) * 0.01  # Number of top genes to consider in List B
M = length(genes2test)  # Number of genes in List A

# Calculate the overlap
top_K_genes = dtf.sexMeta.ageAdj.shared[order(-effectSize.WB)]$gene[1:K]
k = length(intersect(genes2test, top_K_genes))
k
phyper(k - 1, K, N - K, M, lower.tail = FALSE)

top_K_genes = dtf.sexMeta.ageAdj.shared[order(-effectSize.PBMC)]$gene[1:K]
k = length(intersect(genes2test, top_K_genes))
k
phyper(k - 1, K, N - K, M, lower.tail = FALSE)


# -- Framingham, scores vs. age and BMI ----------------------------------------

framingham_score_pheno = fread(paste0(
  dir.results, "Framingham/FH.score_pheno.csv"))
framingham_score_pheno$immature_mature_ratio =
  framingham_score_pheno$immature / framingham_score_pheno$mature

vars = c("immature", "degranulating", "antiprotease", "mature", "total",
         "immature_mature_ratio", "IFN")

df.fh.hc = melt(framingham_score_pheno[healthy == 1],
                id.vars = c("age", "sex", "bmi", "regular_smokers"), measure.vars = vars)
df.fh.hc = as.data.table(df.fh.hc)
df.fh.hc[, bmi_categ := ifelse(bmi > 30, "high", "low")]
df.fh.hc$variable = factor(df.fh.hc$variable, levels = vars)

### one panel per score, against age
p.fh.hc.list.age = list()
for (score in vars) {
  tmp = df.fh.hc[variable == score]
  tmp$facetvar = "Age"
  plottitle = score
  p.fh.hc.list.age[[score]] =
    ggplot(tmp, aes(x = age, y = value, color = sex)) +
    geom_point(size = 0.1, alpha = 0.6, shape = 16) +
    geom_smooth(method = "loess", show.legend = FALSE, size = 0.5) +
    facet_wrap(. ~ facetvar, ncol = 1, scales = "free") +
    theme_expresso(axis_text_size = 9, legend_text_size = 8,
                   plot_background = "grey97", legend_position = c(-0.14, -0.22),
                   legend_justification = c(0.5, 0.5), legend_key_width = 0.1,
                   legend_key_height = 0.1, legend_key_spacing_x = 0.01,
                   legend_key_spacing_y = 0.01) +
    labs(y = "Score", x = NULL, title = cap_first(plottitle)) +
    theme(axis.text       = element_text(margin = margin(0, 0, 0, 0, "line")),
          axis.title      = element_text(margin = margin(0, 0, 0, 0, "line")),
          strip.text.x = element_text(margin = margin(0.02, 0.02, 0.2, 0.02, "line")),
          strip.text.y = element_text(margin = margin(0.1, 0.1, 0.1, 0.1, "line")),
          plot.title      = element_text(margin = margin(0, 0, 0, 0, "line")),
          legend.box = element_blank(), legend.margin = margin(2.75, 2.75, 2.75, 2.75, "pt"),
          legend.key.size = unit(1, "lines")) +
    scale_y_continuous(breaks = seq(0, 10, 1)) +
    scale_color_manual(values = c("F" = "#ffab65", "M" = "#65b9ff")) +
    guides(color = guide_legend(override.aes = list(size = 2), nrow = 1))
}

### the same scores against BMI, y axis suppressed so they sit beside the age panels
p.fh.hc.list.bmi = list()
for (score in vars) {
  tmp = df.fh.hc[variable == score]
  tmp$facetvar = "BMI"
  p.fh.hc.list.bmi[[score]] =
    ggplot(tmp, aes(x = bmi, y = value, color = sex)) +
    geom_point(size = 0.1, alpha = 0.6, shape = 16) +
    geom_smooth(method = "loess", show.legend = FALSE, size = 0.5) +
    facet_wrap(. ~ facetvar, ncol = 1, scales = "free") +
    theme_expresso(axis_text_size = 9, plot_background = "grey97",
                   legend_position = "none", legend_key_spacing_x = 0.01,
                   legend_key_spacing_y = 0.01) +
    labs(y = NULL, x = NULL, title = " ") +
    theme(axis.text       = element_text(margin = margin(0, 0, 0, 0, "line")),
          axis.title      = element_text(margin = margin(0, 0, 0, 0, "line")),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          strip.text.x = element_text(margin = margin(0.02, 0.02, 0.2, 0.02, "line")),
          strip.text.y = element_text(margin = margin(0.1, 0.1, 0.1, 0.1, "line")),
          plot.title      = element_text(margin = margin(0, 0, 0, 0, "line")),
          legend.box = element_blank(), legend.key.size = unit(1, "lines")) +
    scale_y_continuous(breaks = seq(0, 10, 1)) +
    scale_color_manual(values = c("F" = "#ffab65", "M" = "#65b9ff")) +
    guides(color = guide_legend(override.aes = list(size = 2)))
}


# -- Framingham, GLM helpers ---------------------------------------------------

fit_glm_table = function(formula, data, rowlabels, stratify_by = NULL) {
  fit = glm(formula, data = data, family = gaussian())

  # Raw numeric values
  raw_glm = summary(fit)$coefficients[-1, c(1, 4), drop = FALSE]
  colnames(raw_glm) = c("coef", "pval")
  rownames(raw_glm) = rowlabels

  raw_mg    = as.data.frame(summary(margins(fit)))
  raw_mg_tb = raw_mg[, c("factor", "AME", "p")]
  rownames(raw_mg_tb) = raw_mg_tb$factor
  raw_mg_tb$factor = NULL
  colnames(raw_mg_tb) = c("AME", "pval")

  # VIF
  vif_vals = tryCatch(car::vif(fit, type = "predictor"), error = function(e) NULL)

  # Formatted tables
  tb = as.data.frame(raw_glm)
  colnames(tb) = c("regression\ncoefficient", "p value")
  tb[] = lapply(tb, formatC, digits = 1, format = "e")

  mg_tb = as.data.frame(raw_mg_tb)
  colnames(mg_tb) = c("AME", "p(AME)")
  mg_tb[] = lapply(mg_tb, formatC, digits = 1, format = "e")

  # Stratified margins
  strat_raw = NULL
  strat_tb  = NULL
  if (!is.null(stratify_by)) {
    strat_vars   = setdiff(all.vars(formula)[-1], stratify_by)
    strat_levels = sort(unique(data[[stratify_by]]))
    strat_list = lapply(strat_levels, function(lv) {
      mg_s = as.data.frame(summary(margins(
        fit, data = as.data.frame(data[data[[stratify_by]] == lv, ]), variables = strat_vars)))
      raw = mg_s[, c("factor", "AME", "p")]
      rownames(raw) = raw$factor
      raw$factor = NULL
      colnames(raw) = c("AME", "pval")
      raw
    })
    names(strat_list) = strat_levels
    strat_raw = strat_list
    strat_tb = do.call(rbind, lapply(strat_levels, function(lv) {
      df = strat_list[[lv]]
      df[] = lapply(df, formatC, digits = 1, format = "e")
      colnames(df) = c("AME", "p(AME)")
      cbind(group = lv, factor = rownames(df), df)
    }))
    rownames(strat_tb) = NULL
  }

  list(glm = tb,       margins = mg_tb, raw_glm = raw_glm,  raw_margins = raw_mg_tb,
       vif = vif_vals, stratified = strat_tb, raw_stratified = strat_raw)
}

# Circle size scales with FDR (smaller FDR = larger circle).
# A grey outline ring appears only when FDR <= fdr_cutoff.
make_cellfun = function(val_mat, fdr_mat, cf, dot_scale = 0.35, fdr_cutoff = 0.1, text_size = 6) {
  sz_mat  = 1.2 - fdr_mat
  fd2     = fdr_mat; fd2[fd2 > fdr_cutoff] = 1.2
  sz2_mat = 1.2 - fd2

  function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height,
              gp = gpar(lwd = 0.5, col = "grey53", fill = "white", alpha = 0.2))
    v = val_mat[i, j]; f = fdr_mat[i, j]
    if (is.na(v) || is.na(f)) return(invisible(NULL))
    r_base = dot_scale * min(unit.c(width, height))
    grid.circle(x, y, r = sz_mat[i, j] * r_base, gp = gpar(fill = cf(v), col = NA, alpha = 0.9))
    grid.circle(x, y, r = sz2_mat[i, j] * r_base, gp = gpar(fill = NA, col = "grey17", alpha = 1))
    grid.text(sprintf("%.2f", v), x, y, gp = gpar(fontsize = text_size,
                        col = ifelse(abs(v) > 0.15, "white", "grey20")))
  }
}


# -- Framingham, GLM fits ------------------------------------------------------

base_data  = framingham_score_pheno[healthy == 1][!is.na(bmi)]
base_rows  = c("age", "sex", "BMI", "smoking")
score_cols = c("immature", "mature", "degranulating", "antiprotease", "total",
               "IFN", "immature_mature_ratio")

scores = c("immature", "mature", "degranulating", "antiprotease", "total",
           "IFN", "immature_mature_ratio", "age", "bmi")
scores_z = paste0(scores, "_z")

base_data[, (scores_z) := lapply(.SD, function(x) as.numeric(scale(x))), .SDcols = scores]

library(margins)
glm_results = setNames(lapply(scores_z[1:7], function(sc) {
  formula = as.formula(paste(
    sc, "~ age_z + sex + bmi_z + regular_smokers + age_z*sex + age_z*bmi_z + bmi_z*sex"))
  fit_glm_table(formula, base_data, c(base_rows, "age:sex", "age:BMI", "sex:BMI"),
                stratify_by = "sex")
}), scores_z[1:7])


# -- Framingham, AME heatmaps for figure 5 -------------------------------------

ht_opt$TITLE_PADDING = unit(c(0.2, 0.2), "line")

### AME matrices (main effects only)
ame_rows = rownames(glm_results[[1]]$raw_margins)
mat_ame  = sapply(glm_results, function(r) r$raw_margins[ame_rows, "AME"])
mat_ame_pval = sapply(glm_results, function(r) r$raw_margins[ame_rows, "pval"])

colnames(mat_ame) = gsub("_z", "", colnames(mat_ame))
colnames(mat_ame_pval) = gsub("_z", "", colnames(mat_ame_pval))
rownames(mat_ame) = c("age", "BMI", "smoking", "sex")
rownames(mat_ame_pval) = c("age", "BMI", "smoking", "sex")
mat_ame = mat_ame[base_rows, score_cols]
mat_ame_pval = mat_ame_pval[base_rows, score_cols]
colnames(mat_ame) = gsub("immature_mature_ratio", "immature/mature", colnames(mat_ame))
colnames(mat_ame_pval) = gsub("immature_mature_ratio", "immature/mature", colnames(mat_ame_pval))

### stratified AME matrices per level
strat_levels = names(glm_results[[1]]$raw_stratified)
strat_rows   = rownames(glm_results[[1]]$raw_stratified[[1]])

mat_strat_ame  = lapply(strat_levels, function(lv)
  sapply(glm_results, function(r) r$raw_stratified[[lv]][strat_rows, "AME"]))
names(mat_strat_ame) = strat_levels
mat_strat_pval = lapply(strat_levels, function(lv)
  sapply(glm_results, function(r) r$raw_stratified[[lv]][strat_rows, "pval"]))
names(mat_strat_pval) = strat_levels

colnames(mat_strat_ame$F) = gsub("_z", "", colnames(mat_strat_ame$F))
colnames(mat_strat_pval$F) = gsub("_z", "", colnames(mat_strat_pval$F))
rownames(mat_strat_ame$F) = c("age", "BMI", "smoking")
rownames(mat_strat_pval$F) = c("age", "BMI", "smoking")
mat_strat_ame$F = mat_strat_ame$F[c("age", "BMI"), score_cols]
mat_strat_pval$F = mat_strat_pval$F[c("age", "BMI"), score_cols]
colnames(mat_strat_ame$F) = gsub("immature_mature_ratio", "immature/mature",
                                 colnames(mat_strat_ame$F))
colnames(mat_strat_pval$F) = gsub("immature_mature_ratio", "immature/mature",
                                  colnames(mat_strat_pval$F))

colnames(mat_strat_ame$M) = gsub("_z", "", colnames(mat_strat_ame$M))
colnames(mat_strat_pval$M) = gsub("_z", "", colnames(mat_strat_pval$M))
rownames(mat_strat_ame$M) = c("age", "BMI", "smoking")
rownames(mat_strat_pval$M) = c("age", "BMI", "smoking")
mat_strat_ame$M = mat_strat_ame$M[c("age", "BMI"), score_cols]
mat_strat_pval$M = mat_strat_pval$M[c("age", "BMI"), score_cols]
colnames(mat_strat_ame$M) = gsub("immature_mature_ratio", "immature/mature",
                                 colnames(mat_strat_ame$M))
colnames(mat_strat_pval$M) = gsub("immature_mature_ratio", "immature/mature",
                                  colnames(mat_strat_pval$M))

### both sexes
# FDR correction by row (each predictor corrected across scores)
mat_ame_fdr = apply(mat_ame_pval, 1, p.adjust, method = "BH")
mat_ame = t(mat_ame)

# Color scale centered at 0
clamp_ame = 0.3
cf_ame = colorRamp2(c(-clamp_ame, 0, clamp_ame), c("#2A0BD9", "white", "#A60021"))

# Cell function: dot size ~ FDR
ht_ame = Heatmap(
  mat_ame[1:6, ], rect_gp = gpar(type = "none"), col = cf_ame, na_col = "grey",
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_column_names = TRUE, column_names_side = "bottom", column_names_rot = 90,
  column_names_gp   = gpar(fontsize = 9, col = textcolor),
  column_labels = cap_first(colnames(mat_ame)),
  row_labels = cap_first(rownames(mat_ame)[1:6]),
  row_names_gp = gpar(fontsize = 9, col = textcolor), name = "both sexes",
  column_title = " ", column_title_gp = gpar(fontsize = 9, fontface = "plain", fill = "grey87",
                         col = NA), show_heatmap_legend = FALSE, row_names_side = "left",
  cell_fun          = make_cellfun(mat_ame[1:6, ], mat_ame_fdr[1:6, ], cf_ame,
                                   dot_scale = 0.45, fdr_cutoff = 0.05, text_size = 5),
  heatmap_legend_param = list(
    title = "AME", title_gp = gpar(fontsize = 9, col = textcolor),
    labels_gp = gpar(fontsize = 9, col = textcolor), direction = "vertical",
    title_position = "leftcenter-rot", at = seq(-clamp_ame, clamp_ame, length.out = 5)))

### sex stratified, female
mat_ame_fdr.F = apply(mat_strat_pval$F, 1, p.adjust, method = "BH")
mat_ame.F = t(mat_strat_ame$F)

# Cell function: dot size ~ FDR
ht_ame.F = Heatmap(
  mat_ame.F[1:6, ], rect_gp = gpar(type = "none", fill = "white"), col = cf_ame, na_col = "grey",
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_column_names = TRUE, column_names_side = "bottom", column_names_rot = 90,
  column_names_gp   = gpar(fontsize = 9, col = textcolor),
  column_labels = cap_first(colnames(mat_ame.F)),
  row_names_gp = gpar(fontsize = 9, col = textcolor), name = "female",
  column_title = " ", column_title_gp = gpar(fontsize = 9, fontface = "plain", fill = "grey87",
                         col = NA), show_heatmap_legend = FALSE,
  show_row_names = FALSE, cell_fun = make_cellfun(mat_ame.F[1:6, ], mat_ame_fdr.F[1:6, ],
                                   cf_ame, dot_scale = 0.9, fdr_cutoff = 0.05, text_size = 5),
  heatmap_legend_param = list(
    title = "AME", title_gp = gpar(fontsize = 9, col = textcolor),
    labels_gp = gpar(fontsize = 9, col = textcolor), direction = "vertical",
    title_position = "leftcenter-rot", at = seq(-clamp_ame, clamp_ame, length.out = 5)))

### sex stratified, male
mat_ame_fdr.M = apply(mat_strat_pval$M, 1, p.adjust, method = "BH")
mat_ame.M = t(mat_strat_ame$M)

# Cell function: dot size ~ FDR
ht_ame.M = Heatmap(
  mat_ame.M[1:6, ], rect_gp = gpar(type = "none", fill = "white"), col = cf_ame, na_col = "grey",
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_column_names = TRUE, column_names_side = "bottom", column_names_rot = 90,
  column_names_gp   = gpar(fontsize = 9, col = textcolor),
  column_labels = cap_first(colnames(mat_ame.M)),
  row_names_gp = gpar(fontsize = 9, col = textcolor), name = "male",
  column_title = " ", column_title_gp = gpar(fontsize = 9, fontface = "plain", fill = "grey87",
                         col = NA), show_heatmap_legend = FALSE,
  show_row_names = FALSE, cell_fun = make_cellfun(mat_ame.M[1:6, ], mat_ame_fdr.M[1:6, ],
                                   cf_ame, dot_scale = 0.9, fdr_cutoff = 0.05, text_size = 5),
  heatmap_legend_param = list(
    title = "AME", title_gp = gpar(fontsize = 9, col = textcolor),
    labels_gp = gpar(fontsize = 9, col = textcolor), direction = "vertical",
    title_position = "leftcenter-rot", at = seq(-clamp_ame, clamp_ame, length.out = 5)))


# -- Framingham, per-score heatmaps for the supplement -------------------------

col_fun = colorRamp2(c(-0.602, 0, 0.602), c("#2A0BD9", "white", "#A60021"))

### regression coefficients, one narrow heatmap per score
ht_coef_list = lapply(seq_along(scores_z[1:7]), function(k) {
  sc  = scores_z[1:7][k]
  res = glm_results[[sc]]
  coef_mat = matrix(as.numeric(res$raw_glm[, "coef"]), ncol = 1,
                    dimnames = list(rownames(res$glm), sc))
  pval_vec = res$raw_glm[, "pval"]
  anno_name = paste0("pval_", k)
  ht_name   = paste0("ht_", k)

  Heatmap(
    coef_mat, rect_gp = gpar(col = "grey17", lwd = .2), col = col_fun, name = ht_name,
    width = unit(0.7, "cm"), cluster_rows = FALSE,
    row_labels = cap_first(rownames(coef_mat)),
    row_names_side = "left", show_row_names = (k == 1),
    column_title = " ", column_title_gp = gpar(fontsize = 8, col = NA),
    show_column_names = FALSE, cell_fun = function(j, i, x, y, width, height, fill) {
      p  = pval_vec[i]
      bg = if (!is.na(p) && p > 0.05) "white" else fill
      grid.rect(x, y, width, height, gp = gpar(fill = bg, col = "grey17"))
      grid.text(sprintf("%.2f", coef_mat[i, j]), x, y, gp = gpar(fontsize = 8))
    }, right_annotation = do.call(rowAnnotation, c(
      setNames(list(row_anno_text(
        res$glm$`p value`, width = max_text_width(res$glm$`p value`, gp = gpar(fontsize = 8)),
        gp = gpar(fontsize = 8))), anno_name), list(annotation_name_gp = gpar(fontsize = 8),
           annotation_name_side = "top", annotation_name_rot = 0, show_annotation_name = FALSE))),
    row_names_gp = gpar(fontsize = 8), show_heatmap_legend = FALSE)
})
ht_coef_all = Reduce(`+`, ht_coef_list)

### average marginal effects, same layout
row_map = c(age_z = "age", age_cs = "age", sexM = "sex",
            bmi_z = "BMI", bmi_cs = "BMI", regular_smokers = "smoking", smoking = "smoking")
row_order = c("age", "sex", "BMI", "smoking")

ht_ame_list = lapply(seq_along(scores_z[1:7]), function(k) {
  sc  = scores_z[1:7][k]
  res = glm_results[[sc]]

  ame_raw  = res$raw_margins
  old_rows = rownames(ame_raw)
  new_rows = row_map[old_rows]

  ame_mat = matrix(as.numeric(ame_raw[, "AME"]), ncol = 1, dimnames = list(new_rows, sc))
  pval_vec = setNames(as.numeric(ame_raw[, "pval"]), new_rows)

  # Keep only recognised rows, then enforce order
  ame_mat  = ame_mat[row_order[row_order %in% rownames(ame_mat)], , drop = FALSE]
  pval_vec = pval_vec[rownames(ame_mat)]

  # Pad any missing rows with NA so dimensions stay consistent
  missing  = setdiff(row_order, rownames(ame_mat))
  if (length(missing)) {
    pad_mat = matrix(NA_real_, nrow = length(missing), ncol = 1, dimnames = list(missing, sc))
    ame_mat  = rbind(ame_mat, pad_mat)[row_order, , drop = FALSE]
    pval_vec = c(pval_vec, setNames(rep(NA_real_, length(missing)), missing))[row_order]
  }

  pval_fmt = ifelse(is.na(pval_vec), "NA", formatC(pval_vec, digits = 1, format = "e"))
  anno_name = paste0("pval_ame_", k)
  ht_name   = paste0("ht_ame_", k)

  Heatmap(
    ame_mat, rect_gp = gpar(col = "grey17", lwd = .2), col = col_fun, name = ht_name,
    width = unit(0.7, "cm"), cluster_rows = FALSE,
    row_labels = cap_first(rownames(ame_mat)),
    row_names_side = "left", show_row_names = (k == 1),
    column_title = " ", column_title_gp = gpar(fontsize = 8, col = NA),
    show_column_names = FALSE, cell_fun = function(j, i, x, y, width, height, fill) {
      p  = pval_vec[i]
      bg = if (!is.na(p) && p > 0.05) "white" else fill
      grid.rect(x, y, width, height, gp = gpar(fill = bg, col = "grey17"))
      grid.text(sprintf("%.2f", ame_mat[i, j]), x, y, gp = gpar(fontsize = 8,
                          col = if (is.na(ame_mat[i, j])) NA else "black"))
    }, right_annotation = do.call(rowAnnotation, c(
      setNames(list(row_anno_text(
        pval_fmt, width = max_text_width(pval_fmt, gp = gpar(fontsize = 8)),
        gp = gpar(fontsize = 8))), anno_name), list(annotation_name_gp = gpar(fontsize = 8),
           annotation_name_side = "top", annotation_name_rot = 0, show_annotation_name = FALSE))),
    row_names_gp = gpar(fontsize = 8), show_heatmap_legend = FALSE)
})
ht_ame_all = Reduce(`+`, ht_ame_list)


# -- Supplementary figure 2 ----------------------------------------------------

cairo_pdf(file = paste0(dir.fig, "Fig.s2.pdf"), width = 6.4, height = 4.8)
pushViewport(viewport(layout = grid.layout(nrow = 100, ncol = 100)))

### regression coefficients (top half)
pushViewport(viewport(layout.pos.row = 15:55, layout.pos.col = 6:100))
draw(ht_coef_all, newpage = FALSE)
for (k in seq_along(scores_z[1:7])) {
  sc_display = gsub("_z", "", scores_z[1:7][k])
  sc_display = gsub("immature_mature_ratio", "immature/mature", sc_display)
  decorate_column_title(paste0("ht_", k), {
    grid.text(cap_first(sc_display), x = unit(0.5, "npc"), y = unit(0.1, "npc"),
              rot = 90, gp = gpar(fontsize = 8, fontface = "plain"), just = c("left", "center"))
  })
  decorate_annotation(paste0("pval_", k), {
    grid.text("pval", x = unit(0.5, "npc"), y = unit(1, "npc") + unit(0.2, "line"),
              rot = 90, gp = gpar(fontsize = 8), just = c("left", "center"))
  })
}
popViewport()

### AME (bottom half)
pushViewport(viewport(layout.pos.row = 73:98, layout.pos.col = 6:100))
draw(ht_ame_all, newpage = FALSE)
for (k in seq_along(scores_z[1:7])) {
  sc_display = gsub("_z", "", scores_z[1:7][k])
  sc_display = gsub("immature_mature_ratio", "immature/mature", sc_display)
  decorate_column_title(paste0("ht_ame_", k), {
    grid.text(cap_first(sc_display), x = unit(0.5, "npc"), y = unit(0.1, "npc"),
              rot = 90, gp = gpar(fontsize = 8, fontface = "plain"), just = c("left", "center"))
  })
  decorate_annotation(paste0("pval_ame_", k), {
    grid.text("pval", x = unit(0.5, "npc"), y = unit(1, "npc") + unit(0.2, "line"),
              rot = 90, gp = gpar(fontsize = 8), just = c("left", "center"))
  })
}
popViewport()

### panel labels and titles
grid.text(x = unit(0.018, "npc"), y = unit(0.98, "npc"), label = "a",
          gp = gpar(fontsize = 8, fontface = "bold"))
grid.text(x = unit(0.018, "npc"), y = unit(0.4, "npc"), label = "b",
          gp = gpar(fontsize = 8, fontface = "bold"))

grid.text(x = unit(0.05, "npc"), y = unit(0.98, "npc"),
          label = "GLM, regression coefficient, Framingham cohort study",
          gp = gpar(fontsize = 8, fontface = "plain"), just = "left")
grid.text(x = unit(0.05, "npc"), y = unit(0.4, "npc"),
          label = "GLM, average marginal effects, Framingham cohort study",
          gp = gpar(fontsize = 8, fontface = "plain"), just = "left")

dev.off()


# -- Figure 5, sample tables and legend helpers --------------------------------

### single cell sample table
dt.sc = data.table(
  "datasets" = 5, "females" = 17, "males" = 16)

dt.sc = tableGrob(
  t(dt.sc), theme = ttheme_default(
    core = list(
      fg_params = list(fontsize = 10, fontface = "plain", lineheight = 0.8),
      bg_params = list(fill = c("grey97")), padding = unit(c(0.3, 0.3), "line")), colhead = list(
      fg_params = list(fontsize = 10, fontface = "plain"), padding = unit(c(0.3, 0.3), "line")
    ), rowhead = list(
      fg_params = list(fontsize = 10, fontface = "plain"), padding = unit(c(0.3, 0.3), "line"))))

### E-GEAD-397 sample table
dt.ota = data.table(
  " " = c("females", "males"), "neutrophils" = c(58, 20))

dt.ota = tableGrob(
  t(dt.ota), theme = ttheme_default(
    core = list(
      fg_params = list(fontsize = 10, fontface = "plain", lineheight = 0.8),
      bg_params = list(fill = c("grey97")), padding = unit(c(0.3, 0.3), "line")), colhead = list(
      fg_params = list(fontsize = 10, fontface = "plain"), padding = unit(c(0.3, 0.3), "line")
    ), rowhead = list(
      fg_params = list(fontsize = 10, fontface = "plain"), padding = unit(c(0.3, 0.3), "line"))))

### bulk meta-analysis sample table
dt.meta = data.table(
  " " = c("datasets", "females", "males"), "WB" = c("10", "337", "279"),
  "PBMC" = c("8", "204", "157"))

dt.meta = tableGrob(
  t(dt.meta), theme = ttheme_default(
    core = list(
      fg_params = list(fontsize = 10, fontface = "plain", lineheight = 0.8),
      bg_params = list(fill = c("grey97")), padding = unit(c(0.3, 0.3), "line")), colhead = list(
      fg_params = list(fontsize = 10, fontface = "plain"), padding = unit(c(0.3, 0.3), "line")
    ), rowhead = list(
      fg_params = list(fontsize = 10, fontface = "plain"), padding = unit(c(0.3, 0.3), "line"))))

### FDR size legend, three circles
circle_data = data.frame(
  x = c(2, 3.8, 6), y = c(2, 2, 2), xtext = c(2, 4, 8.5), radius = c(0.4, 0.65, 0.9),
  label = c("", "", " low"))

pthreecircle = ggplot(data = circle_data) +
  geom_circle(aes(x0 = x, y0 = y, r = radius), linewidth = 0.1, size = 0.1,
              fill = "grey83", color = "grey83") +
  geom_text(aes(x = xtext, y = y, label = label), size = 3, color = textcolor) +
  theme_void() +
  labs(y = "FDR high", title = NULL) +
  scale_x_continuous(limits = c(1.6, 10)) +
  theme(text = element_text(size = 8, color = textcolor),
        plot.title = element_text(size = 8, hjust = 0.5), axis.title.y = element_text(size = 8))

### FDR outline legend, two circles
df.twocircles = data.table(
  x = rep(1, 2), y = 1:2, color = c("grey88", "grey17"), labeltext = c("FDR > 0.05", "FDR ≤ 0.05"))

p.twocircles = ggplot(df.twocircles, aes(x = x, y = y, label = labeltext, color = color)) +
  geom_point(shape = 21, fill = "grey88", size = 2) +
  geom_text(angle = 0, size = 2.6, hjust = -0.2, color = textcolor) +
  theme(text = element_text(size = textsize, color = textcolor)) +
  scale_color_identity() +
  theme_void() +
  coord_cartesian(clip = "off")

### female / male direction arrow
p.arrow = ggplot() +
  geom_segment(aes(x = 0.1, xend = 0.9, y = 0.5, yend = 0.5), size = .2,
               arrow = arrow(ends = "both", length = unit(.4, "line"), type = "closed")) +
  annotate("text", x = 0.5, y = 0.47, size = 3,
           label = "Enriched in females       Enriched in males") +
  scale_y_continuous(limits = c(0.45, 0.52)) +
  theme_void()


# -- Figure 5 ------------------------------------------------------------------

cairo_pdf(file = paste0(dir.fig, "Figure5.pdf"), width = 8.8, height = 6.5)
pushViewport(viewport(layout = grid.layout(nrow = 118, ncol = 361)))

### panel backgrounds
print(ggplot() + theme_void() +
        theme(plot.background = element_rect(color = "grey97", fill = "grey97")),
      vp = viewport(layout.pos.row = 2:32, layout.pos.col = 3:208))
print(ggplot() + theme_void() +
        theme(plot.background = element_rect(color = "grey97", fill = "grey97")),
      vp = viewport(layout.pos.row = 34:117, layout.pos.col = 3:95))
print(ggplot() + theme_void() +
        theme(plot.background = element_rect(color = "grey97", fill = "grey97"),
              plot.margin = margin(0.2, 0.1, 1, 0.4, "line")),
      vp = viewport(layout.pos.row = 34:117, layout.pos.col = 98:208))

print(ggplot() + theme_void() +
        theme(plot.background = element_rect(color = "grey97", fill = "grey97")),
      vp = viewport(layout.pos.row = 2:117, layout.pos.col = 211:360))

### single cell
pushViewport(viewport(layout.pos.row = 4:32, layout.pos.col = 4:40))
grid.draw(dt.sc)
popViewport()

print(p.es.summary.sex + theme(plot.margin = margin(0.1, 0.04, 0.04, 0.04, "line")),
      vp = viewport(layout.pos.row = 7:29, layout.pos.col = 45:96))

pushViewport(viewport(layout.pos.row = 8:21, layout.pos.col = 97:207))
draw(p.expr.summary.sex, merge_legend = TRUE, heatmap_legend_side = "bottom",
     newpage = FALSE, background = "grey97")
popViewport()

pushViewport(viewport(layout.pos.row = 3:7, layout.pos.col = 5:100))
grid.text(x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"),
          label = expression("Meta-analysis of single cell, male " * italic("vs.") * " female"),
          gp = gpar(fontsize = textsize, fontface = "plain", col = textcolor, lineheight = 0.8))
popViewport()

pushViewport(viewport(layout.pos.row = 27:30, layout.pos.col = 99:128))
col_fun = circlize::colorRamp2(
  seq(-2, 2, length = 13), c(dichromat_pal("DarkRedtoBlue.12")(12)[1:6], "white",
    dichromat_pal("DarkRedtoBlue.12")(12)[7:12]))
draw(Legend(col_fun = col_fun, title = "NES",
            grid_width = unit(0.1, "line"), grid_height = unit(0.1, "line"),
            direction = "horizontal", legend_width = unit(2, "line"), title_position = "lefttop",
            title_gp = gpar(fontsize = 8, fontface = "plain"),
            labels_gp = gpar(fontsize = 8), at = seq(-2, 2, length = 3)))
upViewport(1)

print(pthreecircle + theme(plot.margin = margin(0.1, 0.1, 0.2, 0.04, "line")),
      vp = viewport(layout.pos.row = 26:29, layout.pos.col = 135:185))

print(p.arrow + theme(plot.margin = margin(0.04, 0.04, 0.3, 0.04, "line")),
      vp = viewport(layout.pos.row = 21:25, layout.pos.col = 97:207))

### purified neutrophils
pushViewport(viewport(layout.pos.row = 36:37, layout.pos.col = 11:94))
grid.text(x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"),
          label = "Bulk RNA-seq of neutrophils \n(E-GEAD-397)",
          gp = gpar(fontsize = textsize, fontface = "plain", col = textcolor, lineheight = 0.8))
popViewport()

pushViewport(viewport(layout.pos.row = 38:52, layout.pos.col = 3:94))
grid.draw(dt.ota)
popViewport()

print(p.otz.ageAdj + theme(plot.margin = margin(0.1, 0.1, 0, 0.1, "line")),
      vp = viewport(layout.pos.row = 51:97, layout.pos.col = 4:94))
print(p.sex.bar.otz.ageAdj + theme(plot.margin = margin(0.2, 0.1, 0, 0.1, "line")),
      vp = viewport(layout.pos.row = 93:103, layout.pos.col = 16:92))

pushViewport(viewport(layout.pos.row = 104:117, layout.pos.col = 16:86))
draw(ht.NES.ota.ageAdj, newpage = FALSE, background = "transparent")
decorate_annotation("pval", {
  grid.text("pval", x = unit(0.5, "npc"), y = unit(1, "npc") + unit(.2, "line"),
            gp = gpar(fontsize = 8), just = "bottom")
})
decorate_heatmap_body("NES", {
  grid.text("NES", x = unit(0.5, "npc"), y = unit(1, "npc") + unit(.2, "line"),
            gp = gpar(fontsize = 8), just = "bottom")
})
popViewport()

### bulk transcriptome
pushViewport(viewport(layout.pos.row = 36:37, layout.pos.col = 110:205))
grid.text(x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"),
          label = "Meta-analysis of bulk transcriptome\n",
          gp = gpar(fontsize = textsize, fontface = "plain", col = textcolor, lineheight = 0.8))
popViewport()

pushViewport(viewport(layout.pos.row = 35:52, layout.pos.col = 102:205))
grid.draw(dt.meta)
popViewport()

print(p.sex.PBMC.WB.ageAdj + theme(plot.margin = margin(0.1, 0.1, 0, 0.1, "line")),
      vp = viewport(layout.pos.row = 51:97, layout.pos.col = 110:206))
print(p.sex.WB.ageAdj + theme(plot.margin = margin(0.2, 0.1, 0.1, 0.1, "line")),
      vp = viewport(layout.pos.row = 93:102, layout.pos.col = 120:205))
print(p.sex.PBMC.ageAdj + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "line")),
      vp = viewport(layout.pos.row = 57:87, layout.pos.col = 100:118))

pushViewport(viewport(layout.pos.row = 103:117, layout.pos.col = 99:169))
draw(ht.NES.PBMC.ageAdj, newpage = FALSE, background = "transparent")
decorate_annotation("pval", {
  grid.text("pval", x = unit(0.5, "npc"), y = unit(1, "npc") + unit(.2, "line"),
            gp = gpar(fontsize = 8), just = "bottom")
})
decorate_heatmap_body("NES", {
  grid.text("PBMC\nNES", x = unit(0.5, "npc"), y = unit(1, "npc") + unit(.2, "line"),
            gp = gpar(fontsize = 8), just = "bottom")
})
popViewport()

pushViewport(viewport(layout.pos.row = 103:117, layout.pos.col = 168:206))
draw(ht.NES.WB.ageAdj, newpage = FALSE, background = "transparent")
decorate_annotation("pval", {
  grid.text("pval", x = unit(0.5, "npc"), y = unit(1, "npc") + unit(.2, "line"),
            gp = gpar(fontsize = 8), just = "bottom")
})
decorate_heatmap_body("NES", {
  grid.text("WB\nNES", x = unit(0.5, "npc"), y = unit(1, "npc") + unit(.2, "line"),
            gp = gpar(fontsize = 8), just = "bottom")
})
popViewport()

### framingham
pushViewport(viewport(layout.pos.row = 3:7, layout.pos.col = 228:280))
grid.text(x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"),
          label = paste0("Framingham cohort study, whole blood\n",
                         "822 healthy subjects (433 females, 389 males)"),
          gp = gpar(fontsize = textsize, fontface = "plain", col = textcolor, lineheight = 0.8))
popViewport()

print(p.fh.hc.list.bmi$immature + theme(plot.margin = margin(0.2, 0, 0.2, 0, "line")),
      vp = viewport(layout.pos.row = 8:30, layout.pos.col = 255:284))
print(p.fh.hc.list.age$immature + theme(plot.margin = margin(0.2, 0, 0.2, 0, "line")),
      vp = viewport(layout.pos.row = 8:30, layout.pos.col = 212:253))

print(p.fh.hc.list.bmi$mature + theme(plot.margin = margin(0.2, 0, 0.2, 0, "line")),
      vp = viewport(layout.pos.row = 8:30, layout.pos.col = 331:360))
print(p.fh.hc.list.age$mature + theme(plot.margin = margin(0.2, 0, 0.2, 0, "line")),
      vp = viewport(layout.pos.row = 8:30, layout.pos.col = 289:330))

print(p.fh.hc.list.bmi$degranulating + theme(plot.margin = margin(0.2, 0, 0.2, 0, "line")),
      vp = viewport(layout.pos.row = 31:53, layout.pos.col = 255:284))
print(p.fh.hc.list.age$degranulating + theme(plot.margin = margin(0.2, 0, 0.2, 0, "line")),
      vp = viewport(layout.pos.row = 31:53, layout.pos.col = 212:253))

print(p.fh.hc.list.bmi$antiprotease +
        scale_y_continuous(limits = c(4, 9.9), breaks = seq(5, 9, 2)) +
        theme(plot.margin = margin(0.2, 0, 0.2, 0, "line")),
      vp = viewport(layout.pos.row = 31:53, layout.pos.col = 331:360))
print(p.fh.hc.list.age$antiprotease +
        scale_y_continuous(limits = c(4, 9.9), breaks = seq(5, 9, 2)) +
        theme(plot.margin = margin(0.2, 0, 0.2, 0, "line")),
      vp = viewport(layout.pos.row = 31:53, layout.pos.col = 289:330))

print(p.fh.hc.list.bmi$total + theme(plot.margin = margin(0.2, 0, 0.2, 0, "line")),
      vp = viewport(layout.pos.row = 54:76, layout.pos.col = 255:284))
print(p.fh.hc.list.age$total + theme(plot.margin = margin(0.2, 0, 0.2, 0, "line")),
      vp = viewport(layout.pos.row = 54:76, layout.pos.col = 212:253))

print(p.fh.hc.list.bmi$IFN + theme(plot.margin = margin(0.2, 0, 0.2, 0, "line")),
      vp = viewport(layout.pos.row = 54:76, layout.pos.col = 331:360))
print(p.fh.hc.list.age$IFN + theme(plot.margin = margin(0.2, 0, 0.2, 0, "line")),
      vp = viewport(layout.pos.row = 54:76, layout.pos.col = 289:330))

pushViewport(viewport(layout.pos.row = 81:117, layout.pos.col = 212:330))
ht_opt(heatmap_border = FALSE)
draw(ht_ame + ht_ame.F + ht_ame.M, merge_legend = TRUE,
     heatmap_legend_side = "bottom", newpage = FALSE, background = "grey97")
ht_opt(RESET = TRUE)
for (nm in c("both sexes", "female", "male")) {  # use your actual heatmap names
  decorate_column_title(nm, {
    grid.rect(gp = gpar(fill = "grey91", col = NA))
    grid.text(cap_first(nm), gp = gpar(fontsize = 9, col = textcolor))
  })
}
upViewport(1)

pushViewport(viewport(layout.pos.row = 79:80, layout.pos.col = 225:280))
grid.text(x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"),
          label = "GLM, average marginal effects",
          gp = gpar(fontsize = textsize, fontface = "plain", col = textcolor, lineheight = 0.8))
popViewport()

print(p.twocircles + theme(plot.margin = margin(0.2, 0, 0.2, 0, "line")),
      vp = viewport(layout.pos.row = 107:110, layout.pos.col = 329:338))
print(pthreecircle + theme(plot.margin = margin(0.1, 0.1, 0.2, 0.04, "line")),
      vp = viewport(layout.pos.row = 112:115, layout.pos.col = 310:360))

### panel labels and titles
grid.text(x = unit(0.018, "npc"), y = unit(0.9, "npc"), label = "a",
          gp = gpar(fontsize = 8, fontface = "bold"))
grid.text(x = unit(0.12, "npc"), y = unit(0.955, "npc"), label = "b",
          gp = gpar(fontsize = 8, fontface = "bold"))
grid.text(x = unit(0.28, "npc"), y = unit(0.955, "npc"), label = "c",
          gp = gpar(fontsize = 8, fontface = "bold"))
grid.text(x = unit(0.32, "npc"), y = unit(0.94, "npc"), label = "Expression",
          gp = gpar(fontsize = 10, fontface = "plain"))

grid.text(x = unit(0.012, "npc"), y = unit(0.694, "npc"), label = "d",
          gp = gpar(fontsize = 8, fontface = "bold"))
grid.text(x = unit(0.29, "npc"), y = unit(0.694, "npc"), label = "e",
          gp = gpar(fontsize = 8, fontface = "bold"))

grid.text(x = unit(0.59, "npc"), y = unit(0.98, "npc"), label = "f",
          gp = gpar(fontsize = 8, fontface = "bold"))
grid.text(x = unit(0.59, "npc"), y = unit(0.333, "npc"), label = "g",
          gp = gpar(fontsize = 8, fontface = "bold"))

dev.off()
