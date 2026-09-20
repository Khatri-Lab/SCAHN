# ==============================================================================
# GSE153263 sorted progenitors, and EHR outcomes by ANC quartile
#
# Two unrelated extended figures. Extended figure 1 is a marker heatmap across
# the four sorted populations of the GSE153263 bulk RNA-seq series (eNeP, N1
# without eNeP, N2, mature neutrophils). Extended figure 9 is survival and
# cumulative-incidence curves for the Charlson comorbidity count and for death,
# in a Stanford EHR cohort and a national one, each stratified by quartile of
# absolute neutrophil count, plus a table of the same hazard ratios across ten
# outcomes.
#
# Requires 00_setup.R (packages, figure settings, textsize / panellabelsize /
# textcolor, cap_first / pt2mm, dir.data / dir.results / dir.fig).
#
# Reads:
#   data/GSE153263.expr.rds
#   results/clinical/EHR/stanford_prop_surv_<outcome>.csv       10 outcomes
#   results/clinical/EHR/stanford_prop_iqr_surv_<outcome>.csv   10 outcomes
#   results/clinical/EHR/national_prop_surv_<outcome>.csv        9 outcomes
#
# Writes to figures/post_acceptance/:
#   Fig.e1.pdf   Extended figure 1
#   Fig.e9.pdf   Extended figure 9
# ==============================================================================

SCAHN_SCRIPTS = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/scripts/figures_post_acceptance"
source(file.path(SCAHN_SCRIPTS, "00_setup.R"))


# -- GSE153263 sorted neutrophils ----------------------------------------------

dtf.GSE153263 = readRDS(paste0(dir.data, "GSE153263.expr.rds"))
dtf.GSE153263$group = gsub("MatureNeuts", "Mature Neu", dtf.GSE153263$group)

# Columns are ordered by sorted population, not clustered, so the maturation
# sequence reads left to right.
colorder = c("GSM4637184", "GSM4637182", "GSM4637183", "GSM4637179", "GSM4637180", "GSM4637181",
             "GSM4637175", "GSM4637173", "GSM4637174", "GSM4637178", "GSM4637176", "GSM4637177")

dtf.GSE153263 = dtf.GSE153263[order(match(sampleid, colorder))]

genes = c("ITGA4", "TFRC", "KIT", "ITGAM", "FCGR3B", "CD101", "MME", "CR1",
          "AZU1", "BPI", "CAMP", "CEACAM8", "DEFA3", "DEFA4",
          "ELANE", "LCN2", "LTF", "MMP8", "MPO", "OLFM4",
          "CXCR1", "CXCR2", "CXCR4", "STMN1", "MKI67", "ARG1", "MMP9", "S100A8", "OLR1")


# -- Row labels ----------------------------------------------------------------

cd_map = c(
  ITGA4 = "CD49d", TFRC = "CD71", KIT = "CD117", ITGAM = "CD11b",
  FCGR3B = "CD16b", CD101 = "CD101", MME = "CD10", CR1 = "CD35",
  CEACAM8 = "CD66b", CXCR1 = "CD181", CXCR2 = "CD182", CXCR4 = "CD184")

label_calls = lapply(genes, function(g) {
  prot = cd_map[g]
  if (!is.na(prot) && prot != g) {
    # gene ≠ protein → italic gene + plain protein in parens
    bquote(italic(.(g)) * " (" * .(prot) * ")")
  } else {
    # gene == protein (or no map) → just italic gene
    bquote(italic(.(g)))
  }
})
# turn into an expression vector and name it
expr_labels = do.call("expression", label_calls)
names(expr_labels) = genes


# -- Scaled expression matrix --------------------------------------------------

mat = as.matrix(dtf.GSE153263[,genes, with = FALSE])
mat = scale(mat)
mat = t(mat)
colnames(mat) = dtf.GSE153263$sampleid

col_type = c("eNeP" = "#F6DE20", "N1_wo_eNeP" = "#F0BB46",
             "N2" = "#AFBC65", "Mature Neu" = "#0E5FDB")

col_fun = circlize::colorRamp2(
  seq(-max(abs(mat)), max(abs(mat)), length = 13),
  c(dichromat_pal("DarkRedtoBlue.12")(12)[1:6], "white",
    dichromat_pal("DarkRedtoBlue.12")(12)[7:12]))


# -- Heatmap -------------------------------------------------------------------

ha = HeatmapAnnotation(
  type = factor(dtf.GSE153263$group, levels = names(col_type)), col = list(type = col_type),
  annotation_legend_param = list(
    type = list(title = "type", title_position = "leftcenter",
    title_gp = gpar(fontsize = textsize, col = textcolor),
                labels_gp = gpar(fontsize = textsize, col = textcolor),
                legend_gp = gpar(fontsize = textsize), direction = "horizontal",
                nrow = 1, legend_width = unit(0.2, "line"), legend_height = unit(.01, "line"),
                grid_width = unit(0.2, "line"), grid_height = unit(.01, "line"))
  ), annotation_name_gp = gpar(fontsize = textsize,
                               col = textcolor), show_annotation_name = TRUE, gap = unit(.05,
  "line"),
  show_legend = TRUE, simple_anno_size = unit(0.5, "line"))

pht =
  Heatmap(mat, col = col_fun, na_col = "grey", cluster_columns = FALSE, cluster_rows = TRUE,
          show_row_dend = FALSE, show_column_dend = TRUE, column_dend_reorder = FALSE,
          column_dend_side = "top", row_dend_gp = gpar(lwd = 0.3),
          column_dend_gp = gpar(lwd = 0.3), bottom_annotation = ha,
          row_labels = expr_labels[rownames(mat)], clustering_method_rows = "ward.D2",
          clustering_method_columns = "ward.D2", clustering_distance_rows = "spearman",
          clustering_distance_columns  = "euclidean",
          show_column_names = FALSE, show_row_names = TRUE,
          column_names_side = "top", column_names_rot = 90,
          column_names_gp = gpar(fontsize = textsize, col = textcolor),
          row_names_gp = gpar(fontsize = textsize, fontface = "italic",
                              col = textcolor), row_names_side = "right",
          row_names_max_width = unit(30, "line"), column_title = "", column_title_side = "bottom",
          column_title_gp = gpar(fontsize = textsize, fontface = "plain", col = textcolor),
          show_heatmap_legend = TRUE,
          heatmap_legend_param = list(title = "Scaled expression",
          title_gp = gpar(fontsize = textsize, col = textcolor),
                                      labels_gp = gpar(fontsize = textsize, col = textcolor),
                                      legend_width = unit(.2, "line"),
                                      legend_height = unit(2, "line"),
                                      grid_width = unit(.2, "line"), grid_height = unit(2, "line"),
                                      by_row = FALSE, direction = "vertical",
                                      title_position = "leftcenter-rot", at = c(-2,0,2)))


# -- Extended figure 1 ---------------------------------------------------------

draw_hline_label = function(line_rows, cols, label, gap = 2,
                            text_gap = 0.7,           # label to rule, in mm
                            # Multiple of fontsize between the baselines of a "\n" label.
                            # 1.2 is grid's default; lower it to tighten a two-line label.
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

cairo_pdf(file = paste0(dir.fig, "Fig.e1.pdf"), width = 2.9, height = 3.6)
pushViewport(viewport(layout = grid.layout(nrow = 90, ncol = 86)))
pushViewport(viewport(layout.pos.row = 6:90, layout.pos.col = 2:86))
draw(pht, merge_legend = FALSE, heatmap_legend_side = "right",
     align_heatmap_legend = "heatmap_top",
     annotation_legend_side = "bottom", newpage = FALSE)
upViewport(1)


draw_hline_label(line_rows = 5:5, cols = 5:37, label = "Immature Neu",
                 textsize = textsize, textcol = textcolor)
draw_hline_label(line_rows = 5:5, cols = 40:49, label = "Mature Neu",
                 textsize = textsize, textcol = textcolor)
dev.off()


# -- EHR survival curves -------------------------------------------------------
dir.EHR = paste0(dir.results, "clinical/EHR/")

outcomes.condition = c("chf", "ckd", "copd", "cancer", "dementia", "diabetes", "liver", "pvd")
outcomes = c("charlson", outcomes.condition)

read_EHR = function(file, incidence) {
  dt = fread(paste0(dir.EHR, file, ".csv"))
  if (incidence) dt[, c("surv", "lower", "upper") := .(1 - surv, 1 - lower, 1 - upper)]
  dt[]
}

stanford = sapply(c("death", outcomes), function(o)
  rbind(read_EHR(paste0("stanford_prop_surv_", o), o != "death"),
        read_EHR(paste0("stanford_prop_iqr_surv_", o), o != "death")[strata == "IQR"]),
  simplify = FALSE)

national = sapply(outcomes, function(o) read_EHR(paste0("national_prop_surv_", o), TRUE),
                  simplify = FALSE)


# -- Cox model results ---------------------------------------------------------

hr.stanford.H = c(
  death = "1.98 (1.62 - 2.42)", charlson = "1.18 (1.13 - 1.23)", chf = "1.40 (1.24 - 1.58)",
  ckd = "1.34 (1.20 - 1.48)", copd = "1.10 (1.03 - 1.18)", cancer = "1.06 (0.98 - 1.14)",
  dementia = "1.60 (1.30 - 1.98)", diabetes = "1.79 (1.61 - 1.98)",
  liver = "1.11 (0.78 - 1.59)", pvd = "1.27 (1.13 - 1.43)")
p.stanford.H = c(
  death = "p < 0.001", charlson = "p < 0.001", chf = "p < 0.001", ckd = "p < 0.001",
  copd = "p = 0.003", cancer = "p = 0.17", dementia = "p < 0.001", diabetes = "p < 0.001",
  liver = "p = 0.56", pvd = "p < 0.001")

hr.stanford.I = c(
  death = "1.43 (1.20 - 1.69)", charlson = "1.16 (1.12 - 1.20)", chf = "1.37 (1.24 - 1.52)",
  ckd = "1.26 (1.16 - 1.38)", copd = "1.11 (1.05 - 1.17)", cancer = "1.08 (1.01 - 1.15)",
  dementia = "1.51 (1.27 - 1.79)", diabetes = "1.53 (1.41 - 1.67)",
  liver = "1.08 (0.78 - 1.50)", pvd = "1.13 (1.03 - 1.24)")
p.stanford.I = c(
  death = "p < 0.001", charlson = "p < 0.001", chf = "p < 0.001", ckd = "p < 0.001",
  copd = "p < 0.001", cancer = "p = 0.02", dementia = "p < 0.001", diabetes = "p < 0.001",
  liver = "p = 0.63", pvd = "p = 0.01")

hr.national.H = c(
  charlson = "1.27 (1.21 - 1.33)", chf = "1.57 (1.33 - 1.86)", ckd = "1.39 (1.22 - 1.58)",
  copd = "1.22 (1.14 - 1.30)", cancer = "1.25 (1.12 - 1.41)", dementia = "1.96 (1.50 - 2.56)",
  diabetes = "1.45 (1.32 - 1.59)", liver = "1.40 (0.79 - 2.49)", pvd = "1.34 (1.21 - 1.48)")
p.national.H = c(
  charlson = "p < 0.001", chf = "p < 0.001", ckd = "p < 0.001", copd = "p < 0.001",
  cancer = "p < 0.001", dementia = "p < 0.001", diabetes = "p < 0.001",
  liver = "p = 0.25", pvd = "p < 0.001")


# -- EHR panels ----------------------------------------------------------------

col_quartile = c(HIGHEST.QUARTILE = "#FF7F0E", IQR = "goldenrod2", LOWEST.QUARTILE = "#214DC8")

lab_quartile = c(HIGHEST.QUARTILE = "4th quartile (H)", IQR = "IQR (I)",
                 LOWEST.QUARTILE = "1st quartile")

titles = c(charlson = "Charlson comorbidities", chf = "congestive heart failure",
           ckd = "chronic kidney disease", copd = "COPD", cancer = "cancer",
           dementia = "dementia", diabetes = "diabetes", liver = "liver disease",
           pvd = "peripheral vascular diseases")

pnum = function(p) if (grepl("^p <", p)) 0 else as.numeric(sub("^p = ", "", p))
sigcol = function(p) if (pnum(p) < 0.05) textcolor else "grey67"

EHR_base = function(dt, ptitle, laby, xbreaks) {
  dt = copy(dt)
  dt[, c("time", "surv", "lower", "upper") :=
       .(time / 365, surv * 100, lower * 100, upper * 100)]
  ggplot(dt, aes(x = time, y = surv, color = strata, fill = strata)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
    labs(title = cap_first(ptitle), x = "Time (years)", y = cap_first(laby)) +
    scale_x_continuous(breaks = xbreaks, limits = range(xbreaks)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_color_manual(values = col_quartile, labels = lab_quartile, name = NULL) +
    scale_fill_manual(values = col_quartile, labels = lab_quartile, name = NULL) +
    guides(fill = "none") +
    theme_expresso() + theme(legend.position = "none")
}

EHRplotStanford = function(dt, ptitle, laby, labelx, labely, hjust, vjust, hr1, p1, hr2, p2) {
  EHR_base(dt, ptitle, laby, seq(0, 10, 2)) +
    annotate("text", x = labelx, y = labely, hjust = hjust, vjust = vjust, lineheight = 0.9,
             color = textcolor, size = pt2mm(textsize),
             label = paste0("HR(H): ", hr1, "\n", p1, "\n", "HR(I): ", hr2, "\n", p2))
}

EHRplotStanford2 = function(dt, ptitle, laby, labelx, labely, labely2, hjust, vjust,
                            hr1, p1, hr2, p2,
                            labelcolor1 = textcolor, labelcolor2 = textcolor) {
  hrgrob = function(label, y, col)
    annotation_custom(grob = textGrob(
      label = label, x = unit(labelx, "npc"), y = unit(y, "npc"), just = c(hjust, vjust),
      gp = gpar(col = col, fontsize = textsize, lineheight = 0.8)))
  EHR_base(dt, ptitle, laby, seq(0, 10, 2)) + coord_cartesian(clip = "off") +
    hrgrob(paste0("HR(H): ", hr1, "\n", p1), labely,  labelcolor1) +
    hrgrob(paste0("HR(I): ", hr2, "\n", p2), labely2, labelcolor2)
}

EHRplotNational = function(dt, ptitle, laby, labelx, labely, vjust, hjust, hr1, p1,
                           labelcolor = textcolor) {
  EHR_base(dt, ptitle, laby, seq(0, 5, 1)) +
    annotate("text", x = labelx, y = labely, hjust = hjust, vjust = vjust, lineheight = 0.9,
             color = labelcolor, size = pt2mm(textsize),
             label = paste0("HR(H): ", hr1, "\n", p1))
}

death_plot = EHRplotStanford(stanford$death, "Stanford", "survival rate (%)",
                             labelx = 0.1, labely = 98, hjust = 0, vjust = 1,
                             hr.stanford.H[["death"]], p.stanford.H[["death"]],
                             hr.stanford.I[["death"]], p.stanford.I[["death"]]) +
  scale_y_continuous(breaks = seq(0, 100, 1)) + theme_legend_right

charlson_plot = EHRplotStanford(
  stanford$charlson, "Stanford", "Charlson comorbidities (%)",
  labelx = 0, labely = Inf, hjust = 0, vjust = 1,
  hr.stanford.H[["charlson"]], p.stanford.H[["charlson"]],
  hr.stanford.I[["charlson"]], p.stanford.I[["charlson"]]) + theme_legend_right

charlson_national_plot = EHRplotNational(
  national$charlson, "National", "Charlson comorbidities (%)",
  labelx = 0, labely = Inf, vjust = 1, hjust = 0,
  hr.national.H[["charlson"]], p.national.H[["charlson"]]) + theme_legend_right

plots.stanford = sapply(outcomes.condition, function(o)
  EHRplotStanford2(stanford[[o]], titles[[o]], "percentage",
                   labelx = 0.01, labely = 1, labely2 = 0.83, hjust = 0, vjust = 1,
                   hr.stanford.H[[o]], p.stanford.H[[o]], hr.stanford.I[[o]], p.stanford.I[[o]],
                   labelcolor1 = sigcol(p.stanford.H[[o]]),
                   labelcolor2 = sigcol(p.stanford.I[[o]])), simplify = FALSE)

plots.national = sapply(outcomes.condition, function(o)
  EHRplotNational(national[[o]], titles[[o]], "percentage",
                  labelx = 0, labely = Inf, vjust = 1, hjust = 0,
                  hr.national.H[[o]], p.national.H[[o]],
                  labelcolor = sigcol(p.national.H[[o]])), simplify = FALSE)


# -- EHR hazard ratio table ----------------------------------------------------

table.rows = c("diabetes", "dementia", "chf", "ckd", "pvd", "copd", "cancer", "liver",
               "charlson")
table.labels = c(titles[setdiff(table.rows, "charlson")], charlson = "composite")

# Cell fill by p value: the two shades the key names, white for a p that clears
# neither threshold.
col_pfill = c(`p<0.001` = "#D1C4E9", `p<0.05` = "#EDE7F6")
pfill = function(p) if (pnum(p) < 0.001) col_pfill[["p<0.001"]] else
                    if (pnum(p) < 0.05)  col_pfill[["p<0.05"]]  else "white"

hr.cells = cbind(hr.stanford.H[table.rows], hr.stanford.I[table.rows],
                 hr.national.H[table.rows])
p.cells  = cbind(p.stanford.H[table.rows], p.stanford.I[table.rows],
                 p.national.H[table.rows])

cells = rbind(c("", "", "", ""),
              c("", "4th vs. 1st", "IQR vs. 1st", "4th vs. 1st"),
              cbind(cap_first(table.labels), hr.cells))
fills = cbind("white", rbind("white", "white",
                             apply(p.cells, 2, function(x) sapply(x, pfill))))
# Row labels sit just off the rule to their right, everything else is centred.
hj = cbind(1,    matrix(0.5, nrow(cells), ncol(cells) - 1))
xx = cbind(0.97, matrix(0.5, nrow(cells), ncol(cells) - 1))

hr.table = tableGrob(
  cells, rows = NULL, cols = NULL, theme = ttheme_default(
    core = list(
      fg_params = list(fontsize = textsize, fontface = "plain", col = textcolor,
                       hjust = as.vector(hj), x = as.vector(xx)),
      bg_params = list(fill = as.vector(fills), col = NA),
      padding = unit(c(0.5, 0.35), "line"))))

### the three spanning header cells
hr.table = gtable_add_grob(
  hr.table, textGrob("Hazard ratios (95% CI)", x = 0.97, hjust = 1,
                     gp = gpar(fontsize = textsize, col = textcolor)),
  t = 1, b = 2, l = 1, r = 1, name = "stub")
hr.table = gtable_add_grob(
  hr.table, textGrob("Stanford", gp = gpar(fontsize = textsize, col = textcolor)),
  t = 1, b = 1, l = 2, r = 3, name = "head.stanford")
hr.table = gtable_add_grob(
  hr.table, textGrob("National", gp = gpar(fontsize = textsize, col = textcolor)),
  t = 1, b = 1, l = 4, r = 4, name = "head.national")

### rules

hrule = function(y) segmentsGrob(x0 = 0, x1 = 1, y0 = y, y1 = y,
                                 gp = gpar(lwd = 0.4, col = textcolor))
vrule = function(x) segmentsGrob(x0 = x, x1 = x, y0 = 0, y1 = 1,
                                 gp = gpar(lwd = 0.4, col = textcolor))
rules = list(list(hrule(1), 1,  1,  1, 4),   # above the header
             list(hrule(0), 1,  1,  2, 3),   # under Stanford
             list(hrule(0), 1,  1,  4, 4),   # under National
             list(hrule(0), 2,  2,  1, 4),   # under the sub-headers
             list(hrule(0), 10, 10, 1, 4),   # between the conditions and the composite
             list(hrule(0), 11, 11, 1, 4),   # under the composite
             list(vrule(1), 1,  11, 1, 1),   # right of the row labels
             list(vrule(1), 1,  11, 3, 3))   # between the two cohorts
for (i in seq_along(rules))
  hr.table = gtable_add_grob(hr.table, rules[[i]][[1]], t = rules[[i]][[2]],
                             b = rules[[i]][[3]], l = rules[[i]][[4]], r = rules[[i]][[5]],
                             name = paste0("rule.", i))

df.pkey = data.table(label = names(col_pfill), fill = unname(col_pfill), y = c(2, 1))

p.pvalkey = ggplot(df.pkey, aes(x = 1, y = y, label = label, fill = fill)) +
  geom_label(size = pt2mm(textsize), color = textcolor, linewidth = NA,
             label.r = unit(0, "pt"), label.padding = unit(0.3, "line")) +
  scale_fill_identity() + expand_limits(y = c(0.3, 2.7)) +
  coord_cartesian(clip = "off") + theme_void()


# -- Extended figure 9 ---------------------------------------------------------
cairo_pdf(file = paste0(dir.fig, "Fig.e9.pdf"), width = 6.6, height = 4.1)
pushViewport(viewport(layout = grid.layout(nrow = 62, ncol = 95)))

print(charlson_plot + theme(plot.margin = margin(0.04, 0.1, 0.04, 0.1, "line")),
      vp = viewport(layout.pos.row = 3:31, layout.pos.col = 2:28))
print(charlson_national_plot + theme(plot.margin = margin(0.04, 0.1, 0.04, 0.1, "line")),
      vp = viewport(layout.pos.row = 3:31, layout.pos.col = 59:85))

print(death_plot + theme(plot.margin = margin(0.04, 0.1, 0.04, 0.2, "line")),
      vp = viewport(layout.pos.row = 33:61, layout.pos.col = 59:85))

pushViewport(viewport(layout.pos.row = 36:57, layout.pos.col = 2:54))
pushViewport(viewport(x = 0, width = sum(hr.table$widths), just = "left"))
grid.draw(hr.table)
upViewport(1)
print(p.pvalkey, vp = viewport(x = sum(hr.table$widths), y = 0.8,
                               width = unit(14, "mm"), height = unit(12, "mm"),
                               just = c("left", "centre")))
upViewport(1)

### panel labels
grid.text(x = unit(0.02, "npc"), y = unit(0.98, "npc"), label = "a",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", fontfamily = "Arial",
                    col = textcolor))

grid.text(x = unit(0.02, "npc"), y = unit(0.5, "npc"), label = "b",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", fontfamily = "Arial",
                    col = textcolor))

grid.text(x = unit(0.62, "npc"), y = unit(0.98, "npc"), label = "c",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", fontfamily = "Arial",
                    col = textcolor))

grid.text(x = unit(0.62, "npc"), y = unit(0.5, "npc"), label = "d",
          gp = gpar(fontsize = panellabelsize, fontface = "bold", fontfamily = "Arial",
                    col = textcolor))

dev.off()
