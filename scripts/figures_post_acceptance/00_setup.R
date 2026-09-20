# ==============================================================================
# Shared setup: packages, options, figure settings and palettes
#
# Sourced first by every figure script. Defines the paths (dir.data /
# dir.results / dir.fig), the *size* constants, the post-acceptance style
# overrides and all col_* palettes.
# ==============================================================================


# -- paths ---------------------------------------------------------------------

dir.data    = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/data/"
dir.results = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/results/"
dir.fig     = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/figures/post_acceptance/"


# -- packages ------------------------------------------------------------------

library(expresso)
library(data.table)
library(ComplexHeatmap)
library(ggalluvial)
library(ggpubr)
library(ggrastr)
library(ggsci)
library(ggrepel)
library(ggstance)
library(ggforce)
library(ggtext)
library(ggthemes)
library(gtable)
library(pBrackets)
library(magick)
library(circlize)
library(viridis)
library(gridExtra)
library(scales)
library(fgsea)
library(Hmisc)
library(matrixStats)
library(survival)

ht_opt$TITLE_PADDING = unit(c(0.2, 0.2), "line")


# -- figure settings -----------------------------------------------------------

textsize       = 6
panellabelsize = 6
pointsizen     = 0.1
linesize       = 0.1
gridlwd        = 0.4
boxwidth       = 0.2

textcolor = linecolor = "#252525"
backgroundcol = "white"

pt2mm = function(size) size / ggplot2::.pt

# Sentence case for a displayed label, first letter only, leaving the rest alone so
# gene symbols and abbreviations inside a label survive: "peripheral blood" becomes
# "Peripheral blood", "PBMC" stays "PBMC". Applied where a label is *drawn* -- a
# scale's labels, a Heatmap's column_labels, a geom_text column -- and never to the
# data behind it, since the same strings are filtered and joined on elsewhere.
cap_first = function(x) sub("^(.)", "\\U\\1", x, perl = TRUE)

# The severity qualifier trails the condition in the data ("COVID-19,severe") and
# leads it in the manuscript ("Severe COVID-19"), so a drawn label swaps the two
# and then goes through cap_first(). A string without one of the two qualifiers --
# "bacterial sepsis", "healthy/convalescent", "flu,pregnancy" -- comes back
# sentence-cased and otherwise untouched, so this can stand in for cap_first()
# wherever a key mixes conditions. Drawn labels only, as with cap_first(): the
# comma forms are what the rows are subset, ordered and joined on.
cap_severity = function(x) {
  x = sub("^(.*),severe$",    "Severe \\1",     x)
  x = sub("^(.*),nonsevere$", "Non-severe \\1", x)
  cap_first(x)
}

# -- post-acceptance style overrides -------------------------------------------

# What the editors asked for, applied to every figure in this arm: no grey behind
# panels, facet strips or heatmap titles, no gridlines, one box width (boxwidth)
# for every box plot. expresso defaults the other way and reaches those defaults
# from inside its own helpers, where a call-site edit cannot follow, so the
# constructors are wrapped once here.

theme_nogrey = theme(panel.grid.major.x = element_blank(),
                     panel.grid.major.y = element_blank(),
                     panel.grid.minor   = element_blank(),
                     strip.background   = element_blank())

theme_expresso = function(..., text_size = textsize, axis_text_size = textsize,
                          axis_title_size = textsize, legend_text_size = textsize,
                          grid = "none")
  expresso::theme_expresso(..., text_size = text_size, axis_text_size = axis_text_size,
                           axis_title_size = axis_title_size,
                           legend_text_size = legend_text_size, grid = grid) + theme_nogrey
plot_scatter = function(..., label_size = pt2mm(textsize), grid = "none")
  expresso::plot_scatter(..., label_size = label_size, grid = grid) + theme_nogrey
plot_dotplot = function(..., col_fontsize = textsize, row_fontsize = textsize,
                        legend_textsize = textsize, grid = "none")
  expresso::plot_dotplot(..., col_fontsize = col_fontsize, row_fontsize = row_fontsize,
                         legend_textsize = legend_textsize, grid = grid) + theme_nogrey
theme_expresso_void = function(..., text_size = textsize)
  expresso::theme_expresso_void(..., text_size = text_size) + theme_nogrey
# plot_forest's own argument is called textsize, the same name as the constant, so
# the default is fetched explicitly rather than self-referentially.
#
# Three expresso defaults are overridden, and each is overridden by restating one
# property only, so the rest of what plot_forest() set survives: the study labels
# move to the left with the axes drawn (y_position = "right" gets expresso's side
# back), the plottitle's hjust is restated but not its size or margin, and the
# facet strip's size but not its face or margin. The y guide is repositioned
# rather than the scale replaced, because plot_forest() computes that scale's
# breaks, labels and limits from the data.

plot_forest = function(..., textsize = get("textsize", envir = globalenv()),
                       label_size = pt2mm(textsize), summary_space = 1,
                       y_position = c("left", "right"),
                       title_position = c("left", "center")) {
  y_position     = match.arg(y_position)
  title_position = match.arg(title_position)
  axisline = element_line(color = linecolor, linewidth = linesize)
  y_theme = if (y_position == "left")
    theme(axis.line.y.left = axisline, axis.ticks.y.left = axisline)
  else
    theme(axis.line.y.right = axisline, axis.ticks.y.right = axisline)
  p = expresso::plot_forest(..., textsize = textsize) + theme_nogrey +
    guides(y = guide_axis(position = y_position)) +
    theme(axis.line.x.bottom  = axisline,
          axis.ticks.x.bottom = axisline,
          strip.text.x        = element_text(size = textsize),
          plot.title          = element_text(hjust = if (title_position == "center") 0.5
                                                     else 0)) +
    y_theme
  # The pooled row is drawn "Summary", and capitalised here on the built plot
  # rather than in the data: plot_forest() finds that row by matching
  # rowname == "summary", so a capitalised one is drawn as a study, losing its
  # diamond, its rules and its place at the foot of the panel. Only the exact
  # string is touched -- 04_infection_cancer.R's study rows carry reference
  # numbers, and capitalising those would draw "Ref 14".
  labelscale = p$scales$get_scales("y")
  if (is.character(labelscale$labels))
    labelscale$labels[labelscale$labels == "summary"] = "Summary"
  extra = summary_space - 1
  center = -extra / 2
  if (extra != 0) {
    summaryrow = which(is.na(p$data$lower))
    if (length(summaryrow) != 1L)
      stop("summary_space needs one summary row, found ", length(summaryrow))
    yscale = p$scales$get_scales("y")
    breaks = yscale$breaks
    breaks[summaryrow] = center
    yscale$breaks = breaks
    yscale$limits = c(-0.51 - extra, max(p$data$index) + 0.51)
  }
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomText")) {
      p$layers[[i]]$aes_params$size = label_size
      if (extra != 0) p$layers[[i]]$data$y = center
    }
    if (extra != 0 && inherits(p$layers[[i]]$geom, "GeomHline"))
      p$layers[[i]]$data$yintercept =
        ifelse(p$layers[[i]]$data$yintercept < 0,
               p$layers[[i]]$data$yintercept - extra, p$layers[[i]]$data$yintercept)
    if (extra != 0 && inherits(p$layers[[i]]$geom, "GeomPolygon"))
      p$layers[[i]]$data$y = p$layers[[i]]$data$y + center
  }
  p
}


# -- palettes ------------------------------------------------------------------

# study of origin.
# Keyed by corrected publication year: the raw pid values myin2023 / reyfman2018 /
# xue2022 are renamed to myin2024 / reyfman2019 / xue2023 before this is applied.
col_study = c(
  chan2021 = "#FFF143", chen2020 = "#FDD834", combes2021 = "#000080", deng2021 = "#EE4C97",
  garridotrigo2023 = "#FFE099", gupta2020 = "#A60021",
  habermann2020 = "#2F4F4F", han2020 = "#F3D3E7", hu2022 = "#5E2D30", hu2023 = "#73DAFF",
  jardine2021 = "#FFC0CB", kaiser2024 = "#F47F17", kwok2023 = "grey67", liao2020 = "#2A0BD9",
  maynard2020 = "#40A1FF", montaldo2022 = "#FFAD73", moreno2022 = "#28231D", myin2024 = "#EAFF56",
  qian2020 = "#8CC269", qiu2021 = "#1A5E1F", ren2021 = "#FFFFBF", reyfman2019 = "#FF00FF",
  salcher2022 = "#DDB952", schrepping2020 = "#E0FFFF",
  schupp2020 = "#FD8CC1", sinha2021 = "#D1C4E9", suo2022 = "#1C77A3", tabulasapiens2022 = "#6A1A99",
  wang2021 = "#FFFF00", wang2023 = "#CD92D8", wigerblad2022 = "#E9D097", wilhelm2024 = "#CE5E8A",
  wilk2021 = "#D92632", wu2024 = "#C7E5C9", xue2023 = "#00FF00", yang2022 = "#ABF8FF",
  yang2023 = "#7E57C1", zhang2023 = "#2BAE85", zhu2022 = "#F76E5E", zilionis2019 = "#264EFF")

# clinical group
col_group = c(
  healthy = "#1A6840", convalescent = "#4EAE3C",
  `COVID-19` = "#EFC000", `COVID-19,nonsevere` = "#7b72c5",
  `COVID-19,severe` = "#EFC000", `bacterial sepsis` = "#F97D1C",
  flu = "cyan", `flu,pregnancy` = "green", pregnancy = "lightgreen", `lung cancer` = "#DC3023",
  `GI cancer` = "#FF0097", `kidney cancer` = "#6A1A99",
  `colorectal cancer` = "#CE5E8A", `brain tumor` = "#7C1823",
  `other cancer` = "#EEA2A4", asthma = "#142334", `other lung disease` = "#352A87", CVD = "#0E5FDB",
  diabetes = "#C5BB5C", `autoimmune disease` = "#06A5C7",
  `HSC-T` = "#FFDFB2", `healthy,G-CSF/IFN` = "#C7E5C9FF",
  fetus = "#D6BBCF", `other disease` = "grey87")

# tissue
col_tissue = c(
  `peripheral blood` = "#EDD1D8", `cord blood` = "#FBC831",
  `bone marrow` = "#73DAFF", lung = "#DC3023", BALF = "#BEC936", sputum = "#0E5FDB",
  liver = "#B35C44", colonrectum = "#EE6C00", pancreas = "#FF0097", stomach = "#CE5E8A",
  kidney = "#B0A4E3", brain = "#7C1823", spleen = "#A5D6A6", `other GI organs` = "#FFDFB2",
  `urinary system` = "#352A87", other = "grey77")

# neutrophil subsets
col_celltype = c(
  AZU1 = "#DC3023", LTF = "#F47B00", MMP9 = "#F8D626", S100A4 = "#FFDFB2",
  IL1R2 = "#acf300", MME = "#80C684", TXNIP = "#D1BA58", EGR1 = "#1A5E1F",
  `AP-1` = "#2BAE85", PTGS2 = "#AFBC65", G0S2 = "#E2E7BF", IFN1 = "#2E42B8",
  IFN2 = "#0A73DC", IFN3 = "#264EFF", CXCL = "#E4C6D0", VEGFA = "#73DAFF",
  `CCL3/4` = "#DCE318", `NF-κB` = "#5E34B1", IL1B = "#9E9AC8", IL1RN = "#DCB0F2",
  SLPI = "#7CABB1", HSP = "#4D6D93", CD74 = "#C7E5C9", lowDepth = "grey77")

# signature sets, by maturation stage
col_set = c(
  immature = "#a93434", degranulating = "#F47B00", antiprotease = "#7CABB1", mature = "#1A5E1F",
  IFN = "#2E42B8", other = "grey67")


# -- shared theme fragments ----------------------------------------------------

theme_legend_right = theme(
  legend.position = c(1, 1), legend.justification = c(0, 1),
  legend.direction = "vertical",
  legend.background = element_blank(), legend.key = element_blank(),
  legend.margin = margin(0, 0, 0, 0, "line"), legend.key.height = unit(0.45, "line"))
