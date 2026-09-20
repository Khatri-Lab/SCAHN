# ==============================================================================
# Gene detection rates, pseudobulk expression, and per-cell signature scores
#
# Everything downstream of the integrated object that the figure scripts and the
# meta-analyses need, in four blocks:
#   1. the fraction of cells expressing each gene
#   2. geometric-mean pseudobulk of every gene, at five groupings
#   3. per-cell scores for the SCAHN neutrophil signatures
#   4. per-cell scores for published mouse and human neutrophil gene sets
#      (NeuMap, xie2020, GSE165276, ng2024)
#
# Expensive: loads the full integrated object and densifies its expression
# matrix. Run once, on its own, not as part of the figure scripts.
#
# Reads:
#   data/SCAHN/SCAHN.srt.integrated.rds
#   data/cross_species/mousegenes_tohuman.rds
#   data/cross_species/mousegenes_tohuman.manual.csv
#   data/cross_species/xie2020.genesets.csv
#   data/cross_species/xie2020.genesets.human.csv
#   data/cross_species/neutrotime.humangenes.csv
#   data/cross_species/ng2024.genesets.csv
#   data/cross_species/NeuMAP_mouse_degs.xlsx
#   data/cross_species/hNeuMAP_integrated_hub_markers.xlsx
#   results/annotation/SCAHN.cellmeta.csv
#   results/annotation/SCAHN.subset_categ.markers.csv
#
# Writes:
#   results/annotation/gene.exprpct.rds        read by 04_infection_cancer.R
#   results/annotation/SCAHN.cellscores.rds    read by 02_atlas.R
#   results/cross_species/SCAHN.neumap.scores.rds
#                                              read by 03_cross_species.R
#   results/meta_analysis/dtfgenes.allneu.persample.geommean.rds
#                                              read by 04_infection_cancer.R
#   results/meta_analysis/dtfgenes.subsets_categ.persample.geommean.rds
#                                              read by meta_sc_infection.R
#   results/meta_analysis/dtfgenes.allneu.geommean.rds
#   results/meta_analysis/dtfgenes.subsets.geommean.rds
# ==============================================================================

library(data.table)
library(Seurat)
library(readxl)

dir.data    = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/data/"
dir.results = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/results/"


# -- Integrated object and matching cell metadata ------------------------------

dtf = fread(paste0(dir.results, "annotation/SCAHN.cellmeta.csv"))
srt.sct = readRDS(paste0(dir.data, "SCAHN/SCAHN.srt.integrated.rds"))

# Every matrix below is pulled from srt.sct while every grouping variable comes
# from dtf, so the two must be in the same cell order.
all.equal(dtf$cell.pid, as.character(srt.sct$cell.pid))


# -- Fraction of cells expressing each gene ------------------------------------

counts = srt.sct@assays$RNA@counts
pct_expressed = rowSums(counts > 0) / ncol(counts) * 100
saveRDS(pct_expressed, paste0(dir.results, "annotation/gene.exprpct.rds"))


# -- Pseudobulk ----------------------------------------------------------------

# One row per cell, one column per gene, shifted by 1 so the geometric mean is
# defined on zeros. This is the dense step, and the reason the script is run on
# its own.
mat = t(as.matrix(srt.sct@assays$RNA@data)) + 1

dtf$dummy = "all neutrophils"
dtfgenes = as.data.table(cbind(dtf, mat))

### geometric mean over all neutrophils
dtfgenes.mean = dtfgenes[, lapply(.SD, function(x) {expresso::geom_mean(x)}), by = c("dummy"),
                         .SDcols = setdiff(colnames(dtfgenes), colnames(dtf))]
saveRDS(dtfgenes.mean, paste0(
  dir.results, "meta_analysis/dtfgenes.allneu.geommean.rds"))

### geometric mean per neutrophil subset
dtfgenes.mean = dtfgenes[, lapply(.SD, function(x) {expresso::geom_mean(x)}), by = c("celltype"),
                         .SDcols = setdiff(colnames(dtfgenes), colnames(dtf))]
saveRDS(dtfgenes.mean, paste0(
  dir.results, "meta_analysis/dtfgenes.subsets.geommean.rds"))

### geometric mean per sample, all neutrophils
dtfgenes.mean = dtfgenes[, lapply(.SD, function(x) {expresso::geom_mean(x)}),
                         by = c("sampleid.pid", "pid"),
                         .SDcols = setdiff(colnames(dtfgenes), colnames(dtf))]
saveRDS(dtfgenes.mean, paste0(
  dir.results, "meta_analysis/dtfgenes.allneu.persample.geommean.rds"))

### geometric mean per sample, within coarse subset categories
# The 23 subsets are too fine to meta-analyse per sample, so the granule, IFN
# and cytokine subsets are pooled into three categories.
dtfgenes[, celltype2 := ifelse(grepl("AZU1|LTF|MMP9", celltype), "granule subsets", celltype)]
dtfgenes[, celltype2 := ifelse(grepl("IFN", celltype), "IFN subsets", celltype2)]
dtfgenes[, celltype2 := ifelse(grepl("NF-κB|IL1B|IL1RN|CCL3|CXCL|VEGFA", celltype),
                               "cytokine subsets", celltype2)]

dtfgenes.mean2 = dtfgenes[
  celltype2 %in% c("granule subsets", "IFN subsets", "cytokine subsets")][
  , lapply(.SD, function(x) {expresso::geom_mean(x)}), by = c("sampleid.pid", "celltype2", "pid"),
  .SDcols = setdiff(colnames(dtfgenes), colnames(dtf))]
saveRDS(dtfgenes.mean2, paste0(
  dir.results, "meta_analysis/dtfgenes.subsets_categ.persample.geommean.rds"))


# -- Scores for the SCAHN signatures -------------------------------------------

dtf = fread(paste0(dir.results, "annotation/SCAHN.cellmeta.csv"))

genes.granules = c("PRTN3", "CTSG", "ELANE",
                   "AZU1", "MPO", "DEFA4", "DEFA3", "BPI", "CEACAM6", "CEACAM8",
                   "OLFM4", "LTF", "LCN2", "MMP8", "CAMP", "CRISP3", "TCN1", "HP", "FCN1", "PADI4",
                   "PGLYRP1", "ANXA3", "ARG1", "CD177", "MMP9")

# the ribosomal genes that come out as immature-subset markers
markers = fread(paste0(
  dir.results, "annotation/SCAHN.subset_categ.markers.csv"))
genes.ribosome = markers[avg_log2FC >= 1.5 & cluster == "immature"][
  grep("^RPS|^RPL", gene)]$gene

# The same signature panels as 01_load_data.R, repeated here because this script
# does not source it.
genes.immature.neu = c("AZU1", "BPI", "CAMP", "CEACAM6", "CEACAM8", "CRISP3",
                       "DEFA3", "DEFA4", "LCN2", "LTF", "MMP8", "MPO", "OLFM4")
genes.degranulating.neu = c("ARG1", "ANXA3", "CD177", "MMP9", "PGLYRP1")
genes.mature.neu = c("ALPL", "CMTM2", "DGAT2", "KCNJ15", "LRG1", "MGAM", "MME", "TNFRSF10C")
genes.antiprotease.neu = c("SLPI", "PI3")
genes.total.neu = c("ARG1", "AZU1", "BPI", "CAMP", "CEACAM6", "CEACAM8",
                    "CRISP3", "DEFA3", "DEFA4", "LCN2", "LTF", "MMP8", "MPO",
                    "OLFM4", "ALPL", "ANXA3", "CD177", "CMTM2", "CYP4F3",
                    "DGAT2", "KCNJ15", "LRG1", "MGAM", "MME", "MMP9",
                    "PGLYRP1", "PI3", "SLPI", "TNFRSF10C")
genes.ifn = c("APOL6", "EPSTI1", "GBP1", "GBP5", "HERC5", "HES4", "IFI44",
              "IFI44L", "IFI6", "IFIT1", "IFIT2", "IFIT3", "ISG15", "LY6E",
              "MT1X", "MT2A", "MX1", "OAS3", "OASL", "PARP14", "RSAD2", "STAT1", "XAF1")

genes = intersect(unique(c(genes.granules, genes.ribosome, genes.total.neu,
                           genes.ifn)), rownames(srt.sct@assays$RNA@data))
mat = srt.sct@assays$RNA@data[genes, ]
mat = as.matrix(mat) + 1

dtf$score.granules = expresso::get_gene_scores(mat, genes.granules, "")
dtf$score.ribosome = expresso::get_gene_scores(mat, genes.ribosome, "")

dtf$score.immature.neu = expresso::get_gene_scores(mat, genes.immature.neu, "")
dtf$score.degranulating.neu = expresso::get_gene_scores(
  mat, genes.degranulating.neu, "")
dtf$score.antiprotease.neu = expresso::get_gene_scores(
  mat, genes.antiprotease.neu, "")
dtf$score.mature.neu = expresso::get_gene_scores(mat, genes.mature.neu, "")
dtf$score.total.neu = expresso::get_gene_scores(mat, genes.total.neu, "")
dtf$score.ifn = expresso::get_gene_scores(mat, genes.ifn, "")

saveRDS(dtf[, c("cell.pid", "score.granules", "score.ribosome",
                "score.immature.neu", "score.degranulating.neu",
                "score.antiprotease.neu", "score.mature.neu", "score.total.neu", "score.ifn")],
        paste0(dir.results, "annotation/SCAHN.cellscores.rds"))


# -- Gene sets from published neutrophil studies -------------------------------

### mouse-to-human ortholog table
# Built once with biomaRt and cached; the manual file adds pairs biomaRt missed.
mouse_genes_human = readRDS(paste0(dir.data, "cross_species/mousegenes_tohuman.rds"))
mouse_genes_human = as.data.table(mouse_genes_human)
mouse_genes_human2 = fread(
  paste0(dir.data, "cross_species/mousegenes_tohuman.manual.csv"), header = FALSE)
colnames(mouse_genes_human) = c("mouse", "human")
colnames(mouse_genes_human2) = colnames(mouse_genes_human)
mouse_genes_human = rbind(mouse_genes_human, mouse_genes_human2)


### NeuMap, mouse
mouse_neumap_markers = read_xlsx(paste0(
  dir.data, "cross_species/NeuMAP_mouse_degs.xlsx"), sheet = 1)
mouse_neumap_markers = as.data.table(mouse_neumap_markers)
mouse_neumap_markers = mouse_neumap_markers[!is.na(cluster)]

mouse_neumap_markers$gene_human = mouse_genes_human[
  match(mouse_neumap_markers$gene, mouse_genes_human$mouse)]$human
unique(mouse_neumap_markers[is.na(gene_human)]$gene)
mouse_neumap_markers[!is.na(gene_human)][avg_log2FC >= 1.6][, .N]

genes.mouse_neumap_markers = unique(mouse_neumap_markers[
  !grepl("ENSG|^RPS|^RPL|^MT-|^C1Q", gene_human)][
  !is.na(gene_human) & avg_log2FC > 0.8][
  order(-avg_log2FC), head(.SD, 30), by = cluster]$gene_human)


### NeuMap, human
hneumap = read_xlsx(paste0(
  dir.data, "cross_species/hNeuMAP_integrated_hub_markers.xlsx"), sheet = 2)
hneumap = as.data.table(hneumap)
hneumap$cluster = gsub("hub_hs_", "", hneumap$cluster)
table(hneumap$cluster)

genes.hneumap_markers = unique(hneumap[
  !grepl("ENSG|^MT-|^RPS|^RPL|^IG[HLK]|JCHAIN|^C1Q", gene)][
  avg_log2FC > 0.8][order(-avg_log2FC), head(.SD, 30), by = cluster]$gene)

# Top 30 markers per cluster at log2FC > 1.2, then a second pass at > 0.6 for
# the clusters that pass nothing at the stricter cut.
genesets.NeuMap.mice = list()
for (i in unique(mouse_neumap_markers$cluster)) {
  tmp = mouse_neumap_markers[
    !grepl("ENSG|^RPS|^RPL|^MT-|^C1Q", gene_human)][
    !is.na(gene_human) & avg_log2FC > 1.2][
    order(-avg_log2FC), head(.SD, 30), by = cluster][cluster == i]
  if (nrow(tmp) > 0) {
    genesets.NeuMap.mice[[i]] = unique(tmp$gene_human)
  }
}

for (i in setdiff(unique(mouse_neumap_markers$cluster), names(genesets.NeuMap.mice))) {
  tmp = mouse_neumap_markers[
    !is.na(gene_human) & avg_log2FC > 0.6][
    order(-avg_log2FC), head(.SD, 30), by = cluster][cluster == i]
  if (nrow(tmp) > 0) {
    genesets.NeuMap.mice[[i]] = unique(tmp$gene_human)
  }
}

genesets.NeuMap.human = list()
for (i in unique(hneumap$cluster)) {
  tmp = hneumap[
    !grepl("ENSG|^MT-|^RPS|^RPL|^IG[HLK]|JCHAIN|^C1Q", gene)][
    !is.na(gene) & avg_log2FC > 1.2][
    order(-avg_log2FC), head(.SD, 30), by = cluster][cluster == i]
  if (nrow(tmp) > 0) {
    genesets.NeuMap.human[[i]] = unique(tmp$gene)
  }
}

for (i in setdiff(unique(hneumap$cluster), names(genesets.NeuMap.mice))) {
  tmp = hneumap[
    !grepl("ENSG|^MT-|^RPS|^RPL|^IG[HLK]|JCHAIN|^C1Q", gene)][
    !is.na(gene) & avg_log2FC > 0.6][
    order(-avg_log2FC), head(.SD, 30), by = cluster][cluster == i]
  if (nrow(tmp) > 0) {
    genesets.NeuMap.human[[i]] = unique(tmp$gene)
  }
}


### xie2020, mouse
df.marker = fread(paste0(dir.data, "cross_species/xie2020.genesets.csv"))
df.marker$gene_human = mouse_genes_human[
  match(df.marker$`gene symbol`, mouse_genes_human$mouse)]$human

genesets.xie2020.mice = list()
for (i in unique(df.marker$cluster)) {
  tmp = df.marker[!grepl("ENSG|^RPL|^RPS|^MT-", gene_human)][
    !is.na(gene_human) & avg_logFC > 1][
    order(-avg_logFC), head(.SD, 30), by = cluster][cluster == i]
  if (nrow(tmp) > 0) {
    genesets.xie2020.mice[[i]] = unique(tmp$gene_human)
  }
}

names(which(sapply(genesets.xie2020.mice, length) <= 1))

# second pass for the clusters that came out empty or near-empty
for (i in c(setdiff(unique(df.marker$cluster), names(genesets.xie2020.mice)),
            names(which(sapply(genesets.xie2020.mice, length) <= 2)))) {
  tmp = df.marker[!grepl("ENSG|^RPL|^RPS|^MT-", gene_human)][
    !is.na(gene_human) & avg_logFC > 0.6][
    order(-avg_logFC), head(.SD, 30), by = cluster][cluster == i]
  if (nrow(tmp) > 0) {
    genesets.xie2020.mice[[i]] = unique(tmp$gene_human)
  }
}

genesets.xie2020.mice$GM = NULL


### xie2020, human
df.marker = fread(paste0(dir.data, "cross_species/xie2020.genesets.human.csv"))
df.marker$`gene symbol`

genesets.xie2020.human = list()
for (i in unique(df.marker$cluster)) {
  tmp = df.marker[!grepl("ENSG|^RPL|^RPS|^MT-", `gene symbol`)][
    !is.na(`gene symbol`) & avg_logFC > 1][
    order(-avg_logFC), head(.SD, 30), by = cluster][cluster == i]
  if (nrow(tmp) > 0) {
    genesets.xie2020.human[[i]] = unique(tmp$`gene symbol`)
  }
}

for (i in c(setdiff(unique(df.marker$cluster), names(genesets.xie2020.human)),
            names(which(sapply(genesets.xie2020.human, length) <= 2)))) {
  tmp = df.marker[!grepl("ENSG|^RPL|^RPS|^MT-", `gene symbol`)][
    !is.na(`gene symbol`) & avg_logFC > 0.4][
    order(-avg_logFC), head(.SD, 30), by = cluster][cluster == i]
  if (nrow(tmp) > 0) {
    genesets.xie2020.human[[i]] = unique(tmp$`gene symbol`)
  }
}


### GSE165276 (grieshaber-bouyer 2021)
# The four progenitor stages are listed as mouse symbols and converted here;
# neutrotime is the published pseudotime correlation, already orthologous.
df.neutrotime = fread(paste0(dir.data, "cross_species/neutrotime.humangenes.csv"))

p1_genes = c("Ube2c", "Stmn1", "Tuba1b", "Ptma", "Ube2s", "Tubb5",
             "Chil3", "Cebpe", "Hmgb1", "H2afz", "Camp", "Hmgn2",
             "Hmgb2", "Ngp", "Ltf", "Arhgdib", "Calm2", "mt-Co1")
p2_genes = c("Lcn2", "Lyz2", "Ifitm6", "Wfdc21", "Anxa1", "Mmp8", "Cybb",
             "Dstn", "Ly6c2", "Ly6g", "Cd177", "Serpinb1a", "Prdx5", "Lgals3",
             "AA467197", "Mmp9", "Pglyrp1", "Mgst1")
p3_genes = c("Retnlg", "S100a6", "Prr13", "Tmcc1", "Gm5483")
p4_genes = c("Malat1", "Rps27", "Wfdc17", "Hbb-bs", "Il1b", "Ifitm1", "Hba-a1",
             "Hba-a2", "Dusp1", "Jund")

genesets.GSE165276 = list(
  P1 = setdiff(mouse_genes_human[match(p1_genes, mouse_genes_human$mouse)]$human, NA),
  P2 = setdiff(mouse_genes_human[match(p2_genes, mouse_genes_human$mouse)]$human, NA),
  P3 = setdiff(mouse_genes_human[match(p3_genes, mouse_genes_human$mouse)]$human, NA),
  P4 = setdiff(mouse_genes_human[match(p4_genes, mouse_genes_human$mouse)]$human, NA),
  neutrotime.pos = df.neutrotime[
    spearman >= 0.4 & hsapiens_homolog_orthology_confidence == 1
  ]$hsapiens_homolog_associated_gene_name, neutrotime.neg = df.neutrotime[
    spearman <= -0.4 & hsapiens_homolog_orthology_confidence == 1
  ]$hsapiens_homolog_associated_gene_name)


### ng2024, human tumour neutrophils
df.marker = fread(paste0(dir.data, "cross_species/ng2024.genesets.csv"))

genesets.ng2024 = list()
for (i in unique(df.marker$TumourNeutrophilCluster)) {
  tmp = df.marker[TumourNeutrophilCluster == i]
  if (nrow(tmp) > 0) {
    genesets.ng2024[[i]] = tmp$GeneID
  }
}


# -- Scores for the published gene sets ----------------------------------------

# dtf picked up the SCAHN score columns above, so start again from the metadata
dtf = fread(paste0(dir.results, "annotation/SCAHN.cellmeta.csv"))

allgenes = unique(unlist(c(genesets.NeuMap.mice, genesets.NeuMap.human,
                           genesets.xie2020.mice, genesets.xie2020.human,
                           genesets.GSE165276, genesets.ng2024)))

mat = srt.sct@assays$RNA@data[
  intersect(allgenes, rownames(srt.sct@assays$RNA@data)), ]
mat = as.matrix(mat) + 1

for (i in names(genesets.NeuMap.mice)) {
  dtf[[i]] = expresso::get_gene_scores(mat, genesets.NeuMap.mice[[i]], "")
}
for (i in names(genesets.NeuMap.human)) {
  dtf[[i]] = expresso::get_gene_scores(mat, genesets.NeuMap.human[[i]], "")
}
for (i in names(genesets.xie2020.mice)) {
  dtf[[i]] = expresso::get_gene_scores(mat, genesets.xie2020.mice[[i]], "")
}
for (i in names(genesets.xie2020.human)) {
  dtf[[i]] = expresso::get_gene_scores(mat, genesets.xie2020.human[[i]], "")
}
for (i in names(genesets.GSE165276)) {
  dtf[[i]] = expresso::get_gene_scores(mat, genesets.GSE165276[[i]], "")
}
for (i in names(genesets.ng2024)) {
  dtf[[i]] = expresso::get_gene_scores(mat, genesets.ng2024[[i]], "")
}

saveRDS(dtf[, c("cell.pid", "UMAP_1", "UMAP_2",
                names(genesets.NeuMap.mice), names(genesets.NeuMap.human),
                names(genesets.xie2020.mice), names(genesets.xie2020.human),
                names(genesets.GSE165276), names(genesets.ng2024)), with = FALSE],
        paste0(dir.results, "cross_species/SCAHN.neumap.scores.rds"))
