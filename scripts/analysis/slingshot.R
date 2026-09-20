# ==============================================================================
# Slingshot trajectory and pseudotime over the harmony embedding
#
# Fits principal curves with slingshot on a 50,000-cell subsample of the 25-PC
# harmony embedding, using a coarse regrouping of the neutrophil subsets as the
# cluster labels and AZU1, the most immature subset, as the start cluster. The
# curves are then projected onto every non-lowDepth cell with predict(), and
# separately embedded into UMAP space so the figure script can draw them.
#
# Requires 00_setup.R (data.table, dir.results).
#
# Reads:
#   results/annotation/SCAHN.cellmeta.csv
#   results/integration/harmony/SCAHN.harmony.Z_corr.v2000.PC25.csv
#     written by 02a_run_harmony_SCT.R
#
# Writes:
#   results/annotation/SCAHN.slingshot.res.rds    read by 02_atlas.R
# ==============================================================================

SCAHN_SCRIPTS = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/scripts/figures"
source(file.path(SCAHN_SCRIPTS, "00_setup.R"))  # packages, paths

library(slingshot)


# -- Cells and embedding -------------------------------------------------------

dtf = fread(paste0(dir.results, "annotation/SCAHN.cellmeta.csv"))

emb = fread(paste0(dir.results, "integration/harmony/SCAHN.harmony.Z_corr.v2000.PC25.csv"))
emb = as.matrix(emb)
rownames(emb) = dtf$cell.pid
table(dtf$celltype)

dtf.s = dtf[!celltype %in% c("lowDepth")]


# -- Coarse cluster labels for the fit -----------------------------------------

dtf.s[, celltype2 := ifelse(grepl("IFN", celltype), "IFN", celltype)]
dtf.s[, celltype2 := ifelse(grepl("IL1RN|IL1B|NF-κB", celltype), "NF-κB", celltype2)]
dtf.s[, celltype2 := ifelse(grepl("S100A4|MMP9", celltype), "S100", celltype2)]
dtf.s[, celltype2 := ifelse(grepl("AZU|LTF|MME|SLPI|S100|IFN|NF-κB|CCL|CXCL|VEGFA|G0S2",
                                  celltype2), celltype2, "other")]

# -- Slingshot fit on a subsample ----------------------------------------------

set.seed(42)
row2keep = sample(1:nrow(dtf.s), 50000)
dtf.ss = dtf.s[row2keep, ]

emb.ss = emb[dtf.ss$cell.pid, ]

set.seed(42)
sds = slingshot(emb.ss[, 1:25], clusterLabels = dtf.ss$celltype2, start.clus = "AZU1",
                dist.method = "slingshot", allow.breaks = TRUE)

pst = slingshot::slingPseudotime(sds)
curves = slingshot::slingCurves(sds)


# -- Project onto all cells and embed the curves in UMAP space -----------------

allPseudotime = predict(sds, newdata = emb[dtf.s$cell.pid, 1:25])
pseudotime_matrix = slingPseudotime(allPseudotime)

sds.embedded = embedCurves(sds, as.matrix(dtf.ss[, c("UMAP_1", "UMAP_2")]), stretch = 1,
                           approx_points = 500, shrink = 1, shrink.method = "cosine")

curves_emb = slingshot::slingCurves(sds.embedded)
names(curves_emb) = names(curves)

slingshotlist = list(curves_emb = curves_emb, pseudotime_matrix = pseudotime_matrix)
saveRDS(slingshotlist, paste0(dir.results, "annotation/SCAHN.slingshot.res.rds"))
