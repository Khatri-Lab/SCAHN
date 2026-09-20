# ==============================================================================
# Shared data: cell metadata, sample metadata, proportions, gene signatures
#
# Sourced after 00_setup.R by every figure script except 07_other.R, which uses
# none of these objects.
#
# Requires 00_setup.R (data.table, dir.data / dir.results) and SCAHN_SCRIPTS,
# used to source 00_genesets.R below.
#
# Reads:
#   data/samplemeta.csv
#   data/cellnumber.rds
#   results/annotation/SCAHN.cellmeta.csv
#   results/annotation/SCAHN.subsetproportion.rds
#
# Defines:
#   samplemeta       one row per sample, with sampleid.pid as the join key
#   dtf              one row per cell: cell metadata joined to sample metadata,
#                    restricted to samples with >= 30 cells, plus celltype2
#   dtf.prop         subset proportions, one row per retained sample
#   celltotalnumber  named table of total cells per sample, all cell types, so
#                    it is the denominator for the TNK / myeloid percentages
# and, by sourcing 00_genesets.R:
#   signaturegenes   marker panel shown across the figures
#   genes.*          neutrophil signature gene vectors
#   gene_sets        the same signature vectors as a named list
# ==============================================================================


# -- cell and sample metadata --------------------------------------------------

samplemeta = fread(paste0(dir.data, "samplemeta.csv"))
samplemeta$sampleid.pid = paste0(samplemeta$sampleid, ".", samplemeta$pid)
samplemeta[, source := ifelse(tissue %in% c("bone marrow", "cord blood"), "other", source)]

samplemeta[, group := ifelse(tissue == "cord blood" & pid == "myin2023", "healthy", group)]

dtf = fread(paste0(dir.results, "annotation/SCAHN.cellmeta.csv"))

# expand sample metadata to one row per cell, dropping the three key columns so
# the cbind does not duplicate them
samplemeta.m = samplemeta[match(dtf$sampleid.pid, samplemeta$sampleid.pid)]
samplemeta.m$pid = NULL
samplemeta.m$sampleid = NULL
samplemeta.m$sampleid.pid = NULL
dtf = cbind(dtf, samplemeta.m)

# keep only samples with at least 30 cells
samplestokeep = names(which(table(dtf$sampleid.pid) >= 30))
dtf = dtf[sampleid.pid %in% samplestokeep]

# coarse subset grouping
dtf[, celltype2 := ifelse(celltype %in% c("AZU1", "LTF"), "immature", "mature")]
dtf[, celltype2 := ifelse(celltype %in% c("IFN1", "IFN2", "IFN3"), "IFN", celltype2)]

dtf$pid = gsub("a$|b$|c$|d$", "", dtf$pid)

dtf.prop = readRDS(paste0(dir.results, "annotation/SCAHN.subsetproportion.rds"))
dtf.prop[, group := ifelse(tissue == "cord blood" & pid == "myin2023", "healthy", group)]

celltotalnumber = readRDS(paste0(dir.data, "cellnumber.rds"))

# -- gene panels ---------------------------------------------------------------

# signaturegenes, genes.* and gene_sets.
source(file.path(SCAHN_SCRIPTS, "00_genesets.R"), local = TRUE)
