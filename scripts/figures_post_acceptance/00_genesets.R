# ==============================================================================
# Gene panels: the figure marker panel and the neutrophil signature vectors
#
# Requires nothing. 01_load_data.R sources this, so any script sourcing that
# gets these too.
#
# Defines:
#   signaturegenes   marker panel shown across the figures
#   genes.*          neutrophil signature gene vectors
#   gene_sets        the same signature vectors as a named list
# ==============================================================================


# -- genes shown in the figures ------------------------------------------------

signaturegenes = c("STMN1", "MKI67", "FUT4", "PRTN3", "CTSG", "ELANE",
                   "AZU1", "MPO",  "DEFA4", "DEFA3", "BPI", "CEACAM6", "CEACAM8",
                   "OLFM4", "LTF", "LCN2", "MMP8", "CAMP", "CRISP3", "TCN1", "HP", "FCN1", "PADI4",
                   "PGLYRP1", "ANXA3", "ARG1", "CD177", "MMP9",
                   "NQO2", "CYP4F3", "S100A8", "S100A9", "S100A12", "S100P", "TXN", "TSPO",
                   "PLAC8", "GPI", "OLR1",  "ORM1", "S100A4", "IL1R2", "MME", "TXNIP",
                   "EGR1", "FOS", "FOSB", "JUN", "PTGS2", "G0S2", "MT2A",  "HES4",
                   "LY6E", "ISG15", "IFI44", "OASL", "IFI6", "OAS3",
                   "RSAD2", "HERC5", "MX1", "IFIT1", "IFIT2",  "IFIT5",
                   "DDX58", "XAF1", "IRF7", "EPSTI1", "APOL6", "GBP5", "GBP4", "GBP2",
                   "GBP1", "PARP14", "TNFSF10", "IFITM3", "CD274", "CXCL1", "CXCL2", "CXCL8",
                   "PPIF", "SPP1", "CDKN1A", "CD83", "SQSTM1", "PLAU", "CXCR4",
                   "CSTB", "LGALS3", "HMOX1", "VEGFA", "CCL3", "CCL4",
                   "TNFAIP3", "NFKBIA", "IL1B", "IL1RN", "CD69", "SLPI", "PI3", "SIGLEC10",
                   "HSPA1B", "DNAJB1", "CD74", "HLA-DRA", "ADM",
                   "AQP9", "BCL6", "CR1", "APOBEC3A", "PTEN", "SELL",
                   "TNFRSF10C", "ALPL", "CMTM2", "DGAT2", "KCNJ15", "LRG1", "MGAM",
                   "CSF3R", "NAMPT", "FCGR3B", "CXCR2")


# -- neutrophil signature definitions ------------------------------------------

genes.immature.neu      = c("AZU1", "BPI", "CAMP", "CEACAM6", "CEACAM8",
                            "CRISP3", "DEFA3", "DEFA4", "LCN2", "LTF", "MMP8", "MPO", "OLFM4")
genes.degranulating.neu = c("ARG1", "ANXA3", "CD177", "MMP9", "PGLYRP1")
genes.antiprotease.neu  = c("SLPI", "PI3")
genes.mature.neu = c("ALPL", "CMTM2", "DGAT2", "KCNJ15", "LRG1", "MGAM", "MME", "TNFRSF10C")
genes.ifn               = c("APOL6", "EPSTI1", "GBP1", "GBP5", "HERC5", "HES4",
                            "IFI44", "IFI44L", "IFI6", "IFIT1", "IFIT2",
                            "IFIT3", "ISG15", "LY6E", "MT1X", "MT2A", "MX1",
                            "OAS3", "OASL", "PARP14", "RSAD2", "STAT1", "XAF1")

# all four subset panels, plus CYP4F3
genes.total.neu = unique(c(genes.immature.neu, genes.degranulating.neu,
                           genes.mature.neu, genes.antiprotease.neu, "CYP4F3"))

# the same vectors as a named list
gene_sets = list(
  immature = genes.immature.neu, degranulating = genes.degranulating.neu,
  antiprotease = genes.antiprotease.neu, mature = genes.mature.neu,
  total = genes.total.neu, IFN = genes.ifn)
