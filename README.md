# SCAHN — Single-cell Atlas of Human Neutrophils

Analysis code for *A single-cell atlas of human neutrophils reveals neutrophil
heterogeneity across health and disease*.

> Published in **Nature Immunology** (2026).
> DOI: [10.1038/s41590-026-02671-8](https://doi.org/10.1038/s41590-026-02671-8)

## Abstract

Neutrophils, the most abundant leukocytes in human blood, have multifaceted
roles in health and disease. Here we assembled a single-cell atlas of human
neutrophils (SCAHN) by integrating transcriptome data from 1 million human
neutrophils from 781 samples and 43 independent datasets. Using SCAHN, we
identified multiple neutrophil subsets with distinct functions across tissues
and robust changes in gene expression and subset proportions in infection,
cancer, immune suppression and sex differences at single-cell resolution.
Cross-species reference mapping of mouse and nonhuman primate neutrophils onto
SCAHN revealed both conserved and species-specific features. Using 17,023 bulk
blood transcriptome profiles from 40 datasets, we found that neutrophil-specific
signatures predicted mortality in acutely ill patients at presentation, were
associated with cancer diagnosis and survival, and correlated with metabolic
dysfunction and long-term all-cause mortality. Our study highlights how
neutrophil biology can improve immune risk stratification and develop effective
diagnostic and therapeutic strategies.

## Web app

Explore the atlas online at **<https://scahn.khatrilab.org/>**. Enter a gene or a gene list to see where it is expressed across neutrophil subsets, tissues and conditions:

- **UMAP**: expression on the atlas UMAP, one panel per gene, plus a combined
  score for a gene list
- **Heatmap**: dot plot and heatmap of mean expression and % expressing across
  subsets, tissues, conditions or studies
- **FAQ**: how the numbers behind the plots are computed

Both plot pages can restrict the cells to a tissue, condition or study, and
every figure downloads as PNG and PDF.

## Layout

```
scripts/
  analysis/                downstream analyses: signature scoring and
                           pseudobulk, single-cell meta-analyses (infection,
                           cancer, sex difference), bulk sex-difference
                           analyses, Framingham processing, trajectory. Other
                           analysis code is in figures_post_acceptance/
  figures_post_acceptance/ figure scripts as revised after acceptance
```

## Environment

R 4.3.1.

R dependencies: [expresso](https://github.com/zhengh42/expresso), Seurat
(4.4.0), slingshot (2.10.0), DESeq2 (1.40.2), MetaIntegrator (2.1.3), fgsea
(1.26.0), ComplexHeatmap (2.18.0), circlize (0.4.18), data.table (1.18.4),
Matrix (1.6-1.1), matrixStats (1.5.0), survival (3.8-6), margins (0.3.28),
magick (2.9.1), car, Hmisc, readxl, future, and the ggplot2 (4.0.3) ecosystem
(ggpubr, ggrepel, ggrastr, ggalluvial, ggforce, ggbio, ggfortify, ggstance,
ggtext, ggsci, ggthemes, ggbeeswarm, ggsignif, gridExtra, gtable, grid,
gridtext, scales, viridis, dichromat, pBrackets, reshape2).

## License

MIT — see [LICENSE](LICENSE).
