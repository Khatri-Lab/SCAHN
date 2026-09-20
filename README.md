# SCAHN — Single-cell Atlas of Human Neutrophils

Analysis code for *A single-cell atlas of human neutrophils reveals neutrophil
heterogeneity across health and disease*.

> Published in **Nature Immunology** (2026).
> DOI: [10.1038/s41590-026-02671-8](https://doi.org/10.1038/s41590-026-02671-8)

## Layout

```
scripts/
  integration/             atlas construction: SCTransform, Harmony, SNN,
                           consensus clustering, annotation and the Symphony
                           reference
  analysis/                downstream analyses: signature scoring and
                           pseudobulk, single-cell meta-analyses (infection,
                           cancer, sex difference), bulk sex-difference
                           analyses, Framingham processing, trajectory. Other
                           analysis code is in figures/ or
                           figures_post_acceptance/
  cross_species/           NeuMap ortholog mapping and cross-species Symphony
                           projection
  figures/                 original figure scripts
  figures_post_acceptance/ figure scripts as revised after acceptance
```

### Running the figure scripts

`figures/00_setup.R` (packages, palettes) and `figures/01_load_data.R` (shared
objects and signatures) are sourced by the rest; the numbered scripts `02`–`07`
are otherwise independent of one another. `figures_post_acceptance/` is the same
pipeline with the post-acceptance revisions.

## Environment

R 4.3.1 and Python 3.9.

This code depends on **[expresso](https://github.com/zhengh42/expresso)**, an R
package by the same author that supplies the plotting theme (`theme_expresso()`)
and the meta-analysis helpers used throughout. Install it first:

```r
remotes::install_github("zhengh42/expresso")
```

Other R dependencies: Seurat (4.4.0), harmony (1.2.0), symphony (0.1.2),
slingshot (2.10.0), uwot (0.1.16), BiocNeighbors (1.20.2), DESeq2 (1.40.2),
MetaIntegrator (2.1.3), fgsea (1.26.0), ComplexHeatmap (2.18.0), circlize
(0.4.18), data.table (1.18.4), Matrix (1.6-1.1), matrixStats (1.5.0), survival
(3.8-6), margins (0.3.28), magick (2.9.1), and the ggplot2 (4.0.3) ecosystem
(ggpubr, ggrepel, ggrastr, ggalluvial, ggforce, ggbio, ggfortify, ggstance,
ggtext, ggsci, ggthemes, ggbeeswarm, ggsignif, gridExtra, gtable, grid,
gridtext, scales, viridis, dichromat, pBrackets, reshape2, tibble).

Python dependencies: numpy, pandas, scipy, scikit-learn, leidenalg,
python-igraph, matplotlib, seaborn, joblib, tqdm.

## License

MIT — see [LICENSE](LICENSE).
