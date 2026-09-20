"""
Concordance of the selected clusterings
=======================================
Second of the three clustering-evaluation steps. Takes the
one-resolution-per-pipeline selection made by 05_run_cluster_aggregate.R and
scores how far those clusterings agree: pairwise ARI, NMI and VI per pair, a
per-cell stability score, and per-cluster robustness measured against one
reference clustering. The metrics themselves are in clustering_concordance.py;
this script wires them to the sweep outputs and then prints a few summaries off
the saved matrices.

The three steps run in order, each from its own script:
  1. 05_run_cluster_aggregate.R   the sweep table and the selection
  2. 06a_run_evalcluster.py        this script
  3. 07_plot_cluster_eval.R       the two summary figures

The reference clustering is fixed at v2000.PC25.algo1.res0.8, the mid-sweep
Louvain run. It has to be one of the selected columns, so if the selection in
clusternumber.sub.csv changes, this name has to change with it.

Usage:
  python 06a_run_evalcluster.py

Recorded output for the selection currently on disk:
  === Pairwise Metric Summary ===
    ARI: mean=0.4960, std=0.1296, min=0.3370, max=0.9658
    NMI: mean=0.6214, std=0.1124, min=0.4951, max=0.9616
    VI:  mean=2.3942, std=0.7220, min=0.2257, max=3.3883

  Computing cell-level stability scores...
    Strategy: co-clustering agreement with 30-NN
    Building k-NN index on 200,000 cells (Hamming on cluster labels)...
    Querying all 1,053,128 cells...
    Done in 2562.7s
    Stability: mean=0.9245, std=0.1103, min=0.2833, max=1.0000

Requires:
  numpy pandas, plus whatever clustering_concordance.py needs
  (scipy scikit-learn matplotlib seaborn tqdm joblib)

Reads:
  results/integration/clustering/clusterings.csv
  results/integration/clustering/clusternumber.sub.csv
    both written by 05_run_cluster_aggregate.R

Writes, all into results/integration/clustering/concordance_results/:
  pairwise_ARI.npy, pairwise_NMI.npy, pairwise_VI.npy
  stability_scores.npy, cluster_robustness.csv
  pairwise_{ARI,NMI}.png, pairwise_metric_distributions.png,
    cluster_robustness.png, cell_stability_histogram.png,
    resolution_vs_concordance.png
  cell_stability_scores.csv
    read by 07_plot_cluster_eval.R
"""

import numpy as np
import pandas as pd
from clustering_concordance import evaluate_clustering_concordance

DIR_RESULTS = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/results/"
DIR_CLUSTERING = DIR_RESULTS + "integration/clustering/"
DIR_CONCORD = DIR_CLUSTERING + "concordance_results"

REFERENCE = "v2000.PC25.algo1.res0.8"


# -- The selected clusterings --------------------------------------------------

df = pd.read_csv(f"{DIR_CLUSTERING}clusterings.csv")
nameselect = pd.read_csv(f"{DIR_CLUSTERING}clusternumber.sub.csv")

cols_to_keep = [c for c in nameselect["name"] if c in df.columns]
df = df[cols_to_keep]

clustering_matrix = df.values.astype(np.int32)
col_names = df.columns.tolist()


# -- Concordance ---------------------------------------------------------------

results = evaluate_clustering_concordance(
    clustering_matrix, clustering_names=col_names,
    reference_col=col_names.index(REFERENCE), k_neighbors=30, n_jobs=-1, output_dir=DIR_CONCORD,)

pd.Series(results["cell_stability"], name="stability_score").to_csv(
    f"{DIR_CONCORD}/cell_stability_scores.csv", index=False)


# -- Summaries off the saved matrices ------------------------------------------

ari_mat = np.load(f"{DIR_CONCORD}/pairwise_ARI.npy")
nmi_mat = np.load(f"{DIR_CONCORD}/pairwise_NMI.npy")

n = len(col_names)
i_idx, j_idx = np.triu_indices(n, k=1)

print("ARI - mean:", np.mean(ari_mat[i_idx, j_idx]), "median:", np.median(ari_mat[i_idx, j_idx]))
print("NMI - mean:", np.mean(nmi_mat[i_idx, j_idx]), "median:", np.median(nmi_mat[i_idx, j_idx]))

pairs_df = pd.DataFrame(
    {
        "clust_1": [col_names[i] for i in i_idx],
        "clust_2": [col_names[j] for j in j_idx],
        "ARI": ari_mat[i_idx, j_idx],
        "NMI": nmi_mat[i_idx, j_idx],
    })

# The extremes are what matter here: the worst-agreeing pair bounds how much the
# choice of pipeline can move the answer.
print("=== Bottom 10 ARI ===")
print(pairs_df.nsmallest(10, "ARI").to_string(index=False))
print("=== Top 10 ARI ===")
print(pairs_df.nlargest(10, "ARI").to_string(index=False))
print("=== Bottom 10 NMI ===")
print(pairs_df.nsmallest(10, "NMI").to_string(index=False))


# -- The same, with PC20 dropped -----------------------------------------------

# clusternumber.sub.csv should already exclude PC20, so this is a check that it
# did rather than a filter that is expected to remove anything.
mask = [i for i, name in enumerate(col_names) if "PC20" not in name]

ari_sub = ari_mat[np.ix_(mask, mask)]
nmi_sub = nmi_mat[np.ix_(mask, mask)]

ari_upper = ari_sub[np.triu_indices(len(mask), k=1)]
nmi_upper = nmi_sub[np.triu_indices(len(mask), k=1)]

print(f"Mean ARI (excl. PC20): {ari_upper.mean():.4f}")
print(f"Mean NMI (excl. PC20): {nmi_upper.mean():.4f}")
print(f"Median ARI (excl. PC20): {np.median(ari_upper):.4f}")
print(f"Median NMI (excl. PC20): {np.median(nmi_upper):.4f}")
