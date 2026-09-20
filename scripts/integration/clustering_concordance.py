"""
Clustering Concordance Evaluation for Single-Cell Data
=======================================================
Evaluate stability and concordance across ~100 clusterings of ~1M cells.

Modules:
  1. Pairwise clustering comparison (ARI, NMI, VI) with heatmaps
  2. Cell-level stability scores via approximate co-clustering
  3. Per-cluster robustness scores
  4. Summary visualizations

Requirements:
  pip install numpy pandas scipy scikit-learn matplotlib seaborn tqdm joblib

Input:
  A pandas DataFrame or numpy array of shape (n_cells, n_clusterings),
  where each column is one clustering (integer labels).

Usage:
  See the __main__ block at the bottom for a full worked example with
  synthetic data (adaptable to real data).
"""

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
from joblib import Parallel, delayed
try:
    from tqdm import tqdm
except ImportError:
    # Simple fallback if tqdm is not installed
    def tqdm(iterable, desc="", **kwargs):
        items = list(iterable)
        total = len(items)
        for i, item in enumerate(items):
            if (i + 1) % max(1, total // 10) == 0 or i == total - 1:
                print(f"  {desc}: {i+1}/{total}")
            yield item
import warnings
import time
import os

# ──────────────────────────────────────────────────────────────────────
# 1. PAIRWISE CLUSTERING COMPARISON (ARI, NMI, Variation of Information)
# ──────────────────────────────────────────────────────────────────────

def variation_of_information(labels_a, labels_b):
    """Compute Variation of Information between two clusterings.
    VI(A,B) = H(A|B) + H(B|A). Lower = more similar. 0 = identical."""
    n = len(labels_a)
    # Build contingency
    classes_a = np.unique(labels_a)
    classes_b = np.unique(labels_b)
    # Map to contiguous ints for speed
    map_a = {v: i for i, v in enumerate(classes_a)}
    map_b = {v: i for i, v in enumerate(classes_b)}
    a_mapped = np.array([map_a[x] for x in labels_a])
    b_mapped = np.array([map_b[x] for x in labels_b])

    contingency = sparse.coo_matrix(
        (np.ones(n, dtype=np.int64), (a_mapped, b_mapped)), shape=(len(classes_a), len(classes_b))
    ).toarray()

    # Joint, marginals
    p_joint = contingency / n
    p_a = p_joint.sum(axis=1)
    p_b = p_joint.sum(axis=0)

    # H(A|B) = - sum p(a,b) log(p(a,b)/p(b))
    # H(B|A) = - sum p(a,b) log(p(a,b)/p(a))
    nz = p_joint > 0
    log_ratio_b = np.zeros_like(p_joint)
    log_ratio_a = np.zeros_like(p_joint)
    log_ratio_b[nz] = np.log(p_joint[nz] / np.broadcast_to(p_b[None, :], p_joint.shape)[nz])
    log_ratio_a[nz] = np.log(p_joint[nz] / np.broadcast_to(p_a[:, None], p_joint.shape)[nz])
    h_a_given_b = -np.sum(p_joint[nz] * log_ratio_b[nz])
    h_b_given_a = -np.sum(p_joint[nz] * log_ratio_a[nz])
    return h_a_given_b + h_b_given_a


def _compute_pair_metrics(args):
    """Worker function for parallel pairwise comparison."""
    i, j, labels_i, labels_j = args
    ari = adjusted_rand_score(labels_i, labels_j)
    nmi = normalized_mutual_info_score(labels_i, labels_j, average_method='arithmetic')
    vi = variation_of_information(labels_i, labels_j)
    return i, j, ari, nmi, vi


def pairwise_clustering_comparison(clustering_matrix, n_jobs=-1, max_pairs=None):
    """
    Compute ARI, NMI, and VI for all pairs of clusterings.

    Parameters
    ----------
    clustering_matrix : np.ndarray of shape (n_cells, n_clusterings)
    n_jobs : int, number of parallel workers (-1 = all cores)
    max_pairs : int or None, if set, randomly subsample pairs for speed

    Returns
    -------
    dict with keys 'ARI', 'NMI', 'VI', each a symmetric matrix (n_clusterings x n_clusterings)
    """
    n_cells, n_clust = clustering_matrix.shape
    print(f"Computing pairwise metrics for {n_clust} clusterings on {n_cells:,} cells...")

    # If cells are very numerous, subsample for pairwise metrics
    # (ARI/NMI on 1M cells is fine, just takes a moment per pair)
    all_pairs = list(combinations(range(n_clust), 2))
    if max_pairs and len(all_pairs) > max_pairs:
        rng = np.random.default_rng(42)
        idx = rng.choice(len(all_pairs), max_pairs, replace=False)
        all_pairs = [all_pairs[k] for k in idx]
        print(f"  Subsampled to {max_pairs} pairs")

    tasks = [
        (i, j, clustering_matrix[:, i], clustering_matrix[:, j]) for i, j in all_pairs]

    results = Parallel(n_jobs=n_jobs, verbose=5)(
        delayed(_compute_pair_metrics)(t) for t in tasks)

    # Assemble matrices
    ari_mat = np.eye(n_clust)
    nmi_mat = np.eye(n_clust)
    vi_mat = np.zeros((n_clust, n_clust))

    for i, j, ari, nmi, vi in results:
        ari_mat[i, j] = ari_mat[j, i] = ari
        nmi_mat[i, j] = nmi_mat[j, i] = nmi
        vi_mat[i, j] = vi_mat[j, i] = vi

    return {'ARI': ari_mat, 'NMI': nmi_mat, 'VI': vi_mat}


def plot_pairwise_heatmaps(metrics_dict, clustering_names=None, output_dir='.'):
    """Plot heatmaps for ARI, NMI, VI matrices."""
    for metric_name, mat in metrics_dict.items():
        fig, ax = plt.subplots(figsize=(10, 8))
        cmap = 'YlOrRd' if metric_name != 'VI' else 'YlOrRd_r'
        sns.heatmap(mat, cmap=cmap, ax=ax, square=True, xticklabels=clustering_names or False,
                    yticklabels=clustering_names or False)
        ax.set_title(f'Pairwise {metric_name} across clusterings')
        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, f'pairwise_{metric_name}.png'), dpi=150)
        plt.close(fig)
        print(f"  Saved pairwise_{metric_name}.png")


def summarize_pairwise(metrics_dict):
    """Print summary statistics of pairwise metrics."""
    print("\n=== Pairwise Metric Summary ===")
    for name, mat in metrics_dict.items():
        upper = mat[np.triu_indices_from(mat, k=1)]
        print(f"  {name}: mean={upper.mean():.4f}, std={upper.std():.4f}, "
              f"min={upper.min():.4f}, max={upper.max():.4f}")


# ──────────────────────────────────────────────────────────────────────
# 2. CELL-LEVEL STABILITY via APPROXIMATE CO-CLUSTERING
# ──────────────────────────────────────────────────────────────────────

def cell_stability_scores(clustering_matrix, k_neighbors=30, subsample_size=None,
                          n_jobs=-1, random_state=42):
    """
    Compute a per-cell stability score based on co-clustering with neighbors.

    Instead of building the full N×N co-clustering matrix (infeasible at 1M cells),
    we use a k-NN graph approach:
      1. For each clustering, build a co-cluster indicator for k nearest neighbors
         (based on a "consensus embedding" or simply the clustering label agreement).
      2. The stability score for cell i = average fraction of its k neighbors
         that co-cluster with it across all clusterings.

    The "neighbors" are defined by a co-clustering frequency graph:
    we first compute, for a subsample, how often pairs co-cluster, then use
    that to define neighborhoods.

    Faster alternative (implemented here):
      - For each cell, look at its row in the clustering_matrix.
      - Use a one-hot encoding of cluster labels to find nearest neighbors
        in "clustering agreement space".
      - Then, the stability score = mean agreement with those neighbors.

    Parameters
    ----------
    clustering_matrix : np.ndarray (n_cells, n_clusterings)
    k_neighbors : int
    subsample_size : int or None (if None, uses min(200_000, n_cells))
    n_jobs : int

    Returns
    -------
    stability_scores : np.ndarray of shape (n_cells,), values in [0, 1]
    """
    n_cells, n_clust = clustering_matrix.shape
    rng = np.random.default_rng(random_state)

    if subsample_size is None:
        subsample_size = min(200_000, n_cells)

    print(f"\nComputing cell-level stability scores...")
    print(f"  Strategy: co-clustering agreement with {k_neighbors}-NN")

    # Step 1: Build a compact representation of each cell's clustering profile.
    # For each clustering c, cell i has label clustering_matrix[i, c].
    # Two cells are "similar" if they co-cluster frequently.
    # We represent this as: for each clustering, encode the label as a column,
    # then compute Hamming-like agreement.

    # For efficiency with 1M cells and 100 clusterings, we directly use the
    # integer label matrix and compute agreement as fraction of matching labels.

    # Subsample for building the NN index
    if n_cells > subsample_size:
        sub_idx = rng.choice(n_cells, subsample_size, replace=False)
        sub_idx.sort()
    else:
        sub_idx = np.arange(n_cells)

    sub_matrix = clustering_matrix[sub_idx]

    # Build NN in "clustering profile space" using Hamming distance
    # sklearn's NearestNeighbors supports hamming via BallTree
    print(f"  Building k-NN index on {len(sub_idx):,} cells "
          f"(Hamming distance on cluster labels)...")
    t0 = time.time()
    nn = NearestNeighbors(n_neighbors=k_neighbors + 1, metric='hamming',
                          algorithm='brute', n_jobs=n_jobs)
    nn.fit(sub_matrix.astype(np.float32))

    # Query all cells against the subsampled index
    print(f"  Querying all {n_cells:,} cells...")
    distances, indices = nn.kneighbors(clustering_matrix.astype(np.float32))
    # distances are fraction of mismatching clusterings (Hamming)
    # stability = 1 - mean Hamming distance to neighbors (excluding self)
    # If a cell is in the subsample, its first neighbor is itself (distance 0).
    # For cells not in subsample, all neighbors are genuine.

    # Use columns 1: to skip self-match for subsampled cells,
    # and columns :k for non-subsampled cells.
    # Simplest: just take the mean of the k smallest distances excluding 0.
    # Since we requested k+1 neighbors, take columns 1: for all.
    neighbor_dists = distances[:, 1:]  # (n_cells, k_neighbors)
    stability_scores = 1.0 - neighbor_dists.mean(axis=1)

    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s")
    print(f"  Stability: mean={stability_scores.mean():.4f}, " f"std={stability_scores.std():.4f}, "
          f"min={stability_scores.min():.4f}, max={stability_scores.max():.4f}")

    return stability_scores


def cell_cocluster_frequency_sparse(clustering_matrix, k_neighbors=30,
                                     subsample_size=100_000, n_jobs=-1, random_state=42):
    """
    Build a SPARSE approximate co-clustering frequency matrix.

    For each cell, only store co-clustering frequencies with its k nearest
    neighbors (in clustering-profile space). Returns a sparse matrix.

    This is useful for downstream consensus clustering on a subsample.

    Parameters
    ----------
    clustering_matrix : np.ndarray (n_cells, n_clusterings)
    k_neighbors : int
    subsample_size : int

    Returns
    -------
    cocluster_freq : scipy.sparse.csr_matrix (subsample_size, subsample_size)
    sub_idx : np.ndarray, indices of subsampled cells
    """
    n_cells, n_clust = clustering_matrix.shape
    rng = np.random.default_rng(random_state)

    sub_idx = rng.choice(n_cells, min(subsample_size, n_cells), replace=False)
    sub_idx.sort()
    sub_matrix = clustering_matrix[sub_idx]
    n_sub = len(sub_idx)

    print(f"\nBuilding sparse co-clustering matrix on {n_sub:,} cells...")

    # For each pair of cells in the subsample, co-clustering freq =
    # fraction of clusterings where they share a label.
    # We only store k-NN entries.

    nn = NearestNeighbors(n_neighbors=k_neighbors + 1, metric='hamming',
                          algorithm='brute', n_jobs=n_jobs)
    nn.fit(sub_matrix)
    distances, indices = nn.kneighbors(sub_matrix)

    # Build sparse matrix
    rows = np.repeat(np.arange(n_sub), k_neighbors)
    cols = indices[:, 1:].ravel()
    vals = (1.0 - distances[:, 1:]).ravel()  # co-cluster freq = 1 - hamming

    cocluster_freq = sparse.csr_matrix((vals, (rows, cols)), shape=(n_sub, n_sub))
    # Symmetrize
    cocluster_freq = (cocluster_freq + cocluster_freq.T) / 2.0

    print(f"  Sparse matrix: {cocluster_freq.nnz:,} non-zero entries "
          f"({100 * cocluster_freq.nnz / n_sub**2:.4f}% dense)")

    return cocluster_freq, sub_idx


# ──────────────────────────────────────────────────────────────────────
# 3. PER-CLUSTER ROBUSTNESS
# ──────────────────────────────────────────────────────────────────────

def per_cluster_robustness(clustering_matrix, reference_col=0):
    """
    For each cluster in a reference clustering, compute how robustly
    its members stay together across other clusterings.

    The "cohesion score" for cluster C in reference clustering r:
      For each other clustering c, find the best-matching cluster in c
      (highest Jaccard overlap with C). The cohesion = mean of these
      best Jaccard scores across all other clusterings.

    The "fragmentation score": across other clusterings, how many
    distinct clusters do members of C get split into (entropy-based).

    Parameters
    ----------
    clustering_matrix : np.ndarray (n_cells, n_clusterings)
    reference_col : int, which clustering to use as reference

    Returns
    -------
    pd.DataFrame with columns: cluster, size, mean_jaccard, mean_entropy
    """
    n_cells, n_clust = clustering_matrix.shape
    ref_labels = clustering_matrix[:, reference_col]
    other_cols = [c for c in range(n_clust) if c != reference_col]

    clusters = np.unique(ref_labels)
    print(f"\nComputing per-cluster robustness for {len(clusters)} clusters "
          f"(reference clustering {reference_col})...")

    results = []
    for cl in tqdm(clusters, desc="  Clusters"):
        mask = ref_labels == cl
        size = mask.sum()
        jaccards = []
        entropies = []

        for c in other_cols:
            other_labels = clustering_matrix[:, c]
            # Labels of cells in this cluster under clustering c
            sub_labels = other_labels[mask]
            unique_sub, counts_sub = np.unique(sub_labels, return_counts=True)

            # Best Jaccard: for each unique label in sub_labels, compute
            # Jaccard with the reference cluster
            best_jac = 0
            for lbl, cnt in zip(unique_sub, counts_sub):
                # intersection = cnt
                # union = size + total_with_lbl - cnt
                total_with_lbl = (other_labels == lbl).sum()
                jac = cnt / (size + total_with_lbl - cnt)
                best_jac = max(best_jac, jac)
            jaccards.append(best_jac)

            # Entropy of label distribution within the cluster
            probs = counts_sub / counts_sub.sum()
            entropy = -np.sum(probs * np.log2(probs + 1e-15))
            entropies.append(entropy)

        results.append({
            'cluster': cl,
            'size': size,
            'mean_jaccard': np.mean(jaccards),
            'std_jaccard': np.std(jaccards),
            'mean_entropy': np.mean(entropies),
            'std_entropy': np.std(entropies),
        })

    df = pd.DataFrame(results)
    df = df.sort_values('mean_jaccard', ascending=False).reset_index(drop=True)
    print(f"  Most robust cluster: {df.iloc[0]['cluster']} "
          f"(Jaccard={df.iloc[0]['mean_jaccard']:.3f})")
    print(f"  Least robust cluster: {df.iloc[-1]['cluster']} "
          f"(Jaccard={df.iloc[-1]['mean_jaccard']:.3f})")
    return df


def plot_cluster_robustness(robustness_df, output_dir='.'):
    """Plot per-cluster robustness."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Jaccard
    ax = axes[0]
    ax.barh(range(len(robustness_df)), robustness_df['mean_jaccard'],
            xerr=robustness_df['std_jaccard'], color='steelblue', alpha=0.8)
    ax.set_yticks(range(len(robustness_df)))
    ax.set_yticklabels(robustness_df['cluster'], fontsize=7)
    ax.set_xlabel('Mean best Jaccard')
    ax.set_title('Cluster cohesion (higher = more robust)')
    ax.invert_yaxis()

    # Entropy
    ax = axes[1]
    ax.barh(range(len(robustness_df)), robustness_df['mean_entropy'],
            xerr=robustness_df['std_entropy'], color='coral', alpha=0.8)
    ax.set_yticks(range(len(robustness_df)))
    ax.set_yticklabels(robustness_df['cluster'], fontsize=7)
    ax.set_xlabel('Mean label entropy')
    ax.set_title('Cluster fragmentation (lower = more robust)')
    ax.invert_yaxis()

    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, 'cluster_robustness.png'), dpi=150)
    plt.close(fig)
    print(f"  Saved cluster_robustness.png")


# ──────────────────────────────────────────────────────────────────────
# 4. SUMMARY VISUALIZATIONS
# ──────────────────────────────────────────────────────────────────────

def plot_stability_histogram(stability_scores, output_dir='.'):
    """Histogram of per-cell stability scores."""
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(stability_scores, bins=100, color='steelblue', alpha=0.8, edgecolor='white')
    ax.axvline(stability_scores.mean(), color='red', ls='--',
               label=f'Mean = {stability_scores.mean():.3f}')
    ax.set_xlabel('Cell stability score')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of per-cell clustering stability')
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, 'cell_stability_histogram.png'), dpi=150)
    plt.close(fig)
    print(f"  Saved cell_stability_histogram.png")


def plot_metric_distributions(metrics_dict, output_dir='.'):
    """Violin/strip plots of pairwise metric distributions."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    for ax, (name, mat) in zip(axes, metrics_dict.items()):
        upper = mat[np.triu_indices_from(mat, k=1)]
        ax.hist(upper, bins=50, color='steelblue', alpha=0.8, edgecolor='white')
        ax.axvline(upper.mean(), color='red', ls='--', label=f'Mean = {upper.mean():.3f}')
        ax.set_title(f'{name} distribution')
        ax.set_xlabel(name)
        ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, 'pairwise_metric_distributions.png'), dpi=150)
    plt.close(fig)
    print(f"  Saved pairwise_metric_distributions.png")


def plot_n_clusters_vs_stability(clustering_matrix, metrics_dict, output_dir='.'):
    """Scatter: number of clusters in each clustering vs its mean ARI with others."""
    n_clust = clustering_matrix.shape[1]
    n_clusters_per = [len(np.unique(clustering_matrix[:, c])) for c in range(n_clust)]
    mean_ari_per = [metrics_dict['ARI'][c, :].sum() / (n_clust - 1) for c in range(n_clust)]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(n_clusters_per, mean_ari_per, alpha=0.6, s=30, c='steelblue')
    ax.set_xlabel('Number of clusters')
    ax.set_ylabel('Mean ARI with other clusterings')
    ax.set_title('Resolution vs. concordance')
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, 'resolution_vs_concordance.png'), dpi=150)
    plt.close(fig)
    print(f"  Saved resolution_vs_concordance.png")


# ──────────────────────────────────────────────────────────────────────
# 5. FULL PIPELINE
# ──────────────────────────────────────────────────────────────────────

def evaluate_clustering_concordance(clustering_matrix, clustering_names=None,
                                     reference_col=0, k_neighbors=30,
                                     n_jobs=-1, output_dir='concordance_results'):
    """
    Full pipeline: run all concordance evaluations.

    Parameters
    ----------
    clustering_matrix : np.ndarray of shape (n_cells, n_clusterings)
        Integer cluster labels. Each column is one clustering.
    clustering_names : list of str or None
    reference_col : int, which column to use as reference for per-cluster robustness
    k_neighbors : int, for cell-level stability
    n_jobs : int
    output_dir : str

    Returns
    -------
    dict with all results
    """
    os.makedirs(output_dir, exist_ok=True)
    n_cells, n_clust = clustering_matrix.shape
    print(f"{'='*60}")
    print(f"CLUSTERING CONCORDANCE EVALUATION")
    print(f"  {n_cells:,} cells × {n_clust} clusterings")
    print(f"  Output directory: {output_dir}")
    print(f"{'='*60}")

    # 1. Pairwise metrics
    print(f"\n--- Step 1: Pairwise clustering comparison ---")
    metrics = pairwise_clustering_comparison(clustering_matrix, n_jobs=n_jobs)
    summarize_pairwise(metrics)
    plot_pairwise_heatmaps(metrics, clustering_names, output_dir)
    plot_metric_distributions(metrics, output_dir)
    plot_n_clusters_vs_stability(clustering_matrix, metrics, output_dir)

    # 2. Cell-level stability
    print(f"\n--- Step 2: Cell-level stability ---")
    stability = cell_stability_scores(clustering_matrix, k_neighbors=k_neighbors, n_jobs=n_jobs)
    plot_stability_histogram(stability, output_dir)

    # 3. Per-cluster robustness
    print(f"\n--- Step 3: Per-cluster robustness ---")
    robustness = per_cluster_robustness(clustering_matrix, reference_col=reference_col)
    plot_cluster_robustness(robustness, output_dir)

    # Save numerical results
    np.save(os.path.join(output_dir, 'stability_scores.npy'), stability)
    robustness.to_csv(os.path.join(output_dir, 'cluster_robustness.csv'), index=False)
    for name, mat in metrics.items():
        np.save(os.path.join(output_dir, f'pairwise_{name}.npy'), mat)

    print(f"\n{'='*60}")
    print(f"DONE — all results saved to {output_dir}/")
    print(f"{'='*60}")

    return {
        'pairwise_metrics': metrics,
        'cell_stability': stability,
        'cluster_robustness': robustness,
    }


# ──────────────────────────────────────────────────────────────────────
# DEMO / USAGE EXAMPLE
# ──────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    # -----------------------------------------------------------
    # Generate synthetic data mimicking 100 clusterings of 1M cells.
    # Replace this section with your real data.
    # -----------------------------------------------------------
    print("Generating synthetic clustering data for demo...")
    np.random.seed(42)

    N_CELLS = 50_000         # Demo size; set to 1_000_000 for real data
    N_CLUSTERINGS = 20       # Demo size; set to 100 for real data
    # NOTE: At 1M cells × 100 clusterings, expect:
    #   - Pairwise metrics: ~30 min (4,950 pairs × ~0.4s each)
    #   - Cell stability: ~20 min (k-NN on 1M × 100 features)
    #   - Per-cluster robustness: ~10 min (depends on # clusters)
    # Total: ~1 hour on a 16-core machine. Highly parallelizable.
    N_TRUE_CLUSTERS = 15

    # Ground truth
    true_labels = np.random.randint(0, N_TRUE_CLUSTERS, size=N_CELLS)

    # Generate clusterings with varying noise levels and resolutions
    clustering_matrix = np.zeros((N_CELLS, N_CLUSTERINGS), dtype=np.int32)
    names = []
    for c in range(N_CLUSTERINGS):
        noise_rate = np.random.uniform(0.05, 0.40)  # 5-40% label noise
        n_clusters = np.random.choice([10, 12, 15, 18, 20, 25])
        labels = true_labels.copy()
        # Add noise
        noise_mask = np.random.rand(N_CELLS) < noise_rate
        labels[noise_mask] = np.random.randint(0, n_clusters, size=noise_mask.sum())
        # Possibly merge/split clusters
        if n_clusters < N_TRUE_CLUSTERS:
            labels = labels % n_clusters
        clustering_matrix[:, c] = labels
        names.append(f"run_{c}_k{n_clusters}_n{noise_rate:.0%}")

    print(f"  Generated matrix: {clustering_matrix.shape}")

    # -----------------------------------------------------------
    # Run the full evaluation pipeline
    # -----------------------------------------------------------
    results = evaluate_clustering_concordance(
        clustering_matrix, clustering_names=names, reference_col=0, k_neighbors=30,
        n_jobs=-1, output_dir='concordance_results')

    # -----------------------------------------------------------
    # HOW TO USE WITH YOUR REAL DATA:
    # -----------------------------------------------------------
    #
    # import anndata as ad
    #
    # # Option A: clusterings stored in adata.obs
    # adata = ad.read_h5ad('my_data.h5ad')
    # clustering_cols = [c for c in adata.obs.columns if c.startswith('cluster_')]
    # # Convert to integer matrix
    # from sklearn.preprocessing import LabelEncoder
    # matrix = np.zeros((adata.n_obs, len(clustering_cols)), dtype=np.int32)
    # for i, col in enumerate(clustering_cols):
    #     le = LabelEncoder()
    #     matrix[:, i] = le.fit_transform(adata.obs[col].values)
    #
    # results = evaluate_clustering_concordance(
    #     matrix,
    #     clustering_names=clustering_cols,
    #     reference_col=0,  # pick your "best guess" reference
    #     k_neighbors=30,
    #     n_jobs=-1,
    #     output_dir='concordance_results'
    # )
    #
    # # Add stability scores back to adata
    # adata.obs['clustering_stability'] = results['cell_stability']
    # # Now you can plot stability on UMAP:
    # import scanpy as sc
    # sc.pl.umap(adata, color='clustering_stability', cmap='RdYlGn')
