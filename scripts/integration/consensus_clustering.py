"""
Consensus Clustering for Single-Cell Data (1M+ cells)
======================================================
Given a matrix of (n_cells, n_clusterings), produce a single consensus
clustering that reflects the agreement across all input clusterings.

Strategy (scalable to 1M+ cells):
  1. Subsample cells (e.g., 100-200K)
  2. Build a weighted co-clustering k-NN graph on the subsample
     (edge weight = fraction of clusterings where two cells co-cluster)
  3. Cluster the co-clustering graph using Leiden
  4. Train a k-NN classifier to propagate consensus labels to all cells
  5. Assign confidence scores to each cell

Requirements:
  pip install numpy pandas scipy scikit-learn leidenalg igraph

Usage:
  from consensus_clustering import consensus_clustering
  labels, confidence = consensus_clustering(clustering_matrix)
"""

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.neighbors import NearestNeighbors, KNeighborsClassifier
import time
import os


def _build_cocluster_graph(clustering_matrix_sub, k_neighbors=30, n_jobs=-1):
    """
    Build a weighted k-NN graph where edge weights represent co-clustering
    frequency (1 - Hamming distance).
    """
    n_sub, n_clust = clustering_matrix_sub.shape
    print(f"  Building co-clustering k-NN graph ({n_sub:,} cells, k={k_neighbors})...")
    t0 = time.time()

    nn = NearestNeighbors(
        n_neighbors=k_neighbors + 1, metric='hamming', algorithm='brute', n_jobs=n_jobs)
    nn.fit(clustering_matrix_sub.astype(np.float32))
    distances, indices = nn.kneighbors(clustering_matrix_sub.astype(np.float32))

    # Convert Hamming distance to co-clustering frequency, skip self (col 0)
    rows = np.repeat(np.arange(n_sub), k_neighbors)
    cols = indices[:, 1:].ravel()
    weights = (1.0 - distances[:, 1:]).ravel()

    mask = weights > 0
    rows, cols, weights = rows[mask], cols[mask], weights[mask]

    adj = sparse.csr_matrix((weights, (rows, cols)), shape=(n_sub, n_sub))
    adj = adj.maximum(adj.T)

    elapsed = time.time() - t0
    print(f"  Graph built in {elapsed:.1f}s: {adj.nnz:,} edges")
    return adj


def _cluster_graph_leiden(adjacency, resolution=1.0):
    """Cluster the co-clustering graph using Leiden algorithm."""
    import leidenalg
    import igraph as ig

    print(f"  Running Leiden clustering (resolution={resolution})...")
    t0 = time.time()

    coo = sparse.triu(adjacency).tocoo()
    edges = list(zip(coo.row.tolist(), coo.col.tolist()))
    weights = coo.data.tolist()

    g = ig.Graph(n=adjacency.shape[0], edges=edges, directed=False)
    g.es['weight'] = weights

    partition = leidenalg.find_partition(
        g, leidenalg.RBConfigurationVertexPartition,
        weights='weight', resolution_parameter=resolution, n_iterations=-1, seed=42)

    labels = np.array(partition.membership)
    elapsed = time.time() - t0
    n_clusters = len(np.unique(labels))
    print(f"  Leiden found {n_clusters} consensus clusters in {elapsed:.1f}s")
    return labels


def _propagate_labels(clustering_matrix_full, sub_idx, sub_labels, k_neighbors=30, n_jobs=-1):
    """
    Propagate consensus labels from subsample to all cells using
    k-NN classifier in clustering-profile space.
    Returns full labels and per-cell confidence.
    """
    n_cells = clustering_matrix_full.shape[0]
    n_sub = len(sub_idx)
    print(f"  Propagating labels from {n_sub:,} subsampled to {n_cells:,} total cells...")
    t0 = time.time()

    knn = KNeighborsClassifier(
        n_neighbors=min(k_neighbors, n_sub - 1), metric='hamming',
        algorithm='brute', weights='distance', n_jobs=n_jobs)
    knn.fit(
        clustering_matrix_full[sub_idx].astype(np.float32), sub_labels)

    full_labels = knn.predict(clustering_matrix_full.astype(np.float32))
    proba = knn.predict_proba(clustering_matrix_full.astype(np.float32))
    confidence = proba.max(axis=1)

    # Subsampled cells keep their original labels
    full_labels[sub_idx] = sub_labels

    elapsed = time.time() - t0
    print(f"  Propagation done in {elapsed:.1f}s")
    print(f"  Confidence: mean={confidence.mean():.4f}, "
          f"min={confidence.min():.4f}, max={confidence.max():.4f}")

    return full_labels, confidence


def consensus_clustering(clustering_matrix, subsample_size=200_000,
                         k_neighbors_graph=30, k_neighbors_propagate=30,
                         leiden_resolution=1.0, n_jobs=-1, random_state=42, output_dir=None):
    """
    Compute consensus clustering from multiple clusterings.

    Parameters
    ----------
    clustering_matrix : np.ndarray of shape (n_cells, n_clusterings)
        Integer cluster labels. Each column is one clustering.
    subsample_size : int
        Number of cells to subsample for building the co-clustering graph.
        200K is a good default for 1M cells.
    k_neighbors_graph : int
        Number of neighbors for the co-clustering graph.
    k_neighbors_propagate : int
        Number of neighbors for label propagation.
    leiden_resolution : float
        Resolution for Leiden clustering. Higher = more clusters.
    n_jobs : int
        Parallel workers (-1 = all cores).
    random_state : int
        Random seed.
    output_dir : str or None
        If provided, save results and plots here.

    Returns
    -------
    dict with keys:
        'labels' : np.ndarray (n_cells,) — consensus cluster labels
        'confidence' : np.ndarray (n_cells,) — per-cell confidence [0, 1]
        'sub_idx' : np.ndarray — indices of subsampled cells
        'sub_labels' : np.ndarray — consensus labels on subsample
        'cluster_sizes' : pd.DataFrame
    """
    n_cells, n_clust = clustering_matrix.shape
    rng = np.random.default_rng(random_state)

    print(f"{'='*60}")
    print(f"CONSENSUS CLUSTERING")
    print(f"  {n_cells:,} cells x {n_clust} clusterings")
    print(f"  Subsample: {min(subsample_size, n_cells):,}, "
          f"Leiden resolution: {leiden_resolution}")
    print(f"{'='*60}")

    # Step 1: Subsample
    if n_cells > subsample_size:
        sub_idx = rng.choice(n_cells, subsample_size, replace=False)
        sub_idx.sort()
    else:
        sub_idx = np.arange(n_cells)

    sub_matrix = clustering_matrix[sub_idx]
    print(f"\n--- Step 1: Subsample {len(sub_idx):,} cells ---")

    # Step 2: Build co-clustering graph
    print(f"\n--- Step 2: Build co-clustering graph ---")
    adjacency = _build_cocluster_graph(sub_matrix, k_neighbors=k_neighbors_graph, n_jobs=n_jobs)

    # Step 3: Cluster the graph with Leiden
    print(f"\n--- Step 3: Cluster co-clustering graph ---")
    sub_labels = _cluster_graph_leiden(adjacency, resolution=leiden_resolution)

    # Step 4: Propagate to all cells
    print(f"\n--- Step 4: Propagate labels to all cells ---")
    if n_cells > subsample_size:
        full_labels, confidence = _propagate_labels(
            clustering_matrix, sub_idx, sub_labels, k_neighbors=k_neighbors_propagate, n_jobs=n_jobs
        )
    else:
        full_labels = sub_labels
        confidence = np.ones(n_cells)

    # Summary
    n_consensus = len(np.unique(full_labels))
    print(f"\n{'='*60}")
    print(f"CONSENSUS RESULT: {n_consensus} clusters")
    print(f"  Confidence: mean={confidence.mean():.3f}, " f"median={np.median(confidence):.3f}")
    low_conf = (confidence < 0.5).sum()
    print(f"  Low-confidence cells (<0.5): {low_conf:,} ({100*low_conf/n_cells:.1f}%)")
    print(f"{'='*60}")

    # Cluster sizes
    unique, counts = np.unique(full_labels, return_counts=True)
    size_df = pd.DataFrame({'cluster': unique, 'size': counts,
                            'pct': np.round(100 * counts / n_cells, 2)})
    size_df = size_df.sort_values('size', ascending=False).reset_index(drop=True)
    print(f"\nCluster sizes:")
    print(size_df.to_string(index=False))

    # Save outputs
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        np.save(os.path.join(output_dir, 'consensus_labels.npy'), full_labels)
        np.save(os.path.join(output_dir, 'consensus_confidence.npy'), confidence)
        size_df.to_csv(os.path.join(output_dir, 'consensus_cluster_sizes.csv'), index=False)

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 4))
        ax.hist(confidence, bins=100, color='steelblue', alpha=0.8, edgecolor='white')
        ax.axvline(confidence.mean(), color='red', ls='--', label=f'Mean = {confidence.mean():.3f}')
        ax.set_xlabel('Consensus confidence')
        ax.set_ylabel('Count')
        ax.set_title('Per-cell consensus clustering confidence')
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, 'consensus_confidence_histogram.png'), dpi=150)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 4))
        ax.bar(range(len(size_df)), size_df['size'], color='steelblue', alpha=0.8)
        ax.set_xlabel('Consensus cluster')
        ax.set_ylabel('Number of cells')
        ax.set_title(f'Consensus clustering: {n_consensus} clusters')
        ax.set_xticks(range(len(size_df)))
        ax.set_xticklabels(size_df['cluster'], fontsize=7, rotation=45)
        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, 'consensus_cluster_sizes.png'), dpi=150)
        plt.close(fig)

        print(f"\nResults saved to {output_dir}/")

    return {
        'labels': full_labels,
        'confidence': confidence,
        'sub_idx': sub_idx,
        'sub_labels': sub_labels,
        'cluster_sizes': size_df,
    }


def compare_consensus_to_inputs(clustering_matrix, consensus_labels):
    """Compare the consensus to each input clustering (ARI, NMI)."""
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    n_clust = clustering_matrix.shape[1]
    results = []
    for c in range(n_clust):
        ari = adjusted_rand_score(consensus_labels, clustering_matrix[:, c])
        nmi = normalized_mutual_info_score(consensus_labels, clustering_matrix[:, c],
                                            average_method='arithmetic')
        results.append({'clustering': c, 'ARI': ari, 'NMI': nmi})

    df = pd.DataFrame(results)
    print(f"\nConsensus vs. input clusterings:")
    print(f"  ARI: mean={df['ARI'].mean():.4f}, std={df['ARI'].std():.4f}, "
          f"range=[{df['ARI'].min():.4f}, {df['ARI'].max():.4f}]")
    print(f"  NMI: mean={df['NMI'].mean():.4f}, std={df['NMI'].std():.4f}, "
          f"range=[{df['NMI'].min():.4f}, {df['NMI'].max():.4f}]")
    return df
