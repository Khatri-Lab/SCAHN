"""
Leiden clustering of the shared-nearest-neighbour graph
=======================================================
Partitions the SNN graph at one resolution with Leiden, using the
RBConfiguration objective. The Python half of the clustering sweep;
04a_run_clusteringlouvain.R is the other half, and the two write into the same
directory under names that 05_run_cluster_aggregate.R then collects.

Leiden only labels the nodes it partitions, so the run is restricted to the
largest connected component and the cells outside it are afterwards assigned to
the nearest cluster centroid in the harmony embedding. Without that step those
cells would keep the label -1 and the column would not be comparable to the
Louvain columns, which label every cell.

Usage:
  python 04b_run_clusteringleiden.py <vargenes> <nPCs> <seed> <resolution>
  python 04b_run_clusteringleiden.py v2000 25 42 0.8

  The third argument is a spent seed slot. The Leiden seed is fixed at 42 below
  and the output name no longer records a seed, so the value is read and
  ignored; the slot is kept only so that existing four-argument invocations
  still line up.

Requires:
  numpy pandas scipy scikit-learn python-igraph leidenalg

Reads:
  results/integration/clustering/SCAHN.snn_graph.<vargenes>.PC<nPCs>.mtx
    written by 03_run_SNN.R
  results/integration/harmony/SCAHN.harmony.Z_corr.<vargenes>.PC<nPCs>.csv
    written by 02a_run_harmony_SCT.R, for the centroid assignment above

Writes:
  results/integration/clustering/
    SCAHN.clustering.<vargenes>.PC<nPCs>.leiden.res<resolution>.csv
    read by 05_run_cluster_aggregate.R
"""

import sys

import igraph as ig
import leidenalg
import numpy as np
import pandas as pd
from scipy.io import mmread
from sklearn.neighbors import NearestCentroid

DIR_RESULTS = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/results/"
DIR_HARMONY = DIR_RESULTS + "integration/harmony/"
DIR_CLUSTERING = DIR_RESULTS + "integration/clustering/"


# -- Parameters ----------------------------------------------------------------

args = sys.argv[1:]
vargenes = args[0]      # "v1000" "v2000" "v3000"
nPCs = int(args[1])     # 20 25 30
reso = float(args[3])   # 0.5 0.8 1.2


# -- SNN graph as igraph -------------------------------------------------------

snn_mtx = mmread(
    f"{DIR_CLUSTERING}SCAHN.snn_graph.{vargenes}.PC{nPCs}.mtx").tocsr()

# Every cell is its own nearest neighbour, so the diagonal carries a self-loop
# that would inflate every node's weighted degree.
snn_mtx.setdiag(0)
snn_mtx.eliminate_zeros()

sources, targets = snn_mtx.nonzero()
weights = np.array(snn_mtx[sources, targets]).flatten()
g = ig.Graph(
    n=snn_mtx.shape[0], edges=list(zip(sources.tolist(), targets.tolist())), directed=False,)
g.es["weight"] = weights

# The SNN matrix is symmetric, so nonzero() yielded each edge twice.
g.simplify(combine_edges="max")


# -- Leiden on the largest component -------------------------------------------

components = g.connected_components()
largest = max(components, key=len)
print(f"Total cells: {g.vcount()}")
print(f"Connected components: {len(components)}")
print(f"Largest component size: {len(largest)}")
print(f"Singletons: {sum(1 for s in components.sizes() if s == 1)}")

# -1 marks a cell the partition below does not reach; all of them are resolved
# in the next section.
clusters = np.full(g.vcount(), -1, dtype=int)

subgraph = g.subgraph(largest)
partition = leidenalg.find_partition(
    subgraph, leidenalg.RBConfigurationVertexPartition, resolution_parameter=reso, weights="weight",
    seed=42,)

# partition indexes the subgraph, so map its membership back through `largest`.
for i, node in enumerate(largest):
    clusters[node] = partition.membership[i]

n_clusters = len(set(partition.membership))
n_unassigned = np.sum(clusters == -1)
print(f"Found {n_clusters} clusters (+ {n_unassigned} unassigned nodes)")


# -- Assign the cells left out -------------------------------------------------

# Nearest centroid in the corrected embedding, not in the graph: the cells left
# out are exactly the ones the graph does not connect to the rest.
embedding = pd.read_csv(
    f"{DIR_HARMONY}SCAHN.harmony.Z_corr.{vargenes}.PC{nPCs}.csv").values

assigned_mask = clusters != -1
unassigned_mask = clusters == -1

clf = NearestCentroid()
clf.fit(embedding[assigned_mask], clusters[assigned_mask])
clusters[unassigned_mask] = clf.predict(embedding[unassigned_mask])

print(f"After reassignment: {np.sum(clusters == -1)} unassigned remaining")

np.savetxt(
    f"{DIR_CLUSTERING}SCAHN.clustering.{vargenes}.PC{nPCs}" f".leiden.res{reso}.csv",
    clusters, delimiter=",", fmt="%d", header="cluster", comments="",)
