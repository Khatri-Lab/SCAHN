"""
Consensus clustering across the sweep
=====================================
Collapses the many clusterings of the sweep into one labelling, so that the
annotation work downstream has a single partition to name rather than 81
competing ones. The algorithm is in consensus_clustering.py; this script only
selects the input columns and writes the result.

Unlike 06a_run_evalcluster.py, this does not use the one-resolution-per-pipeline
selection: every PC25 and PC30 clustering goes in, at all three sweep
resolutions, so that agreement is measured over the whole retained sweep. PC20
is excluded, matching the same exclusion in 05_run_cluster_aggregate.R.

The resolution argument is the Leiden resolution used on the co-clustering graph
inside consensus_clustering, and is unrelated to the sweep resolutions in the
column names. The runs kept on disk are 0.1, 0.2, 0.4, 0.6 and 0.8.

Usage:
  python 06b_run_consensuscluster.py --resolution 0.4

Requires:
  numpy pandas, plus whatever consensus_clustering.py needs
  (scipy scikit-learn python-igraph leidenalg)

Reads:
  results/integration/clustering/clusterings.csv
    written by 05_run_cluster_aggregate.R

Writes:
  results/integration/clustering/consensus_results.<resolution>.csv
    one row per cell, consensus_label and confidence; read by SCAHN.anno.tmp.Rmd
"""

import argparse
import re

import numpy as np
import pandas as pd
from consensus_clustering import consensus_clustering

DIR_RESULTS = "/labs/khatrilab/hongzheng/SCAHN/workflow_final/results/"
DIR_CLUSTERING = DIR_RESULTS + "integration/clustering/"


# -- Parameters ----------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument(
    "--resolution", type=float, required=True, help="Leiden resolution for the consensus graph",)
args = parser.parse_args()
resolution = args.resolution


# -- Select the input clusterings ----------------------------------------------

df = pd.read_csv(f"{DIR_CLUSTERING}clusterings.csv")
col_names = df.columns.tolist()
selected = [c for c in col_names if re.search(r"\.PC(25|30)\.", c)]
print(selected)

df_selected = df[selected]
clustering_matrix = df_selected.values.astype(np.int32)


# -- Consensus -----------------------------------------------------------------

result = consensus_clustering(clustering_matrix, leiden_resolution=resolution)

pd.DataFrame(
    {
        "consensus_label": result["labels"],
        "confidence": result["confidence"],
    }).to_csv(f"{DIR_CLUSTERING}consensus_results.{resolution}.csv", index=False)
