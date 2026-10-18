#!/usr/bin/env python3
"""Native Squidpy categorical-adjacency comparison on a frozen Clu subset.

This calls the installed complete Squidpy distribution through its public API.
It does not reconstruct or extract Squidpy's internal permutation function.
The endpoint differs from SPF's continuous-score near-minus-far contrast.
"""
from pathlib import Path
import hashlib
import inspect
import json
import os
from importlib.metadata import version
import sys

_cache = Path(__file__).resolve().parents[1] / "results" / "native_cache"
_cache.mkdir(parents=True,exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_cache / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_cache))

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
import squidpy as sq

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
OUT = BASE / "results" / "native_sampled"
OUT.mkdir(parents=True,exist_ok=True)
N_PERMS = 499
SEED = 20261015
K_NEIGHBORS = 6
paths = [BASE / "results/native_sampled_receiver.tsv",
         BASE / "results/native_sampled_source.tsv"]
receiver = pd.read_csv(paths[0], sep="\t")
source = pd.read_csv(paths[1], sep="\t")
assert len(receiver) == 2500 and len(source) == 419
receiver["source_status"] = "receiver"
source["source_status"] = "Clu_source"
nodes = pd.concat([source, receiver], ignore_index=True)
nodes.index = pd.Index([f"node_{i:04d}" for i in range(len(nodes))], name="cell_id")
coords = nodes[["x", "y"]].to_numpy(float)
nodes["component"] = nodes["component"].astype(str).astype("category")
is_receiver = nodes["source_status"].eq("receiver")

# Fixed before native inference: global median and within-identity medians.
# These are categorical endpoint definitions, not tuned to a test result.
global_threshold = float(receiver["raw_clu_score"].median())
identity_thresholds = receiver.groupby("identity")["raw_clu_score"].median().to_dict()
global_labels = np.repeat("Clu_source", len(nodes)).astype(object)
global_labels[is_receiver] = np.where(
    nodes.loc[is_receiver, "raw_clu_score"].to_numpy() > global_threshold,
    "receiver_high", "receiver_low")
nodes["global_median_label"] = pd.Categorical(global_labels)
identity_labels = np.repeat("Clu_source", len(nodes)).astype(object)
for i in np.flatnonzero(is_receiver.to_numpy()):
    identity = nodes.iloc[i]["identity"]
    state = "high" if nodes.iloc[i]["raw_clu_score"] > identity_thresholds[identity] else "low"
    identity_labels[i] = f"{identity}_{state}"
nodes["identity_median_label"] = pd.Categorical(identity_labels)

# A symmetric union of six-nearest-neighbour graphs built separately within
# each pre-exported tissue component. Graph coordinates and edges stay fixed.
rows, cols = [], []
for component in nodes["component"].cat.categories:
    idx = np.flatnonzero(nodes["component"].eq(component).to_numpy())
    if len(idx) < 2:
        continue
    k = min(K_NEIGHBORS, len(idx) - 1)
    neighbors = NearestNeighbors(n_neighbors=k + 1).fit(coords[idx])
    near = neighbors.kneighbors(coords[idx], return_distance=False)
    for local_i, candidates in enumerate(near):
        selected = candidates[candidates != local_i][:k]
        rows.extend(np.repeat(idx[local_i], len(selected)))
        cols.extend(idx[selected])
graph = sparse.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)),
                          shape=(len(nodes), len(nodes)))
graph = graph.maximum(graph.T).tocsr()
graph.setdiag(0)
graph.eliminate_zeros()
assert (graph != graph.T).nnz == 0
coo = graph.tocoo()
assert np.all(nodes["component"].astype(str).to_numpy()[coo.row] ==
              nodes["component"].astype(str).to_numpy()[coo.col])

adata = ad.AnnData(X=sparse.csr_matrix((len(nodes), 0)), obs=nodes.copy())
adata.obsm["spatial"] = coords
adata.obsp["spatial_connectivities"] = graph
summary = []
for key in ["global_median_label", "identity_median_label"]:
    print(f"Running native Squidpy {sq.__version__}: {key}, {N_PERMS} permutations", flush=True)
    sq.gr.nhood_enrichment(adata, cluster_key=key, library_key="component",
                          connectivity_key="spatial", n_perms=N_PERMS,
                          seed=SEED, n_jobs=1, show_progress_bar=False, copy=False)
    result = adata.uns[f"{key}_nhood_enrichment"]
    labels = adata.obs[key].cat.categories.tolist()
    count = pd.DataFrame(result["count"], index=labels, columns=labels)
    zscore = pd.DataFrame(result["zscore"], index=labels, columns=labels)
    count.to_csv(OUT / f"{key}_native_counts.tsv", sep="\t", index_label="from_label")
    zscore.to_csv(OUT / f"{key}_native_zscores.tsv", sep="\t", index_label="from_label")
    for to_label in labels:
        summary.append({"endpoint": key, "from_label": "Clu_source", "to_label": to_label,
                        "native_edge_count": int(count.loc["Clu_source", to_label]),
                        "native_zscore": float(zscore.loc["Clu_source", to_label]),
                        "n_perms": N_PERMS, "native_p_value": None,
                        "p_value_status": "not returned by native API; not inferred from z"})
    (OUT / f"{key}_native_return_keys.json").write_text(json.dumps(list(result), indent=2))
pd.DataFrame(summary).to_csv(OUT / "source_adjacency_summary.tsv", sep="\t", index=False)
manifest = {
    "python": sys.version, "squidpy": sq.__version__, "squidpy_import_file": sq.__file__,
    "native_function": "squidpy.gr.nhood_enrichment", "signature": str(inspect.signature(sq.gr.nhood_enrichment)),
    "package_versions": {p: version(p) for p in ["squidpy", "anndata", "scanpy", "spatialdata", "numpy", "scipy", "pandas", "numba", "scikit-learn"]},
    "receiver_count": len(receiver), "source_count": len(source), "graph_nodes": len(nodes),
    "graph_directed_edges": graph.nnz, "knn_k": K_NEIGHBORS,
    "component_counts": nodes["component"].value_counts().to_dict(),
    "global_receiver_median_threshold": global_threshold,
    "identity_receiver_median_thresholds": identity_thresholds,
    "threshold_ties": "low: score <= median; high: score > median",
    "library_key": "component", "permutations": N_PERMS, "seed": SEED,
    "input_sha256": {str(p.relative_to(BASE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
    "status": "native_package_execution_completed",
    "p_value": "Native API returns edge counts and null-standardized z scores only. No P value computed.",
    "interpretation": "Illustrative categorical-adjacency endpoint on a fixed receiver subset; not SPF continuous near-far endpoint, not known biological false positives, and not a full-chip benchmark.",
}
(OUT / "run_manifest.json").write_text(json.dumps(manifest, indent=2))
print(pd.DataFrame(summary).to_string(index=False))
