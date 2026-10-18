#!/usr/bin/env python3
"""Native Squidpy on all eligible source/receiver cells in one observed chip."""
from pathlib import Path
import hashlib
import inspect
import json
import os
import sys
import time
from importlib.metadata import version

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
CACHE = BASE / "results" / "native_cache"
CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(CACHE / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(CACHE))

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from sklearn.neighbors import NearestNeighbors
import squidpy as sq

started = time.monotonic()
N_PERMS, SEED, K = 499, 20261015, 6
OUT = BASE / "results" / "native"
OUT.mkdir(parents=True, exist_ok=True)
input_path = BASE / "results" / "native_input.tsv.gz"
if not input_path.exists():
    raise FileNotFoundError("First run: Rscript run.R native-export")
full = pd.read_csv(input_path, sep="\t")
assert full["chip_id"].nunique() == 1

# Use every mapped cell to reconstruct the original occupied 250-um tissue
# mask; do not infer components from the source/receiver subset alone.
grid = np.floor(full[["x_um", "y_um"]].to_numpy() / 250.0).astype(np.int64)
occupied, inverse = np.unique(grid, axis=0, return_inverse=True)
lookup = {tuple(x): i for i, x in enumerate(occupied)}
edge_from, edge_to = [], []
for i, (x, y) in enumerate(occupied):
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            if dx == dy == 0:
                continue
            j = lookup.get((x + dx, y + dy))
            if j is not None:
                edge_from.append(i)
                edge_to.append(j)
mask_graph = sparse.csr_matrix((np.ones(len(edge_from)), (edge_from, edge_to)),
                               shape=(len(occupied), len(occupied)))
n_components, component_labels = connected_components(mask_graph, directed=False)
full["component"] = component_labels[inverse].astype(str)
nodes = full[full["guide_call"].isin(["Neg", "sgrna-clu"])].copy()
nodes["map_row"] = nodes.index
nodes.index = pd.Index([f"map_row_{i}" for i in nodes.index], name="cell_id")
is_receiver = nodes["guide_call"].eq("Neg")
assert int(is_receiver.sum()) == 119664
assert int((~is_receiver).sum()) == 419
assert len(nodes) == 120083
threshold = float(nodes.loc[is_receiver, "raw_clu_score"].median())
labels = np.repeat("Clu_source", len(nodes)).astype(object)
labels[is_receiver] = np.where(nodes.loc[is_receiver, "raw_clu_score"] > threshold,
                               "receiver_high", "receiver_low")
nodes["global_median_label"] = pd.Categorical(labels)
nodes["component"] = nodes["component"].astype("category")
coords = nodes[["x_um", "y_um"]].to_numpy(float)

rows, cols = [], []
for component in nodes["component"].cat.categories:
    idx = np.flatnonzero(nodes["component"].eq(component).to_numpy())
    if len(idx) < 2:
        continue
    k = min(K, len(idx) - 1)
    model = NearestNeighbors(n_neighbors=k + 1, algorithm="kd_tree").fit(coords[idx])
    near = model.kneighbors(coords[idx], return_distance=False)
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
comp = nodes["component"].astype(str).to_numpy()
assert np.all(comp[coo.row] == comp[coo.col])
print(f"Prepared full chip: {len(nodes)} nodes, {graph.nnz} directed edges; "
      f"elapsed={time.monotonic() - started:.1f}s", flush=True)

adata = ad.AnnData(X=sparse.csr_matrix((len(nodes), 0)), obs=nodes)
adata.obsm["spatial"] = coords
adata.obsp["spatial_connectivities"] = graph
sq.gr.nhood_enrichment(adata, cluster_key="global_median_label", library_key="component",
                      connectivity_key="spatial", n_perms=N_PERMS, seed=SEED,
                      n_jobs=1, show_progress_bar=False, copy=False)
result = adata.uns["global_median_label_nhood_enrichment"]
categories = adata.obs["global_median_label"].cat.categories.tolist()
count = pd.DataFrame(result["count"], index=categories, columns=categories)
zscore = pd.DataFrame(result["zscore"], index=categories, columns=categories)
assert int(count.to_numpy().sum()) == graph.nnz
assert np.array_equal(count.to_numpy(), count.to_numpy().T)
count.to_csv(OUT / "global_median_label_native_counts.tsv", sep="\t", index_label="from_label")
zscore.to_csv(OUT / "global_median_label_native_zscores.tsv", sep="\t", index_label="from_label")
summary = pd.DataFrame([{"from_label": "Clu_source", "to_label": label,
                        "native_edge_count": int(count.loc["Clu_source", label]),
                        "native_zscore": float(zscore.loc["Clu_source", label]),
                        "n_perms": N_PERMS, "native_p_value": None,
                        "p_value_status": "not returned by native API; not inferred from z"}
                       for label in categories])
summary.to_csv(OUT / "source_adjacency_summary.tsv", sep="\t", index=False)
manifest = {
    "status": "native_fullchip_execution_completed", "chip": str(full["chip_id"].iloc[0]),
    "input": str(input_path), "input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
    "source_count": 419, "receiver_count": 119664, "nodes": len(nodes),
    "sampling": "all eligible Clu sources and all guide-undetected receivers in the exported chip; no source/receiver downsampling",
    "mask_cells": len(full), "component_mask_grid_um": 250.0, "mask_components": n_components,
    "component_node_counts": nodes["component"].value_counts().to_dict(),
    "graph_k": K, "graph_directed_edges": int(graph.nnz),
    "graph_unique_undirected_edges": int(graph.nnz // 2),
    "label_counts": nodes["global_median_label"].value_counts().to_dict(),
    "receiver_median_threshold": threshold, "threshold_ties": "low: score <= median; high: score > median",
    "native_api": "squidpy.gr.nhood_enrichment", "native_import_path": sq.__file__,
    "signature": str(inspect.signature(sq.gr.nhood_enrichment)),
    "n_perms": N_PERMS, "seed": SEED, "n_jobs": 1, "library_key": "component",
    "native_return_keys": list(result), "p_value": "not returned; not reconstructed from z scores",
    "count_convention": "directed adjacency incidences; diagonal same-label undirected edges count twice, cross-label edge counts once in each symmetric off-diagonal entry",
    "source_source_unique_undirected_edges": int(count.loc["Clu_source", "Clu_source"] // 2),
    "runtime_seconds_excluding_import": time.monotonic() - started,
    "python": sys.version,
    "versions": {p: version(p) for p in ["squidpy", "anndata", "scanpy", "spatialdata", "numpy", "scipy", "pandas", "numba", "scikit-learn"]},
    "scope": "one full observed chip; categorical adjacency under fixed graph, not the continuous-score SPF near-far endpoint or known biological ground truth"
}
(OUT / "run_manifest.json").write_text(json.dumps(manifest, indent=2))
print(summary.to_string(index=False), flush=True)
print(f"FULLCHIP_COMPLETED in {manifest['runtime_seconds_excluding_import']:.1f}s", flush=True)
