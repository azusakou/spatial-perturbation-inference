#!/usr/bin/env python3
"""Native Squidpy categorical adjacency on exported known-truth cases.

Run after the guide-error simulation export. All nodes from the same simulated
focal guide remain in the graph: guide_count is a UMI count, not guide identity
multiplicity. Positive counts define observed sources. This endpoint is distinct
from SPF's continuous near-minus-far endpoint; no superiority claim is implied.
"""
from pathlib import Path
from importlib.metadata import version
import argparse
import hashlib
import inspect
import json
import os
import sys
import time

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("input", type=Path)
parser.add_argument("output", type=Path)
parser.add_argument("--limit", type=int, default=None, help="Smoke check only; omit for all cases")
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
cache = args.output / "native_cache"
cache.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(cache))

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
import squidpy as sq

N_PERMS, K, BASE_SEED, Z_CUTOFF = 199, 6, 20261016, 1.96
start = time.monotonic()
data = pd.read_csv(args.input, sep="\t")
required = {"case", "replicate", "amplitude", "node", "x", "y", "component",
            "observed_source", "true_source", "guide_count", "score"}
if not required.issubset(data.columns):
    raise ValueError(f"Missing columns: {sorted(required-set(data.columns))}")
for col in ("observed_source", "true_source"):
    if data[col].dtype != bool:
        val = data[col].astype(str).str.lower()
        if not val.isin(["true", "false", "1", "0"]).all():
            raise ValueError(f"Invalid Boolean {col}")
        data[col] = val.isin(["true", "1"])
assert np.array_equal(data.observed_source.to_numpy(), data.guide_count.to_numpy() > 0)
case_order = sorted(data["case"].unique())
graph_cache = {}
records, matrix_records = [], []


def make_graph(coords, component):
    key = hashlib.sha256(coords.tobytes() + component.astype("U").tobytes()).hexdigest()
    if key in graph_cache:
        return graph_cache[key]
    rows, cols = [], []
    for label in np.unique(component):
        idx = np.flatnonzero(component == label)
        if len(idx) < 2:
            continue
        k = min(K, len(idx)-1)
        neighbors = NearestNeighbors(n_neighbors=k+1, algorithm="kd_tree").fit(coords[idx])
        near = neighbors.kneighbors(coords[idx], return_distance=False)
        for local_i, candidates in enumerate(near):
            selected = candidates[candidates != local_i][:k]
            rows.extend(np.repeat(idx[local_i], len(selected)))
            cols.extend(idx[selected])
    graph = sparse.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)),
                              shape=(len(coords), len(coords)))
    graph = graph.maximum(graph.T).tocsr()
    graph.setdiag(0)
    graph.eliminate_zeros()
    coo = graph.tocoo()
    assert (graph != graph.T).nnz == 0
    assert np.all(component[coo.row] == component[coo.col])
    graph_cache[key] = graph
    return graph


def wilson(k, n):
    if n == 0:
        return (np.nan, np.nan)
    z = 1.959963984540054
    mid = (k/n + z*z/(2*n))/(1+z*z/n)
    half = z*np.sqrt((k/n)*(1-k/n)/n+z*z/(4*n*n))/(1+z*z/n)
    return mid-half, mid+half


groups = data.groupby(["case", "replicate", "amplitude"], sort=True, observed=True)
for index, ((case, replicate, amplitude), nodes) in enumerate(groups):
    if args.limit is not None and index >= args.limit:
        break
    nodes = nodes.sort_values("node").reset_index(drop=True)
    assert not nodes.node.duplicated().any()
    coords = nodes[["x", "y"]].to_numpy(float)
    component = nodes.component.astype(str).to_numpy()
    source = nodes.observed_source.to_numpy(bool)
    receiver = ~source
    seed = BASE_SEED + int(replicate) + case_order.index(case)*100000
    row = dict(case=case, replicate=int(replicate), amplitude=float(amplitude),
               nodes=len(nodes), observed_sources=int(source.sum()), receivers=int(receiver.sum()),
               hidden_sources_in_receivers=int((nodes.true_source.to_numpy(bool) & receiver).sum()),
               n_perms=N_PERMS, seed=seed, cutoff=Z_CUTOFF,
               decision_rule="abs(native source-to-receiver_high z) > 1.96",
               native_p_value=np.nan, p_value_status="native API does not return empirical P")
    if source.sum() == 0 or receiver.sum() < 2:
        row.update(status="unestimable_missing_source_or_receiver", rejected=np.nan)
        records.append(row)
        continue
    threshold = float(np.median(nodes.loc[receiver, "score"]))
    labels = np.repeat("source", len(nodes)).astype(object)
    labels[receiver] = np.where(nodes.loc[receiver, "score"].to_numpy() > threshold,
                               "receiver_high", "receiver_low")
    cats = ["source", "receiver_high", "receiver_low"]
    if not all(np.any(labels == label) for label in cats):
        row.update(status="unestimable_empty_high_or_low", rejected=np.nan)
        records.append(row)
        continue
    graph = make_graph(coords, component)
    obs = pd.DataFrame({"label": pd.Categorical(labels, categories=cats),
                        "component": pd.Categorical(component)},
                       index=[f"node_{n}" for n in nodes.node])
    adata = ad.AnnData(X=sparse.csr_matrix((len(nodes), 0)), obs=obs)
    adata.obsm["spatial"] = coords
    adata.obsp["spatial_connectivities"] = graph
    sq.gr.nhood_enrichment(adata, cluster_key="label", library_key="component",
                          connectivity_key="spatial", n_perms=N_PERMS,
                          seed=seed, n_jobs=1, show_progress_bar=False, copy=False)
    result = adata.uns["label_nhood_enrichment"]
    counts = np.asarray(result["count"])
    zscore = np.asarray(result["zscore"])
    assert adata.obs.label.cat.categories.tolist() == cats
    # Independently evaluate label adjacency by matrix multiplication. This
    # does not call or reproduce Squidpy's permutation implementation.
    onehot = np.column_stack([labels == label for label in cats]).astype(np.int64)
    independent = onehot.T @ graph.astype(np.int64) @ onehot
    assert np.array_equal(counts, independent)
    assert counts.sum() == graph.nnz and np.array_equal(counts, counts.T)
    high_z, low_z = float(zscore[0, 1]), float(zscore[0, 2])
    estimable = np.isfinite(high_z)
    row.update(status="estimable" if estimable else "unestimable_native_z", graph_directed_edges=int(graph.nnz),
               receiver_median=threshold, source_high_count=int(counts[0, 1]), source_high_z=high_z,
               source_low_count=int(counts[0, 2]), source_low_z=low_z,
               rejected=bool(abs(high_z) > Z_CUTOFF) if estimable else np.nan,
               independent_count_max_error=int(np.max(np.abs(counts-independent))))
    records.append(row)
    for j, from_label in enumerate(cats):
        for k, to_label in enumerate(cats):
            matrix_records.append(dict(case=case, replicate=int(replicate), amplitude=float(amplitude),
                from_label=from_label, to_label=to_label, native_count=int(counts[j,k]),
                independent_count=int(independent[j,k]), native_z=float(zscore[j,k])))
    if (index+1) % 25 == 0:
        pd.DataFrame(records).to_csv(args.output / "native_known_truth_records.partial.tsv", sep="\t", index=False)
        print(f"Completed {index+1} cases in {time.monotonic()-start:.1f}s", flush=True)

records = pd.DataFrame(records)
records.to_csv(args.output / "native_known_truth_records.tsv", sep="\t", index=False, na_rep="NA")
pd.DataFrame(matrix_records).to_csv(args.output / "native_known_truth_matrices.tsv", sep="\t", index=False, na_rep="NA")
summaries = []
for (case, amplitude), group in records.groupby(["case", "amplitude"], sort=True):
    valid = group.status.eq("estimable")
    n = int(valid.sum())
    k = int(group.loc[valid, "rejected"].astype(bool).sum())
    lo, hi = wilson(k, n)
    summaries.append(dict(case=case, amplitude=amplitude, attempted=len(group), estimable=n,
                          unestimable=int((~valid).sum()), rejected=k,
                          gaussian_z_rule_rate=k/n if n else np.nan, wilson_lower=lo, wilson_upper=hi,
                          endpoint="source-to-receiver_high categorical adjacency",
                          interpretation="zero-extrinsic-field rejection" if amplitude==0 else "nonzero-extrinsic-field rejection",
                          decision_rule=f"abs(native z) > {Z_CUTOFF}", native_empirical_p=False))
summary = pd.DataFrame(summaries)
summary.to_csv(args.output / "native_known_truth_summary.tsv", sep="\t", index=False, na_rep="NA")
manifest = dict(status="smoke_completed" if args.limit else "native_known_truth_completed",
    input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(), cases=len(records),
    native_api="squidpy.gr.nhood_enrichment", signature=str(inspect.signature(sq.gr.nhood_enrichment)),
    versions={name: version(name) for name in ["squidpy","anndata","numpy","scipy","pandas","numba","scikit-learn"]},
    python=sys.version, native_import=sq.__file__, n_perms=N_PERMS, knn=K,
    graph="symmetric union of six-nearest-neighbour edges within component, fixed per geometry",
    graph_cache_geometries=len(graph_cache), label_threshold="receiver median per case; high strictly above median",
    source_rule="positive UMI count of the single focal guide; missed true sources become receivers",
    decision_rule=f"abs(native source-high z) > {Z_CUTOFF}; Gaussian z rule, not native empirical P",
    scope="known-truth categorical adjacency; differs from SPF continuous contrast; no method-superiority inference",
    timing_seconds=time.monotonic()-start,
    observed_count_check="all native matrices identical to independently computed onehot.T @ graph @ onehot")
(args.output / "run_manifest.json").write_text(json.dumps(manifest, indent=2))
print(summary.to_string(index=False), flush=True)
print(f"NATIVE_KNOWN_TRUTH_COMPLETED {len(records)} cases in {manifest['timing_seconds']:.1f}s", flush=True)
