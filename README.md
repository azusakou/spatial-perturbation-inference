## Interpretation and evidence boundaries

SPF is a diagnostic workflow for specified receiver scores and spatial references. Real Clu/DBiT tail fractions and BH-transformed fractions are **descriptive reference diagnostics**, without established general type-I-error/FDR control. Guide-selected receivers, unknown source assignment and graph-estimated source specifications do not meet the known-Q premises below. Legacy column names such as `permutation_p` and `fdr_within_chip` are retained for compatibility, not inferential guarantees. This restriction is separate from the documented mouse-paired edgeR and clinical outcome-label tests.

The manuscript and response letters are separate deliverables. The prepared Reviewer Evidence companion supplies supported figures, plotting code, source-data tables and a reviewer index; this numerical package excludes that material. Preparation is not public deposition or journal upload.

## Run

Tested with R 4.4.2; versions in `DEPENDENCIES.tsv`. Core: data.table, RANN, igraph, mgcv, cluster, Matrix, edgeR, limma, statmod. No automatic installation or downloads.

```sh
Rscript run.R all results               # fast statistics, counts, profiles and new checks
Rscript run.R complete results          # additionally regenerate strict Clu / DBiT references
Rscript run.R verify results            # compare completed outputs to current snapshots
Rscript run.R checks results            # Clu count holdouts + new record aggregations only
Rscript run.R full-holdout results/clu   # refit all 20 Clu fold/route/version results from counts
Rscript run.R shape results/shape        # enumerate rigid footprints and real-KP matching
Rscript run.R smoke results/smoke
```

Any mode's second argument is an output directory. `complete` includes the fast new checks, but independent generators remain separate:

```sh
Rscript run.R benchmark results/benchmark 100 199
Rscript run.R guide-error results/guide_error 100 199
Rscript run.R injection results/injection 199
Rscript run.R dropout results/thinning 199
```

The fast `checks` mode preserves all 27 previous checks and adds 16 checks: ten known-assignment record aggregations and six frozen real-Clu reference reaggregations (43 total). Clu fold counts are freshly refit; new outputs containing `frozen` are explicit snapshots/protocols, not cell-level refits.

The 14-scenario benchmark has 50,400 tests and 504 summaries. The common observed/null centre replaces historical null-only centring. The historical formula recovers every previous P; the revised formula changes 126 rejection decisions. Baseline GAM/cluster rejects 1/100 null datasets and detects 74/100 unit-effect datasets; 24/168 zero-effect settings have Wilson lower limits above .05. Source-generator exchangeability and universal calibration are **not established**.

The new guide simulation generates guide counts and four-gene RNA counts, moving hidden sources into observed receivers. Its 57,600 records include 96 unavailable tests (two geometry failures). Attempted/failure counts are retained; rejection rates use the 99 or 100 estimable realizations, with no NA-to-zero replacement. With about 50% detection and direct source effects but no receiver field, primary observed-label rejection is 33/100 versus oracle 3/100; without direct source effects it is 3/100. Apparent detection with a field is not calibrated power. Simulated sensitivity is not an estimate of assay sensitivity. `dropout` is the older fixed-receiver source-thinning example, a different stress test.

## Known-assignment calibration and real-map control audit

From any working directory, pass the package root and output paths explicitly:

```sh
Rscript /path/to/package/optional/independent_bank.R /path/to/package work/knownQ 100 199
Rscript /path/to/package/optional/clu_control_audit.R /path/to/package work/clu_controls
Rscript /path/to/package/optional/clu_control_audit.R /path/to/package work/clu_full --full-only
Rscript /path/to/package/optional/plot_clu_control_audit.R /path/to/package work/clu_controls
```

The independent-bank generator completed six fixed scenarios × 100 realizations, GAM, amplitudes0/1 and two source laws (2,400 records/24 summaries). Each realization draws new199-configuration banks. True sources and the known-Q reference share exactly the same fixed, generator-known anchor/niche law; specifications are not re-estimated from the observed source graph. Geometry/covariates are fixed once per scenario. Under zero field the score transformation and receiver population are source-independent, so the symmetric rank has a conditional exchangeability justification. The complete frozen protocol and conditional argument are in `readRDS("data/computational_checks.rds")$known_assignment$protocol` and `$proof`, and are exported by `checks`. The local scientific-module header's `PROTOCOL_FROZEN.md` refers to that protocol. The optional CLI permits1–100 realizations; the frozen reference budget is199.

Known-Q zero-field rejections are5/6/7/3/1/2 of100, with amplitude1 detection72/69/88/69/83/36. All twelve zero-field Wilson lower limits (including the cluster-only comparator) are≤.05. These finite model frequencies do not validate real guide assignment. They do not replace the historical24/168 flags, receiver contamination33/100, programme-transfer failures or n=2 animal limitation. New geometries/backgrounds/RNG differ from the historical benchmark; changes cannot be attributed solely to bank regeneration.

The real-Clu optional command refits five backgrounds crossed with two source banks in three chips (30 endpoints,2,970 null statistics), with common banks across backgrounds. The raw→identity→identity/library→full-GAM path is nested; the global-linear model is a separate sensitivity. Cluster-reference sequences are.01/.01/.32/.78,.17/.31/.95/.91 and.06/.13/.21/.99; clustering is not uniformly more conservative. The short `--full-only` command actually refits the three full-GAM/cluster endpoints and297 null statistics, restoring.78/.91/.99. Its output directory must differ from the complete30-row mode. Plotting requires the complete30-row result.

Mixed-guide geometry reports another detected single-guide source within100 micrometres for84.6%,75.1%,62.6% of Clu-near receivers. Other single guides include mSafe; multiple-guide neighbours are separate and fractions can overlap. This measures spatial co-exposure, not molecular contamination or another guide's causal effect.

Native Squidpy 1.6.5 / Python 3.12 (existing pinned environment executed; fresh environment installation not tested):

```sh
Rscript run.R native-export results
python optional/native_squidpy.py
python optional/native_sampled.py
python optional/native_known_truth.py results/guide_error/native_known_truth_input.tsv.gz results/native_truth
```

The last command regenerates 600 cases, 199 label permutations each. Native categorical adjacency z scores use the stated Gaussian cutoff; they are not empirical molecular-field P values. Perfect-call zero-field rejection was 87/100, all negative adjacency depletion under clustered source geometry. Native and continuous SPF endpoints differ; no superiority claim follows. The full-chip and sampled scripts read the default results/ exports and write results/native and results/native_sampled; known-truth input/output paths are explicit.

## Raw reconstruction and independent mouse folds

Obtain the five GSE274058 10x directories, each with matrix/features/barcodes triplets. Seurat 5.4.0, Harmony 1.2.4 and additional optional dependencies are recorded.

```sh
Rscript optional/reference.R preprocess /path/to/GSE274058_extracted work/reference
Rscript optional/reference.R export work/reference/perturbseq_qc_harmony.rds work/export
Rscript run.R reference work/refit work/export/reference_counts.rds
Rscript optional/reference.R identities work/reference/perturbseq_qc_harmony.rds work/reference
Rscript optional/reference.R markers work/reference/perturbseq_qc_harmony.rds work/reference
Rscript optional/full_holdout.R /path/to/GSE274058_extracted work/full_holdout
```

The joint five-mouse replay exactly reproduced 41,140 QC cells, 24,295 genes, features, coordinates, matches, weights and markers. That numerical replay is separate from the last command, which actually rebuilds preprocessing in four training mice per fold. Training HVG/scaling/PCA/Harmony, matching, DE, mSafe calibration and weights exclude the query. Query matching uses fixed training PCA projection because no unseen-batch Harmony transform is defined; the auxiliary route also matches training cells in PCA. Both routes are reported, with unavailable folds retained. Five raw folds generated 340 rows for 17 non-mSafe guides × two versions × two routes. Primary Clu no-IEG AUCs are .530864/.359375, with three unavailable mSafe folds; stable generalization is not established.

The compact Clu+mSafe counts refit 20 rows exactly. All-guide metrics and scored weights are a **frozen snapshot**, not a compact all-guide raw reconstruction. The optional raw command regenerates all guides. This distinction saves approximately 13 MB; no large processed Seurat objects are required in the release.

Raw Clu/DBiT readers remain in `optional/`. HCC raw geometry defaults to correct hex spacing; pixel/legacy branches are explicit:

```sh
Rscript optional/external_features.R hcc PAIRED_H5AD EXTERNAL_VISIUM_DIRECTORY work/hcc hex
Rscript optional/threshold_sensitivity.R mi MI_H5AD_DIRECTORY work/mi_windows
Rscript optional/threshold_sensitivity.R skin SKIN_10X_DIRECTORY work/skin_windows
```

The threshold script accepts `mi RAW_DIRECTORY OUTPUT_DIRECTORY` or `skin RAW_DIRECTORY OUTPUT_DIRECTORY`; inspect its usage for the input layout. Clinical, atlas, injury and window fast statistics start from documented derived features/endpoints. The package does not claim their full matrices are bundled.



## Eight compact inputs

| File | Boundary and use |
|---|---|
| `clu_cells.rds` | All 229,775 cells, guide labels, coordinates, identities, programme scores, native GEF libraries and components; refits strict spatial backgrounds/nulls. |
| `reference_counts.rds` | Fresh matched full-gene pseudobulks/heldout cells; refits original weights and historical shared-embedding weight holdouts. |
| `perturbmap_spots.tsv.gz` | All 6,363 spots, native hex indices, image coordinates, 42 raw targeted counts, whole libraries and original annotations; regenerates lesion/animal and shape controls. |
| `dbit_scores.rds` | Four focal sections, source masks, coordinates, scores and library proxies; refits background and section references. |
| `hallmark_2026_1.tsv.gz` | Exact eight Hallmark memberships, MSigDB 2026.1.Hs. |
| `endpoints.rds` | Corrected hex HCC, MI, atlas inputs; original exploratory screening; new benchmark records and validation snapshots. |
| `computational_checks.rds` | Clu+mSafe fold counts; all-guide frozen snapshot; complete guide/native records; shape/HCC/window/injection results; complete known-Q records/protocol and frozen30-row Clu control audit. Fast aggregation is distinguished from generation in filenames. |
| `clu_control_cluster_bank.tsv` | 297 original99-draw cluster configurations across three chips, with chip-local row indices; reusable real-map audit input, not new biological data. |

`DATA_DICTIONARY.tsv` identifies objects and provenance. `VALIDATION.md` states actual executions and limits. There are only two focal Perturb-map animals: KP_1/KP_2 share mouse 4; KP_4 is mouse 7. Fourteen annotations are dependent subsamples. Shape enumeration reveals insufficient matching support, not a valid matched null. Real-KP/composition sensitivities retain direction but do not establish causality. Unknown cross-section clone identity and two animals preclude population significance or mediation inference.

HCC uses x=array_col/2, y=sqrt(3)*array_row/2, with section median nearest-neighbour normalization; units are hex spacing, not micrometres. Pixel sensitivity was executed. Eleven fixed geometry/threshold/window specifications retain all empty bins and coverage changes. The 62-curve ordinal atlas has k=2, silhouette .488826, 48/14 shapes; profile comparisons use B1–B5 rather than a shared physical radius. Technology holdout uses Stereo-seq (10) versus Visium (52), excluding the whole technology. Clinical SPF models exclude proxy-source scores; no sign or penalty was tuned to improve poor held-out results.

For new data:

```sh
Rscript optional/apply.R receivers.tsv sources.tsv output_prefix 499 20261015
```

Receivers: x,y,component,niche,identity,score, optional log_lib; sources: x,y,component,niche. Coordinates in micrometres. Supply justified receiver eligibility, candidate-pool exchangeability and biological units. The CLI cannot establish these assumptions from coordinates or replace biological replication.

Data and code rights are separate; see `LICENSE_NOTICE.md`. No public deposition identifier is invented.
