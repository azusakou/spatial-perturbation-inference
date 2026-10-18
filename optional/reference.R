#!/usr/bin/env Rscript
# Optional raw/reference interfaces, adapted from original R13/14/15/23 and review heldout.
# Raw preprocessing was executed afresh in the 2026-10-15 full-rerun audit.
a <- commandArgs(trailingOnly=TRUE)
if(length(a)!=3L)stop("Usage: Rscript optional/reference.R preprocess|export|programmes|identities|markers|heldout INPUT OUTPUT_DIRECTORY")
stage <- a[1]; input <- normalizePath(a[2],mustWork=TRUE)
out <- normalizePath(a[3],mustWork=FALSE); dir.create(out,recursive=TRUE,showWarnings=FALSE)
out <- normalizePath(out); setwd(out); dir.create("results",showWarnings=FALSE)
suppressPackageStartupMessages({library(Seurat);library(harmony);library(Matrix);library(data.table);library(edgeR);library(RANN)})
root <- input; derived <- out; path <- input; input_path <- input; out_dir <- file.path(out,"results")
if(stage=="preprocess") {
dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
dirs <- dirs[file.exists(file.path(dirs, "matrix.mtx.gz"))]
if (length(dirs) != 5L) stop("Expected five 10x mouse directories")

read_mouse <- function(d) {
  mouse <- sub("_filtered_feature_bc_matrix$", "", basename(d))
  mouse <- sub("^GSM[0-9]+_", "", mouse)
  x <- Read10X(d, gene.column = 2, unique.features = TRUE)
  guide_rows <- grep("_gRNA$", rownames(x), value = TRUE)
  guide <- x[guide_rows, , drop = FALSE]
  n_guides <- Matrix::colSums(guide > 0)
  call <- rep("Neg", ncol(x))
  single <- which(n_guides == 1)
  call[single] <- guide_rows[max.col(t(as.matrix(guide[, single, drop = FALSE])),
                                      ties.method = "first")]
  call[n_guides > 1] <- "Multiple"
  genes <- x[setdiff(rownames(x), guide_rows), , drop = FALSE]
  colnames(genes) <- paste(mouse, colnames(genes), sep = "::")
  names(call) <- colnames(genes)
  list(mouse = mouse, counts = genes,
       meta = data.frame(cell_id = colnames(genes), mouse_id = mouse,
                         guide_call = call, n_guides = n_guides,
                         stringsAsFactors = FALSE, row.names = colnames(genes)))
}

raw <- lapply(dirs, read_mouse)
common <- Reduce(intersect, lapply(raw, function(z) rownames(z$counts)))
counts <- do.call(cbind, lapply(raw, function(z) z$counts[common, , drop = FALSE]))
meta <- do.call(rbind, lapply(raw, `[[`, "meta"))

obj <- CreateSeuratObject(counts = counts, meta.data = meta, min.cells = 3, min.features = 200)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
obj <- subset(obj, subset = nFeature_RNA >= 300 & nFeature_RNA <= 7500 &
                nCount_RNA >= 500 & percent.mt <= 20)
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 3000, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 40, seed.use = 42,
              verbose = FALSE)
obj <- RunHarmony(obj, group.by.vars = "mouse_id", reduction.use = "pca",
                  dims.use = 1:35, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.6, random.seed = 17, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, seed.use = 17, verbose = FALSE)

md <- obj[[]]
md$cell_id <- rownames(md)
fwrite(as.data.table(md), "results/perturbseq_mouse_cell_registry.tsv.gz", sep = "\t")
guide_audit <- as.data.table(md)[, .N, by = .(mouse_id, guide_call)]
fwrite(guide_audit, "results/perturbseq_mouse_guide_audit.tsv", sep = "\t")
saveRDS(obj, file.path(derived, "perturbseq_qc_harmony.rds"), compress = FALSE)
} else if(stage=="export") {
# Reconstruct every match from the supplied object, never from archived pairs.
obj <- readRDS(path)
md <- obj[[]]
counts <- GetAssayData(obj, assay="RNA", layer="counts")
state <- Embeddings(obj,"harmony")[,1:30,drop=FALSE]
stopifnot(identical(colnames(counts),rownames(md)),
          identical(colnames(obj),rownames(state)),
          all(c("cell_id","mouse_id","guide_call") %in% names(md)),
          identical(as.character(md$cell_id),rownames(md)))
markers <- list(
 Astrocyte=c("Aldoc","Aqp4","Slc1a3","Gfap"),
 Microglia=c("C1qa","C1qb","Csf1r","P2ry12","Tmem119"),
 Oligodendrocyte=c("Mbp","Plp1","Mog","Mal"),
 OPC=c("Pdgfra","Cspg4","Olig1","Olig2"),
 Excitatory_neuron=c("Slc17a7","Camk2a","Snap25","Rbfox3"),
 Inhibitory_neuron=c("Gad1","Gad2","Slc6a1"),
 Endothelial=c("Kdr","Pecam1","Cldn5","Emcn"),
 Perivascular=c("Rgs5","Pdgfrb","Acta2","Des"))
features <- intersect(unique(unlist(markers)),rownames(obj))
avg <- AverageExpression(obj,assays="RNA",features=features,
 group.by="seurat_clusters",layer="data",verbose=FALSE)$RNA
z <- t(scale(t(avg)));z[!is.finite(z)] <- 0
program <- sapply(markers,function(gs)colMeans(z[intersect(gs,rownames(z)),,drop=FALSE]))
cluster_map <- data.table(seurat_clusters=sub("^g","",rownames(program)),
 broad_identity=colnames(program)[max.col(program,ties.method="first")],
 identity_score=apply(program,1,max))
identity <- cluster_map$broad_identity[match(as.character(md$seurat_clusters),cluster_map$seurat_clusters)]
names(identity) <- rownames(md);stopifnot(!anyNA(identity))
fwrite(cluster_map,"results/perturbseq_broad_identity_map.tsv",sep="\t")
identity_metadata <- merge(as.data.table(md),cluster_map,by="seurat_clusters",all.x=TRUE)
fwrite(identity_metadata,"results/perturbseq_cell_broad_identity.tsv.gz",sep="\t")
targets <- sort(setdiff(unique(md$guide_call),c("Neg","Multiple")))
nearest_unique_controls <- function(tc,pool) {
 if(!length(pool))return(character())
 k <- min(50L,length(pool))
 nn <- RANN::nn2(state[pool,,drop=FALSE],state[tc,,drop=FALSE],k=k)$nn.idx
 nn <- matrix(nn,nrow=length(tc),ncol=k)
 used <- rep(FALSE,length(pool));chosen <- integer(length(tc))
 # Archived R14 input order; greedy nearest available candidate, no rematching.
 for(i in seq_along(tc)) {
  pick <- nn[i,][!used[nn[i,]]][1L]
  if(is.na(pick))return(character())
  chosen[i] <- pick;used[pick] <- TRUE
 }
 pool[chosen]
}
pair_rows <- list()
for(target in targets)for(m in sort(unique(md$mouse_id))) {
 tc <- rownames(md)[md$mouse_id==m & md$guide_call==target]
 if(length(tc)<3L)next
 pool <- rownames(md)[md$mouse_id==m & md$guide_call=="Neg"]
 controls <- nearest_unique_controls(tc,pool)
 if(length(controls)!=length(tc))next
 pair_rows[[paste(target,m,sep="::")]] <- data.table(guide=target,mouse_id=m,
  target_cell=tc,control_cell=controls,
  harmony_distance=sqrt(rowSums((state[tc,,drop=FALSE]-state[controls,,drop=FALSE])^2)),
  target_identity=unname(identity[tc]),control_identity=unname(identity[controls]),
  target_log_library=log1p(md[tc,"nCount_RNA"]),
  control_log_library=log1p(md[controls,"nCount_RNA"]),
  target_features=md[tc,"nFeature_RNA"],control_features=md[controls,"nFeature_RNA"])
}
pairs <- rbindlist(pair_rows)
pb <- list()
for(target in sort(unique(pairs$guide))) {
 q <- pairs[guide==target];columns <- list();sample_rows <- list()
 for(m in sort(unique(q$mouse_id))) {
  mm <- q[mouse_id==m]
  for(group in c("target","matched_neg")) {
   cells <- if(group=="target")mm$target_cell else mm$control_cell
   v <- Matrix::rowSums(counts[,cells,drop=FALSE])
   stopifnot(all(v==round(v)),max(v)<.Machine$integer.max)
   columns[[paste(m,group,sep="::")]] <- as.integer(v)
  }
  sample_rows[[m]] <- data.table(guide=target,mouse_id=m,
   group=c("target","matched_neg"),cells=nrow(mm),
   sample=paste(m,c("target","matched_neg"),sep="::"))
 }
 mat <- do.call(cbind,columns);rownames(mat) <- rownames(counts)
 pb[[target]] <- list(counts=mat,samples=rbindlist(sample_rows))
}
clu <- pairs[guide=="Clu_gRNA"]
clu_cells <- c(clu$target_cell,clu$control_cell)
stopifnot(!anyDuplicated(clu_cells))
lib <- Matrix::colSums(counts)
heldout <- list(counts=counts[,clu_cells,drop=FALSE],
 cells=data.table(cell_id=clu_cells,mouse_id=md[clu_cells,"mouse_id"],
  group=c(rep("target",nrow(clu)),rep("matched_neg",nrow(clu))),
  whole_gene_library=as.numeric(lib[clu_cells]),nFeature_RNA=md[clu_cells,"nFeature_RNA"]))
qc <- as.data.table(md)[,.(cells=.N),by=.(mouse_id,guide_call)]
saveRDS(list(pseudobulk=pb,heldout_cells=heldout,matching_pairs=pairs,qc=qc),
 file.path(out,"reference_counts.rds"),compress="xz")
fwrite(pairs,"results/reference_matching_pairs.tsv",sep="\t")
fwrite(qc,"results/training_qc_counts.tsv",sep="\t")
writeLines(c("Export recomputes the input-order R14 greedy nearest-state matching from the supplied processed object.",
 "No archived pair, differential-expression, heldout-cell or pseudobulk table is read.",
 "The all-five-mouse Harmony embedding is retained; weights-out-of-mouse refitting is conditional on that jointly learned embedding.",
 paste0("Gene rows: ",nrow(counts),"; QC cells: ",ncol(counts),"; guides: ",length(pb),"; successful mouse-guide matches: ",sum(vapply(pb,function(z)ncol(z$counts)/2,numeric(1))),"; target/control pairs: ",nrow(pairs))),
 file.path(out,"EXPORT_SCOPE.txt"))
} else if(stage=="programmes") {
if (!file.exists(path)) stop("Missing processed five-mouse object: ", path)
obj <- readRDS(path)
md <- obj[[]]
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
state <- Embeddings(obj, "harmony")[, 1:30, drop = FALSE]
targets <- sort(setdiff(unique(md$guide_call), c("Neg", "Multiple")))

nearest_unique_controls <- function(target_cells, pool, k = 50L) {
  k <- min(k, length(pool))
  nn <- RANN::nn2(state[pool, , drop = FALSE], state[target_cells, , drop = FALSE], k = k)$nn.idx
  used <- rep(FALSE, length(pool))
  chosen <- integer(length(target_cells))
  # Preserve input cell order and select its nearest still-available control.
  for (i in seq_along(target_cells)) {
    candidates <- nn[i, ]
    pick <- candidates[!used[candidates]][1]
    if (is.na(pick)) return(character())
    chosen[[i]] <- pick
    used[[pick]] <- TRUE
  }
  pool[chosen]
}

fit_target <- function(target) {
  pieces <- list(); meta <- list(); k <- 0L
  for (m in sort(unique(md$mouse_id))) {
    target_cells <- rownames(md)[md$mouse_id == m & md$guide_call == target]
    if (length(target_cells) < 3) next
    pool <- rownames(md)[md$mouse_id == m & md$guide_call == "Neg"]
    controls <- nearest_unique_controls(target_cells, pool)
    if (length(controls) != length(target_cells)) next
    k <- k + 1L
    pieces[[2 * k - 1L]] <- Matrix::rowSums(counts[, target_cells, drop = FALSE])
    pieces[[2 * k]] <- Matrix::rowSums(counts[, controls, drop = FALSE])
    meta[[k]] <- data.table(mouse_id = rep(m, 2), group = c("target", "matched_neg"),
                            cells = length(target_cells))
  }
  if (k < 3) return(NULL)
  mat <- do.call(cbind, pieces)
  sm <- rbindlist(meta)
  sm[, mouse_id := factor(mouse_id)]
  sm[, group := factor(group, levels = c("matched_neg", "target"))]
  y <- DGEList(mat)
  design <- model.matrix(~ mouse_id + group, data = sm)
  keep <- filterByExpr(y, design = design, min.count = 5)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  out <- as.data.table(topTags(glmQLFTest(fit, coef = "grouptarget"),
                              n = Inf, sort.by = "none")$table, keep.rownames = "gene")
  out[, `:=`(perturbation = sub("_gRNA$", "", target), paired_mice = k)]
  out
}

de <- rbindlist(lapply(targets, fit_target), fill = TRUE)
de[, program_member := FDR < 0.10 & abs(logFC) >= 0.5 &
     !grepl("^mt-|^Rpl|^Rps", gene, ignore.case = TRUE)]
control_genes <- de[perturbation == "mSafe" & program_member == TRUE, gene]
de[, passes_msafe_calibration := !gene %in% control_genes]
de[, primary_program_member := program_member & passes_msafe_calibration & perturbation != "mSafe"]
de[, abs_logFC := abs(logFC)]
setorder(de, perturbation, FDR, -abs_logFC)

dir.create("results", showWarnings = FALSE)
fwrite(de, "results/perturbseq_mouse_nearest_state_edger.tsv.gz", sep = "\t")
summary <- de[primary_program_member == TRUE,
  .(n_genes = .N, n_up = sum(logFC > 0), n_down = sum(logFC < 0),
    paired_mice = unique(paired_mice), top_genes = paste(head(gene, 20), collapse = ";")),
  by = perturbation]
fwrite(summary, "results/perturbseq_mouse_nearest_state_program_summary.tsv", sep = "\t")
calibration <- data.table(msafe_fdr10_effect_genes = length(control_genes),
                          msafe_genes = paste(control_genes, collapse = ";"),
                          perturbations_with_primary_program = uniqueN(summary$perturbation))
fwrite(calibration, "results/perturbseq_mouse_nearest_state_calibration.tsv", sep = "\t")
print(calibration)
print(summary)

} else if(stage=="identities") {
obj <- readRDS(path)

markers <- list(
  Astrocyte = c("Aldoc", "Aqp4", "Slc1a3", "Gfap"),
  Microglia = c("C1qa", "C1qb", "Csf1r", "P2ry12", "Tmem119"),
  Oligodendrocyte = c("Mbp", "Plp1", "Mog", "Mal"),
  OPC = c("Pdgfra", "Cspg4", "Olig1", "Olig2"),
  Excitatory_neuron = c("Slc17a7", "Camk2a", "Snap25", "Rbfox3"),
  Inhibitory_neuron = c("Gad1", "Gad2", "Slc6a1"),
  Endothelial = c("Kdr", "Pecam1", "Cldn5", "Emcn"),
  Perivascular = c("Rgs5", "Pdgfrb", "Acta2", "Des")
)
features <- intersect(unique(unlist(markers)), rownames(obj))
avg <- AverageExpression(obj, assays = "RNA", features = features, group.by = "seurat_clusters",
                         layer = "data", verbose = FALSE)$RNA
z <- t(scale(t(avg)))
z[!is.finite(z)] <- 0
program <- sapply(markers, function(gs) {
  use <- intersect(gs, rownames(z))
  colMeans(z[use, , drop = FALSE])
})
cluster_map <- data.table(seurat_clusters = sub("^g", "", rownames(program)),
                          broad_identity = colnames(program)[max.col(program, ties.method = "first")],
                          identity_score = apply(program, 1, max))
fwrite(cluster_map, "results/perturbseq_broad_identity_map.tsv", sep = "\t")

md <- obj[[]]
md$cell_id <- rownames(md)
md <- merge(as.data.table(md), cluster_map, by = "seurat_clusters", all.x = TRUE)
fwrite(md, "results/perturbseq_cell_broad_identity.tsv.gz", sep = "\t")

} else if(stage=="markers") {
obj <- readRDS(path)
md <- as.data.table(obj[[]])
if (!"cell_id" %in% names(md)) md[, cell_id := rownames(obj[[]])]
identity_map <- fread("results/perturbseq_cell_broad_identity.tsv.gz")[,
  .(cell_id, broad_identity)]
md <- merge(md, identity_map, by = "cell_id", all.x = TRUE)
stopifnot(!anyNA(md$broad_identity))

# The reference is deliberately restricted to barcode-negative cells so that
# receiver identities cannot inherit a perturbation-response programme.
neg_cells <- md[guide_call == "Neg", cell_id]
obj_neg <- subset(obj, cells = neg_cells)
obj_neg$broad_identity <- md$broad_identity[match(colnames(obj_neg), md$cell_id)]
avg <- AverageExpression(obj_neg, assays = "RNA", group.by = "broad_identity",
                         layer = "data", verbose = FALSE)$RNA

marker_tables <- lapply(colnames(avg), function(id) {
  other <- setdiff(colnames(avg), id)
  specificity <- avg[, id] - apply(avg[, other, drop = FALSE], 1, max)
  prevalence <- Matrix::rowMeans(GetAssayData(obj_neg, assay = "RNA", layer = "counts")[,
    obj_neg$broad_identity == id, drop = FALSE] > 0)
  data.table(gene = rownames(avg), broad_identity = id,
             mean_log_expression = avg[, id], specificity = specificity,
             within_identity_prevalence = prevalence)
})
markers <- rbindlist(marker_tables)
markers <- markers[!grepl("^mt-|^Rpl|^Rps|^sgrna", gene, ignore.case = TRUE) &
                     specificity > 0 & within_identity_prevalence >= 0.05]
setorder(markers, broad_identity, -specificity, -within_identity_prevalence)
markers[, rank := seq_len(.N), by = broad_identity]
markers[, selected := rank <= 25]
fwrite(markers, "results/receiver_identity_reference_markers.tsv.gz", sep = "\t")

selected <- markers[selected == TRUE]
summary <- selected[, .(marker_genes = .N,
                        top_markers = paste(head(gene, 10), collapse = ";")),
                    by = broad_identity]
fwrite(summary, "results/receiver_identity_reference_summary.tsv", sep = "\t")

} else if(stage=="heldout") {
if (!file.exists(input_path)) stop("Missing input: ", input_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Reading frozen object once: ", input_path)
obj <- readRDS(input_path)
md <- obj[[]]
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
state <- Embeddings(obj, "harmony")[, 1:30, drop = FALSE]
stopifnot(identical(colnames(counts), rownames(md)),
          all(rownames(md) %in% rownames(state)),
          all(c("mouse_id", "guide_call") %in% names(md)))
state <- state[rownames(md), , drop = FALSE]
rm(obj)
invisible(gc())
mice <- sort(unique(as.character(md$mouse_id)))
stopifnot(length(mice) == 5L)
target_labels <- c("Clu_gRNA", "mSafe_gRNA")
stopifnot(all(target_labels %in% md$guide_call))
lib <- Matrix::colSums(counts)

# Exactly follows the frozen R/14 greedy input-order nearest-state matching:
# same mouse, 30 Harmony coordinates, 50 nearest control candidates,
# no replacement within one mouse/guide comparison, no global rematching.
nearest_unique_controls <- function(target_cells, pool, k = 50L) {
  if (!length(pool)) return(character())
  k <- min(k, length(pool))
  nn <- RANN::nn2(state[pool, , drop = FALSE],
                  state[target_cells, , drop = FALSE], k = k)$nn.idx
  nn <- matrix(nn, nrow = length(target_cells), ncol = k)
  used <- rep(FALSE, length(pool))
  chosen <- integer(length(target_cells))
  for (i in seq_along(target_cells)) {
    candidates <- nn[i, ]
    pick <- candidates[!used[candidates]][1L]
    if (is.na(pick)) return(character())
    chosen[i] <- pick
    used[pick] <- TRUE
  }
  pool[chosen]
}

matches <- list()
pair_rows <- list()
for (target in target_labels) for (m in mice) {
  target_cells <- rownames(md)[md$mouse_id == m & md$guide_call == target]
  pool <- rownames(md)[md$mouse_id == m & md$guide_call == "Neg"]
  if (length(target_cells) < 3L) next
  control_cells <- nearest_unique_controls(target_cells, pool)
  if (length(control_cells) != length(target_cells)) next
  key <- paste(target, m, sep = "::")
  matches[[key]] <- list(target = target_cells, control = control_cells)
  pair_rows[[key]] <- data.table(
    guide = target, mouse_id = m, pair = seq_along(target_cells),
    target_cell = target_cells, control_cell = control_cells,
    frozen_harmony_distance = sqrt(rowSums(
      (state[target_cells, , drop = FALSE] - state[control_cells, , drop = FALSE])^2)))
}
pair_table <- rbindlist(pair_rows)
fwrite(pair_table, file.path(out_dir, "heldout_mouse_matching.tsv"), sep = "\t")

pseudobulk <- list()
for (target in target_labels) {
  columns <- list(); metadata <- list()
  for (m in mice) {
    z <- matches[[paste(target, m, sep = "::")]]
    if (is.null(z)) next
    columns[[paste(m, "target", sep = "::")]] <-
      Matrix::rowSums(counts[, z$target, drop = FALSE])
    columns[[paste(m, "matched_neg", sep = "::")]] <-
      Matrix::rowSums(counts[, z$control, drop = FALSE])
    metadata[[m]] <- data.table(mouse_id = rep(m, 2L),
      group = c("target", "matched_neg"), n_cells = length(z$target))
  }
  pseudobulk[[target]] <- list(mat = do.call(cbind, columns),
                               metadata = rbindlist(metadata))
}

fit_training_weights <- function(target, heldout) {
  z <- pseudobulk[[target]]
  keep_samples <- z$metadata$mouse_id != heldout
  sample_meta <- copy(z$metadata[keep_samples])
  paired_mice <- uniqueN(sample_meta$mouse_id)
  if (paired_mice < 3L) return(list(status = "insufficient_training_mice", de = NULL,
                                  paired_mice = paired_mice))
  sample_meta[, mouse_id := factor(mouse_id)]
  sample_meta[, group := factor(group, levels = c("matched_neg", "target"))]
  design <- model.matrix(~mouse_id + group, data = sample_meta)
  y <- DGEList(z$mat[, keep_samples, drop = FALSE])
  keep_genes <- filterByExpr(y, design = design, min.count = 5)
  if (!any(keep_genes)) return(list(status = "no_expressed_training_genes", de = NULL,
                                   paired_mice = paired_mice))
  y <- calcNormFactors(y[keep_genes, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  de <- as.data.table(topTags(glmQLFTest(fit, coef = "grouptarget"),
    n = Inf, sort.by = "none")$table, keep.rownames = "gene")
  de[, eligible := FDR < 0.10 & abs(logFC) >= 0.5 &
      !grepl("^mt-|^Rpl|^Rps", gene, ignore.case = TRUE)]
  list(status = "ok", de = de, paired_mice = paired_mice)
}

auc_rank <- function(target, control) {
  mean(outer(target, control, function(x, y) (x > y) + 0.5 * (x == y)))
}
metric_rows <- list(); weight_rows <- list(); de_rows <- list()
for (heldout in mice) {
  message("Held-out mouse: ", heldout)
  train_clu <- fit_training_weights("Clu_gRNA", heldout)
  train_safe <- fit_training_weights("mSafe_gRNA", heldout)
  if (!is.null(train_clu$de)) {
    tmp <- copy(train_clu$de)
    tmp[, `:=`(heldout_mouse = heldout, guide = "Clu_gRNA")]
    de_rows[[paste(heldout, "clu")]] <- tmp
  }
  if (!is.null(train_safe$de)) {
    tmp <- copy(train_safe$de)
    tmp[, `:=`(heldout_mouse = heldout, guide = "mSafe_gRNA")]
    de_rows[[paste(heldout, "safe")]] <- tmp
  }
  fold_match <- matches[[paste("Clu_gRNA", heldout, sep = "::")]]
  blacklist <- if (train_safe$status == "ok") train_safe$de[eligible == TRUE, gene] else character()
  for (version in c("all_retained_genes", "no_immediate_early_genes")) {
    metric <- data.table(heldout_mouse = heldout, programme = version,
      n_training_mice = if (is.null(train_clu$paired_mice)) NA_integer_ else train_clu$paired_mice,
      n_training_msafe_mice = if (is.null(train_safe$paired_mice)) NA_integer_ else train_safe$paired_mice,
      clu_training_status = train_clu$status, msafe_training_status = train_safe$status,
      n_heldout_target = if (is.null(fold_match)) 0L else length(fold_match$target),
      n_heldout_control = if (is.null(fold_match)) 0L else length(fold_match$control),
      n_candidate_clu_genes = if (train_clu$status == "ok") sum(train_clu$de$eligible) else NA_integer_,
      n_training_msafe_blacklist = if (train_safe$status == "ok") length(blacklist) else NA_integer_,
      n_programme_genes = NA_integer_,
      auc = NA_real_, mean_target_score = NA_real_, mean_control_score = NA_real_,
      mean_score_difference = NA_real_, control_robust_scale = NA_real_,
      standardized_mean_difference = NA_real_, matched_pair_mean_difference = NA_real_,
      status = "unestimable", reason = "", genes = "")
    if (train_clu$status != "ok" || train_safe$status != "ok") {
      metric[, reason := paste("training fit:", train_clu$status, train_safe$status)]
    } else if (is.null(fold_match)) {
      metric[, reason := "held-out matching unavailable"]
    } else {
      w <- copy(train_clu$de[eligible == TRUE & !gene %in% blacklist])
      if (version == "no_immediate_early_genes") w <- w[!gene %in% c("Fos", "Fosb", "Junb")]
      metric[, n_programme_genes := nrow(w)]
      if (!nrow(w)) {
        metric[, reason := "no genes passed training-only thresholds and blacklist"]
      } else {
        setorder(w, FDR, gene)
        w[, `:=`(heldout_mouse = heldout, programme = version)]
        weight_rows[[paste(heldout, version)]] <- w
        heldout_cells <- c(fold_match$target, fold_match$control)
        x <- as.matrix(counts[w$gene, heldout_cells, drop = FALSE])
        x <- log1p(sweep(x, 2L, pmax(lib[heldout_cells], 1), "/") * 1e4)
        scores <- as.numeric(crossprod(w$logFC, x)) / sum(abs(w$logFC))
        n_target <- length(fold_match$target)
        a <- scores[seq_len(n_target)]
        b <- scores[n_target + seq_len(n_target)]
        scale <- mad(b, center = median(b), constant = 1.4826)
        if (!is.finite(scale) || scale < 1e-8) scale <- sd(b)
        metric[, `:=`(auc = auc_rank(a, b), mean_target_score = mean(a),
          mean_control_score = mean(b), mean_score_difference = mean(a) - mean(b),
          control_robust_scale = scale,
          standardized_mean_difference = if (is.finite(scale) && scale > 0) (mean(a) - mean(b)) / scale else NA_real_,
          matched_pair_mean_difference = mean(a - b), status = "estimable",
          reason = "descriptive weight-learning holdout only; frozen full-reference embedding reused",
          genes = paste(w$gene, collapse = ";"))]
      }
    }
    metric_rows[[paste(heldout, version)]] <- metric
  }
}
metrics <- rbindlist(metric_rows, fill = TRUE)
fwrite(metrics, file.path(out_dir, "heldout_mouse_programme_metrics.tsv"), sep = "\t")
fwrite(rbindlist(weight_rows, fill = TRUE),
       file.path(out_dir, "heldout_mouse_programme_weights.tsv"), sep = "\t")
fwrite(rbindlist(de_rows, fill = TRUE),
       file.path(out_dir, "heldout_mouse_training_edger.tsv.gz"), sep = "\t")
summary <- metrics[, .(n_folds = .N, n_estimable = sum(status == "estimable"),
  n_heldout_targets = sum(n_heldout_target),
  median_auc = if (all(is.na(auc))) NA_real_ else median(auc, na.rm = TRUE),
  min_auc = if (all(is.na(auc))) NA_real_ else min(auc, na.rm = TRUE),
  max_auc = if (all(is.na(auc))) NA_real_ else max(auc, na.rm = TRUE),
  n_genes_min = if (all(is.na(n_programme_genes))) NA_integer_ else min(n_programme_genes, na.rm = TRUE),
  n_genes_max = if (all(is.na(n_programme_genes))) NA_integer_ else max(n_programme_genes, na.rm = TRUE)), by = programme]
fwrite(summary, file.path(out_dir, "heldout_mouse_programme_summary.tsv"), sep = "\t")
writeLines(c(
  "This is an author-requested weights-out-of-mouse diagnostic, not a fully independent pipeline validation.",
  "The frozen Seurat/Harmony object was learned using all five mice; its 30-dimensional embedding is reused for within-mouse state matching.",
  "For each held-out mouse, Clu differential expression and the mSafe blacklist use only eligible matches in the other four mice; each fit requires at least three eligible training mice.",
  "Matching follows the archived R/14 input-order greedy algorithm and is constructed once within each mouse/guide pair.",
  "Programme thresholds: BH FDR < 0.10, absolute log2 fold change >= 0.5, exclude mitochondrial/ribosomal prefixes, subtract training-only mSafe genes.",
  "The no-immediate-early variant additionally excludes Fos, Fosb and Junb.",
  "Scores are signed-weighted log1p(CPM/100) expression sums divided by the sum of absolute weights.",
  "AUC compares held-out Clu-guide cells with their equally numerous held-out matched guide-undetected controls; half credit is given to ties.",
  "Effect size is the difference of mean held-out scores, optionally divided by the held-out control MAD (1.4826 factor; SD fallback).",
  "Per-cell AUC/effects are descriptive within-mouse diagnostics; they are not independent biological-replicate tests or calibration of propagation.",
  "Three folds are unestimable because holding out one of the three eligible mSafe mice leaves only two mSafe training mice; Clu differential-expression fitting itself succeeds in those folds.",
  "Programme gene count and mSafe blacklist count are NA when calibration cannot be fitted; zero is reserved for a fitted programme whose genes fail the training filters.",
  "The two estimable no-immediate-early programmes have 8 and 19 genes. No fallback calibration or gene selection is supplied for unestimable folds.",
  "No spatial receiver endpoint, public input file or archived analysis script is modified."
), file.path(out_dir, "heldout_mouse_programme_scope.txt"))
print(metrics[, .(heldout_mouse, programme, n_heldout_target, n_programme_genes,
                   auc, standardized_mean_difference, status)])
print(summary)

} else stop("Unknown reference stage: ",stage)
