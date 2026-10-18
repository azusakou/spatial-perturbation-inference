#!/usr/bin/env Rscript
# Reconstruct the full Clu cell/score/identity input from native GEF files and
# the bundled independently learned programme/marker weights.
a <- commandArgs(trailingOnly=TRUE)
if(length(a)!=2L)stop("Usage: Rscript optional/build_clu.R GEF_DIRECTORY OUTPUT_RDS")
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
base <- normalizePath(file.path(dirname(script),"..")); raw_dir <- normalizePath(a[1]); output <- normalizePath(a[2],mustWork=FALSE)
suppressPackageStartupMessages({library(rhdf5);library(Matrix);library(data.table)})
source(file.path(base,"R/clu.R")); d <- readRDS(file.path(base,"data/endpoints.rds"))
marker_ref <- d$marker_weights
clu <- d$programme_weights[perturbation=="Clu" & !gene %in% c("Fos","Fosb","Junb"),.(gene,weight=logFC)]
extract_chip <- function(path) {
  chip <- sub("[.]adjusted[.]cellbin[.]gef$", "", basename(path))
  chip <- sub("^GSM[0-9]+_", "", chip)
  message("Identity-adjusted geometry: ", chip)
  cell <- as.data.table(h5read(path, "/cellBin/cell"))
  gene <- as.data.table(h5read(path, "/cellBin/gene"))
  resolution_nm <- as.numeric(h5readAttributes(path, "/")$resolution)
  coords <- as.matrix(cell[, .(x, y)]) * resolution_nm / 1000
  cell[, component := component_grid(x, y)]

  guides <- gene[grepl("^sgrna[-_]", geneName, ignore.case = TRUE)]
  guide_names <- tolower(sub("^sgrna_", "sgrna-", guides$geneName))
  hits <- rbindlist(lapply(seq_len(nrow(guides)), function(i) {
    g <- guides[i]
    z <- as.data.table(h5read(path, "/cellBin/geneExp", start = as.numeric(g$offset) + 1,
                              count = as.numeric(g$cellCount)))
    z[count > 0, .(cellID, guide = guide_names[[i]])]
  }))
  calls <- hits[, .(n_guides = uniqueN(guide),
                    guide_call = if (uniqueN(guide) == 1) unique(guide) else "Multiple"),
                by = cellID]
  cell[, `:=`(n_guides = 0L, guide_call = "Neg")]
  ii <- match(calls$cellID, cell$id)
  cell$n_guides[ii] <- calls$n_guides; cell$guide_call[ii] <- calls$guide_call
  source <- which(cell$guide_call == "sgrna-clu")
  neg <- which(cell$guide_call == "Neg")
  stopifnot(length(source) >= 5)

  union_genes <- unique(c(marker_ref$gene, clu$gene))
  sm <- read_sparse_genes(path, union_genes, cell, gene)
  lib <- pmax(as.numeric(cell$expCount), 1)
  sm <- sm %*% Diagonal(x = 1 / lib); sm@x <- log1p(1e4 * sm@x)

  ids <- sort(unique(as.character(marker_ref$broad_identity)))
  marker_genes <- intersect(marker_ref$gene, rownames(sm))
  ref <- marker_ref[gene %in% marker_genes]
  mu <- Matrix::rowMeans(sm[marker_genes, neg, drop = FALSE])
  second <- Matrix::rowMeans(sm[marker_genes, neg, drop = FALSE]^2)
  sdev <- sqrt(pmax(second - mu^2, 1e-6))
  W <- sparseMatrix(i = match(ref$gene, marker_genes),
                    j = match(ref$broad_identity, ids), x = pmax(ref$specificity, 0),
                    dims = c(length(marker_genes), length(ids)),
                    dimnames = list(marker_genes, ids))
  W <- W %*% Diagonal(x = 1 / pmax(Matrix::colSums(W), 1e-8))
  Wz <- Diagonal(x = 1 / sdev) %*% W
  identity_score <- as.matrix(crossprod(Wz, sm[marker_genes, , drop = FALSE]))
  offset <- as.numeric(crossprod(Wz, mu))
  identity_score <- sweep(identity_score, 1, offset, "-")
  ord <- apply(identity_score, 2, order, decreasing = TRUE)
  identity <- ids[ord[1, ]]
  confidence <- identity_score[cbind(ord[1, ], seq_len(ncol(identity_score)))] -
    identity_score[cbind(ord[2, ], seq_len(ncol(identity_score)))]
  identity[!is.finite(confidence) | confidence <= 0] <- "Unresolved"

  w <- clu[gene %in% rownames(sm)]
  weight <- setNames(w$weight, w$gene)[rownames(sm)]
  weight[is.na(weight)] <- 0
  raw_score <- as.numeric(crossprod(weight, sm)) / sum(abs(weight))
  center <- median(raw_score[neg]); spread <- mad(raw_score[neg], constant = 1.4826)
  if (!is.finite(spread) || spread == 0) spread <- sd(raw_score[neg])
  raw_score <- (raw_score - center) / spread

  data.table(chip_id=chip,x_um=coords[,1],y_um=coords[,2],
    identity=factor(identity),guide_call=factor(cell$guide_call),raw_clu_score=raw_score,
    library_size=lib,component=as.integer(cell$component))
}
files <- sort(list.files(raw_dir,pattern="[.]gef$",full.names=TRUE))
stopifnot(length(files)==3L)
z <- rbindlist(lapply(files,extract_chip)); saveRDS(z,output,compress="xz",version=2)
message("Saved ",nrow(z)," cells to ",output)
