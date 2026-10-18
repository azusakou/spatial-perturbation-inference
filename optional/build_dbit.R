#!/usr/bin/env Rscript
# Rebuild DBiT score-level inputs from exact GEO filenames and frozen Hallmarks.
a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=2L)stop("Usage: Rscript optional/build_dbit.R GSE319277_EXTRACTED OUTPUT_RDS")
script<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1]);ROOT<-normalizePath(file.path(dirname(script),".."))
raw_dir<-normalizePath(a[1]);output<-normalizePath(a[2],mustWork=FALSE)
suppressPackageStartupMessages({library(data.table);library(mgcv)})
frozen<-readRDS(file.path(ROOT,"data/dbit_scores.rds"))
hallmark<-fread(file.path(ROOT,"data/hallmark_2026_1.tsv.gz"));frozen$atg4a$gene_sets<-split(hallmark$gene_symbol,hallmark$gs_name)
make_atg4a <- function() {
gene_sets<-frozen$atg4a$gene_sets;all_genes<-sort(unique(unlist(gene_sets)))
read_spot_by_gene <- function(path, genes) {
  con <- gzfile(path, "rt")
  header <- strsplit(readLines(con, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
  close(con)
  wanted <- paste0(genes, "__exon")
  idx <- match(wanted, header)
  keep <- which(!is.na(idx))
  x <- fread(path, select = c(1L, idx[keep]), showProgress = FALSE)
  barcodes <- x[[1L]]
  m <- as.matrix(x[, -1L])
  storage.mode(m) <- "numeric"
  colnames(m) <- genes[keep]
  rownames(m) <- barcodes
  m
}

read_gene_by_spot <- function(path, genes) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  header <- readLines(con, n = 1L, warn = FALSE)
  barcodes <- gsub('"', "", strsplit(header, "\t", fixed = TRUE)[[1L]], fixed = TRUE)
  rows <- list(); found <- character()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break
    tab <- regexpr("\t", line, fixed = TRUE)[[1L]]
    if (tab < 1L) next
    gene <- gsub('"', "", substr(line, 1L, tab - 1L), fixed = TRUE)
    if (!gene %chin% genes) next
    values <- as.numeric(strsplit(substr(line, tab + 1L, nchar(line)), "\t",
                                  fixed = TRUE)[[1L]])
    if (length(values) == length(barcodes)) {
      rows[[length(rows) + 1L]] <- values
      found <- c(found, gene)
    }
  }
  m <- t(do.call(rbind, rows))
  rownames(m) <- barcodes
  colnames(m) <- found
  m
}

read_position <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con))
  x <- strsplit(readLines(con, n = 1L, warn = FALSE), ",", fixed = TRUE)[[1L]]
  x[nzchar(x)]
}

read_guides <- function(path) {
  x <- fread(path, header = FALSE, fill = TRUE,
             col.names = c("barcode", "guide", "count", "extra"))
  x[, barcode := sub("^bc", "", barcode)]
  x[, target := sub("^sg[0-9]+_", "", guide)]
  x[!is.na(barcode) & !is.na(target) & !is.na(count)]
}

make_section <- function(sample) {
  if (sample == "Bru.LC.6") {
    mat <- read_spot_by_gene(file.path(raw_dir, "GSM9514752_Bru.LC.6_allRNA.tsv.gz"),
                             all_genes)
    position <- read_position(file.path(raw_dir, "GSM9514752_position_Bru.LC.6.txt.gz"))
    guides <- read_guides(file.path(raw_dir,
      "GSM9514752_sgCounts_Bru.LC.6_grnalib.txt.gz"))
  } else {
    mat <- read_gene_by_spot(file.path(raw_dir, "GSM9514756_Bru.LC.8.tsv.gz"), all_genes)
    position <- read_position(file.path(raw_dir, "GSM9514756_position_Bru.LC.8.txt.gz"))
    guides <- read_guides(file.path(raw_dir,
      "GSM9514756_sgCounts_Bru.LC.8_grnalib.txt.gz"))
  }
  position <- intersect(position, rownames(mat))
  mat <- mat[position, , drop = FALSE]
  lib_proxy <- pmax(rowSums(mat), 1)
  norm <- log1p(1e4 * mat / lib_proxy)
  dt <- data.table(barcode = rownames(mat))
  xy <- tstrsplit(dt$barcode, "x", fixed = TRUE)
  dt[, `:=`(x = as.numeric(xy[[1L]]) * 50, y = as.numeric(xy[[2L]]) * 50,
            lib_proxy = lib_proxy, any_guide = barcode %chin% guides$barcode,
            atg4a = barcode %chin% guides[target == "ATG4A", barcode], sample = sample)]
  for (nm in names(gene_sets)) {
    g <- intersect(gene_sets[[nm]], colnames(norm))
    dt[[nm]] <- rowMeans(norm[, g, drop = FALSE])
  }
  list(dt = dt, guides = guides,
       gene_coverage = data.table(sample, programme = names(gene_sets),
         genes_total = lengths(gene_sets), genes_observed = vapply(gene_sets, function(g)
           length(intersect(g, colnames(norm))), integer(1))))
}

s6<-make_section("Bru.LC.6");s8<-make_section("Bru.LC.8")
list(s6=s6$dt,s8=s8$dt,programme_names=names(gene_sets),gene_sets=gene_sets,coverage=rbind(s6$gene_coverage,s8$gene_coverage))
}
make_tgfbr2 <- function() {
target_genes <- c("Serpine1", "Ctgf", "Smad7", "Acta2", "Tagln", "Timp1",
                  "Col1a1", "Col1a2", "Col3a1", "Fn1", "Tnc")
housekeeping <- c("Actb", "Gapdh", "B2m", "Rplp0", "Rpl13a", "Rps18", "Rps27a",
                  "Hprt", "Ppia", "Tbp", "Gusb", "Pgk1", "Ywhaz", "Ubc")
genes <- unique(c(target_genes, housekeeping))

read_gene_by_spot <- function(path, genes) {
  con <- gzfile(path, "rt"); on.exit(close(con))
  barcodes <- gsub('"', "", strsplit(readLines(con, n = 1L, warn = FALSE),
    "\t", fixed = TRUE)[[1L]], fixed = TRUE)
  barcodes <- barcodes[nzchar(barcodes)]
  rows <- list(); found <- character()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break
    tab <- regexpr("\t", line, fixed = TRUE)[[1L]]
    if (tab < 1L) next
    gene <- gsub('"', "", substr(line, 1L, tab - 1L), fixed = TRUE)
    gene <- sub("^X(?=[0-9])", "", gene, perl = TRUE)
    if (!gene %chin% genes) next
    v <- as.numeric(strsplit(substr(line, tab + 1L, nchar(line)), "\t", fixed = TRUE)[[1L]])
    if (length(v) == length(barcodes)) { rows[[length(rows) + 1L]] <- v; found <- c(found, gene) }
  }
  m <- t(do.call(rbind, rows)); rownames(m) <- barcodes; colnames(m) <- found; m
}

read_spot_by_gene <- function(path, genes) {
  con <- gzfile(path, "rt")
  h <- strsplit(readLines(con, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
  close(con)
  idx <- match(genes, h); keep <- which(!is.na(idx))
  x <- fread(path, select = c(1L, idx[keep]), showProgress = FALSE)
  m <- as.matrix(x[, -1L]); storage.mode(m) <- "numeric"
  rownames(m) <- x[[1L]]; colnames(m) <- genes[keep]; m
}

read_position <- function(path) {
  z <- gzfile(path); on.exit(close(z))
  x <- strsplit(readLines(z, n = 1L, warn = FALSE), ",", fixed = TRUE)[[1L]]
  x[nzchar(x)]
}

read_guides <- function(path) {
  x <- fread(path, header = FALSE, fill = TRUE,
             col.names = c("barcode", "guide", "count", "extra"))
  x[, barcode := sub("^bc", "", barcode)]
  x[!is.na(barcode) & !is.na(guide)]
}

make_section <- function(sample) {
  if (sample == "mTSG.LV.2") {
    m <- read_gene_by_spot(file.path(raw_dir, "GSM9514753_mTSG.LV.2.tsv.gz"), genes)
    pos <- read_position(file.path(raw_dir, "GSM9514753_position_mTSG.LV.2.txt.gz"))
    guides <- read_guides(file.path(raw_dir, "GSM9514753_sgCounts_mTSG.LV.2.txt.gz"))
  } else {
    m <- read_spot_by_gene(file.path(raw_dir, "GSM9514754_mTSG.LV.3.updated.tsv.gz"), genes)
    pos <- read_position(file.path(raw_dir, "GSM9514754_position_mTSG.LV.3.txt.gz"))
    guides <- read_guides(file.path(raw_dir, "GSM9514754_sgCounts_mTSG.LV.3.txt.gz"))
  }
  pos <- intersect(pos, rownames(m)); m <- m[pos, , drop = FALSE]
  lib <- pmax(rowSums(m), 1); norm <- log1p(1e4 * m / lib)
  observed_targets <- intersect(target_genes, colnames(norm))
  dt <- data.table(barcode = rownames(m), score = rowMeans(norm[, observed_targets, drop = FALSE]),
                   lib_proxy = lib)
  xy <- tstrsplit(dt$barcode, "x", fixed = TRUE)
  dt[, `:=`(x = as.numeric(xy[[1L]]) * 50, y = as.numeric(xy[[2L]]) * 50,
            any_guide = barcode %chin% guides$barcode,
            tgfbr2 = barcode %chin% guides[grepl("^Tgfbr2", guide), barcode], sample = sample)]
  recv <- !dt$any_guide
  kxy <- min(40L, max(10L, floor(sum(recv) / 40)))
  fit <- gam(score ~ s(x, y, k = kxy) + s(log1p(lib_proxy), k = 5),
             data = dt[recv], method = "REML")
  dt[, adjusted := score - as.numeric(predict(fit, newdata = dt))]
  list(dt = dt, observed_targets = observed_targets)
}

a<-make_section("mTSG.LV.2");b<-make_section("mTSG.LV.3")
a$dt[,adjusted:=NULL];b$dt[,adjusted:=NULL];list(t2=a,t3=b)
}
frozen$atg4a<-make_atg4a();frozen$atg4a$gene_sets<-NULL;frozen$tgfbr2<-make_tgfbr2();saveRDS(frozen,output,compress="xz",version=2)
