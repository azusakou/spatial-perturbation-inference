#!/usr/bin/env Rscript
# Optional raw five-fold validation; the heldout mouse never fits preprocessing.
a <- commandArgs(trailingOnly=TRUE)
if(!length(a)%in%c(2L,3L))stop("Usage: Rscript optional/full_holdout.R RAW_10X_DIRECTORY OUTPUT_DIRECTORY [HELDOUT_MOUSE]")
raw_root <- normalizePath(a[1]);out <- normalizePath(a[2],mustWork=FALSE)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
release <- normalizePath(file.path(dirname(script),".."))
if(length(a)==2L) {
 # Separate workers prevent multi-fold Seurat allocations accumulating in R.
 dd <- list.dirs(raw_root,recursive=FALSE,full.names=TRUE)
 dd <- dd[file.exists(file.path(dd,"matrix.mtx.gz"))]
 mice <- sub("^GSM[0-9]+_","",sub("_filtered_feature_bc_matrix$","",basename(dd)))
 if(length(mice)!=5L)stop("Expected five mouse directories")
 for(m in mice) {
  rc <- system2(file.path(R.home("bin"),"Rscript"),shQuote(c(normalizePath(script),raw_root,out,m)))
  if(rc!=0L)stop("Fold worker failed: ",m," (exit ",rc,")")
 }
 suppressPackageStartupMessages(library(data.table))
 input <- list(folds=setNames(lapply(mice,function(m)readRDS(file.path(out,m,"compact_fold_counts.rds"))),mice),
  original_primary_programmes=c("Cfap410","Clu","Dpp6","Fasn","Flcn","Gfap","Rbfox","Rraga","Stk39","Tbk1","Trem2"),
  scope=readLines(file.path(out,"FROZEN_PROTOCOL.txt")))
 for(n in c("fold_preprocessing_checks.tsv","fold_qc.tsv","stage_timings.tsv"))
  fwrite(rbindlist(lapply(mice,function(m)fread(file.path(out,m,n)))),file.path(out,n),sep="\t")
 saveRDS(input,file.path(out,"full_holdout_counts.rds"),compress="xz")
 source(file.path(release,"R/full_holdout.R"));run_full_holdout(input,file.path(out,"results"))
 writeLines(capture.output(sessionInfo()),file.path(out,"refit_sessionInfo.txt"))
 quit(status=0L)
}
suppressPackageStartupMessages({library(Seurat);library(harmony);library(Matrix);library(data.table);library(RANN);library(digest);library(matrixStats)})
dirs <- list.dirs(raw_root,recursive=FALSE,full.names=TRUE)
dirs <- dirs[file.exists(file.path(dirs,"matrix.mtx.gz"))]
if(length(dirs)!=5L)stop("Expected five mouse 10x directories")
files <- unlist(lapply(dirs,function(d)file.path(d,c("matrix.mtx.gz","features.tsv.gz","barcodes.tsv.gz"))))
fwrite(data.table(path=files,bytes=file.size(files),sha256=vapply(files,digest,character(1),file=TRUE,algo="sha256")),file.path(out,"raw_input_sha256.tsv"),sep="\t")
read_mouse <- function(d) {
 m <- sub("^GSM[0-9]+_","",sub("_filtered_feature_bc_matrix$","",basename(d)))
 x <- Read10X(d,gene.column=2,unique.features=TRUE)
 gr <- grep("_gRNA$",rownames(x),value=TRUE);g <- x[gr,,drop=FALSE];ng <- Matrix::colSums(g>0)
 call <- rep("Neg",ncol(x));single <- which(ng==1)
 call[single] <- gr[max.col(t(as.matrix(g[,single,drop=FALSE])),ties.method="first")];call[ng>1] <- "Multiple"
 x <- x[setdiff(rownames(x),gr),,drop=FALSE];colnames(x) <- paste(m,colnames(x),sep="::")
 list(mouse=m,counts=x,meta=data.frame(cell_id=colnames(x),mouse_id=m,guide_call=call,n_guides=ng,row.names=colnames(x)))
}
raw <- lapply(dirs,read_mouse);names(raw) <- vapply(raw,`[[`,character(1),"mouse")
all_guides <- sort(setdiff(unique(unlist(lapply(raw,function(z)z$meta$guide_call))),c("Neg","Multiple")))
primary <- c("Cfap410","Clu","Dpp6","Fasn","Flcn","Gfap","Rbfox","Rraga","Stk39","Tbk1","Trem2")
scope <- c(
 "Five heldout mice; four training mice rebuild gene filtering, cell QC, LogNormalize(10000), vst HVG(3000), scaling, PCA(40), Harmony(35), neighbors(30), clustering(0.6) and UMAP(30).",
 "No heldout count, cell, guide label, embedding or response is used to fit any training preprocessing, matching, DE, mSafe blacklist, weight, marker or threshold.",
 "Primary training_harmony route uses fold-specific training Harmony for within-training-mouse matching; query matching uses training PCA projection because the archived Harmony fit has no unseen-mouse prediction operator.",
 "Auxiliary training_pca route uses the same frozen training PCA projection for both training and query matching, separating preprocessing leakage from choice of matching representation.",
 "The query uses training-retained genes and fixed cell QC. Query normalized HVGs are centered/scaled with training-only gene means/sample SDs, clipped above 10 as in Seurat ScaleData, and multiplied by training PCA loadings; no query scaling/PCA/Harmony/marker fit occurs.",
 "Matching: same mouse and guide; at least three targets; input-order greedy nearest available of 50 guide-undetected controls, no replacement within comparison.",
 "edgeR paired mouse+group fits require >=3 training mice; filterByExpr min.count=5, robust dispersion/QL; FDR<0.10, |log2FC|>=0.5, exclude mt-/Rpl/Rps prefixes and subtract training-only mSafe eligible genes.",
 "If mSafe has <3 eligible training mice, the calibrated programme is unestimable; no fallback blacklist, threshold relaxation or no-blacklist result replaces it.",
 "All 17 non-mSafe guides are tabulated; original_primary_programme labels the historical 11 solely for reporting. The no-IEG version removes Fos/Fosb/Junb without retuning.",
 "Query AUC/effects are within-mouse descriptive transfer diagnostics; they are not new biological replicates, propagation calibration, or tests of spatial receiver fields.",
 "This is a prespecified revision-time extension of the original joint-Harmony pipeline, not an exact replay of the original joint representation. Query biological independence from the spatial cohort remains unknown.")
writeLines(scope,file.path(out,"FROZEN_PROTOCOL.txt"))
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
fwrite(data.table(file=c(script,file.path(release,"R/full_holdout.R")),sha256=vapply(c(script,file.path(release,"R/full_holdout.R")),digest,character(1),file=TRUE,algo="sha256")),file.path(out,"frozen_code_sha256.tsv"),sep="\t")
match_pairs <- function(md,state,targets) {
 rows <- list()
 for(guide in targets)for(m in sort(unique(md$mouse_id))) {
  tc <- rownames(md)[md$mouse_id==m&md$guide_call==guide];pool <- rownames(md)[md$mouse_id==m&md$guide_call=="Neg"]
  if(length(tc)<3L||!length(pool))next
  k <- min(50L,length(pool));nn <- matrix(RANN::nn2(state[pool,,drop=FALSE],state[tc,,drop=FALSE],k=k)$nn.idx,nrow=length(tc),ncol=k)
  used <- logical(length(pool));chosen <- integer(length(tc))
  for(i in seq_along(tc)) {p <- nn[i,][!used[nn[i,]]][1L];if(is.na(p))break;chosen[i] <- p;used[p] <- TRUE}
  if(any(chosen==0L))next
  cc <- pool[chosen]
  rows[[paste(guide,m)]] <- data.table(guide=guide,mouse_id=m,pair=seq_along(tc),target_cell=tc,control_cell=cc,distance=sqrt(rowSums((state[tc,,drop=FALSE]-state[cc,,drop=FALSE])^2)))
 }
 rbindlist(rows)
}
make_pb <- function(counts,pairs) {
 ans <- setNames(vector("list",length(all_guides)),all_guides)
 for(guide in all_guides) {
  guide_name <- guide;q <- pairs[guide==guide_name];cols <- list();sm <- list()
  for(m in sort(unique(q$mouse_id))) {
   z <- q[mouse_id==m]
   cols[[paste(m,"target",sep="::")]] <- as.integer(Matrix::rowSums(counts[,z$target_cell,drop=FALSE]))
   cols[[paste(m,"matched_neg",sep="::")]] <- as.integer(Matrix::rowSums(counts[,z$control_cell,drop=FALSE]))
   sm[[m]] <- data.table(mouse_id=m,group=c("target","matched_neg"),cells=nrow(z))
  }
  if(length(cols)){mat <- do.call(cbind,cols);rownames(mat) <- rownames(counts);ans[[guide]] <- list(counts=mat,samples=rbindlist(sm))}
 }
 ans
}
folds <- list();checks <- list();qc <- list();timing <- list()
if(!a[3]%in%names(raw))stop("Unknown heldout mouse: ",a[3])
for(heldout in a[3]) {
 t0 <- proc.time()[["elapsed"]];message("START heldout ",heldout," ",Sys.time())
 train <- raw[setdiff(names(raw),heldout)]
 training_mice <- names(train);qraw <- raw[[heldout]]
 genes <- Reduce(intersect,lapply(train,function(z)rownames(z$counts)))
 ct <- do.call(cbind,lapply(train,function(z)z$counts[genes,,drop=FALSE]));md <- do.call(rbind,lapply(train,`[[`,"meta"));rownames(md) <- md$cell_id
 rm(raw,train);gc()
 obj <- CreateSeuratObject(ct,meta.data=md,min.cells=3,min.features=200)
 rm(ct,md);gc()
 obj[["percent.mt"]] <- PercentageFeatureSet(obj,pattern="^mt-")
 obj <- subset(obj,subset=nFeature_RNA>=300&nFeature_RNA<=7500&nCount_RNA>=500&percent.mt<=20)
 obj <- NormalizeData(obj,verbose=FALSE);obj <- FindVariableFeatures(obj,nfeatures=3000,verbose=FALSE)
 obj <- ScaleData(obj,features=VariableFeatures(obj),verbose=FALSE)
 obj <- RunPCA(obj,features=VariableFeatures(obj),npcs=40,seed.use=42,verbose=FALSE)
 obj <- RunHarmony(obj,group.by.vars="mouse_id",reduction.use="pca",dims.use=1:35,verbose=FALSE)
 obj <- FindNeighbors(obj,reduction="harmony",dims=1:30,verbose=FALSE)
 obj <- FindClusters(obj,resolution=0.6,random.seed=17,verbose=FALSE)
 obj <- RunUMAP(obj,reduction="harmony",dims=1:30,seed.use=17,verbose=FALSE)
 md <- obj[[]];md$cell_id <- rownames(md)
 ct <- GetAssayData(obj,layer="counts");genes <- rownames(ct);hvg <- rownames(Loadings(obj,"pca"))
 # Gene/cell blocks retain the complete calculation without large dense copies.
 train_norm <- GetAssayData(obj,layer="data")
 center <- sdev <- setNames(numeric(length(hvg)),hvg)
 for(ii in split(seq_along(hvg),ceiling(seq_along(hvg)/100L))) {
  xx <- as.matrix(train_norm[hvg[ii],,drop=FALSE]);center[ii] <- rowMeans(xx)
  sdev[ii] <- matrixStats::rowSds(xx,center=center[ii])
 }
 train_loadings <- Loadings(obj,"pca")
 train_projection_error <- 0
 for(ii in split(seq_len(ncol(train_norm)),ceiling(seq_len(ncol(train_norm))/1000L))) {
  xx <- as.matrix(train_norm[hvg,ii,drop=FALSE])
  projected <- t(pmin(sweep(sweep(xx,1,center,"-"),1,sdev,"/"),10))%*%train_loadings
  train_projection_error <- max(train_projection_error,max(abs(projected-Embeddings(obj,"pca")[ii,,drop=FALSE])))
 }
 if(!is.finite(train_projection_error)||train_projection_error>1e-8)stop("Training PCA projection mismatch: ",train_projection_error)
 train_states <- list(training_harmony=Embeddings(obj,"harmony")[,1:30,drop=FALSE],training_pca=Embeddings(obj,"pca")[,1:30,drop=FALSE])
 snapshot <- list(training_pca=Embeddings(obj,"pca"),training_harmony=Embeddings(obj,"harmony"),training_umap=Embeddings(obj,"umap"),training_pca_loadings=train_loadings,training_scaling=list(center=center,sd=sdev))
 rm(obj,train_norm,xx,projected);gc()
 qcount <- qraw$counts[genes,,drop=FALSE]
 query <- CreateSeuratObject(qcount,meta.data=qraw$meta,min.cells=0,min.features=200)
 rm(qraw,qcount);gc()
 query[["percent.mt"]] <- PercentageFeatureSet(query,pattern="^mt-")
 query <- subset(query,subset=nFeature_RNA>=300&nFeature_RNA<=7500&nCount_RNA>=500&percent.mt<=20)
 qmd <- query[[]];qmd$cell_id <- rownames(qmd);qcount <- GetAssayData(query,layer="counts")
 qlib <- Matrix::colSums(qcount)
 qstate <- matrix(NA_real_,nrow=ncol(qcount),ncol=ncol(train_loadings),dimnames=list(colnames(qcount),colnames(train_loadings)))
 for(ii in split(seq_len(ncol(qcount)),ceiling(seq_len(ncol(qcount))/1000L))) {
  xx <- log1p(sweep(as.matrix(qcount[hvg,ii,drop=FALSE]),2,pmax(qlib[ii],1),"/")*1e4)
  qstate[ii,] <- t(pmin(sweep(sweep(xx,1,center,"-"),1,sdev,"/"),10))%*%train_loadings
 }
 rm(xx);gc()
 query_pairs <- match_pairs(qmd,qstate[,1:30,drop=FALSE],all_guides)
 routes <- lapply(train_states,function(st){p <- match_pairs(md,st,all_guides);list(pseudobulk=make_pb(ct,p),matching_pairs=p)})
 query_cells <- unique(c(query_pairs$target_cell,query_pairs$control_cell))
 folds[[heldout]] <- list(heldout=heldout,training_mice=training_mice,routes=routes,query_pairs=query_pairs,
  query_counts=qcount[,query_cells,drop=FALSE],query_library=Matrix::colSums(qcount)[query_cells])
 fold_dir <- file.path(out,heldout);dir.create(fold_dir,showWarnings=FALSE)
 fwrite(as.data.table(md),file.path(fold_dir,"training_metadata.tsv.gz"),sep="\t")
 fwrite(as.data.table(qmd),file.path(fold_dir,"query_metadata.tsv.gz"),sep="\t")
 fwrite(data.table(gene=hvg,training_mean=center,training_sd=sdev),file.path(fold_dir,"training_hvg_scaling.tsv"),sep="\t")
 fwrite(query_pairs,file.path(fold_dir,"query_matching_pairs.tsv"),sep="\t")
 for(route in names(routes))fwrite(routes[[route]]$matching_pairs,file.path(fold_dir,paste0(route,"_matching_pairs.tsv")),sep="\t")
 snapshot$query_pca <- qstate;saveRDS(snapshot,file.path(fold_dir,"preprocessing_snapshot.rds"),compress="xz")
 fwrite(data.table(file=c(script,file.path(release,"R/full_holdout.R")),sha256=vapply(c(script,file.path(release,"R/full_holdout.R")),digest,character(1),file=TRUE,algo="sha256")),file.path(fold_dir,"executed_code_sha256.tsv"),sep="\t")
 checks[[heldout]] <- data.table(heldout_mouse=heldout,training_mice=paste(training_mice,collapse=";"),n_training_cells=ncol(ct),n_query_cells=ncol(qcount),n_training_genes=nrow(ct),n_hvg=length(hvg),n_training_clusters=uniqueN(md$seurat_clusters),train_projection_max_abs_error=train_projection_error,query_excluded_from_training=all(md$mouse_id!=heldout),n_training_harmony_pairs=nrow(routes$training_harmony$matching_pairs),n_training_pca_pairs=nrow(routes$training_pca$matching_pairs),n_query_pairs=nrow(query_pairs))
 qc[[paste(heldout,"train")]] <- as.data.table(md)[,.(cells=.N),by=.(mouse_id,guide_call)][,`:=`(heldout_mouse=heldout,partition="training")]
 qc[[paste(heldout,"query")]] <- as.data.table(qmd)[,.(cells=.N),by=.(mouse_id,guide_call)][,`:=`(heldout_mouse=heldout,partition="query")]
 timing[[heldout]] <- data.table(heldout_mouse=heldout,seconds=proc.time()[["elapsed"]]-t0)
 fwrite(rbindlist(checks),file.path(out,"fold_preprocessing_checks.tsv"),sep="\t")
 fwrite(rbindlist(qc),file.path(out,"fold_qc.tsv"),sep="\t")
 fwrite(rbindlist(timing),file.path(out,"stage_timings.tsv"),sep="\t")
 fwrite(checks[[heldout]],file.path(fold_dir,"fold_preprocessing_checks.tsv"),sep="\t")
 fwrite(rbindlist(qc),file.path(fold_dir,"fold_qc.tsv"),sep="\t")
 fwrite(timing[[heldout]],file.path(fold_dir,"stage_timings.tsv"),sep="\t")
 saveRDS(folds[[heldout]],file.path(fold_dir,"compact_fold_counts.rds"),compress="xz")
 rm(query,ct,md,qcount,qmd,qstate,train_states,routes,snapshot);gc()
 message("DONE heldout ",heldout," ",Sys.time())
}
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
