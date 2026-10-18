# Foldwise raw-preprocessing holdout: compact count refits, no joint embedding.
run_full_holdout <- function(input,outdir) {
 suppressPackageStartupMessages({library(data.table);library(edgeR);library(Matrix)})
 dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
 write <- function(x,n)fwrite(x,file.path(outdir,n),sep="\t",na="NA")
 fit <- function(z) {
  nm <- if(is.null(z))0L else uniqueN(z$samples$mouse_id)
  if(nm<3L)return(list(status="insufficient_training_mice",mice=nm,de=NULL))
  sm <- copy(z$samples);sm[,mouse_id:=factor(mouse_id)]
  sm[,group:=factor(group,levels=c("matched_neg","target"))]
  design <- model.matrix(~mouse_id+group,sm);y <- DGEList(z$counts)
  keep <- filterByExpr(y,design=design,min.count=5)
  if(!any(keep))return(list(status="no_expressed_training_genes",mice=nm,de=NULL))
  y <- calcNormFactors(y[keep,,keep.lib.sizes=FALSE])
  y <- estimateDisp(y,design,robust=TRUE)
  q <- glmQLFit(y,design,robust=TRUE)
  d <- as.data.table(topTags(glmQLFTest(q,coef="grouptarget"),n=Inf,sort.by="none")$table,keep.rownames="gene")
  d[,eligible:=FDR<0.10&abs(logFC)>=0.5&!grepl("^mt-|^Rpl|^Rps",gene,ignore.case=TRUE)]
  list(status="ok",mice=nm,de=d)
 }
 metrics <- weights <- de_rows <- scores <- list(); k <- 0L
 for(fold in input$folds)for(route in names(fold$routes)) {
  fits <- lapply(fold$routes[[route]]$pseudobulk,fit)
  safe <- fits[["mSafe_gRNA"]]
  blacklist <- if(safe$status=="ok")safe$de[eligible==TRUE,gene] else character()
  for(guide in names(fits)) {
   a <- fits[[guide]]
   if(!is.null(a$de))de_rows[[paste(fold$heldout,route,guide)]] <- copy(a$de)[,`:=`(heldout_mouse=fold$heldout,matching_route=route,guide=guide)]
   if(guide=="mSafe_gRNA")next
   guide_name <- guide
   mt <- fold$query_pairs[guide==guide_name]
   for(version in c("all_retained_genes","no_immediate_early_genes")) {
    k <- k+1L
    r <- data.table(heldout_mouse=fold$heldout,matching_route=route,
      perturbation=sub("_gRNA$","",guide),programme=version,
      original_primary_programme=sub("_gRNA$","",guide)%in%input$original_primary_programmes,
      n_training_mice=a$mice,n_training_msafe_mice=safe$mice,
      target_training_status=a$status,msafe_training_status=safe$status,
      n_heldout_target=nrow(mt),n_heldout_control=nrow(mt),
      n_candidate_genes=if(a$status=="ok")sum(a$de$eligible) else NA_integer_,
      n_msafe_blacklist=if(safe$status=="ok")length(blacklist) else NA_integer_,
      n_programme_genes=NA_integer_,auc=NA_real_,mean_target_score=NA_real_,mean_control_score=NA_real_,
      mean_score_difference=NA_real_,control_robust_scale=NA_real_,standardized_mean_difference=NA_real_,
      status="unestimable",reason="",genes="")
    if(a$status!="ok"||safe$status!="ok")r[,reason:=paste("training fit:",a$status,safe$status)] else {
     w <- copy(a$de[eligible==TRUE&!gene%in%blacklist])
     if(version=="no_immediate_early_genes")w <- w[!gene%in%c("Fos","Fosb","Junb")]
     r[,n_programme_genes:=nrow(w)]
     if(!nrow(w))r[,reason:="no genes passed training-only thresholds and mSafe subtraction"] else if(!nrow(mt))r[,reason:="heldout target has fewer than three QC cells or matching failed"] else {
      setorder(w,FDR,gene)
      weights[[paste(fold$heldout,route,guide,version)]] <- copy(w)[,`:=`(heldout_mouse=fold$heldout,matching_route=route,guide=guide,programme=version)]
      cells <- c(mt$target_cell,mt$control_cell);n <- nrow(mt)
      lib <- fold$query_library[cells]
      x <- log1p(sweep(as.matrix(fold$query_counts[w$gene,cells,drop=FALSE]),2,pmax(lib,1),"/")*1e4)
      score <- as.numeric(crossprod(w$logFC,x))/sum(abs(w$logFC))
      av <- score[seq_len(n)];bv <- score[n+seq_len(n)]
      sc <- mad(bv,center=median(bv),constant=1.4826)
      if(!is.finite(sc)||sc<1e-8)sc <- sd(bv)
      r[,`:=`(auc=mean(outer(av,bv,function(x,y)(x>y)+0.5*(x==y))),
       mean_target_score=mean(av),mean_control_score=mean(bv),mean_score_difference=mean(av)-mean(bv),
       control_robust_scale=sc,standardized_mean_difference=if(is.finite(sc)&&sc>0)(mean(av)-mean(bv))/sc else NA_real_,
       status="estimable",reason="strict training-only weights; unseen mouse matched in frozen training-PCA projection",genes=paste(w$gene,collapse=";"))]
      scores[[paste(fold$heldout,route,guide,version)]] <- data.table(heldout_mouse=fold$heldout,matching_route=route,guide=guide,programme=version,
       cell_id=cells,group=rep(c("target","matched_neg"),each=n),pair=rep(seq_len(n),2),score=score)
     }
    }
    metrics[[k]] <- r
   }
  }
 }
 m <- rbindlist(metrics,fill=TRUE)
 write(m,"full_holdout_metrics.tsv")
 write(rbindlist(weights,fill=TRUE),"full_holdout_weights.tsv")
 write(rbindlist(de_rows,fill=TRUE),"full_holdout_training_edger.tsv.gz")
 write(rbindlist(scores,fill=TRUE),"full_holdout_scores.tsv.gz")
 summary <- m[,.(n_folds=.N,n_estimable=sum(status=="estimable"),
  min_auc=if(all(is.na(auc)))NA_real_ else min(auc,na.rm=TRUE),
  median_auc=if(all(is.na(auc)))NA_real_ else median(auc,na.rm=TRUE),
  max_auc=if(all(is.na(auc)))NA_real_ else max(auc,na.rm=TRUE)),
  by=.(matching_route,perturbation,programme,original_primary_programme)]
 write(summary,"full_holdout_summary.tsv")
 write(m[,.(n_folds=.N),by=.(matching_route,status,reason)],"full_holdout_status_counts.tsv")
 writeLines(input$scope,file.path(outdir,"FULL_HOLDOUT_SCOPE.txt"))
 invisible(m)
}
