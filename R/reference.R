# Count-level reference fits after matching from the fully reconstructed joint five-mouse Harmony model.
run_reference <- function(input,outdir) {
suppressPackageStartupMessages({library(data.table);library(edgeR);library(Matrix)})
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
out <- function(x,name)fwrite(x,file.path(outdir,name),sep="\t",na="NA")
pb <- input$pseudobulk;clu_cells <- input$heldout_cells;pair <- input$matching_pairs
fit_pb <- function(target,heldout=NULL) {
 z <- pb[[target]]
 use <- if(is.null(heldout))rep(TRUE,nrow(z$samples)) else z$samples$mouse_id!=heldout
 sm <- copy(z$samples[use]);nm <- uniqueN(sm$mouse_id)
 if(nm<3L)return(list(status="insufficient_training_mice",de=NULL,mice=nm))
 sm[,mouse_id:=factor(mouse_id)];sm[,group:=factor(group,levels=c("matched_neg","target"))]
 design <- model.matrix(~mouse_id+group,data=sm)
 y <- DGEList(z$counts[,use,drop=FALSE])
 keep <- filterByExpr(y,design=design,min.count=5)
 if(!any(keep))return(list(status="no_expressed_training_genes",de=NULL,mice=nm))
 y <- calcNormFactors(y[keep,,keep.lib.sizes=FALSE])
 y <- estimateDisp(y,design,robust=TRUE)
 fit <- glmQLFit(y,design,robust=TRUE)
 de <- as.data.table(topTags(glmQLFTest(fit,coef="grouptarget"),n=Inf,sort.by="none")$table,keep.rownames="gene")
 de[,program_member:=FDR<0.10&abs(logFC)>=0.5&!grepl("^mt-|^Rpl|^Rps",gene,ignore.case=TRUE)]
 list(status="ok",de=de,mice=nm)
}
de <- rbindlist(lapply(names(pb),function(target) {
 fit <- fit_pb(target)
 if(fit$status!="ok")return(NULL)
 fit$de[,`:=`(perturbation=sub("_gRNA$","",target),paired_mice=fit$mice)]
 fit$de
}))
blacklist <- de[perturbation=="mSafe"&program_member==TRUE,gene]
de[,passes_msafe_calibration:=!gene%in%blacklist]
de[,primary_program_member:=program_member&passes_msafe_calibration&perturbation!="mSafe"]
de[,abs_logFC:=abs(logFC)]
setorder(de,perturbation,FDR,-abs_logFC)
out(de,"nearest_state_edger.tsv.gz")
summary <- de[primary_program_member==TRUE,.(n_genes=.N,n_up=sum(logFC>0),n_down=sum(logFC<0),
 paired_mice=unique(paired_mice),top_genes=paste(head(gene,20),collapse=";")),by=perturbation]
out(summary,"nearest_state_programme_summary.tsv")
out(data.table(gene=blacklist),"mSafe_blacklist.tsv")
out(de[primary_program_member==TRUE,.(perturbation,gene,weight=logFC,
 no_immediate_early_member=!gene%in%c("Fos","Fosb","Junb"))],"programme_weights.tsv")

mice <- sort(unique(clu_cells$cells$mouse_id))
auc <- function(a,b)mean(outer(a,b,function(x,y)(x>y)+0.5*(x==y)))
rows <- list();weights <- list();training <- list();i <- 0L
for(m in mice) {
 a <- fit_pb("Clu_gRNA",m);safe <- fit_pb("mSafe_gRNA",m)
 if(!is.null(a$de))training[[paste(m,"clu")]] <- copy(a$de)[,`:=`(heldout_mouse=m,guide="Clu_gRNA")]
 if(!is.null(safe$de))training[[paste(m,"msafe")]] <- copy(safe$de)[,`:=`(heldout_mouse=m,guide="mSafe_gRNA")]
 match <- pair[guide=="Clu_gRNA"&mouse_id==m]
 blacklist_m <- if(safe$status=="ok")safe$de[program_member==TRUE,gene] else character()
 for(version in c("all_retained_genes","no_immediate_early_genes")) {
  i <- i+1L
  r <- data.table(heldout_mouse=m,programme=version,n_training_mice=a$mice,
    n_training_msafe_mice=safe$mice,clu_training_status=a$status,msafe_training_status=safe$status,
    n_heldout_target=nrow(match),n_heldout_control=nrow(match),
    n_candidate_clu_genes=if(a$status=="ok")sum(a$de$program_member) else NA_integer_,
    n_training_msafe_blacklist=if(safe$status=="ok")length(blacklist_m) else NA_integer_,
    n_programme_genes=NA_integer_,auc=NA_real_,mean_target_score=NA_real_,mean_control_score=NA_real_,
    mean_score_difference=NA_real_,control_robust_scale=NA_real_,standardized_mean_difference=NA_real_,
    matched_pair_mean_difference=NA_real_,status="unestimable",reason="",genes="")
  if(a$status!="ok"||safe$status!="ok") {
   r[,reason:=paste("training fit:",a$status,safe$status)]
  } else {
   w <- copy(a$de[program_member==TRUE&!gene%in%blacklist_m])
   if(version=="no_immediate_early_genes")w <- w[!gene%in%c("Fos","Fosb","Junb")]
   r[,n_programme_genes:=nrow(w)]
   if(!nrow(w))r[,reason:="no genes passed training-only thresholds and blacklist"] else {
    setorder(w,FDR,gene)
    weights[[paste(m,version)]] <- copy(w)[,`:=`(heldout_mouse=m,programme=version)]
    cells <- c(match$target_cell,match$control_cell)
    n <- nrow(match);lib <- setNames(clu_cells$cells$whole_gene_library,clu_cells$cells$cell_id)[cells]
    x <- as.matrix(clu_cells$counts[w$gene,cells,drop=FALSE])
    x <- log1p(sweep(x,2L,pmax(lib,1),"/")*1e4)
    score <- as.numeric(crossprod(w$logFC,x))/sum(abs(w$logFC))
    av <- score[seq_len(n)];bv <- score[n+seq_len(n)]
    sc <- mad(bv,center=median(bv),constant=1.4826)
    if(!is.finite(sc)||sc<1e-8)sc <- sd(bv)
    r[,`:=`(auc=auc(av,bv),mean_target_score=mean(av),mean_control_score=mean(bv),
     mean_score_difference=mean(av)-mean(bv),control_robust_scale=sc,
     standardized_mean_difference=if(is.finite(sc)&&sc>0)(mean(av)-mean(bv))/sc else NA_real_,
     matched_pair_mean_difference=mean(av-bv),status="estimable",
     reason="descriptive weight-learning holdout only; frozen full-reference embedding reused",
     genes=paste(w$gene,collapse=";"))]
   }
  }
  rows[[i]] <- r
 }
}
out(rbindlist(rows),"heldout_Clu_metrics.tsv")
out(rbindlist(weights),"heldout_Clu_weights.tsv")
out(rbindlist(training),"heldout_training_edger.tsv.gz")
writeLines(c(
 "This compact module starts after joint five-mouse Harmony matching. optional/reference.R preprocess/export reconstructs that boundary from public raw matrices; the complete raw replay was executed.",
 "R/14 nearest-state programmes feed R/16 spatial projection; R/13 cluster-matched programmes are a separate historical analysis.",
 "edgeR::filterByExpr min.count=5, TMM calcNormFactors, estimateDisp robust=TRUE, glmQLFit robust=TRUE, grouptarget QL F-test.",
 "Design: ~mouse_id+group; at least three eligible mice, at least three target cells per mouse-guide match.",
 "Programme filters: FDR<0.10, |log2FC|>=0.5, exclude ^mt-|^Rpl|^Rps (ignore.case=TRUE); subtract mSafe eligible genes; exclude mSafe as perturbation.",
 "The no-immediate-early variant additionally removes Fos/Fosb/Junb.",
 "Heldout fitting uses other animals only for Clu weights and mSafe blacklist; the frozen matching/embedding is retained.",
 "Whole-gene library lengths are used in log1p(10000*count/library) cell scoring.",
 "AUC and score differences are descriptive within heldout mice, not propagation calibration or independent biological-replicate tests.",
 paste0("R ",getRversion(),"; edgeR ",packageVersion("edgeR"),"; limma ",packageVersion("limma"),
  "; statmod ",packageVersion("statmod"),"; Matrix ",packageVersion("Matrix"))),
 file.path(outdir,"MODEL_AND_SCOPE.txt"))
print(summary)
print(rbindlist(rows)[,.(heldout_mouse,programme,n_programme_genes,auc,status)])

invisible(summary)
}
