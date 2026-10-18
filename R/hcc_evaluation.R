# Fixed-set revision diagnostic; LPO ranks a held-out pair within one fit. No LPO Brier or binomial CI.
suppressPackageStartupMessages(library(data.table))
run_hcc_evaluation <- function(d,out) {
dir.create(out,recursive=TRUE,showWarnings=FALSE)
z<-d$hcc_internal;y<-z$outcome;n<-length(y);grid<-c(.01,.1,1,10,100)
sets<-list('Direct killing only'=c('delta_apoptosis','post_apoptosis_q95'),'Pretreatment state'=paste0(c('interferon','cytotoxic','fibrosis','hypoxia'),'_Pre'),'SPF geometry only'=paste0('field_slope_',c('interferon','cytotoxic','fibrosis','hypoxia')),'Direct + SPF'=c('delta_apoptosis','post_apoptosis_q95',paste0('field_slope_',c('interferon','cytotoxic','fibrosis','hypoxia'))))
fit<-function(X,yy,lam){mu<-colMeans(X);sdv<-apply(X,2,sd);sdv[!is.finite(sdv)|sdv==0]<-1;D<-cbind(1,sweep(sweep(X,2,mu,'-'),2,sdv,'/'));f<-function(b){eta<-drop(D%*%b);sum(log1p(exp(-abs(eta)))+pmax(eta,0)-yy*eta)+lam*sum(b[-1]^2)/2};g<-function(b){v<-drop(crossprod(D,plogis(D%*%b)-yy));v[-1]<-v[-1]+lam*b[-1];v};o<-optim(numeric(ncol(D)),f,g,method='BFGS',control=list(maxit=500,reltol=1e-10));stopifnot(o$convergence==0);list(b=o$par,mu=mu,sd=sdv)}
predict1<-function(f,X)drop(plogis(cbind(1,sweep(sweep(X,2,f$mu,'-'),2,f$sd,'/'))%*%f$b))
loss<-function(yy,p)-sum(yy*log(pmax(p,1e-8))+(1-yy)*log(pmax(1-p,1e-8)))
auc<-function(yy,p)mean(outer(p[yy==1],p[yy==0],function(a,b)(a>b)+.5*(a==b)))
rows<-pairs<-list();maxerr<-0
for(nm in names(sets)){
 X<-as.matrix(z[,sets[[nm]],with=FALSE]);cache<-new.env(parent=emptyenv())
 getfit<-function(ids,lam){key<-paste(paste(ids,collapse=','),lam,sep=':');if(!exists(key,cache,inherits=FALSE))assign(key,fit(X[ids,,drop=FALSE],y[ids],lam),cache);get(key,cache,inherits=FALSE)}
 choose<-function(ids){ll<-vapply(grid,function(la){p<-vapply(ids,function(j)predict1(getfit(setdiff(ids,j),la),X[j,,drop=FALSE]),numeric(1));loss(y[ids],p)},numeric(1));grid[which.min(ll)]}
 for(i in seq_len(n)){
  ids<-setdiff(seq_len(n),i);la<-choose(ids);pr<-predict1(getfit(ids,la),X[i,,drop=FALSE]);old<-d$hcc_internal_predictions_expected[patient_id==z$patient_id[i]&model==nm]$probability;maxerr<-max(maxerr,abs(old-pr));rows[[length(rows)+1L]]<-data.table(model=nm,patient_id=z$patient_id[i],y=y[i],prediction=pr,lambda=la)
 }
 for(i in which(y==1))for(j in which(y==0)){
  ids<-setdiff(seq_len(n),c(i,j));la<-choose(ids);pr<-predict1(getfit(ids,la),X[c(i,j),,drop=FALSE]);pairs[[length(pairs)+1L]]<-data.table(model=nm,responder=z$patient_id[i],nonresponder=z$patient_id[j],train_responders=sum(y[ids]),train_n=length(ids),lambda=la,responder_prediction=pr[1],nonresponder_prediction=pr[2],rank=(pr[1]>pr[2])+.5*(pr[1]==pr[2]))
 }
 cat(nm,'done\n');flush.console()
}
r<-rbindlist(rows);p<-rbindlist(pairs);s<-merge(r[,.(pooled_nested_LOPO_auc=auc(y,prediction),pooled_nested_LOPO_brier=mean((prediction-y)^2)),by=model],p[,.(nested_LPO_auc=mean(rank),outer_pairs=.N),by=model],by='model')
s<-rbind(s,data.table(model='Intercept only',pooled_nested_LOPO_auc=auc(y,(sum(y)-y)/(n-1)),pooled_nested_LOPO_brier=mean(((sum(y)-y)/(n-1)-y)^2),nested_LPO_auc=.5,outer_pairs=sum(y==1)*sum(y==0)))
ex<-d$hcc_external_predictions_expected;external<-ex[,.(auc=auc(outcome,probability),brier=mean((outcome-probability)^2),log_loss=loss(outcome,probability)/.N),by=model]
p0<-mean(y);external<-rbind(external,data.table(model='Training prevalence intercept',auc=.5,brier=mean((p0-d$hcc_external$outcome)^2),log_loss=loss(d$hcc_external$outcome,rep(p0,nrow(d$hcc_external)))/nrow(d$hcc_external)))
features<-c('post_apoptosis_q95',paste0('field_slope_',c('interferon','cytotoxic','fibrosis','hypoxia')))
shift<-rbindlist(lapply(features,function(f){v<-z[[f]];w<-d$hcc_external[[f]];data.table(feature=f,training_min=min(v),training_max=max(v),external_min=min(w),external_max=max(w),external_outside_training_range=sum(w<min(v)|w>max(v)),external_n=length(w))}))
stopifnot(maxerr<1e-10,nrow(p)==length(sets)*sum(y==1)*sum(y==0),all(p$train_responders==sum(y)-1L),all(p$train_n==n-2L))
fwrite(r,file.path(out,'reproduced_LOPO.tsv'),sep='\t');fwrite(p,file.path(out,'nested_LPO_pairs.tsv'),sep='\t');fwrite(s,file.path(out,'CV_comparison.tsv'),sep='\t');fwrite(external,file.path(out,'external_baselines.tsv'),sep='\t');fwrite(shift,file.path(out,'covariate_shift.tsv'),sep='\t');writeLines(c(paste('LOPO maximum error',maxerr),capture.output(sessionInfo())),file.path(out,'verification.txt'));print(s);print(external);print(shift)

invisible(list(CV_comparison=s,pair_records=p,LOPO_records=r,external_baselines=external,covariate_shift=shift))
}
# Portable full diagnostics: retain n11 mixed-organ primary; exclude only the
# documented non-liver post-treatment sample for an anatomy-defined sensitivity.
run_hcc_diagnostics <- function(package_root,out,mode="diagnostics") {
 stopifnot(mode %in% c("primary","diagnostics"))
 package_root<-normalizePath(package_root);d<-readRDS(file.path(package_root,"data/endpoints.rds"))
 ans<-run_hcc_evaluation(d,file.path(out,"primary_n11"))
 if(mode=="diagnostics") {
  stopifnot(nrow(d$hcc_internal)==11L,sum(d$hcc_internal$patient_id==85567)==1L)
  # Anatomy was specified from sample metadata before this sensitivity was run.
  # Source: revision_code/reproduce/results/hcc_visium_sample_audit.tsv.
  provenance<-data.table(patient_id=85567,pre_sample="cytassist_84_pre",pre_organ="Liver",pre_site="Primary",pre_spots=1155L,post_sample="cytassist_84_post",post_organ="Lymph_node",post_site="Metastatis",post_spots=3089L,reason="Post-treatment lymph-node specimen; exclude by anatomy for liver-only sensitivity")
  z<-d;z$hcc_internal<-copy(d$hcc_internal[patient_id!=85567])
  od<-file.path(out,"liver_only_n10");dir.create(od,recursive=TRUE,showWarnings=FALSE)
  fwrite(provenance,file.path(od,"metadata_exclusion.tsv"),sep="\t")
  e<-new.env(parent=globalenv());source(file.path(package_root,"R/clinical.R"),local=e)
  e$run_clinical(z,od)
  z$hcc_internal_predictions_expected<-fread(file.path(od,"hcc_nested_loo_predictions.tsv"))
  z$hcc_external_predictions_expected<-fread(file.path(od,"hcc_external_locked_predictions.tsv"))
  # The comparison uses newly fit n10 predictions, not stale n11 expectations.
  ans$liver_only<-run_hcc_evaluation(z,file.path(od,"evaluation"))
  writeLines(c("Anatomy-defined n10 sensitivity; n11 retained as primary mixed-organ descriptive benchmark.","Exclude patient85567 by post-treatment organ metadata, not outcome or prediction performance.","The clinical feature sets, ridge grid and nested tuning rule are unchanged.","LPO ranks only the two cases evaluated by one fit; dependent pairs are not independent subjects and no LPO Brier is calculated."),file.path(od,"protocol.txt"))
 }
 invisible(ans)
}
if(sys.nframe()==0L) {
 a<-commandArgs(TRUE)
 if(length(a)<2L||length(a)>3L)stop("Usage: Rscript hcc_evaluation.R PACKAGE_ROOT OUTPUT [primary|diagnostics]; or ENDPOINTS_RDS OUTPUT")
 if(dir.exists(a[1]))run_hcc_diagnostics(a[1],a[2],if(length(a)==3L)a[3]else"diagnostics") else {
  stopifnot(length(a)==2L);run_hcc_evaluation(readRDS(a[1]),a[2])
 }
}
