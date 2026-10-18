# Refit compact Clu holdouts and reaggregate complete new benchmark records.
# All-guide holdouts / HCC / shape geometries below are labelled frozen outputs.
run_computational_checks <- function(d,out) {
 dir.create(out,recursive=TRUE,showWarnings=FALSE)
 checks <- list()
 check <- function(label,x,y,tol=1e-10) {
  x<-as.numeric(x);y<-as.numeric(y)
  ok<-length(x)==length(y)&&identical(is.na(x),is.na(y))
  err<-if(ok&&any(is.finite(x)&is.finite(y)))max(abs(x[is.finite(x)&is.finite(y)]-y[is.finite(x)&is.finite(y)]))else 0
  ok<-ok&&all(is.finite(x)==is.finite(y))&&err<=tol
  checks[[length(checks)+1L]]<<-data.table(check=label,max_absolute_difference=err,pass=ok)
 }
 compare <- function(label,a,b,keys,cols) {
  a<-copy(as.data.table(a));b<-copy(as.data.table(b));setorderv(a,keys);setorderv(b,keys)
  stopifnot(nrow(a)==nrow(b));for(k in keys)stopifnot(identical(as.character(a[[k]]),as.character(b[[k]])))
  for(n in cols)check(paste(label,n),a[[n]],b[[n]])
 }
 source('R/full_holdout.R');run_full_holdout(d$full_holdout_clu_counts,file.path(out,'full_holdout'))
 fm<-fread(file.path(out,'full_holdout/full_holdout_metrics.tsv'))
 compare('Clu count refit',fm,d$full_holdout_all_guides$metrics[perturbation=='Clu'],
  c('heldout_mouse','matching_route','programme'),c('n_programme_genes','auc','mean_score_difference','standardized_mean_difference'))
 for(n in c('metrics','weights','programme_summary','fold_preprocessing'))fwrite(d$full_holdout_all_guides[[n]],file.path(out,paste0('all_guide_frozen_',n,'.tsv')),sep='\t')
 g<-d$guide_error$replicates
 gs<-g[,{
  ok<-is.finite(p);n<-sum(ok);k<-sum(p[ok]<=.05);ci<-wilson(k,n)
  list(replicates=.N,estimable=n,rejections=k,rate=k/n,lower=ci[1],upper=ci[2],mean_effect=mean(effect[ok]),
   mean_contamination=mean(contaminated_receivers/n_receiver),legacy_rejections=sum(p_null_only<=.05,na.rm=TRUE))
 },by=.(case,direct_effect,field_amplitude,receiver_definition,adjustment,null_model)]
 compare('guide records',gs,d$guide_error$summary,c('case','direct_effect','field_amplitude','receiver_definition','adjustment','null_model'),
  c('replicates','estimable','rejections','rate','lower','upper','mean_effect','mean_contamination','legacy_rejections'))
 fwrite(gs,file.path(out,'guide_error_summary.tsv'),sep='\t')
 det<-d$guide_error$detection
 ds<-det[,lapply(.SD,mean),by=case,.SDcols=c('sensitivity','specificity','receiver_contamination','true_sources','detected_true','false_negative','false_positive','observed_sources')]
 check('guide FN accounting',det$false_negative,det$true_sources-det$detected_true)
 check('guide receiver contamination',det$receiver_contamination,det$contaminating_sources/det$receivers)
 fwrite(ds,file.path(out,'guide_detection_summary.tsv'),sep='\t')
 n<-d$native$native_known_truth_records
 ns<-n[,{ok<-is.finite(source_high_z);nn<-sum(ok);k<-sum(abs(source_high_z[ok])>1.96);ci<-wilson(k,nn)
  list(attempted=.N,estimable=nn,rejected=k,gaussian_z_rule_rate=k/nn,wilson_lower=ci[1],wilson_upper=ci[2],
       positive_rejections=sum(source_high_z>1.96,na.rm=TRUE),negative_rejections=sum(source_high_z< -1.96,na.rm=TRUE))},by=.(case,amplitude)]
 compare('native z rule',ns,d$native$native_known_truth_summary,c('case','amplitude'),c('attempted','estimable','rejected','gaussian_z_rule_rate','wilson_lower','wilson_upper'))
 fwrite(ns,file.path(out,'native_known_truth_summary.tsv'),sep='\t')
 for(kind in c('mi','skin')) {
  z<-copy(d$windows[[kind]]);if(kind=='skin')z[,is_control:=time=='Unwounded']
  tests<-z[,{
   obs<-mean(near_far_effect[!is_control])-mean(near_far_effect[is_control])
   nul<-combn(.N,sum(is_control),function(ii)mean(near_far_effect[-ii])-mean(near_far_effect[ii]))
   list(observed_injury_minus_control=obs,exact_assignments=length(nul),exact_p_two_sided=mean(abs(nul)>=abs(obs)-1e-12))
  },by=.(program,distance_scale)]
  compare(paste(kind,'window tests'),tests,d$windows[[paste0(kind,'_expected')]],c('program','distance_scale'),
   c('observed_injury_minus_control','exact_assignments','exact_p_two_sided'))
  fwrite(tests,file.path(out,paste0(kind,'window_tests.tsv')),sep='\t')
 }
 for(nm in names(d$hcc))fwrite(d$hcc[[nm]],file.path(out,paste0('hcc_frozen_',nm,'.tsv')),sep='\t')
 for(nm in names(d$shape))fwrite(d$shape[[nm]],file.path(out,paste0('shape_frozen_',nm,'.tsv')),sep='\t')
 fwrite(d$windows$clu,file.path(out,'clu_frozen_window_tests.tsv'),sep='\t')
 for(nm in names(d$injection))fwrite(d$injection[[nm]],file.path(out,paste0('injection_frozen_',nm,'.tsv')),sep='\t')
 kq<-d$known_assignment$replicates
 ks<-kq[,{
  good<-is.finite(p);nn<-sum(good);kk<-sum(p[good]<=.05);ci<-wilson(kk,nn)
  list(attempted_replicates=.N,estimable=nn,unavailable=.N-nn,rejections=kk,
       rate=if(nn)kk/nn else NA_real_,lower=ci[1],upper=ci[2],minimum_valid=min(valid),
       maximum_valid=max(valid),mean_effect=mean(effect[good]))
 },by=.(scenario,amplitude,adjustment,null_model)]
 compare('known assignment records',ks,d$known_assignment$summary,
  c('scenario','amplitude','adjustment','null_model'),
  c('attempted_replicates','estimable','unavailable','rejections','rate','lower','upper',
    'minimum_valid','maximum_valid','mean_effect'))
 fwrite(ks,file.path(out,'known_assignment_summary.tsv'),sep='\t')
 fwrite(kq,file.path(out,'known_assignment_frozen_records.tsv'),sep='\t')
 writeLines(d$known_assignment$protocol,file.path(out,'known_assignment_frozen_protocol.md'))
 writeLines(d$known_assignment$proof,file.path(out,'known_assignment_conditional_argument.md'))
 cr<-copy(d$clu_control$reference);cs<-copy(d$clu_control$summary)
 ca<-cr[,{
  key<-cs[chip_id==.BY$chip_id&reference==.BY$reference&background==.BY$background]
  stopifnot(nrow(key)==1L)
  val<-statistic[is.finite(statistic)];obs<-key$observed;centre<-mean(c(obs,val))
  list(null_mean=mean(val),null_sd=sd(val),excess=obs-mean(val),common_centre=centre,
       tail_fraction=(1+sum(abs(val-centre)>=abs(obs-centre)))/(length(val)+1),valid_nulls=length(val))
 },by=.(chip_id,reference,background)]
 compare('Clu frozen reference reaggregation',ca,cs,c('chip_id','reference','background'),
  c('null_mean','null_sd','excess','common_centre','tail_fraction','valid_nulls'))
 fwrite(ca,file.path(out,'clu_control_reference_reaggregation.tsv'),sep='\t')
 for(nm in c('summary','reference','mixed_receiver','mixed_source','source_banks'))
  fwrite(d$clu_control[[nm]],file.path(out,paste0('clu_control_frozen_',nm,'.tsv')),sep='\t')
 writeLines(d$clu_control$protocol_text,file.path(out,'clu_control_frozen_protocol.md'))
 writeLines(c('full_holdout/: fresh edgeR/score refit from compact Clu+mSafe counts.',
  'guide_error/native/window/known_assignment summaries: reaggregated individual generated records or frozen unit endpoints.',
  'clu_control_reference_reaggregation: recalculates rank fractions/effects from frozen 2970 null statistics and observed contrasts; does not refit cell-level models.',
  '*_frozen_*: exports completed-generator snapshots or protocols; not raw regeneration by this mode.',
  'Real source-reference fractions and BH-adjusted fractions are descriptive diagnostics; no general type-I/FDR guarantee is supplied.',
  'optional/independent_bank.R generates new banks/realizations under a specified known synthetic Q; it does not validate real source assignment.',
  'optional/clu_control_audit.R refits 30 real-map control combinations; --full-only refits the three full-GAM cluster endpoints.',
  'Rscript run.R benchmark/guide-error/injection/dropout/shape and optional raw/native interfaces execute the corresponding generators.'),
  file.path(out,'EXECUTION_SCOPE.txt'))
 v<-rbindlist(checks);fwrite(v,file.path(out,'verification.tsv'),sep='\t')
 if(!all(v$pass))stop('New computational-check mismatch')
 message('PASS: ',nrow(v),' new refit/aggregation checks; frozen output scope is explicit.')
 invisible(v)
}
