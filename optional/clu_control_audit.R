#!/usr/bin/env Rscript
# Real-map ablations; the inference is conditional on scores and source banks.
suppressPackageStartupMessages({library(data.table);library(RANN);library(mgcv)})
a<-commandArgs(TRUE)
if(length(a)<2L||length(a)>3L)stop('Usage: Rscript optional/clu_control_audit.R MINIMAL_PACKAGE_ROOT OUT [--full-only]')
package_root<-normalizePath(a[1]);out<-a[2];full_only<-length(a)==3L&&identical(a[3],'--full-only')
if(length(a)==3L&&!full_only)stop('Only optional flag is --full-only')
input<-normalizePath(file.path(package_root,'data','clu_cells.rds'))
bank_input<-normalizePath(file.path(package_root,'data','clu_control_cluster_bank.tsv'))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
zall<-readRDS(input); stopifnot(nrow(zall)==229775L)
banklog<-fread(bank_input)
check<-list();add<-function(name,ok,value=NA_real_) {
 check[[length(check)+1L]]<<-data.table(check=name,pass=isTRUE(ok),value=value)
 if(!isTRUE(ok))stop('Check failed: ',name)
}
close_num<-function(x,y,tol=1e-10)all((is.infinite(x)&is.infinite(y)&sign(x)==sign(y))|(!is.na(x)&!is.na(y)&abs(x-y)<=tol))
nearest<-function(coords,src,cc) {
 d<-rep(Inf,nrow(coords))
 for(k in intersect(unique(cc),unique(cc[src]))) {
  q<-which(cc==k);s<-src[cc[src]==k]
  d[q]<-RANN::nn2(coords[s,,drop=FALSE],coords[q,,drop=FALSE],k=1)$nn.dists[,1]
 };d
}
tail_fraction<-function(obs,null) {
 null<-null[is.finite(null)];if(!is.finite(obs)||!length(null))return(NA_real_)
 centre<-mean(c(obs,null));(1+sum(abs(null-centre)>=abs(obs-centre)))/(length(null)+1)
}
median_with_inf<-function(d)if(length(d))median(d)else NA_real_
median_finite<-function(d)if(any(is.finite(d)))median(d[is.finite(d)])else NA_real_
summaries<-references<-mixed<-source_mixed<-banks<-point_audits<-list()
backgrounds<-if(full_only)'full_bam' else c('raw','identity','identity_library_smooth','global_linear','full_bam')
references_used<-if(full_only)'cluster' else c('count','cluster')
chips<-unique(as.character(zall$chip_id))
for(chip in chips) {
 cat(format(Sys.time()),chip,'start\n');flush.console()
 z<-copy(zall[as.character(chip_id)==chip]);coords<-as.matrix(z[,.(x_um,y_um)])
 src<-which(z$guide_call=='sgrna-clu');neg<-which(z$guide_call=='Neg');cc<-z$component
 bg<-data.table(score=z$raw_clu_score[neg],identity=factor(as.character(z$identity[neg])),
  log_lib=log1p(z$library_size[neg]),x_um=coords[neg,1],y_um=coords[neg,2])
 score<-if(full_only)matrix(residuals(bam(score~identity+s(log_lib,k=5)+s(x_um,y_um,k=30),
  data=bg,method='fREML',discrete=TRUE,nthreads=2)),ncol=1L,dimnames=list(NULL,'full_bam')) else
  cbind(raw=bg$score,identity=residuals(lm(score~identity,data=bg)),
  identity_library_smooth=residuals(bam(score~identity+s(log_lib,k=5),data=bg,
   method='fREML',discrete=TRUE,nthreads=2)),
  global_linear=residuals(lm(score~identity+log_lib+x_um+y_um,data=bg)),
  full_bam=residuals(bam(score~identity+s(log_lib,k=5)+s(x_um,y_um,k=30),data=bg,
   method='fREML',discrete=TRUE,nthreads=2)))
 stopifnot(identical(colnames(score),backgrounds),all(is.finite(score)))
 ds<-nearest(coords,src,cc);obsnear<-which(ds[neg]<=100);obsfar<-which(ds[neg]>400&ds[neg]<=800)
 add(paste(chip,'observed_support'),length(obsnear)>=30&&length(obsfar)>=30)
 observed<-vapply(seq_len(ncol(score)),function(j)mean(score[obsnear,j])-mean(score[obsfar,j]),numeric(1))
 names(observed)<-backgrounds
 obs_weight<-numeric(length(neg));obs_weight[obsnear]<-1/length(obsnear);obs_weight[obsfar]<- -1/length(obsfar)
 add(paste(chip,'observed_direct_sum_check'),close_num(observed,drop(crossprod(obs_weight,score)),1e-12),
  max(abs(observed-drop(crossprod(obs_weight,score)))))
 old<-banklog[chip_id==chip][order(draw)]
 add(paste(chip,'bank_99_draws'),identical(old$draw,seq_len(99L)))
 clusterbank<-lapply(old$pseudo_indices,function(s)as.integer(strsplit(s,',',fixed=TRUE)[[1]]))
 set.seed(20261016L+sum(utf8ToInt(chip)))
 source_cc<-unique(cc[src]);countbank<-lapply(seq_len(99L),function(b)unlist(lapply(source_cc,function(k) {
  pool<-neg[cc[neg]==k];pool[sample.int(length(pool),sum(cc[src]==k),replace=FALSE)]
 }),use.names=FALSE))
 for(kind in references_used) {
  bank<-if(kind=='cluster')clusterbank else countbank
  null<-matrix(NA_real_,99L,length(backgrounds),dimnames=list(NULL,backgrounds))
  mask<-vector('list',99L)
  for(b in seq_len(99L)) {
   pseudo<-bank[[b]]
   add(paste(chip,kind,b,'source_count'),length(pseudo)==length(src))
   add(paste(chip,kind,b,'unique_eligible'),!anyDuplicated(pseudo)&&all(pseudo%in%neg))
   add(paste(chip,kind,b,'component_counts'),identical(as.integer(table(factor(cc[pseudo],levels=sort(unique(cc))))),
    as.integer(table(factor(cc[src],levels=sort(unique(cc)))))))
   pd<-nearest(coords,pseudo,cc);use<-!neg%in%pseudo
   ni<-which(use&pd[neg]<=100);fi<-which(use&pd[neg]>400&pd[neg]<=800)
   add(paste(chip,kind,b,'no_source_receivers'),!any(neg[c(ni,fi)]%in%pseudo))
   if(length(ni)>=30&&length(fi)>=30) {
    null[b,]<-vapply(seq_along(backgrounds),function(j)mean(score[ni,j])-mean(score[fi,j]),numeric(1))
    independently<-colSums(score[ni,,drop=FALSE])/length(ni)-colSums(score[fi,,drop=FALSE])/length(fi)
    add(paste(chip,kind,b,'null_direct_sum_check'),close_num(null[b,],independently,1e-12),max(abs(null[b,]-independently)))
   }
   if(kind=='cluster')add(paste(chip,kind,b,'original_bank_row_and_mask'),
    old$base_receivers[b]==length(neg)&&old$excluded_pseudo_sources[b]==length(pseudo)&&
    old$near_receivers[b]==length(ni)&&old$far_receivers[b]==length(fi)&&
    old$effective_receivers[b]==sum(use)&&old$pseudo_source_spots[b]==length(pseudo))
   mask[[b]]<-data.table(chip_id=chip,reference=kind,draw=b,pseudo_sources=length(pseudo),
    near_receivers=length(ni),far_receivers=length(fi),receiver_source_overlap=0L,
    pseudo_indices=paste(sort(pseudo),collapse=','))
  }
  banktable<-rbindlist(mask);banks[[length(banks)+1L]]<-banktable
  for(j in seq_along(backgrounds)) {
   v<-null[,j];v<-v[is.finite(v)];centre<-mean(c(observed[j],v));p<-tail_fraction(observed[j],v)
   # Equivalent independent rank against a common pooled centre; no expected input used.
   dev<-abs(c(observed[j],v)-centre)
   rankp<-sum(dev>=dev[1])/length(dev)
   add(paste(chip,kind,backgrounds[j],'pooled_rank_check'),identical(as.numeric(p),as.numeric(rankp)))
   summaries[[length(summaries)+1L]]<-data.table(chip_id=chip,reference=kind,background=backgrounds[j],
    source_cells=length(src),receiver_cells=length(neg),near_receivers=length(obsnear),far_receivers=length(obsfar),
    observed=unname(observed[j]),null_mean=mean(v),null_sd=sd(v),excess=unname(observed[j])-mean(v),
    common_centre=centre,tail_fraction=p,valid_nulls=length(v),attempted_nulls=99L)
  }
  references[[length(references)+1L]]<-data.table(chip_id=chip,reference=kind,
   draw=rep(seq_len(99L),length(backgrounds)),background=rep(backgrounds,each=99L),statistic=as.vector(null))
  cat(format(Sys.time()),chip,kind,'99 complete\n');flush.console()
 }
 other<-which(startsWith(as.character(z$guide_call),'sgrna-')&z$guide_call!='sgrna-clu')
 multi<-which(z$guide_call=='Multiple');do<-nearest(coords,other,cc);dm<-nearest(coords,multi,cc);de<-pmin(do,dm)
 for(w in c('near','far')) {
  wi<-if(w=='near')neg[obsnear]else neg[obsfar]
  mixed[[length(mixed)+1L]]<-data.table(chip_id=chip,window=w,receiver_cells=length(wi),
   other_single_sources=length(other),multiple_sources=length(multi),
   other_single_within100=sum(do[wi]<=100),multiple_within100=sum(dm[wi]<=100),
   either_within100=sum(de[wi]<=100),other_single_closer_than_clu=sum(do[wi]<ds[wi]),
   either_closer_than_clu=sum(de[wi]<ds[wi]),
   fraction_other_single_within100=mean(do[wi]<=100),fraction_multiple_within100=mean(dm[wi]<=100),
   fraction_either_within100=mean(de[wi]<=100),fraction_other_single_closer_than_clu=mean(do[wi]<ds[wi]),
   fraction_either_closer_than_clu=mean(de[wi]<ds[wi]),
   median_clu_distance=median_with_inf(ds[wi]),median_other_single_distance=median_with_inf(do[wi]),
   median_multiple_distance=median_with_inf(dm[wi]),median_either_distance=median_with_inf(de[wi]),
   finite_median_other_single_distance=median_finite(do[wi]),finite_median_multiple_distance=median_finite(dm[wi]),
   no_same_component_other_single=sum(!is.finite(do[wi])),no_same_component_multiple=sum(!is.finite(dm[wi])))
 }
 src_other_clu<-rep(Inf,length(src))
 for(i in seq_along(src)) {
  candidates<-setdiff(src[cc[src]==cc[src[i]]],src[i])
  if(length(candidates))src_other_clu[i]<-sqrt(min(rowSums(sweep(coords[candidates,,drop=FALSE],2,coords[src[i],])^2)))
 }
 source_mixed[[length(source_mixed)+1L]]<-data.table(chip_id=chip,clu_sources=length(src),
  other_single_sources=length(other),multiple_sources=length(multi),
  other_single_within100=sum(do[src]<=100),multiple_within100=sum(dm[src]<=100),either_within100=sum(de[src]<=100),
  other_clu_within100=sum(src_other_clu<=100),fraction_other_single_within100=mean(do[src]<=100),
  fraction_multiple_within100=mean(dm[src]<=100),fraction_either_within100=mean(de[src]<=100),
  fraction_other_clu_within100=mean(src_other_clu<=100),median_other_single_distance=median_with_inf(do[src]),
  median_multiple_distance=median_with_inf(dm[src]),median_other_clu_distance=median_with_inf(src_other_clu))
 # Brute-force check independently samples 100 queries from each full chip.
 set.seed(20261017L+sum(utf8ToInt(chip)));qtest<-sample.int(nrow(z),100L)
 for(label in c('clu','other_single','multiple')) {
  ss<-switch(label,clu=src,other_single=other,multiple=multi)
  actual<-switch(label,clu=ds,other_single=do,multiple=dm)
  brute<-vapply(qtest,function(q){s<-ss[cc[ss]==cc[q]];if(!length(s))Inf else
   sqrt(min((coords[s,1]-coords[q,1])^2+(coords[s,2]-coords[q,2])^2))},numeric(1))
  add(paste(chip,label,'100_bruteforce_distances'),close_num(actual[qtest],brute,1e-10),
   max(abs(actual[qtest][is.finite(brute)]-brute[is.finite(brute)]),0))
 }
 point_audits[[chip]]<-list(scores=score,neg_indices=neg,clu_distance=ds[neg],
  other_single_distance=do[neg],multiple_distance=dm[neg],guide_calls=as.character(z$guide_call))
 cat(format(Sys.time()),chip,'mixed complete\n');flush.console()
}
res<-rbindlist(summaries);ref<-rbindlist(references);mix<-rbindlist(mixed);smix<-rbindlist(source_mixed);bank<-rbindlist(banks)
expected_endpoints<-length(chips)*length(references_used)*length(backgrounds)
add('all_factorial_endpoints',nrow(res)==expected_endpoints)
add('all_finite_reference_statistics',nrow(ref)==expected_endpoints*99L&&all(is.finite(ref$statistic)))
add('six_receiver_windows',nrow(mix)==6L)
add('three_source_coverage_rows',nrow(smix)==3L)
# Post-fit recovery check: these three documented values are never used in a fit.
old_p<-c(C02943C3=.78,B03018A2=.91,A03599E2=.99)
old_excess<-c(C02943C3=.00781146545416495,B03018A2=-.00244112946817131,A03599E2=.00106592445967114)
cur<-res[reference=='cluster'&background=='full_bam']
for(i in seq_len(nrow(cur))) {
 chip<-cur$chip_id[i]
 add(paste(chip,'documented_tail_recovery'),cur$tail_fraction[i]==old_p[chip])
 # A tolerance of 1e-9 covers only printed archival precision, not fit choices.
 add(paste(chip,'documented_effect_recovery'),abs(cur$excess[i]-old_excess[chip])<1e-9,
  abs(cur$excess[i]-old_excess[chip]))
}
fwrite(res,file.path(out,'clu_control_factorial.tsv'),sep='\t')
fwrite(ref,file.path(out,'clu_control_reference.tsv'),sep='\t')
fwrite(mix,file.path(out,'clu_mixed_receiver_exposure.tsv'),sep='\t')
fwrite(smix,file.path(out,'clu_mixed_source_coverage.tsv'),sep='\t')
fwrite(bank,file.path(out,'clu_source_banks.tsv'),sep='\t')
fwrite(rbindlist(check),file.path(out,'machine_checks.tsv'),sep='\t')
saveRDS(list(summary=res,reference=ref,mixed_receiver=mix,mixed_source=smix,source_banks=bank,
 protocol='PROTOCOL_FROZEN.md',seed=20261016L,backgrounds=backgrounds),file.path(out,'clu_control_audit.rds'),compress='xz')
saveRDS(point_audits,file.path(out,'independent_point_audit.rds'),compress='xz')
writeLines(c(capture.output(sessionInfo()),'Input MD5:',paste(names(tools::md5sum(c(input,bank_input))),
 tools::md5sum(c(input,bank_input)))),file.path(out,'session_and_input_md5.txt'))
cat(format(Sys.time()),'all complete; all',length(check),'checks passed\n');flush.console()
