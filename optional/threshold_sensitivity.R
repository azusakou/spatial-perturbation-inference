#!/usr/bin/env Rscript
# Actual near/far-window sensitivity from public MI/skin matrices. No installs.
# Keeps the original source q=.95, whole-gene libraries and assay score formula.
a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=3L||!a[1]%in%c('mi','skin'))stop('Usage: Rscript optional/threshold_sensitivity.R mi|skin RAW_DIRECTORY OUTPUT_DIRECTORY')
suppressPackageStartupMessages({library(data.table);library(Matrix);library(RANN)})
root<-normalizePath(a[2]);out<-normalizePath(a[3],mustWork=FALSE);dir.create(out,recursive=TRUE,showWarnings=FALSE)
SCALES<-c(.75,1,1.25)

window_endpoints<-function(z,source,programs,metadata){
 receiver<-setdiff(seq_len(nrow(z)),source)
 d<-nn2(as.matrix(z[source,.(x_um,y_um)]),as.matrix(z[receiver,.(x_um,y_um)]),k=1)$nn.dists[,1]
 effects<-profiles<-list();j<-0L
 for(w in SCALES)for(pr in programs){
  j<-j+1L;y<-z[[pr]][receiver];near<-d<=200*w;far<-d>400*w&d<=800*w
  estimable<-any(near)&any(far);value<-if(estimable)mean(y[near])-mean(y[far])else NA_real_
  effects[[j]]<-cbind(metadata,data.table(program=pr,distance_scale=w,source_spots=length(source),receiver_spots=length(receiver),
   near_upper=200*w,far_lower=400*w,far_upper=800*w,near_receivers=sum(near),far_receivers=sum(far),
   near_far_effect=value,source_score=mean(z[[pr]][source]),status=if(estimable)'estimable'else'unestimable_empty_near_or_far'))
  cuts<-c(-Inf,c(100,200,400,800)*w,Inf)
  labs<-c(paste0('0-',100*w),paste0(100*w,'-',200*w),paste0(200*w,'-',400*w),paste0(400*w,'-',800*w),paste0('>',800*w))
  pp<-data.table(distance_bin=cut(d,cuts,labels=labs),score=y)[,.(mean_score=mean(score),receivers=.N),by=distance_bin]
  profiles[[j]]<-cbind(metadata,data.table(program=pr,distance_scale=w),pp)
 }
 list(effect=rbindlist(effects),profile=rbindlist(profiles))
}

label_tests<-function(units,expected_n,expected_controls){
 rbindlist(lapply(unique(units$distance_scale),function(w)rbindlist(lapply(unique(units$program),function(pr){
  z<-units[distance_scale==w&program==pr];n<-nrow(z);nc<-sum(z$is_control)
  stopifnot(n==expected_n,nc==expected_controls)
  ok<-all(is.finite(z$near_far_effect));obs<-p<-NA_real_
  if(ok){obs<-mean(z$near_far_effect[!z$is_control])-mean(z$near_far_effect[z$is_control]);
   nul<-combn(n,nc,function(ii)mean(z$near_far_effect[-ii])-mean(z$near_far_effect[ii]));p<-mean(abs(nul)>=abs(obs)-1e-12)}
  data.table(program=pr,distance_scale=w,units=n,controls=nc,estimable_units=sum(is.finite(z$near_far_effect)),
   observed_injury_minus_control=obs,exact_assignments=choose(n,nc),exact_p_two_sided=p,
   status=if(ok)'estimable'else'unestimable_missing_unit_endpoint',
   interpretation='conditional outcome-label reference; non-significance does not establish no field')
 }))))
}

run_mi<-function(){
 suppressPackageStartupMessages(library(rhdf5))
 programs<-list(fibrotic=c('COL1A1','COL1A2','COL3A1','FN1','POSTN','TIMP1','TNC'),
  inflammatory=c('LYZ','S100A8','S100A9','IL1B','CCL2','TYROBP'),vascular=c('PECAM1','VWF','KDR','EMCN','ENG','RAMP2'),
  regenerative=c('MKI67','TOP2A','MYC','IGF1','NRG1','AREG'))
 files<-sort(list.files(root,'^Visium_.*[.]h5ad$',full.names=TRUE));stopifnot(length(files)==27L)
 result<-cache<-audit<-vector('list',length(files))
 for(k in seq_along(files)){
  path<-files[k];sid<-sub('[.]h5ad$','',basename(path));message('MI ',k,'/',length(files),' ',sid)
  patient<-sub('.*_(P[0-9]+)(?:_rep[0-9]+)?$','\\1',sid,perl=TRUE)
  zone<-sub('^Visium_','',sub('_P[0-9]+(?:_rep[0-9]+)?$','',sid,perl=TRUE))
  primary_zone<-if(zone=='control')'Control'else if(grepl('IZ',zone))'Infarct'else if(grepl('GT',zone))'Granulation'else if(grepl('FZ',zone))'Fibrotic'else if(grepl('BZ',zone))'Border'else'Remote'
  xy<-h5read(path,'/obsm/X_spatial');if(nrow(xy)==2L)xy<-t(xy);n<-nrow(xy)
  cm<-as.numeric(h5read(path,'/obs/Cardiomyocyte'));my<-as.numeric(h5read(path,'/obs/Myeloid'))
  si<-as.numeric(scale(-cm))+as.numeric(scale(my));source<-which(si>=quantile(si,.95,na.rm=TRUE));receiver<-setdiff(seq_len(n),source)
  nn<-nn2(xy,xy,k=2)$nn.dists[,2];grid<-median(nn[nn>0]);xy<-xy/grid*100
  genes<-h5read(path,'/raw/var/_index');ids<-which(toupper(genes)%in%toupper(unique(unlist(programs))))
  ip<-as.numeric(h5read(path,'/raw/X/indptr'));stopifnot(length(ip)==n+1L)
  sel<-matrix(0,length(ids),n,dimnames=list(toupper(genes[ids]),NULL));lib<-numeric(n)
  # CSR streaming holds only 500 spots of full counts, preserving the whole-
  # transcriptome denominator and exact selected count matrix.
  for(start in seq.int(1L,n,by=500L)){
   end<-min(n,start+499L);begin<-ip[start];last<-ip[end+1L]
   if(last>begin){
    dd<-as.numeric(h5read(path,'/raw/X/data',start=begin+1,count=last-begin));gi<-as.integer(h5read(path,'/raw/X/indices',start=begin+1,count=last-begin))+1L
    stopifnot(all(dd>=0),all(abs(dd-round(dd))<1e-7))
    cum<-c(0,cumsum(dd));lib[start:end]<-diff(cum[ip[start:(end+1L)]-begin+1L]);keep<-which(gi%in%ids)
    if(length(keep))sel[,start:end]<-as.matrix(sparseMatrix(i=match(gi[keep],ids),j=findInterval(keep-1L,ip[start:(end+1L)]-begin,all.inside=TRUE),x=dd[keep],dims=c(length(ids),end-start+1L)))
   }
  }
  supplied<-as.numeric(h5read(path,'/obs/n_counts'));err<-max(abs(lib-supplied));stopifnot(err<1e-7)
  lognorm<-log1p(sweep(sel,2,pmax(lib,1),'/')*1e4)
  score<-sapply(programs,function(gs){v<-colMeans(lognorm[intersect(toupper(gs),rownames(lognorm)),,drop=FALSE]);sc<-mad(v[receiver],constant=1.4826);if(!is.finite(sc)||sc==0)sc<-sd(v[receiver]);(v-median(v[receiver]))/sc})
  z<-cbind(data.table(x_um=xy[,1],y_um=xy[,2]),as.data.table(score));metadata<-data.table(sample_id=sid,patient,zone,primary_zone)
  result[[k]]<-window_endpoints(z,source,names(programs),metadata)
  cache[[k]]<-list(metadata=metadata,source=source,score_positions=z)
  audit[[k]]<-data.table(sample_id=sid,raw_genes=length(genes),selected_genes=length(ids),spots=n,whole_library_max_error=err,coordinate_median_nn=grid)
 }
 ef<-rbindlist(lapply(result,`[[`,'effect'));pf<-rbindlist(lapply(result,`[[`,'profile'))
 units<-ef[,.(near_far_effect=if(all(is.finite(near_far_effect)))mean(near_far_effect)else NA_real_,is_control=all(primary_zone=='Control'),sections=.N),by=.(patient,program,distance_scale)]
 tests<-label_tests(units,20L,4L)
 fwrite(ef,file.path(out,'mi_window_section_endpoints.tsv'),sep='\t');fwrite(pf,file.path(out,'mi_window_profiles.tsv'),sep='\t');fwrite(units,file.path(out,'mi_window_patient_endpoints.tsv'),sep='\t');fwrite(tests,file.path(out,'mi_window_tests.tsv'),sep='\t');fwrite(rbindlist(audit),file.path(out,'mi_normalization_audit.tsv'),sep='\t')
 saveRDS(cache,file.path(out,'mi_score_position_cache.rds'),compress='xz');invisible(tests)
}

run_skin<-function(){
 suppressPackageStartupMessages(library(jsonlite))
 programs<-list(inflammation_source=c('S100a8','S100a9','Lcn2','Il1b','Cxcl2','Ccl2'),fibrosis=c('Col1a1','Col1a2','Col3a1','Fn1','Ctgf','Timp1','Postn'),regeneration=c('Krt14','Krt17','Krt6a','Krt6b','Mki67','Top2a'),myeloid=c('Ptprc','Lyz2','Tyrobp','Cd68','Csf1r','Ccl2'))
 files<-sort(Sys.glob(file.path(root,'GSM*_matrix.mtx.gz')));stopifnot(length(files)==7L)
 good<-vapply(files,function(f)system2('gzip',c('-t',shQuote(f)),stdout=FALSE,stderr=FALSE)==0L,logical(1))
 fwrite(data.table(file=basename(files),gzip_valid=good),file.path(out,'skin_matrix_integrity.tsv'),sep='\t')
 files<-files[good];stopifnot(length(files)==5L)
 result<-cache<-audit<-vector('list',length(files))
 for(k in seq_along(files)){
  path<-files[k];stem<-sub('_matrix[.]mtx[.]gz$','',basename(path));sid<-sub('^GSM[0-9]+_','',stem);message('Skin ',k,'/',length(files),' ',sid)
  ff<-function(s)file.path(root,paste0(stem,s));counts<-readMM(gzfile(path));features<-fread(ff('_features.tsv.gz'),header=FALSE);barcodes<-fread(ff('_barcodes.tsv.gz'),header=FALSE)[[1]]
  rownames(counts)<-make.unique(features[[2]]);colnames(counts)<-barcodes
  posfile<-ff('_tissue_positions_list.csv.gz');if(!file.exists(posfile))posfile<-ff('_tissue_positions.csv.gz')
  pos<-fread(posfile,header=FALSE,col.names=c('barcode','in_tissue','array_row','array_col','px_row','px_col'))[in_tissue==1]
  common<-intersect(pos$barcode,colnames(counts));pos<-pos[match(common,barcode)];counts<-counts[,common,drop=FALSE]
  sf<-fromJSON(paste(readLines(gzfile(ff('_scalefactors_json.json.gz')),warn=FALSE),collapse=''))
  z<-data.table(x_um=pos$px_col*55/sf$spot_diameter_fullres,y_um=pos$px_row*55/sf$spot_diameter_fullres);lib<-pmax(colSums(counts),1)
  for(pr in names(programs)){
   gs<-intersect(programs[[pr]],rownames(counts));stopifnot(length(gs)>=3L)
   # Skin's frozen score is log of the normalized gene-set mean. MI uses the
   # mean of gene-wise log values. These assay definitions are not interchanged.
   v<-log1p(1e4*as.numeric(colMeans(counts[gs,,drop=FALSE]%*%Diagonal(x=1/lib))));sc<-mad(v,constant=1.4826);if(!is.finite(sc)||sc==0)sc<-sd(v)
   z[[pr]]<-(v-median(v))/sc
  }
  tt<-if(grepl('^UW',sid))'Unwounded'else if(grepl('^D12',sid))'Day 12'else if(grepl('^D21',sid))'Day 21'else'Day 28'
  source<-which(z$inflammation_source>=quantile(z$inflammation_source,.95,na.rm=TRUE));metadata<-data.table(sample_id=sid,time=tt)
  result[[k]]<-window_endpoints(z,source,c('fibrosis','regeneration','myeloid'),metadata);cache[[k]]<-list(metadata=metadata,source=source,score_positions=z)
  audit[[k]]<-data.table(sample_id=sid,genes=nrow(counts),spots=nrow(z),whole_library_min=min(lib),whole_library_max=max(lib),spot_diameter_fullres=sf$spot_diameter_fullres)
 }
 ef<-rbindlist(lapply(result,`[[`,'effect'));pf<-rbindlist(lapply(result,`[[`,'profile'));units<-copy(ef);units[,is_control:=time=='Unwounded']
 tests<-label_tests(units,5L,1L)
 fwrite(ef,file.path(out,'skin_window_section_endpoints.tsv'),sep='\t');fwrite(pf,file.path(out,'skin_window_profiles.tsv'),sep='\t');fwrite(tests,file.path(out,'skin_window_tests.tsv'),sep='\t');fwrite(rbindlist(audit),file.path(out,'skin_normalization_audit.tsv'),sep='\t')
 saveRDS(cache,file.path(out,'skin_score_position_cache.rds'),compress='xz');invisible(tests)
}
if(a[1]=='mi')run_mi()else run_skin()
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
