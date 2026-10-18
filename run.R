#!/usr/bin/env Rscript
# One portable entry point. No package installation or network access is performed.
script <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
CALLER_WD <- getwd()
ROOT <- normalizePath(dirname(script)); setwd(ROOT)
a <- commandArgs(trailingOnly=TRUE)
usage <- paste(
 'Rscript run.R MODE [OUTPUT_DIRECTORY] [PARAMETER1] [PARAMETER2]',
 'Modes: stats | reference | clu | perturbmap | shape | diagnostics | dbit | profiles | checks | full-holdout | guide-error | benchmark | injection | dropout | native-export | smoke | verify | all | complete',
 'clu: B=99; benchmark: replicates=100 B=199; injection/dropout: B=199',
 'reference: optional third argument overrides data/reference_counts.rds (relative to caller directory)',
 'all = fast core; complete additionally runs full strict Clu and DBiT relocations. Simulations/native Python are separate.',sep='\n')
if(!length(a)){cat(usage,'\n');quit(status=0)}
mode <- a[1]
out <- if(length(a)>1){v<-a[2];if(!grepl('^(/|[A-Za-z]:)',v))v<-file.path(CALLER_WD,v);normalizePath(v,mustWork=FALSE)}else file.path(ROOT,'results')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages(library(data.table))
read_inputs <- function()readRDS(file.path(ROOT,'data/endpoints.rds'))
read_cells <- function()readRDS(file.path(ROOT,'data/clu_cells.rds'))
param <- function(i,default)if(length(a)>=i)as.integer(a[i])else default
if(mode%in%c('stats','all','complete')) {
 source('R/summary.R');d<-read_inputs();run_summary(d,out)
 if(file.exists('R/clinical.R')){source('R/clinical.R');run_clinical(d,out)}
}
if(mode%in%c('perturbmap','all','complete')) {
 if(!file.exists('R/perturbmap.R'))stop('Perturb-map module is not installed')
 source('R/perturbmap.R');run_perturbmap(file.path(ROOT,'data/perturbmap_spots.tsv.gz'),file.path(out,'perturbmap'))
}
if(mode=='shape'){source('R/shape_controls.R');run_shape_controls(file.path(ROOT,'data/perturbmap_spots.tsv.gz'),out)}
if(mode=='diagnostics') {
 source('R/clu.R');source('R/clu_exposure_diagnostics.R');source('R/score_influence.R')
 run_clu_exposure_diagnostics(read_cells(),fread('data/clu_control_cluster_bank.tsv'),file.path(out,'clu'))
 run_score_influence(file.path(ROOT,'data/perturbmap_spots.tsv.gz'),file.path(out,'perturbmap'))
 source('R/hcc_evaluation.R');run_hcc_diagnostics(ROOT,file.path(out,'clinical'))
}

if(mode%in%c('checks','all','complete')){source('R/summary.R');source('R/computational_checks.R');run_computational_checks(readRDS('data/computational_checks.rds'),file.path(out,'computational_checks'))}
if(mode=='full-holdout'){source('R/full_holdout.R');run_full_holdout(readRDS('data/computational_checks.rds')$full_holdout_clu_counts,out)}
if(mode=='guide-error'){source('R/spf_core.R');source('R/guide_error.R');run_guide_error(out,param(3,100L),param(4,199L))}
if(mode%in%c('profiles','all','complete')){d<-read_inputs();pm<-file.path(out,'perturbmap/perturb_map_response_field_profiles.tsv');if(file.exists(pm))d$perturb_map_profiles<-fread(pm);source('R/atlas.R');d$atlas<-run_atlas(d,out);source('R/profiles.R');run_profiles(d,out)}
if(mode%in%c('reference','all','complete')){
 reference_path <- file.path(ROOT,'data/reference_counts.rds')
 if(mode=='reference'&&length(a)>2L){reference_path<-a[3];if(!grepl('^(/|[A-Za-z]:)',reference_path))reference_path<-file.path(CALLER_WD,reference_path)}
 source('R/reference.R');run_reference(readRDS(reference_path),file.path(out,'reference'))
}
if(mode%in%c('clu','complete')){source('R/clu.R');run_strict_clu(read_cells(),out,param(3,99L))}
if(mode%in%c('dbit','complete')){source('R/dbit.R');run_dbit(readRDS('data/dbit_scores.rds'),file.path(out,'dbit'))}
if(mode%in%c('benchmark','injection','dropout','smoke')) {
 source('R/spf_core.R')
 if(mode=='benchmark'){source('R/benchmark.R');run_benchmark(out,param(3,100L),param(4,199L))}
 if(mode=='injection'){source('R/injection.R');run_injection(read_cells(),out,param(3,199L))}
 if(mode=='dropout'){source('R/injection.R');run_dropout(read_inputs()$real_injection_design,out,param(3,199L))}
 if(mode=='smoke') {
  r <- data.frame(x=c(0,50,500,700),y=0,component=c(1,1,2,2),niche=c(0,0,1,1),identity=c('A','B','A','B'),score=1:4)
  s <- data.frame(x=c(0,500),y=0,component=c(1,2),niche=c(0,1))
  stopifnot(identical(as.numeric(spf_distance(r,s,1L)),c(0,50,Inf,Inf)))
  set.seed(17);spec<-data.frame(component=c(1,2),size=c(1,1),niche=c(0,1))
  for(i in 1:10){ii<-spf_sample_sources(r,spec,'cluster_anchor_niche');stopifnot(all(sort(r$component[ii])==1:2))}
  stopifnot(all(is.na(spf_contrast(rep(Inf,4)))),is.na(spf_pvalue(NA,1:99)['p']))
  cells<-read_cells();d<-read_inputs()
  stopifnot(nrow(cells)==229775L,length(unique(cells$chip_id))==3L,
    sum(cells$guide_call=='sgrna-clu')==508L,all(cells$library_size>0),
    nrow(d$synthetic_replicates)==50400L,
    all(c('mouse 1','mouse 2; mouse 3','mouse 4')%in%as.character(d$animal_crosswalk$animal_or_block)))
  source('R/summary.R');run_summary(d,out)
  message('PASS: component exclusion, source-pool matching, unavailable contrasts, full cell counts, animal provenance, and endpoint summaries.')
 }
}
if(mode=='verify'){source('R/summary.R');verify_results(read_inputs(),out)}
if(mode=='native-export') {
 z<-read_cells()[chip_id=='B03018A2'];f<-file.path(out,'native_input.tsv.gz')
 fwrite(z[,.(chip_id,x_um,y_um,identity,guide_call,raw_clu_score)],f,sep='\t')
 d<-read_inputs()$real_injection_design
 for(nm in c('receiver','source'))fwrite(d[[nm]][,c('x','y','component','identity','raw_clu_score')],file.path(out,paste0('native_sampled_',nm,'.tsv')),sep='\t')
 message('Exported full-chip and frozen sampled native inputs: ',f)
}
if(!mode%in%c('stats','reference','clu','perturbmap','shape','diagnostics','dbit','profiles','checks','full-holdout','guide-error','benchmark','injection','dropout','native-export','smoke','verify','all','complete'))stop(usage)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
