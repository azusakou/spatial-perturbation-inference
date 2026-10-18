# Revision-time same-estimand benchmark with a newly drawn bank per realization.
# See PROTOCOL_FROZEN.md. Core functions are read, never modified here.
suppressPackageStartupMessages({library(data.table);library(mgcv);library(RANN);library(igraph)})

independent_bank_plan <- function(n_rep=100L) {
  RNGkind("L'Ecuyer-CMRG");set.seed(2026101605L)
  next_state <- .Random.seed
  take <- function() {
    z<-next_state;next_state<<-parallel::nextRNGStream(next_state);z
  }
  scenarios<-c('baseline','disconnected','singleton','source_niche',
               'receiver_composition','dropout_50pct')
  setNames(lapply(scenarios,function(s)list(geometry=take(),replicates=lapply(seq_len(n_rep),function(i)
    setNames(lapply(seq_len(4L),function(k)take()),c('source','background','cluster_component','cluster_anchor_niche'))))),scenarios)
}

run_independent_bank <- function(out,n_rep=100L,B=199L,scenarios=NULL,plan=NULL) {
  if(is.null(plan))plan<-independent_bank_plan(100L)
  if(is.null(scenarios))scenarios<-names(plan)
  stopifnot(n_rep<=length(plan[[1]]$replicates),B==199L)
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  saveRDS(plan,file.path(out,'random_stream_plan.rds'))
  use_stream<-function(z)assign('.Random.seed',z,envir=.GlobalEnv)
  all_rows<-list();null_records<-list();k<-0L;geometry_rows<-list();started<-proc.time()[3]
  for(ss in seq_along(scenarios)) {
    scenario<-scenarios[ss];use_stream(plan[[scenario]]$geometry)
    grid<-expand.grid(x=seq(40,2360,length.out=28),y=seq(40,2360,length.out=28))
    receiver<-transform(grid,x=x+runif(nrow(grid),-25,25),y=y+runif(nrow(grid),-25,25))
    candidate<-expand.grid(x=seq(40,2360,length.out=21),y=seq(40,2360,length.out=21))
    receiver$component<-1L;candidate$component<-1L
    if(scenario=='disconnected') {
      receiver<-receiver[receiver$x<1000|receiver$x>1400,]
      candidate<-candidate[candidate$x<1000|candidate$x>1400,]
      receiver$component<-1L+(receiver$x>1400);candidate$component<-1L+(candidate$x>1400)
    }
    receiver$niche<-as.integer(receiver$x>=1200);candidate$niche<-as.integer(candidate$x>=1200)
    p_identity<-if(scenario=='receiver_composition')plogis((receiver$x-1200)/120)else rep(.5,nrow(receiver))
    receiver$identity<-ifelse(runif(nrow(receiver))<p_identity,'A','B')
    receiver$log_lib<-8+.2*(receiver$identity=='A')+.0001*receiver$x+rnorm(nrow(receiver),0,.3)
    group_size<-if(scenario=='singleton')1L else 3L
    groups<-12L/group_size
    spec<-data.frame(component=rep(unique(candidate$component),length.out=groups),size=group_size,
                     niche=if(scenario=='source_niche')rep(1L,groups)else rep(0:1,length.out=groups))
    if(scenario=='disconnected')spec$niche<-spec$component-1L
    cache<-spf_sampler_cache(candidate,spec)
    geometry_rows[[scenario]]<-list(receiver=receiver,candidate=candidate,spec=spec)
    make_bank<-function(method) vapply(seq_len(B),function(b) {
      ids<-spf_sample_sources(candidate,spec,method,cache)
      dd<-spf_distance(receiver,candidate,ids)
      # No overlap exists in this separate-grid model; retain the core rule explicitly.
      rkey<-paste(receiver$component,sprintf('%.17g',receiver$x),sprintf('%.17g',receiver$y))
      skey<-paste(candidate$component[ids],sprintf('%.17g',candidate$x[ids]),sprintf('%.17g',candidate$y[ids]))
      dd[rkey%in%skey]<-Inf
      spf_contrast(dd)
    },numeric(nrow(receiver)))
    for(rr in seq_len(n_rep)) {
      streams<-plan[[scenario]]$replicates[[rr]]
      use_stream(streams$source)
      source<-spf_sample_sources(candidate,spec,'cluster_anchor_niche',cache)
      distance<-spf_distance(receiver,candidate,source)
      contrast<-spf_contrast(distance)
      near_count<-sum(is.finite(distance)&distance<=100)
      far_count<-sum(is.finite(distance)&distance>400&distance<=800)
      kernel<-exp(-distance/150)
      use_stream(streams$background)
      a<-rnorm(4,0,.35)
      background<-sin(receiver$x/390)+.7*cos(receiver$y/430)+
        a[1]*sin((receiver$x+receiver$y)/500)+a[2]*cos((receiver$x-receiver$y)/600)+
        a[3]*sin(receiver$x/650)+a[4]*cos(receiver$y/650)+
        .8*(receiver$identity=='A')+.2*(receiver$log_lib-8)
      base<-background+rnorm(nrow(receiver),0,.8)
      mask<-if(scenario=='dropout_50pct')runif(nrow(receiver))>.5 else rep(TRUE,nrow(receiver))
      banks<-list()
      for(method in c('cluster_component','cluster_anchor_niche')) {
        use_stream(streams[[method]]);banks[[method]]<-make_bank(method)
      }
      for(amplitude in c(0,1)) {
        score<-(base+amplitude*kernel)*mask
        residual<-spf_adjust(score,receiver,'spatial')
        observed<-sum(contrast*residual)
        for(method in names(banks)) {
          null<-as.numeric(crossprod(residual,banks[[method]]));ans<-spf_pvalue(observed,null)
          k<-k+1L;all_rows[[k]]<-data.table(scenario=scenario,replicate=rr,amplitude=amplitude,
            adjustment='spatial',null_model=method,observed=observed,p=unname(ans['p']),
            effect=unname(ans['effect']),valid=as.integer(ans['valid']),attempted=B,
            near_receivers=near_count,far_receivers=far_count,n_receiver=nrow(receiver),
            n_source=length(source),source_indices=paste(source,collapse=','),
            injected_near_far=amplitude*sum(contrast*kernel*mask))
          null_records[[k]]<-null
        }
      }
      if(rr%%10L==0L) {
        fwrite(rbindlist(all_rows),file.path(out,'independent_bank_replicates.tsv'),sep='\t')
        message(scenario,' ',rr,'/',n_rep,'; elapsed ',round(proc.time()[3]-started,1),' s')
      }
    }
  }
  x<-rbindlist(all_rows)
  summary<-x[,{
    good<-is.finite(p);n<-sum(good);rejected<-sum(p[good]<=.05);ci<-spf_wilson(rejected,n)
    list(attempted_replicates=.N,estimable=n,unavailable=.N-n,rejections=rejected,
         rate=if(n)rejected/n else NA_real_,lower=ci[1],upper=ci[2],
         minimum_valid=min(valid),maximum_valid=max(valid),mean_effect=mean(effect[good]))
  },by=.(scenario,amplitude,adjustment,null_model)]
  fwrite(x,file.path(out,'independent_bank_replicates.tsv'),sep='\t')
  fwrite(summary,file.path(out,'independent_bank_summary.tsv'),sep='\t')
  saveRDS(list(keys=x[,.(scenario,replicate,amplitude,null_model)],
               null_statistics=do.call(rbind,null_records)),
          file.path(out,'null_statistics.rds'),compress='xz')
  saveRDS(geometry_rows,file.path(out,'scenario_geometries.rds'))
  writeLines(c(paste('elapsed_seconds',proc.time()[3]-started),capture.output(sessionInfo())),
             file.path(out,'runtime_session.txt'))
  invisible(summary)
}
