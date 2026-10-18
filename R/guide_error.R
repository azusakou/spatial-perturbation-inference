# Repeated count-level guide detection with explicit receiver reassignment.
# Effects are on gene-count means; this is a model stress test, not assay validation.
run_guide_error <- function(out, n_rep=100L, B=199L) {
  dir.create(out, recursive=TRUE, showWarnings=FALSE)
  cfg <- data.frame(case=c('perfect','drop25','drop50','drop75','drop50_fp01','drop50_niche'),
    dropout=c(0,.25,.5,.75,.5,.5), false_positive=c(0,0,0,0,.01,0),
    niche=c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE))
  fwrite(cfg,file.path(out,'guide_error_design.tsv'),sep='\t')
  records <- list(); detection <- list(); native <- list(); nr <- nd <- nn <- 0L
  make_bank <- function(receiver, source, B) {
    ans <- list()
    if(nrow(source)<3L) return(setNames(rep(list(matrix(NA_real_,nrow(receiver),B)),6L),
      as.vector(t(outer(c('count_component','cluster_component','cluster_anchor_niche'),
                       c('retain_pseudo','exclude_pseudo'),paste,sep=':')))))
    spec <- spf_cluster_spec(source,100); cache <- spf_sampler_cache(receiver,spec)
    for (method in c('count_component','cluster_component','cluster_anchor_niche')) {
      keep <- drop <- matrix(NA_real_,nrow(receiver),B)
      for (b in seq_len(B)) {
        ids <- spf_sample_sources(receiver,spec,method,cache)
        d <- spf_distance(receiver,receiver,ids)
        keep[,b] <- spf_contrast(d)
        d[ids] <- Inf # matched exclusion of a pseudo-source from its receivers
        drop[,b] <- spf_contrast(d)
      }
      ans[[paste(method,'retain_pseudo',sep=':')]] <- keep
      ans[[paste(method,'exclude_pseudo',sep=':')]] <- drop
    }
    ans
  }
  for (rr in seq_len(n_rep)) {
    set.seed(20261016L+rr*1000L)
    grid <- expand.grid(x=seq(40,1400,length.out=18),y=seq(40,1400,length.out=18))
    cell <- transform(grid,x=x+runif(nrow(grid),-5,5),y=y+runif(nrow(grid),-5,5))
    cell$component <- 1L; cell$niche <- as.integer(cell$x>=720)
    cell$identity <- ifelse(runif(nrow(cell))<plogis((cell$x-720)/350),'A','B')
    cell$log_lib <- 8.5+.15*(cell$identity=='A')+rnorm(nrow(cell),0,.25)
    lib <- round(exp(cell$log_lib)); u_count <- matrix(runif(nrow(cell)*4),nrow(cell),4)
    background <- sin(cell$x/390)+.7*cos(cell$y/430)+.8*(cell$identity=='A')+
      .2*(cell$log_lib-8.5)+rnorm(1,0,.2)*sin((cell$x+cell$y)/500)
    for (ss in seq_len(nrow(cfg))) {
      par <- cfg[ss,]; set.seed(20261016L+rr*1000L+ss*100L)
      spec <- data.frame(component=1L,size=4L,niche=if(par$niche)rep(1L,6) else rep(0:1,3))
      truth <- spf_sample_sources(cell,spec,'cluster_anchor_niche')
      is_source <- seq_len(nrow(cell)) %in% truth
      dose <- exp(.2*(cell$log_lib-mean(cell$log_lib)))
      lambda <- if(par$dropout==0)rep(Inf,nrow(cell)) else -log(par$dropout)*dose
      guide <- integer(nrow(cell))
      guide[is_source] <- if(par$dropout==0)1L else rpois(sum(is_source),lambda[is_source])
      guide[!is_source] <- rpois(sum(!is_source),-log(1-par$false_positive)*dose[!is_source])
      observed <- which(guide>0)
      d_true <- spf_distance(cell,cell,truth); kernel <- exp(-d_true/150)
      nd <- nd+1L
      detection[[nd]] <- data.frame(case=par$case,replicate=rr,true_sources=sum(is_source),
        detected_true=sum(guide[is_source]>0),false_negative=sum(guide[is_source]==0),
        false_positive=sum(guide[!is_source]>0),observed_sources=length(observed),
        receivers=sum(guide==0),contaminating_sources=sum(is_source & guide==0),
        sensitivity=mean(guide[is_source]>0),specificity=mean(guide[!is_source]==0),
        receiver_contamination=sum(is_source & guide==0)/sum(guide==0))
      banks <- list(); views <- list()
      for (view in c('oracle','observed')) {
        src <- if(view=='oracle')truth else observed
        recv <- setdiff(seq_len(nrow(cell)),src)
        receiver <- cell[recv,]; source <- cell[src,]
        views[[view]] <- list(src=src,recv=recv,receiver=receiver,
          v=spf_contrast(spf_distance(receiver,source,seq_len(nrow(source)))))
        banks[[view]] <- make_bank(receiver,source,B)
      }
      # The extra direct=0 null isolates contamination from a cell-autonomous effect.
      effects <- rbind(data.frame(direct=2,amplitude=c(0,.5,1)),data.frame(direct=0,amplitude=0))
      for (ee in seq_len(nrow(effects))) {
        direct <- effects$direct[ee]; amplitude <- effects$amplitude[ee]
        mu <- lib*.0007*exp(.4*background+.35*direct*is_source+.35*amplitude*kernel)
        counts <- vapply(seq_len(4),function(g)qnbinom(u_count[,g],mu=mu,size=5),numeric(nrow(cell)))
        score <- rowMeans(log1p(1e4*counts/lib))
        if(rr==1L && ss==1L && ee==1L) saveRDS(list(cell=cell,library=lib,
          counts=counts,score=score,true_source=is_source,guide=guide,mean_count=mu,
          background=background,d_true=d_true),file.path(out,'count_model_example.rds'))
        for (view in names(views)) {
          z <- views[[view]]; y <- score[z$recv]
          for (adjust in c('none','spatial')) {
            r <- spf_adjust(y,z$receiver,adjust); obs <- sum(z$v*r)
            for (method in names(banks[[view]])) {
              nul <- as.numeric(crossprod(r,banks[[view]][[method]]))
              zz <- spf_pvalue(obs,nul); nul <- nul[is.finite(nul)]
              legacy <- if(is.finite(obs)&&length(nul)>=20L)
                (1+sum(abs(nul-mean(nul))>=abs(obs-mean(nul))))/(length(nul)+1) else NA_real_
              nr <- nr+1L; records[[nr]] <- data.frame(case=par$case,replicate=rr,
                direct_effect=direct,field_amplitude=amplitude,receiver_definition=view,
                adjustment=adjust,null_model=method,p=zz['p'],p_null_only=legacy,
                effect=zz['effect'],valid=zz['valid'],n_source=length(z$src),n_receiver=length(z$recv),
                contaminated_receivers=sum(is_source[z$recv]),
                latent_field_contrast=amplitude*sum(z$v*kernel[z$recv]))
          }
        }
        }
        # Native comparison uses fixed, prespecified subsets of the same generated truth.
        if(par$case %in% c('perfect','drop50','drop50_niche') && direct==2 && amplitude %in% c(0,1)) {
          nn <- nn+1L
          native[[nn]] <- data.table(case=par$case,replicate=rr,amplitude=amplitude,
            node=seq_len(nrow(cell)),x=cell$x,y=cell$y,component=cell$component,
            observed_source=guide>0,true_source=is_source,guide_count=guide,score=score)
        }
      }
    }
    if(rr%%10L==0L) {
      fwrite(rbindlist(records),file.path(out,'guide_error_replicates.tsv.gz'),sep='\t')
      message('Count-level guide error: completed ',rr,'/',n_rep)
    }
  }
  x <- rbindlist(records); det <- rbindlist(detection)
  fwrite(x,file.path(out,'guide_error_replicates.tsv.gz'),sep='\t')
  fwrite(det,file.path(out,'guide_detection_replicates.tsv'),sep='\t')
  summary <- x[, {
    ok<-is.finite(p); n<-sum(ok); k<-sum(p[ok]<=.05); ci<-spf_wilson(k,n)
    list(replicates=.N,estimable=n,rejections=k,rate=if(n)k/n else NA_real_,
      lower=ci[1],upper=ci[2],mean_effect=mean(effect[ok]),
      mean_contamination=mean(contaminated_receivers/n_receiver),
      legacy_rejections=sum(p_null_only<=.05,na.rm=TRUE))
  },by=.(case,direct_effect,field_amplitude,receiver_definition,adjustment,null_model)]
  fwrite(summary,file.path(out,'guide_error_summary.tsv'),sep='\t')
  fwrite(rbindlist(native),file.path(out,'native_known_truth_input.tsv.gz'),sep='\t')
  fwrite(det[,lapply(.SD,mean),by=case,.SDcols=c('sensitivity','specificity','receiver_contamination',
    'true_sources','detected_true','false_negative','false_positive','observed_sources')],
    file.path(out,'guide_detection_summary.tsv'),sep='\t')
  saveRDS(list(design=cfg,replicates=x,detection=det,summary=summary,
    seed=20261016L,replications=n_rep,null_configurations=B),file.path(out,'guide_error_results.rds'),compress='xz')
  invisible(summary)
}
