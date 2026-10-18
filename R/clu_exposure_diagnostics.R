# Fixed-score sensitivity and mSafe-placement falsification; descriptive only.
# Requires the existing clu.R helpers; no new source assignment validity is implied.
run_clu_exposure_diagnostics <- function(cells, bank, out) {
  dir.create(out, recursive=TRUE, showWarnings=FALSE)
  result <- refs <- masks <- checks <- list()
  check <- function(name, ok) {
    checks[[length(checks)+1L]] <<- data.table(check=name, pass=isTRUE(ok))
    if (!isTRUE(ok)) stop(name)
  }
  for (chip in unique(as.character(cells$chip_id))) {
    z <- cells[as.character(chip_id)==chip]
    xy <- as.matrix(z[,.(x_um,y_um)]); cc <- z$component
    neg <- which(z$guide_call=='Neg')
    bg <- data.table(score=z$raw_clu_score[neg], identity=factor(z$identity[neg]),
      log_lib=log1p(z$library_size[neg]),x_um=xy[neg,1],y_um=xy[neg,2])
    fit <- mgcv::bam(score~identity+s(log_lib,k=5)+s(x_um,y_um,k=30),
      data=bg,method='fREML',discrete=TRUE,nthreads=2)
    scores <- cbind(raw=bg$score, frozen_GAM=residuals(fit))
    other <- which((startsWith(as.character(z$guide_call),'sgrna-') & z$guide_call!='sgrna-clu') | z$guide_call=='Multiple')
    od <- nearest_component(xy,other,cc)$distance[neg]
    keep_isolated <- od>100
    source_sets <- list(Clu=which(z$guide_call=='sgrna-clu'), mSafe=which(z$guide_call=='sgrna-msafe'))
    for (target in names(source_sets)) {
      src <- source_sets[[target]]
      d <- nearest_component(xy,src,cc)$distance[neg]
      if (target=='Clu') {
        old <- bank[chip_id==chip][order(draw)]
        check(paste(chip,'99 fixed Clu bank draws'),identical(old$draw,1:99))
        draws <- lapply(old$pseudo_indices,function(s)as.integer(strsplit(s,',',fixed=TRUE)[[1]]))
      } else {
        spec <- source_cluster_spec(xy,src,cc)
        set.seed(20261017L+sum(utf8ToInt(chip)))
        draws <- lapply(1:99,function(b)sample_cluster_preserving_sources(xy,neg,cc,spec))
      }
      groups <- if(target=='Clu') c('all_Neg','exclude_detected_other100') else 'mSafe_placement'
      null <- array(NA_real_,c(99,2,length(groups)))
      obs <- matrix(NA_real_,2,length(groups)); support <- list()
      contrast <- function(ds, use) {
        near <- use & ds<=100; far <- use & ds>400 & ds<=800
        v <- if(sum(near)>=30 && sum(far)>=30) colMeans(scores[near,,drop=FALSE])-colMeans(scores[far,,drop=FALSE]) else c(NA_real_,NA_real_)
        list(value=v,n=sum(near),f=sum(far))
      }
      for (g in seq_along(groups)) {
        use <- if(groups[g]=='exclude_detected_other100') keep_isolated else rep(TRUE,length(neg))
        support[[g]] <- contrast(d,use);obs[,g] <- support[[g]]$value
      }
      for (b in 1:99) {
        ps <- draws[[b]]
        check(paste(chip,target,b,'source membership'),length(ps)==length(src) && !anyDuplicated(ps) && all(ps%in%neg))
        check(paste(chip,target,b,'component counts'),identical(as.integer(table(factor(cc[ps],levels=sort(unique(cc))))),as.integer(table(factor(cc[src],levels=sort(unique(cc)))))))
        pd <- nearest_component(xy,ps,cc)$distance[neg]
        for (g in seq_along(groups)) {
          use <- !(neg%in%ps) & if(groups[g]=='exclude_detected_other100')keep_isolated else TRUE
          check(paste(chip,target,b,groups[g],'source receiver exclusion'),!any(neg[use]%in%ps))
          q <- contrast(pd,use);null[b,,g] <- q$value
          masks[[length(masks)+1L]] <- data.table(chip_id=chip,target=target,population=groups[g],draw=b,near_receivers=q$n,far_receivers=q$f,pseudo_source_receiver_overlap=0L)
        }
      }
      for (g in seq_along(groups)) for (j in 1:2) {
        v <- null[,j,g];v <- v[is.finite(v)];o <- obs[j,g]
        centre <- mean(c(o,v));p <- if(is.finite(o)&&length(v))sum(abs(c(o,v)-centre)>=abs(o-centre))/(1+length(v)) else NA_real_
        result[[length(result)+1L]] <- data.table(chip_id=chip,target=target,population=groups[g],score=colnames(scores)[j],source_cells=length(src),base_receivers=length(neg),retained_receivers=if(groups[g]=='exclude_detected_other100')sum(keep_isolated) else length(neg),near_receivers=support[[g]]$n,far_receivers=support[[g]]$f,observed=o,reference_mean=if(length(v))mean(v) else NA_real_,reference_excess=if(length(v))o-mean(v) else NA_real_,P_ref=p,finite_draws=length(v),attempted_draws=99L)
        refs[[length(refs)+1L]] <- data.table(chip_id=chip,target=target,population=groups[g],score=colnames(scores)[j],draw=1:99,statistic=null[,j,g])
      }
    }
  }
  result <- rbindlist(result); refs <- rbindlist(refs); masks <- rbindlist(masks)
  primary <- result[target=='Clu' & population=='all_Neg' & score=='frozen_GAM']
  expected <- c(C02943C3=.78,B03018A2=.91,A03599E2=.99)
  check('unchanged three primary fractions',all(primary$P_ref==expected[primary$chip_id]))
  check('all18 prespecified endpoints retained',nrow(result)==18L)
  check('all1782 reference statistics retained',nrow(refs)==1782L)
  for (name in c('result','refs','masks')) fwrite(get(name),file.path(out,paste0('clu_exposure_',name,'.tsv')),sep='\t',na='NA')
  fwrite(rbindlist(checks),file.path(out,'clu_exposure_checks.tsv'),sep='\t')
  invisible(result)
}
