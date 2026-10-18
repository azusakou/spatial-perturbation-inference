# Synthetic known-truth generator; expensive full run is optional.
run_benchmark <- function(out, n_rep=100L, B=199L) {
seed <- 20261015L; set.seed(seed)
scenarios <- data.frame(
  scenario = c("baseline", "source_6", "source_24", "singleton", "cluster_6",
               "boundary", "disconnected", "radius_75", "radius_300", "radius_600",
               "receiver_composition", "source_niche", "annotation_25pct", "dropout_50pct"),
  n_source = c(12,6,24,12,12,12,12,12,12,12,12,12,12,12),
  cluster_size = c(3,3,3,1,6,3,3,3,3,3,3,3,3,3),
  radius = c(150,150,150,150,150,150,150,75,300,600,150,150,150,150),
  stringsAsFactors = FALSE)
fwrite(scenarios, file.path(out,"synthetic_scenarios.tsv"), sep = "\t")
all_rows <- list(); nr <- 0L
for (ss in seq_len(nrow(scenarios))) {
  par <- scenarios[ss, ]; set.seed(seed + ss * 10000L)
  grid <- expand.grid(x = seq(40,2360,length.out=28), y = seq(40,2360,length.out=28))
  receiver <- transform(grid, x=x+runif(nrow(grid),-25,25), y=y+runif(nrow(grid),-25,25))
  candidate <- expand.grid(x=seq(40,2360,length.out=21),y=seq(40,2360,length.out=21))
  receiver$component <- 1L; candidate$component <- 1L
  if (par$scenario == "disconnected") {
    receiver <- receiver[receiver$x<1000 | receiver$x>1400, ]
    candidate <- candidate[candidate$x<1000 | candidate$x>1400, ]
    receiver$component <- 1L+(receiver$x>1400); candidate$component <- 1L+(candidate$x>1400)
  }
  receiver$niche <- as.integer(receiver$x>=1200); candidate$niche <- as.integer(candidate$x>=1200)
  p_id <- if(par$scenario=="receiver_composition") plogis((receiver$x-1200)/120) else rep(.5,nrow(receiver))
  receiver$identity <- ifelse(runif(nrow(receiver))<p_id,"A","B")
  receiver$log_lib <- 8+.2*(receiver$identity=="A")+.0001*receiver$x+rnorm(nrow(receiver),0,.3)
  if (par$scenario=="annotation_25pct") {
    err <- runif(nrow(receiver))<.25
    receiver$identity_true <- receiver$identity
    receiver$identity[err] <- ifelse(receiver$identity[err]=="A","B","A")
  } else receiver$identity_true <- receiver$identity
  nc <- par$n_source / par$cluster_size
  comp <- rep(unique(candidate$component),length.out=nc)
  spec <- data.frame(component=comp,size=par$cluster_size,
                     niche=if(par$scenario=="source_niche") rep(1L,nc) else rep(0:1,length.out=nc))
  # In disconnected tissue the component determines the coarse x niche.
  if(par$scenario=="disconnected") spec$niche <- spec$component-1L
  cache <- spf_sampler_cache(candidate,spec)
  banks <- spf_null_bank(receiver,candidate,spec,B)
  valid_bank <- vapply(banks,function(m)sum(colSums(is.finite(m))==nrow(m)),integer(1))
  message(sprintf("Synthetic %s: n=%d, valid null configurations %s",par$scenario,nrow(receiver),paste(valid_bank,collapse=",")))
  for (rr in seq_len(n_rep)) {
    source <- spf_sample_sources(candidate,spec,"cluster_anchor_niche",cache)
    d <- spf_distance(receiver,candidate,source); contrast <- spf_contrast(d)
    # A new low-frequency background is drawn independently of source placement.
    a <- rnorm(4,0,.35)
    bg <- sin(receiver$x/390)+.7*cos(receiver$y/430)+
      a[1]*sin((receiver$x+receiver$y)/500)+a[2]*cos((receiver$x-receiver$y)/600)+
      a[3]*sin(receiver$x/650)+a[4]*cos(receiver$y/650)+
      .8*(receiver$identity_true=="A")+.2*(receiver$log_lib-8)
    base_score <- bg+rnorm(nrow(receiver),0,.8)
    kernel <- exp(-d/par$radius)
    if(par$scenario=="boundary") {
      # The vertical anatomical boundary attenuates cross-boundary injection.
      nn <- RANN::nn2(as.matrix(candidate[source,c("x","y")]),as.matrix(receiver[,c("x","y")]),k=1)$nn.idx[,1]
      crossed <- (receiver$x>=1200)!=(candidate$x[source[nn]]>=1200)
      kernel[crossed] <- .2*kernel[crossed]
    }
    for (amplitude in c(0,.5,1)) {
      score <- base_score+amplitude*kernel
      if(par$scenario=="dropout_50pct") {
        # Score-level masking is a proxy, not a mechanistic gene-count model.
        score <- score*(runif(length(score))>.5)
      }
      for (adjust in c("none","linear","spatial")) {
        residual <- spf_adjust(score,receiver,adjust)
        obs <- sum(contrast*residual)
        for (null in names(banks)) {
          values <- as.numeric(crossprod(residual,banks[[null]]))
          zz <- spf_pvalue(obs,values)
          legacy_p <- (1+sum(abs(values-mean(values))>=abs(obs-mean(values))))/(B+1)
          nr<-nr+1L; all_rows[[nr]]<-data.frame(scenario=par$scenario,replicate=rr,
            amplitude=amplitude,adjustment=adjust,null_model=null,p=zz["p"],
            p_null_only=legacy_p,
            null_adjusted_effect=zz["effect"],valid_permutations=zz["valid"],
            injected_near_far=amplitude*sum(contrast*kernel),n_receiver=nrow(receiver))
        }
      }
    }
  }
  # Checkpoint after each scenario.
  fwrite(rbindlist(all_rows),file.path(out,"synthetic_replicates.tsv"),sep="\t")
}
replicates <- rbindlist(all_rows)
summary <- replicates[,{
  ok<-is.finite(p); n<-sum(ok); k<-sum(p[ok]<=.05); ci<-spf_wilson(k,n)
  list(n_replicates=n,rejections=k,rejection_rate=if(n)k/n else NA_real_,
       mc_lower=ci[1],mc_upper=ci[2],mean_effect=mean(null_adjusted_effect[ok]),
       mean_injected_near_far=mean(injected_near_far[ok]),
       min_valid_permutations=min(valid_permutations))
},by=.(scenario,amplitude,adjustment,null_model)]
fwrite(summary,file.path(out,"synthetic_summary.tsv"),sep="\t")


invisible(summary)
}
