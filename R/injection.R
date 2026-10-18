# Single-chip deterministic sensitivities; these do not estimate biological power.
run_injection <- function(cells, out, B=199L) {
seed <- 20261015L
# Real-coordinate injection sensitivity: one exported Clu chip, not three-chip
# re-analysis. Keep all observed sources, sample receivers/candidate controls.
set.seed(seed+900000L)
z <- as.data.frame(cells[cells$chip_id=="B03018A2",])
z$identity <- as.character(z$identity)
z$guide_call <- as.character(z$guide_call)
z$chip_id <- as.character(z$chip_id)
z$library_size <- NULL
z$component <- NULL
stopifnot(length(unique(z$chip_id))==1L)
names(z)[match(c("x_um","y_um"),names(z))] <- c("x","y")
# Reconstruct the published occupied-grid connected components at 250 um.
bins <- unique(data.frame(bx=floor(z$x/250),by=floor(z$y/250)))
adj <- as.matrix(dist(bins)); ed<-which(adj>0 & adj<=sqrt(2)+1e-8,arr.ind=TRUE)
g <- igraph::make_empty_graph(nrow(bins),directed=FALSE)
if(nrow(ed))g<-igraph::add_edges(g,as.vector(t(ed)))
lookup <- setNames(igraph::components(g)$membership,paste(bins$bx,bins$by))
z$component <- unname(lookup[paste(floor(z$x/250),floor(z$y/250))])
z$niche <- as.integer(z$x>median(z$x))+2L*as.integer(z$y>median(z$y))
controls <- which(z$guide_call=="Neg")
receiver <- z[sample(controls,min(2500L,length(controls))),]
candidate <- z[sample(controls,min(8000L,length(controls))),]
actual <- z[z$guide_call=="sgrna-clu",]
spec <- spf_cluster_spec(actual,100)
saveRDS(list(receiver=receiver,candidate=candidate,source=actual,spec=spec,seed=seed),
        file.path(out,"real_coordinate_worked_example.rds"))
fwrite(receiver[,c("x","y","component","niche","identity","raw_clu_score")],
       file.path(out,"worked_example_receivers.tsv"),sep="\t")
fwrite(actual[,c("x","y","component","niche","identity")],
       file.path(out,"worked_example_sources.tsv"),sep="\t")
banks <- spf_null_bank(receiver,candidate,spec,B)
d <- spf_distance(receiver,actual,seq_len(nrow(actual))); contrast<-spf_contrast(d)
raw <- receiver$raw_clu_score
rows <- list(); nr<-0L
for(radius in c(75,150,300,600)) {
  kernel<-exp(-d/radius)
  for(adjust in c("none","linear","spatial")) {
    fixed_fit<-spf_adjust(raw,receiver,adjust,return_fitted=TRUE)
    for(amplitude in c(0,.5,1,2)) for(mode in c("fixed","refit")) {
      score<-raw+amplitude*kernel
      residual<-if(mode=="fixed")score-fixed_fit else spf_adjust(score,receiver,adjust)
      obs<-sum(contrast*residual)
      for(null in names(banks)) {
        zz<-spf_pvalue(obs,as.numeric(crossprod(residual,banks[[null]])))
        nr<-nr+1L; rows[[nr]]<-data.frame(chip=unique(z$chip_id),radius=radius,
          amplitude=amplitude,adjustment=adjust,background_mode=mode,null_model=null,
          p=zz["p"],null_adjusted_effect=zz["effect"],valid_permutations=zz["valid"],
          observed_near_far=obs,injected_near_far=amplitude*sum(contrast*kernel))
      }
    }
  }
}
semi <- rbindlist(rows)
fwrite(semi,file.path(out,"real_coordinate_injection.tsv"),sep="\t")
fwrite(data.table(chip=unique(z$chip_id),total_cells=nrow(z),control_cells=length(controls),
                 sampled_receivers=nrow(receiver),candidate_controls=nrow(candidate),
                 source_cells=nrow(actual),source_clusters=nrow(spec),largest_cluster=max(spec$size),
                 components=length(unique(z$component)),near_receivers=sum(d<=100),
                 far_receivers=sum(d>400 & d<=800),seed=seed,permutations=B),
       file.path(out,"real_coordinate_design.tsv"),sep="\t")

invisible(semi)
}

run_dropout <- function(design, outdir, B=199L) {
z <- design
receiver<-z$receiver;candidate<-z$candidate;truth<-z$source
set.seed(20261015L+700000L);order_source<-sample.int(nrow(truth))
d_true<-spf_distance(receiver,truth,seq_len(nrow(truth)))
kernel<-exp(-d_true/150);score<-receiver$raw_clu_score+kernel
resids<-list()
for(bg in c('none','linear','spatial'))for(mode in c('fixed','refit')){
 fit<-if(mode=='fixed')spf_adjust(receiver$raw_clu_score,receiver,bg,TRUE)else spf_adjust(score,receiver,bg,TRUE)
 resids[[paste(bg,mode,sep=':')]]<-score-fit
}
rows<-list();idx<-0L
for(drop in c(0,.25,.5)) {
 observed<-truth[order_source[seq_len(floor(nrow(truth)*(1-drop)))],]
 spec<-spf_cluster_spec(observed,100)
 set.seed(20261015L+710000L+as.integer(drop*100))
 bank<-spf_null_bank(receiver,candidate,spec,B)
 d<-spf_distance(receiver,observed,seq_len(nrow(observed)));v<-spf_contrast(d)
 for(key in names(resids))for(null in names(bank)){
  bgmode<-strsplit(key,':',fixed=TRUE)[[1]];resid<-resids[[key]]
  out<-spf_pvalue(sum(v*resid),as.numeric(crossprod(resid,bank[[null]])))
  idx<-idx+1L;rows[[idx]]<-data.frame(dropout_fraction=drop,observed_sources=nrow(observed),
   true_sources=nrow(truth),observed_clusters=nrow(spec),largest_observed_cluster=max(spec$size),
   near_receivers=sum(d<=100),far_receivers=sum(d>400&d<=800),
   adjustment=bgmode[1],background_mode=bgmode[2],null_model=null,
   null_adjusted_effect=out['effect'],p=out['p'],valid_permutations=out['valid'],
   known_injected_contrast=sum(v*kernel))
 }
 message('Completed source thinning ',drop)
}
a<-rbindlist(rows);fwrite(a,file.path(outdir,'guide_dropout_sensitivity.tsv'),sep='\t')

invisible(a)
}
