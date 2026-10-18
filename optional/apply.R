#!/usr/bin/env Rscript
# Exploratory new-data interface; all coordinates and windows are in um.
# Rscript optional/apply.R receivers.tsv sources.tsv output_prefix [B=499] [seed=20261015]
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<3L)stop("Usage: Rscript optional/apply.R receivers.tsv sources.tsv output_prefix [B] [seed]")
suppressPackageStartupMessages(library(data.table))
script<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
source(file.path(dirname(script),"../R/spf_core.R"))
receiver <- as.data.frame(fread(args[1])); source <- as.data.frame(fread(args[2]))
required <- c("x","y","component","niche")
stopifnot(all(c(required,"identity","score") %in% names(receiver)),all(required %in% names(source)))
stopifnot(all(is.finite(receiver$score)),all(is.finite(receiver$x)),all(is.finite(receiver$y)),
          all(is.finite(source$x)),all(is.finite(source$y)))
# This coordinate-only interface cannot resolve distinct units at one address.
# Enforce the documented source/receiver separation before fitting any score.
coordinate_key <- function(z) paste(z$component,sprintf("%.17g",z$x),
                                    sprintf("%.17g",z$y),sep=":")
receiver_key <- coordinate_key(receiver); source_key <- coordinate_key(source)
if(anyDuplicated(receiver_key) || anyDuplicated(source_key))
  stop("Duplicate component/x/y coordinates are ambiguous in this coordinate-only interface; resolve the measurement-unit identities and supply one row per unit before analysis.")
overlap <- sum(receiver_key %in% source_key)
if(overlap)
  stop(sprintf("%d supplied receiver row(s) coincide with observed source component/x/y coordinates. Observed source units cannot be receivers; remove them from the receiver input before analysis. Distinct measurements at identical coordinates require an explicitly resolved input design.",overlap))
B <- if(length(args)>=4L)as.integer(args[4])else 499L
seed <- if(length(args)>=5L)as.integer(args[5])else 20261015L
stopifnot(B>=99L,nrow(source)>=1L,nrow(receiver)>nrow(source))
set.seed(seed)
spec <- spf_cluster_spec(source,threshold=100)
bank <- spf_null_bank(receiver,receiver,spec,B)
d <- spf_distance(receiver,source,seq_len(nrow(source)))
v <- spf_contrast(d)
if(anyNA(v))stop("Insufficient near/far receiver support for the specified windows")
ans <- list(); k<-0L
for(adjustment in c("none","linear","spatial")) {
  score <- spf_adjust(receiver$score,receiver,adjustment)
  for(null in names(bank)) {
    zz <- spf_pvalue(sum(v*score),as.numeric(crossprod(score,bank[[null]])))
    k<-k+1L;ans[[k]]<-data.frame(adjustment=adjustment,null_model=null,
                                p=zz["p"],effect=zz["effect"],valid_permutations=zz["valid"])
  }
}
dir.create(dirname(args[3]),recursive=TRUE,showWarnings=FALSE)
fwrite(rbindlist(ans),paste0(args[3],"_tests.tsv"),sep="\t")
fwrite(data.table(x=receiver$x,y=receiver$y,distance=d,score=receiver$score),
       paste0(args[3],"_distances.tsv"),sep="\t")
writeLines(c(paste("seed",seed),paste("B",B),"near <= 100 um; far > 400 and <= 800 um",
             "source-cluster threshold 100 um; minimum 5 receivers per window",
             "niche labels supplied by the user; source cluster anchors are matched",
             "This is an exploratory interface; sample-level replication is not supplied by cell counts.",
             capture.output(sessionInfo())),paste0(args[3],"_manifest.txt"))
