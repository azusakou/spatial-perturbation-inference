#!/usr/bin/env Rscript
# Known-assignment synthetic check; no claim of real-data source exchangeability.
# Rscript optional/independent_bank.R PACKAGE_ROOT OUTPUT [REPLICATES=100] [B=199]
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<2L)stop('Usage: independent_bank.R PACKAGE_ROOT OUTPUT [REPLICATES=100] [B=199]')
root<-normalizePath(args[1],mustWork=TRUE)
out<-args[2]
n<-if(length(args)>=3L)as.integer(args[3])else 100L
B<-if(length(args)>=4L)as.integer(args[4])else 199L
stopifnot(is.finite(n),n>=1L,n<=100L,B==199L)
source(file.path(root,'R','spf_core.R'))
source(file.path(root,'R','independent_bank.R'))
run_independent_bank(out,n_rep=n,B=B)
