#!/usr/bin/env Rscript
# Optional raw-feature interfaces. Requires the original public count matrices.
# They do not download data or install software. No plots are made.
a<-commandArgs(trailingOnly=TRUE)
if(!length(a)||!a[1]%in%c("mi","hcc"))stop("Usage: Rscript optional/external_features.R mi MI_H5AD_DIRECTORY OUTPUT_DIRECTORY; or hcc PAIRED_H5AD EXTERNAL_VISIUM_DIRECTORY OUTPUT_DIRECTORY [hex|pixel|legacy] [training_compatible|legacy_available] [assay|common_unique]")
suppressPackageStartupMessages({library(rhdf5);library(Matrix);library(data.table);library(RANN)})
run_mi_raw <- function(root,out) {
files <- sort(list.files(root, pattern = "^Visium_.*[.]h5ad$", full.names = TRUE))
stopifnot(length(files) == 27L)

programs <- list(
  fibrotic = c("COL1A1", "COL1A2", "COL3A1", "FN1", "POSTN", "TIMP1", "TNC"),
  inflammatory = c("LYZ", "S100A8", "S100A9", "IL1B", "CCL2", "TYROBP"),
  vascular = c("PECAM1", "VWF", "KDR", "EMCN", "ENG", "RAMP2"),
  regenerative = c("MKI67", "TOP2A", "MYC", "IGF1", "NRG1", "AREG")
)

read_selected <- function(path, wanted) {
  genes <- h5read(path, "/raw/var/_index")
  data <- as.numeric(h5read(path, "/raw/X/data"))
  indices <- as.integer(h5read(path, "/raw/X/indices"))
  indptr <- as.integer(h5read(path, "/raw/X/indptr"))
  n_spots <- length(indptr) - 1L
  # AnnData raw/X is stored as spot-by-gene CSR in these files.
  mat <- sparseMatrix(i = rep.int(seq_len(n_spots), diff(indptr)), j = indices + 1L,
                      x = data, dims = c(n_spots, length(genes)))
  keep <- which(toupper(genes) %in% toupper(wanted))
  out <- t(mat[, keep, drop = FALSE])
  rownames(out) <- toupper(genes[keep])
  full_lib <- as.numeric(Matrix::rowSums(mat))
  supplied_lib <- as.numeric(h5read(path, "/obs/n_counts"))
  stopifnot(length(full_lib) == length(supplied_lib),
            all(is.finite(full_lib)), all(full_lib >= 0),
            max(abs(full_lib-supplied_lib)) < 1e-7,
            all(abs(data-round(data)) < 1e-7))
  attr(out, "whole_transcriptome_library") <- full_lib
  attr(out, "normalization_audit") <- data.table(
    sample_id = sub("[.]h5ad$", "", basename(path)),
    raw_genes = length(genes), selected_genes = nrow(out), spots = ncol(out),
    totalcount_vs_obs_n_counts_max_error = max(abs(full_lib-supplied_lib)),
    fraction_programme_panel_counts_median = median(colSums(out)/pmax(full_lib,1)))
  out
}

process_section <- function(path) {
  sample_id <- sub("[.]h5ad$", "", basename(path))
  patient <- sub(".*_(P[0-9]+)(?:_rep[0-9]+)?$", "\\1", sample_id, perl = TRUE)
  zone <- sub("^Visium_", "", sub("_P[0-9]+(?:_rep[0-9]+)?$", "", sample_id, perl = TRUE))
  primary_zone <- if (zone == "control") "Control" else if (grepl("IZ", zone)) "Infarct" else
    if (grepl("GT", zone)) "Granulation" else if (grepl("FZ", zone)) "Fibrotic" else
      if (grepl("BZ", zone)) "Border" else "Remote"
  xy <- h5read(path, "/obsm/X_spatial")
  if (nrow(xy) == 2) xy <- t(xy)
  n <- nrow(xy)
  cm <- as.numeric(h5read(path, "/obs/Cardiomyocyte"))
  my <- as.numeric(h5read(path, "/obs/Myeloid"))
  source_index <- as.numeric(scale(-cm)) + as.numeric(scale(my))
  source <- which(source_index >= quantile(source_index, 0.95, na.rm = TRUE))
  receiver <- setdiff(seq_len(n), source)
  base_nn <- nn2(xy, xy, k = 2)$nn.dists[, 2]
  grid <- median(base_nn[base_nn > 0], na.rm = TRUE)
  coords_um <- xy / grid * 100
  near <- nn2(coords_um[source, , drop = FALSE], coords_um[receiver, , drop = FALSE], k = 1)
  d <- near$nn.dists[, 1]

  mat <- read_selected(path, unique(unlist(programs)))
  normalization_audit <- attr(mat, "normalization_audit")
  lib <- pmax(attr(mat, "whole_transcriptome_library"), 1)
  mat <- mat %*% Diagonal(x = 1 / lib); mat@x <- log1p(1e4 * mat@x)
  score <- sapply(programs, function(gs) {
    use <- intersect(toupper(gs), rownames(mat))
    v <- as.numeric(colMeans(mat[use, , drop = FALSE]))
    center <- median(v[receiver]); spread <- mad(v[receiver], constant = 1.4826)
    if (!is.finite(spread) || spread == 0) spread <- sd(v[receiver])
    (v - center) / spread
  })
  bins <- cut(d, c(-Inf, 100, 200, 400, 800, Inf),
              labels = c("0-100", "100-200", "200-400", "400-800", ">800"))
  prof <- rbindlist(lapply(colnames(score), function(pr) {
    data.table(distance_bin = bins, value = score[receiver, pr])[
      , .(mean_score = mean(value), receivers = .N), by = distance_bin][,
        `:=`(sample_id = sample_id, patient = patient, zone = zone,
             primary_zone = primary_zone, program = pr)]
  }))
  eff <- rbindlist(lapply(colnames(score), function(pr) {
    sc <- score[receiver, pr]
    data.table(sample_id, patient, zone, primary_zone, program = pr,
               source_spots = length(source),
               near_far_effect = mean(sc[d <= 200], na.rm = TRUE) -
                 mean(sc[d > 400 & d <= 800], na.rm = TRUE),
               source_score = mean(score[source, pr]))
  }))
  list(profile = prof, effect = eff, normalization_audit = normalization_audit)
}

res <- lapply(files, process_section)
normalization_audit <- rbindlist(lapply(res, `[[`, "normalization_audit"))
fwrite(normalization_audit, file.path(out,"mi_total_count_normalization_audit.tsv"), sep="\t")
profile <- rbindlist(lapply(res, `[[`, "profile"))
effect <- rbindlist(lapply(res, `[[`, "effect"))
profile[, primary_zone := factor(primary_zone,
  levels = c("Control", "Remote", "Border", "Fibrotic", "Granulation", "Infarct"))]
effect[, primary_zone := factor(primary_zone, levels = levels(profile$primary_zone))]
fwrite(profile, file.path(out,"mi_injury_proxy_field_profiles.tsv"), sep = "\t")
fwrite(effect, file.path(out,"mi_injury_proxy_field_effects.tsv"), sep = "\t")

patient <- effect[, .(near_far_effect = mean(near_far_effect),
                      is_control = all(primary_zone == "Control")), by = .(patient, program)]
tests <- rbindlist(lapply(unique(patient$program), function(pr) {
  z <- patient[program == pr]
  n_control <- sum(z$is_control)
  observed <- mean(z[is_control == FALSE, near_far_effect]) -
    mean(z[is_control == TRUE, near_far_effect])
  null <- combn(nrow(z), n_control, function(control) {
    mean(z[-control, near_far_effect]) - mean(z[control, near_far_effect])
  })
  data.table(program = pr, patients = nrow(z), controls = n_control,
             observed_injury_minus_control = observed,
             exact_assignments = length(null),
             exact_p_two_sided = mean(abs(null) >= abs(observed) - 1e-12))
}))
fwrite(patient, file.path(out,"mi_injury_patient_field_embedding.tsv"), sep = "\t")
fwrite(tests, file.path(out,"mi_injury_exact_tests.tsv"), sep = "\t")


invisible(tests)
}
run_hcc_raw <- function(paired_path,external_root,out,coordinate_mode="hex",panel_mode="training_compatible",library_mode="assay") {
stopifnot(panel_mode %in% c("training_compatible","legacy_available"),library_mode %in% c("assay","common_unique"))
# Explicit denominator sensitivity: intersect symbols unique in every input assay.
# This is a sensitivity to library composition, not a batch calibration procedure.
if(library_mode=="common_unique") {
  training_symbols<-as.character(h5read(paired_path,"/var/_index"))
  universe_dirs<-list.dirs(external_root,recursive=FALSE,full.names=TRUE)
  universe_dirs<-universe_dirs[file.exists(file.path(universe_dirs,"filtered_feature_bc_matrix.h5"))]
  stopifnot(length(universe_dirs)>0L,!anyNA(training_symbols))
  external_symbols<-lapply(universe_dirs,function(dir)toupper(h5read(file.path(dir,"filtered_feature_bc_matrix.h5"),"/matrix/features/name")))
  stopifnot(all(vapply(external_symbols,function(g)!anyNA(g),logical(1))))
  all_symbols<-c(list(training_symbols),external_symbols)
  ambiguous_symbols<-Reduce(union,lapply(all_symbols,function(g)unique(g[duplicated(g)|duplicated(g,fromLast=TRUE)])))
  shared_symbols<-Reduce(intersect,all_symbols)
  common_genes<-setdiff(shared_symbols,ambiguous_symbols)
  stopifnot(length(common_genes)>0L)
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  fwrite(data.table(gene=common_genes),file.path(out,"fixed_common_library_genes.tsv"),sep="\t")
  fwrite(data.table(gene=intersect(shared_symbols,ambiguous_symbols),reason="Duplicated symbol in one or more raw assays; no unambiguous cross-assay match"),file.path(out,"excluded_ambiguous_shared_symbols.tsv"),sep="\t")
}
programs <- list(
  apoptosis=c('BAX','BAK1','CASP3','CASP7','CASP8','FAS','BBC3','PMAIP1'),
  interferon=c('CXCL9','CXCL10','CXCL11','STAT1','IRF1','GBP1','ISG15'),
  cytotoxic=c('CD3D','CD3E','CD8A','NKG7','GNLY','GZMB','PRF1'),
  fibrosis=c('COL1A1','COL1A2','COL3A1','DCN','LUM','ACTA2','TAGLN','POSTN'),
  hypoxia=c('CA9','VEGFA','SLC2A1','PGK1','ENO1'))
# Freeze the 35 genes measured by the training assay for both cohorts.
# The original available-gene branch is retained explicitly for audit only.
if(panel_mode=="legacy_available") {
 programs$cytotoxic<-c('CD3D','CD3E','CD8A','CD8B','NKG7','GNLY','GZMB','PRF1')
 programs$hypoxia<-c('CA9','VEGFA','SLC2A1','LDHA','PGK1','ENO1','BNIP3')
}
read_paired<-function(path) {
  genes <- h5read(path,'/var/_index')
  if(panel_mode=="training_compatible" && !all(unique(unlist(programs)) %in% genes))stop("The fixed training-compatible programme requires all 35 genes; missing genes cannot be silently omitted.")
  ids <- match(intersect(unique(unlist(programs)),genes),genes)
  ip <- as.numeric(h5read(path,'/X/indptr')); n <- length(ip)-1L
  selected <- matrix(0,length(ids),n,dimnames=list(genes[ids],NULL));lib<-numeric(n)
  for(start in seq.int(1L,n,by=1000L)) {
    end <- min(n,start+999L);begin<-ip[start];last<-ip[end+1L]
    d <- as.numeric(h5read(path,'/X/data',start=begin+1,count=last-begin))
    gi <- as.integer(h5read(path,'/X/indices',start=begin+1,count=last-begin))+1L
    stopifnot(all(d>=0),all(abs(d-round(d))<1e-6))
    library_counts<-if(library_mode=="common_unique")d*(genes[gi] %in% common_genes)else d
    cum <- c(0,cumsum(library_counts));lib[start:end] <- diff(cum[ip[start:(end+1L)]-begin+1L])
    keep <- which(gi %in% ids)
    if(length(keep)) {
      ci <- findInterval(keep-1L,ip[start:(end+1L)]-begin,all.inside=TRUE)
      selected[,start:end] <- as.matrix(sparseMatrix(i=match(gi[keep],ids),j=ci,x=d[keep],dims=c(length(ids),end-start+1L)))
    }
    if(start%%10000L==1L) message('Paired count cache spots ',start,'-',end,'/',n)
  }
  lognorm <- log1p(sweep(selected,2,pmax(lib,1),'/')*1e4)
  scores <- sapply(programs,function(gs)colMeans(lognorm[intersect(gs,rownames(lognorm)),,drop=FALSE]))

  read_obs <- function(nm) {
    v<-h5read(path,paste0("/obs/",nm))
    if(is.list(v)&&all(c("categories","codes")%in%names(v)))return(as.character(v$categories)[as.integer(v$codes)+1L])
    as.character(v)
  }
  md<-data.table(spot_id=as.character(h5read(path,"/obs/index")))
  for(nm in c("patient_id","Response","diagnosis","sample_id","array_row","array_col"))md[[nm]]<-read_obs(nm)
  md[,`:=`(array_row=as.numeric(array_row),array_col=as.numeric(array_col))]
  pixel <- h5read(path,"/obsm/spatial")
  if(nrow(pixel)==2L)pixel<-t(pixel)
  stopifnot(nrow(pixel)==nrow(md),ncol(pixel)==2L)
  md[,`:=`(px_col=as.numeric(pixel[,1]),px_row=as.numeric(pixel[,2]))]
  cbind(md,data.table(total_count=lib),as.data.table(scores))
}
# Geometry and thresholds are declared explicitly; outcomes never select a mode.
field_features <- function(z,sample_id,patient_id,response,mode="hex",source_quantile=.98,bin_scale=1) {
  stopifnot(mode%in%c("hex","pixel","legacy"),source_quantile>0,source_quantile<1,bin_scale>0)
  scores <- as.matrix(z[,names(programs),with=FALSE]);zs<-scale(scores)
  source <- zs[,"apoptosis"]>=quantile(zs[,"apoptosis"],source_quantile)
  receiver <- !source
  xy <- switch(mode,hex=as.matrix(z[,.(array_col/2,sqrt(3)*array_row/2)]),
    pixel=as.matrix(z[,.(px_col,px_row)]),legacy=as.matrix(z[,.(array_col,array_row)]))
  stopifnot(all(is.finite(xy)))
  near <- nn2(xy,xy,k=2)$nn.dists[,2];grid<-median(near[near>0])
  d <- nn2(xy[source,,drop=FALSE],xy,k=1)$nn.dists[,1]/grid
  # Snap mathematical bin-boundary ties, not image-coordinate deviations.
  d_bin<-d
  if(mode=="hex")for(bound in c(1,2,4,8)*bin_scale)d_bin[abs(d_bin-bound)<1e-10]<-bound
  feat <- data.table(sample_id,patient_id,Response=response,post_apoptosis_q95=quantile(scores[,"apoptosis"],.95))
  profiles<-list();counts<-data.table(sample_id,patient_id,spots=nrow(z),sources=sum(source),receivers=sum(receiver),coordinate_mode=mode,source_quantile,bin_scale,median_nn_spacing=grid)
  for(pr in setdiff(names(programs),"apoptosis")) {
    feat[[paste0("field_slope_",pr)]] <- coef(lm(zs[receiver,pr]~log1p(d[receiver])))[2]
    bins<-cut(d_bin[receiver],c(-Inf,c(1,2,4,8)*bin_scale,Inf),labels=c("0-1","1-2","2-4","4-8",">8"))
    profiles[[pr]]<-data.table(sample_id,patient_id,Response=response,program=pr,distance_bin=bins,score=zs[receiver,pr])[,.(mean_score=mean(score),spots=.N),by=.(sample_id,patient_id,Response,program,distance_bin)]
  }
  list(feature=feat,profile=rbindlist(profiles),counts=counts,
       distance=data.table(sample_id,spot_index=seq_len(nrow(z)),source,receiver,distance=d))
}

md<-read_paired(paired_path)
dirs <- list.dirs(external_root,recursive=FALSE,full.names=TRUE);dirs<-dirs[file.exists(file.path(dirs,'filtered_feature_bc_matrix.h5'))]
external_scores <- lapply(dirs,function(dir){
  p<-file.path(dir,'filtered_feature_bc_matrix.h5');shape<-as.integer(h5read(p,'/matrix/shape'))
  counts<-sparseMatrix(i=as.integer(h5read(p,'/matrix/indices'))+1L,j=rep.int(seq_len(shape[2]),diff(as.integer(h5read(p,'/matrix/indptr')))),x=as.numeric(h5read(p,'/matrix/data')),dims=shape)
  rownames(counts)<-toupper(h5read(p,'/matrix/features/name'));colnames(counts)<-h5read(p,'/matrix/barcodes')
  if(panel_mode=="training_compatible" && !all(unique(unlist(programs)) %in% rownames(counts)))stop("The fixed training-compatible programme requires all 35 genes in each external assay.")
  pos<-fread(file.path(dir,'spatial/tissue_positions_list.csv'),header=FALSE,col.names=c('barcode','in_tissue','array_row','array_col','px_row','px_col'))[in_tissue==1]
  common<-intersect(pos$barcode,colnames(counts));pos<-pos[match(common,barcode)];counts<-counts[,common,drop=FALSE]
  stopifnot(all(counts@x>=0),all(abs(counts@x-round(counts@x))<1e-7))
  lib<-pmax(if(library_mode=="common_unique")colSums(counts[rownames(counts) %in% common_genes,,drop=FALSE])else colSums(counts),1);gs<-intersect(unique(unlist(programs)),rownames(counts));sm<-counts[gs,,drop=FALSE]%*%Diagonal(x=1/lib);sm@x<-log1p(1e4*sm@x)
  scores<-sapply(programs,function(v)as.numeric(colMeans(sm[intersect(v,rownames(sm)),,drop=FALSE])))
  z<-cbind(pos,as.data.table(scores));sid<-basename(dir)
  z[,`:=`(sample_id=sid,patient_id=sid,Response=if(grepl('NR$',sid))'Non_Responder'else'Responder')]
  z
})
# Reuse these compact scores/coordinates for prespecified sensitivity branches.
process_scores <- function(raw,out,mode="hex",source_quantile=.98,bin_scale=1) {
  md<-raw$paired;post<-md[diagnosis=="Post"]
  stopifnot(uniqueN(post$sample_id)==11L)
  internal_res<-lapply(unique(post$sample_id),function(sid){z<-post[sample_id==sid];field_features(z,sid,z$patient_id[1],z$Response[1],mode,source_quantile,bin_scale)})
  internal<-rbindlist(lapply(internal_res,`[[`,"feature"))
  means<-md[,lapply(.SD,mean),by=.(patient_id,Response,diagnosis),.SDcols=names(programs)]
  wide<-dcast(means,patient_id+Response~diagnosis,value.var=names(programs),sep="_")
  wide[,delta_apoptosis:=apoptosis_Post-apoptosis_Pre]
  internal<-merge(wide,internal,by=c("patient_id","Response"));internal[,outcome:=as.integer(Response=="Responder")]
  external_res<-lapply(raw$external,function(z)field_features(z,z$sample_id[1],z$patient_id[1],z$Response[1],mode,source_quantile,bin_scale))
  external<-rbindlist(lapply(external_res,`[[`,"feature"));external[,outcome:=as.integer(Response=="Responder")]
  vals<-list(hcc_patient_embedding_features=internal,hcc_proxy_field_patient_distance=rbindlist(lapply(internal_res,`[[`,"profile")),
    hcc_external_spf_features=external,hcc_external_spf_distance_profiles=rbindlist(lapply(external_res,`[[`,"profile")),
    hcc_geometry_audit=rbindlist(c(lapply(internal_res,`[[`,"counts"),lapply(external_res,`[[`,"counts"))))
  for(nm in names(vals))fwrite(vals[[nm]],file.path(out,paste0(nm,".tsv")),sep="\t")
  invisible(c(vals,list(distances=rbindlist(c(lapply(internal_res,`[[`,"distance"),lapply(external_res,`[[`,"distance"))))))
}
raw<-list(paired=md,external=external_scores)
ans<-lapply(coordinate_mode,function(mode){od<-if(length(coordinate_mode)>1L)file.path(out,mode)else out;dir.create(od,recursive=TRUE,showWarnings=FALSE);process_scores(raw,od,mode)})
names(ans)<-coordinate_mode
invisible(list(results=ans,spot_scores=raw,process_scores=process_scores))
}

if(a[1]=="mi") {stopifnot(length(a)==3L);out<-normalizePath(a[3],mustWork=FALSE);dir.create(out,recursive=TRUE,showWarnings=FALSE);run_mi_raw(normalizePath(a[2]),out)}
if(a[1]=="hcc") {stopifnot(length(a)%in%c(4L,5L,6L,7L));out<-normalizePath(a[4],mustWork=FALSE);dir.create(out,recursive=TRUE,showWarnings=FALSE);run_hcc_raw(normalizePath(a[2]),normalizePath(a[3]),out,if(length(a)>=5L)a[5]else"hex",if(length(a)>=6L)a[6]else"training_compatible",if(length(a)>=7L)a[7]else"assay")}
