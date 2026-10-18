# Perturb-map: exact spot-level reconstruction followed by descriptive animal summaries.
# mouse4 contributes KP_1 and KP_2; mouse7 contributes KP_4. Repeated sections
# never increase the number of independent animals. Clone independence is unknown.
run_perturbmap <- function(spot_path,outdir,raw=NA_character_) {
 suppressPackageStartupMessages({library(data.table);library(RANN)})
 dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
 here <- dirname(spot_path)
 endpoint_path <- file.path(outdir,"sensitivity_lesion_endpoints.tsv.gz")
 out <- function(d,name)fwrite(d,file.path(outdir,name),sep="\t",na="NA")
sections <- c("KP_1","KP_2","KP_4")
mouse <- c(KP_1="mouse4",KP_2="mouse4",KP_3="mouse5",KP_4="mouse7")
programs <- list(
  primary_stromal=c("Serpine1","Ctgf","Col1a1","Col1a2","Col3a1","Fn1","Tnc","Timp1","Tgfb1","Spp1"),
  downstream_target=c("Serpine1","Ctgf","Smad7","Skil","Pmepa1","Tagln","Acta2"),
  extracellular_matrix=c("Col1a1","Col1a2","Col3a1","Fn1","Tnc","Timp1"),
  compact_stromal=c("Serpine1","Ctgf","Col1a1","Col1a2","Fn1","Timp1"),
  t_cell=c("Cd3d","Cd3e","Cd3g","Trac","Cd8a","Cd8b1","Lck","Lat"))
mechanism_programs <- list(tgfb_ligand=c("Tgfb1","Tgfb2","Tgfb3"),
  contact_ligand=c("Jag1","Jag2","Dll1","Dll4"),
  inflammatory_ligand=c("Ccl2","Ccl5","Cxcl9","Cxcl10","Il6"),
  matrix_mediator=c("Col1a1","Col1a2","Col3a1","Fn1","Tnc"),
  vascular_mediator=c("Vegfa","Vegfb","Angpt1","Angpt2"),
  mechanical_stress=c("Cyr61","Ctgf","Ankrd1","Amotl2"),
  tgfb_receptor=c("Tgfbr1","Tgfbr2","Tgfbr3"),
  tgfb_target=c("Serpine1","Ctgf","Smad7","Acta2","Tagln","Timp1","Col1a1","Col1a2","Col3a1","Fn1","Tnc"))
score_programs <- c(programs,mechanism_programs)
grid <- CJ(program=names(programs),near_um=c(100,150,200,300),
  far_start_um=c(400,500),far_end_um=c(800,Inf))
if(is.finite(nchar(raw)) && !is.na(raw)) {
  suppressPackageStartupMessages({library(Matrix);library(rhdf5);library(RANN);library(jsonlite)})
  stopifnot(dir.exists(raw))
  rawfiles <- list(); spotpieces <- list(); rawpieces <- list(); score_defs <- list()
  for(sec in names(mouse)) {
    annfile <- Sys.glob(file.path(raw,paste0("GSM*_",sec,"_spot_annotation.csv.gz")))
    stopifnot(length(annfile)==1L)
    stem <- sub("_spot_annotation[.]csv[.]gz$","",basename(annfile))
    posfile <- file.path(raw,paste0(stem,"_tissue_positions_list.csv.gz"))
    h5file <- file.path(raw,paste0(stem,"_filtered_feature_bc_matrix.h5"))
    sffile <- file.path(raw,paste0(stem,"_scalefactors_json.json.gz"))
    rawfiles[[sec]] <- data.table(section_id=sec,file=c(annfile,posfile,h5file,sffile))
    z <- merge(fread(annfile),fread(posfile,header=FALSE,col.names=c("barcode","in_tissue",
      "array_row","array_col","px_row","px_col")),by="barcode",all.x=TRUE)
    shape <- as.integer(h5read(h5file,"/matrix/shape"))
    x <- sparseMatrix(i=as.integer(h5read(h5file,"/matrix/indices"))+1L,
      j=rep.int(seq_len(shape[2]),diff(as.integer(h5read(h5file,"/matrix/indptr")))),
      x=as.numeric(h5read(h5file,"/matrix/data")),dims=shape)
    rownames(x) <- h5read(h5file,"/matrix/features/name")
    colnames(x) <- h5read(h5file,"/matrix/barcodes")
    common <- intersect(z$barcode,colnames(x)); z <- z[match(common,barcode)];x <- x[,common,drop=FALSE]
    lib <- pmax(Matrix::colSums(x),1)
    sf <- fromJSON(paste(readLines(gzfile(sffile),warn=FALSE),collapse=""))
    xy <- as.matrix(z[,.(px_col,px_row)])*55/sf$spot_diameter_fullres
    reference <- is.na(z$phenotypes)|z$phenotypes=="periphery"
    z[,lesion_class:=fcase(grepl("^Tgfbr2",phenotypes),"Tgfbr2",
      grepl("^KP_",phenotypes),"KP control",grepl("^Ifngr2",phenotypes),"Ifngr2",
      grepl("^Jak2",phenotypes),"Jak2",default=NA_character_)]
    receiver <- which(is.na(z$lesion_class))
    scores <- sapply(names(score_programs),function(nm) {
      genes <- score_programs[[nm]]
      genes <- intersect(genes,rownames(x));stopifnot(length(genes)>0L)
      if(nm%in%names(programs))stopifnot(length(genes)>=4L)
      score_defs[[length(score_defs)+1L]] <<- data.table(section_id=sec,program=nm,
        requested_genes=paste(score_programs[[nm]],collapse=";"),available_genes=paste(genes,collapse=";"))
      v <- log1p(1e4*as.numeric(Matrix::colMeans(x[genes,,drop=FALSE]%*%Diagonal(x=1/lib))))
      sc <- mad(v[reference],constant=1.4826,na.rm=TRUE)
      if(!is.finite(sc)||sc==0)sc <- sd(v[reference],na.rm=TRUE)
      (v-median(v[reference],na.rm=TRUE))/sc
    })
    tgfbr2_expression <- as.numeric(scale(log1p(1e4*as.numeric(x["Tgfbr2",]/lib))))
    keep <- c("barcode","phenotypes","kmeans","in_tissue","array_row","array_col","px_row","px_col","lesion_class")
    spotpieces[[sec]] <- cbind(data.table(section_id=sec,mouse=mouse[[sec]]),z[,..keep],
      data.table(x_um=xy[,1],y_um=xy[,2],whole_library_size=lib,
        score_reference_receiver=reference,eligible_receiver=is.na(z$lesion_class),tgfbr2_expression=tgfbr2_expression),
      as.data.table(scores))
    gene_union <- intersect(unique(unlist(score_programs,use.names=FALSE)),rownames(x))
    raw_counts <- as.data.table(t(as.matrix(x[gene_union,,drop=FALSE])))
    setnames(raw_counts,paste0("gene__",gene_union))
    rawpieces[[sec]] <- cbind(data.table(section_id=sec,mouse=mouse[[sec]]),z[,..keep],
      data.table(x_um=xy[,1],y_um=xy[,2],whole_library_size=lib,
        score_reference_receiver=reference,eligible_receiver=is.na(z$lesion_class)),raw_counts)
  }
  spots <- rbindlist(spotpieces)
  dir.create(dirname(spot_path),recursive=TRUE,showWarnings=FALSE)
  fwrite(spots,file.path(outdir,"raw_matrix_score_reference.tsv.gz"),sep="\t",compress="gzip",na="NA")
  fwrite(rbindlist(rawpieces,fill=TRUE),spot_path,sep="\t",compress="gzip",na="NA")
  out(rbindlist(rawfiles),"raw_inputs.tsv")
  out(rbindlist(score_defs),"score_gene_definitions.tsv")
} else {
  stopifnot(file.exists(spot_path))
  spots <- fread(spot_path,na.strings=c("","NA"))
}
if(any(grepl("^gene__",names(spots)))) {
 # The compact input stores actual count rows, not subset-library totals or
 # precomputed scores. Missing original rows remain absent/NA, never zero-filled.
 compact_defs <- list()
 for(sec in names(mouse)) {
  ii <- which(spots$section_id==sec);lib <- pmax(spots$whole_library_size[ii],1)
  reference <- spots$score_reference_receiver[ii]
  for(nm in names(score_programs)) {
   genes <- intersect(score_programs[[nm]],sub("^gene__","",grep("^gene__",names(spots),value=TRUE)))
   genes <- genes[vapply(genes,function(g)all(!is.na(spots[[paste0("gene__",g)]][ii])),logical(1))]
   stopifnot(length(genes)>0L)
   gene_columns <- paste0("gene__",genes)
   v <- log1p(1e4*rowMeans(as.matrix(spots[ii,..gene_columns]))/lib)
   sc <- mad(v[reference],constant=1.4826,na.rm=TRUE)
   if(!is.finite(sc)||sc==0)sc <- sd(v[reference],na.rm=TRUE)
   set(spots,i=ii,j=nm,value=(v-median(v[reference],na.rm=TRUE))/sc)
   compact_defs[[length(compact_defs)+1L]] <- data.table(section_id=sec,program=nm,
     requested_genes=paste(score_programs[[nm]],collapse=";"),available_genes=paste(genes,collapse=";"))
  }
  set(spots,i=ii,j="tgfbr2_expression",value=as.numeric(scale(log1p(1e4*spots$gene__Tgfbr2[ii]/lib))))
 }
 out(rbindlist(compact_defs),"score_gene_definitions.tsv")
}
stopifnot(nrow(spots)==6363L,all(spots$eligible_receiver==is.na(spots$lesion_class)),
  all(spots$mouse==unname(mouse[spots$section_id])))
pieces <- list();features <- list();geometry <- list();i <- 0L
for(sec in sections) {
  z <- spots[section_id==sec];xy <- as.matrix(z[,.(x_um,y_um)]);receiver <- which(z$eligible_receiver)
  neighbours <- nn2(xy,xy,k=8)$nn.dists
  spacing <- median(neighbours[,2])
  edge <- rowSums(neighbours[,2:8,drop=FALSE]<1.25*spacing)<6
  stopifnot(any(edge))
  lesions <- unique(z[lesion_class%chin%c("Tgfbr2","KP control"),.(lesion=phenotypes,lesion_class)])
  for(j in seq_len(nrow(lesions))) {
    src <- which(z$phenotypes==lesions$lesion[j]);stopifnot(length(src)>=5L)
    d <- nn2(xy[src,,drop=FALSE],xy[receiver,,drop=FALSE],k=1)$nn.dists[,1]
    near <- receiver[d<=200];far <- receiver[d>400&d<=800]
    key <- data.table(section_id=sec,mouse=mouse[[sec]],lesion=lesions$lesion[j],
      lesion_class=lesions$lesion_class[j],source_spots=length(src))
    features[[length(features)+1L]] <- cbind(key,data.table(
      source_ligand=mean(z$tgfb_ligand[src]),source_contact=mean(z$contact_ligand[src]),
      source_inflammatory=mean(z$inflammatory_ligand[src]),source_matrix=mean(z$matrix_mediator[src]),
      source_vascular=mean(z$vascular_mediator[src]),source_mechanical=mean(z$mechanical_stress[src]),
      source_tgfbr2=mean(z$tgfbr2_expression[src]),near_receptor=mean(z$tgfb_receptor[near]),
      target_near=mean(z$tgfb_target[near]),target_far=mean(z$tgfb_target[far]),
      target_effect=mean(z$tgfb_target[near])-mean(z$tgfb_target[far]),
      tcell_effect=mean(z$t_cell[near])-mean(z$t_cell[far]),
      near_receivers=length(near),far_receivers=length(far)))
    geometry[[length(geometry)+1L]] <- cbind(key,data.table(
      near_far_effect=mean(z$primary_stromal[near])-mean(z$primary_stromal[far]),
      log_source_size=log1p(length(src)),
      mask_edge_distance_um=mean(nn2(xy[edge,,drop=FALSE],xy[src,,drop=FALSE],k=1)$nn.dists),
      near_normal_fraction=mean(z$kmeans[near]=="normal"),near_spots=length(near),
      G=as.integer(lesions$lesion_class[j]=="Tgfbr2")))
    for(k in seq_len(nrow(grid))) {
      pr <- grid$program[k];nr <- d<=grid$near_um[k]
      fr <- d>grid$far_start_um[k]&d<=grid$far_end_um[k]
      stopifnot(sum(nr)>0L,sum(fr)>0L)
      i <- i+1L
      pieces[[i]] <- cbind(data.table(section_id=sec,mouse=mouse[[sec]],lesion=lesions$lesion[j],
        lesion_class=lesions$lesion_class[j],source_spots=length(src)),grid[k],
        data.table(near_far_effect=mean(z[[pr]][receiver][nr])-mean(z[[pr]][receiver][fr]),
          near_receivers=sum(nr),far_receivers=sum(fr)))
    }
  }
}
endpoints <- rbindlist(pieces)
fwrite(endpoints,endpoint_path,sep="\t",compress="gzip")
fixed <- endpoints[program%chin%c("primary_stromal","t_cell")&near_um==200&far_start_um==400&far_end_um==800,
  .(section_id,lesion,lesion_class,program=ifelse(program=="primary_stromal","TGF-beta/stromal","T-cell"),
    source_spots,near_far_effect,near_receivers,far_receivers,mouse)]
out(fixed,"focal_endpoints.tsv")
out(rbindlist(features),"focal_mechanism_features.tsv")
out(rbindlist(geometry),"focal_geometry_covariates.tsv")
stopifnot(nrow(endpoints)==1120L,uniqueN(endpoints$lesion)==14L)
# The original primary near/far values must be exactly recovered up to rounding.
if(file.exists(file.path(here,"input/focal_endpoints.tsv"))) {
 old <- fread(file.path(here,"input/focal_endpoints.tsv"))
 for(pr in c("primary_stromal","t_cell")) {
  a <- endpoints[program==pr & near_um==200 & far_start_um==400 & far_end_um==800]
  b <- old[program==if(pr=="primary_stromal")"TGF-beta/stromal" else "T-cell"]
  stopifnot(max(abs(a$near_far_effect-b$near_far_effect[match(a$lesion,b$lesion)]))<1e-12)
 }
}
byvars <- c("program","near_um","far_start_um","far_end_um")
section_effect <- endpoints[,
  .(effect=mean(near_far_effect[lesion_class=="Tgfbr2"])-mean(near_far_effect[lesion_class=="KP control"])),
  by=c(byvars,"section_id","mouse")]
animal_effect <- section_effect[,.(effect=mean(effect),informative_sections=.N),by=c(byvars,"mouse")]
animal_summary <- animal_effect[,.(animal_equal_effect=mean(effect),mouse4_effect=effect[mouse=="mouse4"],
  mouse7_effect=effect[mouse=="mouse7"],positive_animals=sum(effect>0),negative_animals=sum(effect<0)),by=byvars]
section_summary <- section_effect[,.(section_equal_effect=mean(effect)),by=byvars]
ans <- merge(animal_summary,section_summary,by=byvars)
signs <- as.matrix(expand.grid(mouse4=c(-1,1),mouse7=c(-1,1)))
null <- rbindlist(lapply(seq_len(nrow(ans)),function(k) {
  vals <- signs%*%c(ans$mouse4_effect[k],ans$mouse7_effect[k])/2
  cbind(ans[k,..byvars],data.table(sign_id=1:4,sign_mouse4=signs[,1],sign_mouse7=signs[,2],statistic=as.numeric(vals)))
}))
ans[,sign_flip_p_two_sided:=vapply(seq_len(.N),function(k) {
  nul <- signs%*%c(mouse4_effect[k],mouse7_effect[k])/2
  mean(abs(nul)>=abs(animal_equal_effect[k])-1e-12)
},numeric(1))]
ans[,inference:="descriptive animal contrasts; 4 signs only illustrate resolution under symmetric independent animal null"]
out(section_effect,"sensitivity_section_effects.tsv")
out(animal_effect,"sensitivity_animal_effects.tsv")
out(ans,"sensitivity_grid_animal_summary.tsv")
out(null,"sensitivity_grid_all_animal_signs.tsv")
directions <- ans[,.(specifications=.N,positive_equal_animal=sum(animal_equal_effect>0),
  both_animals_positive=sum(positive_animals==2),both_animals_negative=sum(negative_animals==2),
  positive_mouse4=sum(mouse4_effect>0),positive_mouse7=sum(mouse7_effect>0),
  animal_equal_min=min(animal_equal_effect),animal_equal_max=max(animal_equal_effect),
  mouse4_min=min(mouse4_effect),mouse4_max=max(mouse4_effect),
  mouse7_min=min(mouse7_effect),mouse7_max=max(mouse7_effect),
  sign_flip_two_sided_min=min(sign_flip_p_two_sided),sign_flip_two_sided_max=max(sign_flip_p_two_sided)),by=program]
out(directions,"sensitivity_direction_summary.tsv")
print(directions)


# Recreate all 40 lesion/programme curves, including the non-focal animal.
profiles <- list();pp <- 0L
for(sec in unique(spots$section_id)) {
 z<-spots[section_id==sec];recv<-which(z$eligible_receiver)
 lesions<-unique(z[!is.na(lesion_class),.(lesion=phenotypes,lesion_class)])
 for(j in seq_len(nrow(lesions))) {
  src<-which(z$phenotypes==lesions$lesion[j]);if(length(src)<5L)next
  d<-nn2(as.matrix(z[src,.(x_um,y_um)]),as.matrix(z[recv,.(x_um,y_um)]),k=1)$nn.dists[,1]
  bins<-cut(d,c(-Inf,100,200,400,800,Inf),labels=c("0-100","100-200","200-400","400-800",">800"))
  for(pr in c("primary_stromal","t_cell")) {
   pp<-pp+1L
   q<-data.table(section_id=sec,lesion=lesions$lesion[j],lesion_class=lesions$lesion_class[j],
    program=if(pr=="primary_stromal")"TGF-beta/stromal" else "T-cell",distance_bin=bins,score=z[[pr]][recv])
   profiles[[pp]]<-q[,.(mean_score=mean(score),se=sd(score)/sqrt(.N),receivers=.N),
    by=.(section_id,lesion,lesion_class,program,distance_bin)]
  }
 }
}
out(rbindlist(profiles),"perturb_map_response_field_profiles.tsv")

summarise_perturbmap(outdir)
invisible(ans)
}

summarise_perturbmap <- function(input_dir) {
 outdir <- input_dir
 write_tsv <- function(x,name)write.table(x,file.path(outdir,name),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
 read_tsv <- function(path)read.delim(path,check.names=FALSE,stringsAsFactors=FALSE)
sections <- c("KP_1", "KP_2", "KP_4")
mouse <- c(KP_1="mouse4", KP_2="mouse4", KP_4="mouse7")
focal <- function(x) x[x$section_id %in% sections &
  x$lesion_class %in% c("Tgfbr2", "KP control"), , drop=FALSE]
read_focal <- function(name) {
  path <- file.path(input_dir,name)
  if(!file.exists(path))stop("Missing reconstructed endpoint: ",path)
  focal(read_tsv(path))
}
e <- read_focal("focal_endpoints.tsv")
f <- read_focal("focal_mechanism_features.tsv")
geo <- read_focal("focal_geometry_covariates.tsv")
e$mouse <- unname(mouse[e$section_id]); f$mouse <- unname(mouse[f$section_id])
geo$mouse <- unname(mouse[geo$section_id])

pnames <- unique(e$program)
z <- e[e$program == pnames[1], , drop=FALSE]
stopifnot(nrow(z)==14L, sum(z$lesion_class=="Tgfbr2")==4L)
G <- as.integer(z$lesion_class=="Tgfbr2")
y <- sapply(pnames, function(pr) {
  q <- e[e$program==pr, , drop=FALSE]
  stopifnot(setequal(q$lesion,z$lesion))
  q$near_far_effect[match(z$lesion,q$lesion)]
})
section_effect <- function(v, g=G) sapply(sections, function(sec) {
  ii <- z$section_id==sec
  mean(v[ii & g==1])-mean(v[ii & g==0])
})
animal_effect <- function(d) c(mouse4=mean(d[c("KP_1", "KP_2")]), mouse7=d[["KP_4"]])
observed_section <- sapply(seq_along(pnames), function(k) section_effect(y[,k]))
colnames(observed_section) <- pnames
observed_animal <- apply(observed_section, 2, animal_effect)
section_table <- do.call(rbind, lapply(seq_along(pnames),function(k)
  data.frame(program=pnames[k], section_id=sections, mouse=unname(mouse[sections]),
    Tgfbr2_lesion_annotations=sapply(sections,function(s)sum(z$section_id==s & G==1)),
    KP_lesion_annotations=sapply(sections,function(s)sum(z$section_id==s & G==0)),
    mutant_minus_KP=observed_section[,k], row.names=NULL)))
write_tsv(section_table,"section_effects.tsv")
animal_table <- do.call(rbind,lapply(seq_along(pnames),function(k)
  data.frame(program=pnames[k],mouse=rownames(observed_animal),
    informative_sections=c(2,1),mutant_lesion_annotations=c(2,2),
    KP_lesion_annotations=c(6,4),mutant_minus_KP=observed_animal[,k],
    independent_clone_count="unresolved",row.names=NULL)))
write_tsv(animal_table,"animal_effects.tsv")

# All four whole-animal sign vectors; each vector is shared across outcomes.
signs <- as.matrix(expand.grid(mouse4=c(-1,1),mouse7=c(-1,1)))
sf <- signs %*% observed_animal / 2
sf_out <- data.frame(sign_id=seq_len(nrow(signs)),signs,check.names=FALSE)
for(k in seq_along(pnames)) sf_out[[pnames[k]]] <- sf[,k]
write_tsv(sf_out,"animal_sign_flip_all_four.tsv")
summary <- do.call(rbind,lapply(seq_along(pnames),function(k) {
  obs <- mean(observed_animal[,k]); expected <- if(pnames[k]=="T-cell") -1 else 1
  data.frame(program=pnames[k],supporting_animals=2L,
    section_equal_descriptive_effect=mean(observed_section[,k]),
    animal_equal_descriptive_effect=obs,mouse4_effect=observed_animal[1,k],
    mouse7_effect=observed_animal[2,k],
    animals_in_expected_direction=sum(expected*observed_animal[,k]>0),
    sign_vectors=4L,sign_flip_p_two_sided=sum(abs(sf[,k])>=abs(obs)-1e-12)/4,
    sign_flip_p_expected_one_sided=sum(expected*sf[,k]>=expected*obs-1e-12)/4,
    interpretation="descriptive; sign-flip probability needs independent symmetric zero-centred animal contrasts",
    row.names=NULL)
}))
write_tsv(summary,"animal_summary.tsv")
selection <- do.call(rbind,lapply(seq_along(pnames),function(k)
  data.frame(program=pnames[k],retained_mouse4_section=c("KP_1","KP_2"),
    retained_mouse7_section="KP_4",supporting_animals=2L,
    animal_equal_effect=(observed_section[1:2,k]+observed_section[3,k])/2,
    interpretation="one section per animal; descriptive sensitivity",row.names=NULL)))
write_tsv(selection,"one_section_per_animal_sensitivity.tsv")
deletion <- do.call(rbind,lapply(seq_along(pnames),function(k)
  data.frame(program=pnames[k],omitted_animal=c("mouse7","mouse4"),
    retained_animal=c("mouse4","mouse7"),supporting_animals=1L,
    observed_effect=observed_animal[,k],
    interpretation="influence only; one animal cannot estimate between-animal uncertainty",row.names=NULL)))
write_tsv(deletion,"leave_one_animal_out.tsv")

# Regression products are numerical descriptions, with no P, CI, bootstrap,
# cluster-robust SE, or mixed-model inference from just two independent animals.
f$G <- as.integer(f$lesion_class=="Tgfbr2")
med <- c(tgfb_soluble="source_ligand",contact="source_contact",
  inflammatory="source_inflammatory",matrix="source_matrix",vascular="source_vascular",
  mechanical="source_mechanical")
fit_chain <- function(dd,mediator,weighting="equal_lesion") {
  section_term <- if(length(unique(dd$section_id))>1L)"factor(section_id)" else character()
  w <- if(weighting=="equal_animal")1/as.numeric(table(dd$mouse)[dd$mouse]) else rep(1,nrow(dd))
  a <- lm(reformulate(c("G",section_term),mediator),data=dd,weights=w)
  b <- lm(reformulate(c(mediator,"G","near_receptor",section_term),"target_effect"),data=dd,weights=w)
  total <- lm(reformulate(c("G",section_term),"target_effect"),data=dd,weights=w)
  aa <- unname(coef(a)["G"]); bb <- unname(coef(b)[mediator])
  data.frame(lesion_annotations=nrow(dd),animals=length(unique(dd$mouse)),
    sections=length(unique(dd$section_id)),a=aa,b=bb,product_ab=aa*bb,
    total_target_coefficient=unname(coef(total)["G"]),
    b_residual_df=df.residual(b),b_rank=b$rank,model_columns=ncol(model.matrix(b)),
    inference="descriptive OLS/WLS coefficients only; dependence and causal identification unresolved")
}
fitsets <- list(all_animals=f,mouse4_only=f[f$mouse=="mouse4",],
  mouse7_only=f[f$mouse=="mouse7",],one_section_mouse4_KP1=f[f$section_id%in%c("KP_1","KP_4"),],
  one_section_mouse4_KP2=f[f$section_id%in%c("KP_2","KP_4"),])
chains <- do.call(rbind,lapply(names(fitsets),function(fitname)
  do.call(rbind,lapply(names(med),function(mname) {
    dd <- fitsets[[fitname]]
    do.call(rbind,lapply(if(fitname=="all_animals")c("equal_lesion","equal_animal") else "equal_lesion",
      function(weighting)cbind(data.frame(subset=fitname,mediator=mname,weighting=weighting),
        fit_chain(dd,med[[mname]],weighting))))
  }))))
write_tsv(chains,"mediation_coefficients_descriptive.tsv")

# Covariate-adjusted coefficients: equal-animal WLS and animal deletion are
# descriptive. With mouse7 only, the full model has six data and five
# parameters; one residual degree of freedom is not animal replication.
geo$G <- as.integer(geo$lesion_class=="Tgfbr2")
covars <- list(section_only=character(),source_size="log_source_size",
  mask_edge="mask_edge_distance_um",composition="near_normal_fraction",
  size_edge=c("log_source_size","mask_edge_distance_um"),
  size_edge_composition=c("log_source_size","mask_edge_distance_um","near_normal_fraction"))
geometry <- do.call(rbind,lapply(c("all_animals","mouse4_only","mouse7_only"),function(subset) {
  dd <- if(subset=="all_animals")geo else geo[geo$mouse==if(subset=="mouse4_only")"mouse4" else "mouse7",]
  do.call(rbind,lapply(names(covars),function(model) {
    section_term <- if(length(unique(dd$section_id))>1L)"factor(section_id)" else character()
    do.call(rbind,lapply(if(subset=="all_animals")c("equal_lesion","equal_animal") else "equal_lesion",
      function(weighting) {
        w <- if(weighting=="equal_animal")1/as.numeric(table(dd$mouse)[dd$mouse]) else rep(1,nrow(dd))
        fit <- lm(reformulate(c("G",section_term,covars[[model]]),"near_far_effect"),data=dd,weights=w)
        data.frame(subset=subset,model=model,weighting=weighting,
          animals=length(unique(dd$mouse)),lesion_annotations=nrow(dd),
          coefficient=unname(coef(fit)["G"]),rank=fit$rank,
          model_columns=ncol(model.matrix(fit)),residual_df=df.residual(fit),
          inference="descriptive coefficient only; no animal-valid P or CI")
      }))
  }))
}))
write_tsv(geometry,"geometry_coefficients_descriptive.tsv")

write_tsv(data.frame(
  proposal=c("independent section label allocations","equal-animal statistic with same 240 labels",
    "pool labels within mouse: choose(8,2)*choose(6,2)=420","whole-animal sign flipping",
    "lesion-stratified bootstrap","animal bootstrap with n=2","random-animal-intercept or cluster SE with n=2"),
  disposition=c("historical reference only","descriptive statistic; reference probabilities not dependence robust",
    "not adopted: assumes observational lesions exchangeable and independent of clone duplication/section structure",
    "report all four signs only as resolution illustration under explicit symmetry assumptions",
    "withdraw CI/sign probability as biological uncertainty; unknown clone dependence not reproduced",
    "not adopted as reliable CI: only three distinct resampled pairs and two independent units",
    "not adopted as reliable population inference: too few clusters; no evidence recovered by model choice")),
  "method_disposition.tsv")
print(summary,row.names=FALSE)
print(chains[chains$mediator=="tgfb_soluble",],row.names=FALSE)

invisible(summary)
}
