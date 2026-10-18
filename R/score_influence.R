# Descriptive frozen-window diagnostics only. No new animal P value is computed.
# source('score_influence.R'); run_score_influence('perturbmap_spots.tsv.gz','out')
run_score_influence <- function(input_path, out) {
suppressPackageStartupMessages({library(data.table);library(RANN)})
dir.create(out,recursive=TRUE,showWarnings=FALSE)
z <- fread(input_path)
stopifnot(nrow(z)==6363L,all(z$eligible_receiver==is.na(z$lesion_class)),
  all(c('score_reference_receiver','whole_library_size','section_id','mouse','x_um','y_um')%in%names(z)))
genes <- list(primary_stromal=c('Serpine1','Ctgf','Col1a1','Col1a2','Col3a1','Fn1','Tnc','Timp1','Tgfb1','Spp1'),matrix_mediator=c('Col1a1','Col1a2','Col3a1','Fn1','Tnc'),tgfb_target=c('Serpine1','Ctgf','Smad7','Acta2','Tagln','Timp1','Col1a1','Col1a2','Col3a1','Fn1','Tnc'),tgfb_ligand=c('Tgfb1','Tgfb2','Tgfb3'))
avail <- sub('^gene__','',grep('^gene__',names(z),value=TRUE))
requested<-genes; genes <- lapply(genes,intersect,avail)
stopifnot(length(genes$primary_stromal)==9L,length(genes$matrix_mediator)==5L,length(genes$tgfb_target)==10L)
definitions <- c(genes,list(target_excluding_matrix=setdiff(genes$tgfb_target,genes$matrix_mediator),primary_excluding_matrix=setdiff(genes$primary_stromal,genes$matrix_mediator)))
for(g in genes$primary_stromal)definitions[[paste0('drop_',g)]] <- setdiff(genes$primary_stromal,g)
for(g in genes$primary_stromal)definitions[[paste0('single_',g)]] <- g
scale_ref <- function(v,ref) {s <- mad(v[ref],constant=1.4826);if(!is.finite(s)||s==0)s<-sd(v[ref]);(v-median(v[ref]))/s}
rows <- list();composition <- list();features<-list()
for(sec in c('KP_1','KP_2','KP_4')) {
 x<-z[section_id==sec];rec<-which(x$eligible_receiver);ref<-x$score_reference_receiver;xy<-as.matrix(x[,.(x_um,y_um)])
 stopifnot(!anyNA(x[,paste0('gene__',unique(unlist(genes))),with=FALSE]),all(x$whole_library_size>=1))
 score <- sapply(definitions,function(gs)scale_ref(log1p(1e4*rowMeans(as.matrix(x[,paste0('gene__',gs),with=FALSE]))/x$whole_library_size),ref))
 logged<-log1p(1e4*as.matrix(x[,paste0('gene__',genes$primary_stromal),with=FALSE])/x$whole_library_size)
 score<-cbind(score,primary_mean_log=scale_ref(rowMeans(logged),ref),primary_mean_gene_z=scale_ref(rowMeans(apply(logged,2,scale_ref,ref=ref)),ref))
 normcounts<-1e4*as.matrix(x[,paste0('gene__',genes$primary_stromal),with=FALSE])/x$whole_library_size
 fraction<-colSums(normcounts[ref,,drop=FALSE])/sum(normcounts[ref,,drop=FALSE])
 composition[[sec]]<-data.table(section_id=sec,gene=genes$primary_stromal,receiver_count_share=fraction)
 for(les in unique(x[lesion_class%chin%c('Tgfbr2','KP control'),phenotypes])) {
  src<-which(x$phenotypes==les);d<-nn2(xy[src,,drop=FALSE],xy[rec,,drop=FALSE],k=1)$nn.dists[,1]
  near<-rec[d<=200];far<-rec[d>400&d<=800]
  v<-colMeans(score[near,,drop=FALSE])-colMeans(score[far,,drop=FALSE])
  rows[[length(rows)+1L]]<-data.table(section_id=sec,mouse=x$mouse[1],lesion=les,lesion_class=x$lesion_class[src[1]],definition=names(v),effect=v)
  features[[length(features)+1L]]<-data.table(section_id=sec,mouse=x$mouse[1],lesion=les,G=as.integer(x$lesion_class[src[1]]=='Tgfbr2'),source_ligand=mean(score[src,'tgfb_ligand']),source_matrix=mean(score[src,'matrix_mediator']),receiver_target=v['tgfb_target'],receiver_target_no_matrix=v['target_excluding_matrix'])
 }
}
end<-rbindlist(rows)
se<-end[,.(effect=mean(effect[lesion_class=='Tgfbr2'])-mean(effect[lesion_class=='KP control'])),by=.(definition,section_id,mouse)]
ae<-se[,.(effect=mean(effect)),by=.(definition,mouse)]
summary<-ae[,.(mouse4=effect[mouse=='mouse4'],mouse7=effect[mouse=='mouse7'],equal_animal=mean(effect)),by=definition]
fwrite(end,file.path(out,'score_influence_lesions.tsv'),sep='\t')
fwrite(ae,file.path(out,'score_influence_animals.tsv'),sep='\t')
fwrite(summary,file.path(out,'score_influence_summary.tsv'),sep='\t')
fwrite(rbindlist(composition),file.path(out,'primary_gene_count_share.tsv'),sep='\t')
fwrite(rbindlist(features),file.path(out,'nonoverlap_mechanism_features.tsv'),sep='\t')
overlap<-rbindlist(lapply(names(genes),function(a)rbindlist(lapply(names(genes),function(b)data.table(first=a,second=b,n_first=length(genes[[a]]),n_second=length(genes[[b]]),n_overlap=length(intersect(genes[[a]],genes[[b]])),jaccard=length(intersect(genes[[a]],genes[[b]]))/length(union(genes[[a]],genes[[b]])),shared=paste(intersect(genes[[a]],genes[[b]]),collapse=';'))))))
fwrite(overlap,file.path(out,'programme_gene_overlap.tsv'),sep='\t')
memberships<-rbindlist(lapply(names(definitions),function(n)data.table(definition=n,genes=paste(definitions[[n]],collapse=';'),n_genes=length(definitions[[n]]))))
fwrite(memberships,file.path(out,'score_definitions.tsv'),sep='\t')
checks<-data.table(check=c('fourteen_lesions','four_mutant_ten_KP','all_lesion_scores_finite','all_summary_values_finite','primary_mouse4_reproduced','primary_mouse7_reproduced','primary_equal_animal_reproduced','all_nine_LOGO_definitions_present'),pass=c(uniqueN(end$lesion)==14L,uniqueN(end[lesion_class=='Tgfbr2',lesion])==4L&&uniqueN(end[lesion_class=='KP control',lesion])==10L,all(is.finite(end$effect)),all(is.finite(as.matrix(summary[,.(mouse4,mouse7,equal_animal)]))),abs(summary[definition=='primary_stromal',mouse4]-0.9706527750)<1e-9,abs(summary[definition=='primary_stromal',mouse7]-1.6819931353)<1e-9,abs(summary[definition=='primary_stromal',equal_animal]-1.3263229552)<1e-9,summary[grepl('^drop_',definition),.N]==9L))
fwrite(checks,file.path(out,'score_influence_checks.tsv'),sep='\t')
writeLines(c('Descriptive exploratory diagnostics; the original primary score is unchanged.',
 'Input is compact raw counts, with recorded whole-gene library sizes and receiver masks.',
 'Near <=200 um; far >400 and <=800 um; all named sources excluded from eligible receivers.',
 'Lesion contrasts are averaged within section/class, then sections within mouse, then equally across mouse4/mouse7.',
 'Original score: log1p(1e4*mean(available raw gene counts)/whole_library_size), median/MAD receiver scaling with original SD fallback.',
 'All leave-one-gene-out scores use this original formula and rescale on the same reference receivers.',
 'Single-gene rows are gene-level descriptive endpoints, not independent animal replicates or multiple significant genes.',
 'mean_log and mean_gene_z are explicit alternative scoring diagnostics, not corrections or replacement primary endpoints.',
 'receiver_count_share is a sum-of-normalized-counts share; it is not a share of the nonlinear score contrast or variance.',
 'Disjoint sets remove the five measured matrix-mediator genes. They do not remove mixture confounding or establish independent pathways.',
 paste0('R ',getRversion(),'; data.table ',packageVersion('data.table'),'; RANN ',packageVersion('RANN'))),file.path(out,'SCORE_INFLUENCE_PROTOCOL.txt'))
stopifnot(all(checks$pass))
invisible(summary)
}
if(sys.nframe()==0L) {
 args<-commandArgs(trailingOnly=TRUE)
 if(length(args)!=2L)stop('Usage: Rscript score_influence.R perturbmap_spots.tsv.gz output_directory')
 print(run_score_influence(args[1],args[2]))
}
