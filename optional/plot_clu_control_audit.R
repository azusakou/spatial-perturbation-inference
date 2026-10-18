#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
a<-commandArgs(TRUE);if(length(a)!=2L)stop('Usage: Rscript optional/plot_clu_control_audit.R MINIMAL_PACKAGE_ROOT OUT')
package_root<-normalizePath(a[1]);root<-normalizePath(a[2]);out<-root
stopifnot(file.exists(file.path(package_root,'data','clu_cells.rds')))
x<-readRDS(file.path(out,'clu_control_audit.rds'));d<-copy(x$summary);m<-copy(x$mixed_receiver)
if(nrow(d)!=30L)stop('Plot requires full five-background/two-reference output, not --full-only')
chips<-c('C02943C3','B03018A2','A03599E2')
bg<-c('raw','identity','identity_library_smooth','global_linear','full_bam')
bglabels<-c('Raw','Identity','Identity +\nlib smooth','Global\nlinear','Full GAM')
d[,chip_id:=factor(chip_id,levels=chips)];d[,background:=factor(background,levels=bg,labels=bglabels)]
d[,reference:=factor(reference,levels=c('count','cluster'),labels=c('Count','Cluster'))]
th<-theme_classic(base_size=7.5,base_family='Helvetica')+theme(strip.background=element_blank(),
 strip.text=element_text(face='bold',size=7.5),legend.title=element_blank(),
 axis.title=element_text(size=7.5),axis.text=element_text(size=6.5),
 axis.line=element_line(linewidth=.3),axis.ticks=element_line(linewidth=.3),
 plot.title=element_text(size=8,face='bold'),legend.text=element_text(size=7),plot.margin=margin(4,6,4,4))
p1<-ggplot(d,aes(background,excess,colour=reference,shape=reference))+
 geom_hline(yintercept=0,linewidth=.3,colour='#777777')+
 geom_point(position=position_dodge(width=.35),size=1.8,stroke=.55)+
 facet_wrap(~chip_id,nrow=1)+scale_colour_manual(values=c(Count='#436EAA',Cluster='#D9822B'))+
 scale_shape_manual(values=c(Count=16,Cluster=17))+
 labs(title='a  Observed minus source-reference mean',x=NULL,y='Clu score contrast excess')+th+
 theme(legend.position='top')
d[,row:=factor(paste(chip_id,reference,sep='  |  '),
 levels=rev(as.vector(t(outer(chips,c('Count','Cluster'),paste,sep='  |  ')))))]
p2<-ggplot(d,aes(background,row,fill=tail_fraction))+
 geom_tile(colour='white',linewidth=.6)+
 geom_text(aes(label=sprintf('%.2f',tail_fraction)),size=2.5,colour='black')+
 scale_fill_gradient(low='#DFE9F3',high='#FFFFFF',limits=c(0,1),guide='none')+
 labs(title='b  Conditional common-centre tail fractions; 99 configurations each',x=NULL,y=NULL)+th+
 theme(axis.line=element_blank(),axis.ticks=element_blank())
ml<-melt(m,id.vars=c('chip_id','window','receiver_cells'),measure.vars=c(
 'fraction_other_single_within100','fraction_multiple_within100','fraction_either_within100'),
 variable.name='exposure',value.name='fraction')
ml[,chip_id:=factor(chip_id,levels=chips)];ml[,window:=factor(window,levels=c('near','far'),labels=c('Clu near','Clu far'))]
ml[,exposure:=factor(exposure,levels=c('fraction_other_single_within100','fraction_multiple_within100','fraction_either_within100'),
 labels=c('Other single guide','Multiple guides','Either detected class'))]
p3<-ggplot(ml,aes(window,fraction,colour=exposure,shape=exposure))+
 geom_point(position=position_dodge(width=.45),size=1.8,stroke=.55)+facet_wrap(~chip_id,nrow=1)+
 scale_colour_manual(values=c('#436EAA','#A83B32','#3A8E8C'))+scale_shape_manual(values=c(16,17,15))+
 scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1),labels=c('0','50','100'))+
 labs(title='c  Receiver proximity to other detected guide sources',x=NULL,y='Within 100 µm (%)')+th+
 theme(legend.position='bottom')
p<-p1/p2/p3+plot_layout(heights=c(1.05,1.05,1))
stem<-file.path(root,'clu_real_control_and_mixed_exposure')
pdf(paste0(stem,'.pdf'),width=183/25.4,height=190/25.4,family='Helvetica',useDingbats=FALSE);print(p);dev.off()
svglite::svglite(paste0(stem,'.svg'),width=183/25.4,height=190/25.4);print(p);dev.off()
ragg::agg_tiff(paste0(stem,'.tiff'),width=183,height=190,units='mm',res=600,compression='lzw');print(p);dev.off()
ragg::agg_png(paste0(stem,'.png'),width=183,height=190,units='mm',res=180);print(p);dev.off()
writeLines(c('Real-map controls and mixed detected-guide exposure.',
 'a, All 30 adjusted near–far score contrasts minus their respective source-reference mean, from three complete chips, five fixed backgrounds and two 99-configuration source banks. Count and cluster banks are each held fixed across background models; no point is an independent biological replication. Global linear is an alternative branch, not a sequential addition to the preceding smooth model. b, Two-sided common-centre conditional tail fractions with pseudo-source sites excluded from that configuration\'s receiver contrast. No cross-panel multiplicity or confirmatory comparison is claimed. c, Fraction of observed guide-undetected receivers near (<=100 µm) and far (>400–800 µm) from Clu sources that are also within 100 µm of another detected single guide, a multiple-guide cell or either. Other single guides include mSafe and exclude Clu. The multiple-guide and single-guide proximity sets can overlap; their fractions are not additive. Detected geometric proximity does not establish a non-cell-autonomous effect or identify which guide produced a molecular score. These diagnostics neither validate source exchangeability nor resolve imperfect guide detection or programme transfer.'),file.path(root,'FIGURE_CAPTION.txt'))
cat('Figure exports complete\n')
