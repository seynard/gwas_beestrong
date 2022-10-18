#module load system/R-3.6.1
library(data.table)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(ashr)
library(gridExtra)
library(ggpubr)
library(viridis)
library(RColorBrewer)
library(gplots)
library(reshape2)
library(qvalue)
args<-commandArgs(TRUE)
dir_out<-'results'
population='Mellifera_ncorse'
annot<-fread('data/proteins_48_403979.csv')
colnames(annot)<-gsub('#','',colnames(annot))
colnames(annot)<-gsub(' ','_',colnames(annot))
annot$chr<-NA
annot$chr[annot$Name=='linkage group LG1']<-1
annot$chr[annot$Name=='linkage group LG2']<-2
annot$chr[annot$Name=='linkage group LG3']<-3
annot$chr[annot$Name=='linkage group LG4']<-4
annot$chr[annot$Name=='linkage group LG5']<-5
annot$chr[annot$Name=='linkage group LG6']<-6
annot$chr[annot$Name=='linkage group LG7']<-7
annot$chr[annot$Name=='linkage group LG8']<-8
annot$chr[annot$Name=='linkage group LG9']<-9
annot$chr[annot$Name=='linkage group LG10']<-10
annot$chr[annot$Name=='linkage group LG11']<-11
annot$chr[annot$Name=='linkage group LG12']<-12
annot$chr[annot$Name=='linkage group LG13']<-13
annot$chr[annot$Name=='linkage group LG14']<-14
annot$chr[annot$Name=='linkage group LG15']<-15
annot$chr[annot$Name=='linkage group LG16']<-16

pp.ggplot = function(pval) {
    dat = as_tibble(pval)
    nmax=nrow(dat)
    dat %<>% mutate(stat = -log10(value)) %>%
        mutate(rank = dense_rank(value),
               expect = -log10(rank/(nmax+1)),
               hival = -log10(qbeta(0.975,rank,nmax-rank)),
               loval = -log10(qbeta(0.025,rank,nmax-rank)) )
    dat %<>% arrange(expect)
    x = dat$expect
    x = c(x,rev(x))
    y = c(dat$loval,rev(dat$hival))
    poly.dat=tibble(x,y)
    chisq <- qchisq(1 - pval, 1)
	lambda <- median(chisq) / qchisq(0.5, 1)
    p = dat %>% ggplot(aes(x=expect,y=stat)) +
        geom_polygon(data=poly.dat,aes(x=x,y=y),fill='blue', alpha=0.1) +
        geom_path(data=poly.dat,aes(x=x,y=y), alpha=0.1) +
        geom_point() +
        geom_abline(intercept=0,slope=1,colour='blue') +
        xlab("Expected -log(pvalue)") +
        ylab("Observed -log(pvalue)")+ 
		annotate(geom = "text",x = -Inf,y = Inf,hjust = -0.15,vjust = 1 + 0.15 * 3,label = sprintf("λ = %.2f", lambda),size = 8) +
		theme(axis.ticks = element_line(size = 0.5),panel.grid = element_blank())
    p
}

sum_group<-function(dat){
	dat<-subset(dat,dat$nb>0)
	dat$group<-NA
	x<-1
	n<-1
	while(x<=nrow(dat)){
	if(x==1){dat$group[x]<-n
	}else{
		if((dat$start_ps[x]-1)==dat$stop_ps[x-1] & dat$chr[x]==dat$chr[x-1]){dat$group[x]<-dat$group[x-1]}else{dat$group[x]<-n}
	}
	x<-x+1
	n<-n+1}
	dat<-dat%>%group_by(group)%>%summarise(start_ps=min(start_ps),stop_ps=max(stop_ps),nb=sum(nb),chr=unique(chr))
}

phenotype_id<-'logit_taux_reop_inf'
grm_id='freq'
gwas_id='egs'
input_file=paste0(dir_out,'/gemma_',population,'_',phenotype_id,'_',grm_id,'_lmm_',gwas_id,'_cov.assoc.txt')
prfx=paste0(population,'_',phenotype_id,'_',grm_id,'_',gwas_id)
gwas=fread(input_file)
a<-colsplit(gwas$rs,':',c('CHROM','POS'))
gwas$chr<-a$CHROM
gwas$ps<-a$POS
pv=as.data.frame(gwas%>%select(p_wald))
gwas.ppplot<-pp.ggplot(pv[,1])+theme_bw()
png(paste0(dir_out,'/pplot_plot_wcgalp.png'),width=300,height=300)
gwas.ppplot
dev.off()

fit.ash=ash(gwas$beta,gwas$se,mixcompdist='uniform')
gwas=as_tibble(cbind(gwas,fit.ash$result))
qseuil<-max(gwas$p_wald[gwas$qvalue<0.1],na.rm=T)
sseuil<-max(gwas$p_wald[gwas$svalue<0.1],na.rm=T)
max_log<-round(max(-log10(gwas$p_wald[!is.na(gwas$p_wald)])))
col_chr<-rep(c('deepskyblue','orange'),8)

pval_seg<-list()
qval_seg<-list()
sval_seg<-list()
for(x in 1:length(unique(gwas$chr))){
size<-1000000 #1Mb
a<-1
b<-size
minps<-a
maxps<-b
pval<-NA
qval<-NA
sval<-NA
while(b<max(gwas$ps[gwas$chr==unique(gwas$chr)[x]])){
	dg<-subset(gwas,gwas$chr==unique(gwas$chr)[x] & gwas$ps>=a &gwas$ps<b)
	if(qseuil>0 | sseuil>0){pval<-c(pval,nrow(dg[dg$p_wald<qseuil | dg$p_wald<sseuil,]))
	}else{pval<-c(pval,0)}
	qval<-c(qval,nrow(subset(dg,dg$qvalue<0.1)))
	sval<-c(sval,nrow(subset(dg,dg$qvalue<0.1)))
	a<-a+size
	b<-b+size
	minps<-c(minps,a)
	maxps<-c(maxps,b)
	}
pval_seg[[x]]<-data.frame('start_ps'=minps[-length(minps)],'stop_ps'=maxps[-length(maxps)],'nb'=pval[-1],'chr'=unique(gwas$chr)[x])
qval_seg[[x]]<-data.frame('start_ps'=minps[-length(minps)],'stop_ps'=maxps[-length(maxps)],'nb'=qval[-1],'chr'=unique(gwas$chr)[x])
sval_seg[[x]]<-data.frame('start_ps'=minps[-length(minps)],'stop_ps'=maxps[-length(maxps)],'nb'=sval[-1],'chr'=unique(gwas$chr)[x])
}
pval_seg<-do.call(rbind,pval_seg)
pval_seg<-sum_group(pval_seg)
qval_seg<-do.call(rbind,qval_seg)
qval_seg<-sum_group(qval_seg)
sval_seg<-do.call(rbind,sval_seg)
sval_seg<-sum_group(sval_seg)
print('regions done')

g1<-ggplot()+
	geom_rect(aes(xmin=start_ps,xmax=stop_ps,ymin=0,ymax=max_log,alpha=nb/20),col='red',fill='red',data=pval_seg)+
	geom_point(aes(x=ps,y=-log10(p_wald),col=as.factor(chr)),data=gwas)+
	facet_grid(.~chr,scales='free', space='free')+
	scale_colour_manual(values=col_chr,name='')+
	geom_hline(yintercept=-log10(sseuil),col='purple')+
	geom_hline(yintercept=-log10(qseuil),col='red')+
	scale_x_continuous(labels = function(x) format(x, scientific = TRUE))+
	xlab('position in (bp)')+
	ylab('-log10(p-value)')+
	theme_bw()+theme(legend.position='none',axis.text.x=element_text(angle=90))+ggtitle('pvalue')
g2<-ggplot()+
	geom_rect(aes(xmin=start_ps,xmax=stop_ps,ymin=1,ymax=0,alpha=nb/20),col='red',fill='red',data=qval_seg)+
	geom_point(aes(x=ps,y=qvalue,col=as.factor(chr)),data=gwas)+
	geom_hline(yintercept=0.1,col='red')+
	ylim(1,0)+
	facet_grid(.~chr,scales='free', space='free')+
	scale_colour_manual(values=col_chr,name='')+
	scale_x_continuous(labels = function(x) format(x, scientific = TRUE))+
	xlab('position in (bp)')+
	theme_bw()+theme(legend.position='none',axis.text.x=element_text(angle=90))+ggtitle('qvalue')
g3<-ggplot()+
	geom_rect(aes(xmin=start_ps,xmax=stop_ps,ymin=1,ymax=0,alpha=nb/20),col='red',fill='red',data=sval_seg)+
	geom_point(aes(x=ps,y=svalue,col=as.factor(chr)),data=gwas)+
	geom_hline(yintercept=0.1,col='red')+
	ylim(1,0)+
	facet_grid(.~chr,scales='free', space='free')+
	scale_colour_manual(values=col_chr,name='')+
	scale_x_continuous(labels = function(x) format(x, scientific = TRUE))+
	xlab('position in (bp)')+
	theme_bw()+theme(legend.position='none',axis.text.x=element_text(angle=90))+ggtitle('svalue')
png(paste0(dir_out,'/manhattan_plot_wcgalp.png'),width=1500,height=500)
grid.arrange(g1,g2,g3,nrow=3)
dev.off()

sign<-subset(gwas,gwas$p_wald<qseuil & gwas$p_wald<sseuil & gwas$qvalue<0.1 & gwas$svalue<0.1)
seg_annot<-data.frame('Name'=NA,'Accession'=NA,'Start'=NA,'Stop'=NA,'Strand'=NA,'GeneID'=NA,'Locus'=NA,'Locus_tag'=NA,'Protein_product'=NA,'Length'=NA,'Protein_Name'=NA,'chr'=NA,'type'=NA,'rs'=NA)
for(x in 1:nrow(sign)){
	annot_g<-subset(annot,annot$chr==sign$chr[x] & annot$Start<=sign$ps[x] & annot$Stop>=sign$ps[x])
	if(nrow(annot_g)>=1){
		annot_g$type<-'in'
		annot_g$rs<-sign$rs[x]
	}else if(nrow(annot_g)==0){
		annot_g<-subset(annot,annot$chr==sign$chr[x] & annot$Start<=sign$ps[x])
		diff<-which.min(abs(sign$ps[x]-annot_g$Stop))
		annot_g1<-annot_g[diff]
		annot_g<-subset(annot,annot$chr==sign$chr[x] & annot$Stop>=sign$ps[x])
		diff<-which.min(abs(sign$ps[x]-annot_g$Start))
		annot_g2<-annot_g[diff]
		annot_g<-rbind(annot_g1,annot_g2)
		annot_g$type<-'out'
		annot_g$rs<-sign$rs[x]
	}
seg_annot<-rbind(seg_annot,annot_g)
}
seg_annot<-unique(seg_annot)
seg_annot<-subset(seg_annot,!is.na(seg_annot$type))
print('annotation done')
write.table(seg_annot,paste0(dir_out,'/annot_wcgalp.txt'),col.names=T,row.names=F,quote=F)


