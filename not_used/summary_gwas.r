#module load system/R-3.6.1
library(data.table)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ashr)
library(reshape2)
library(gridExtra)
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

pop<-c('Ligustica_Carnica','Mellifera_ncorse')
pheno<-c('chargevarroaracine4','eb_smr','logit_taux_reop','logit_taux_reop_inf','logit_varroainfestation','smr_brut','varroadepthmitoracine4','varroaphoretic')
grm<-'freq'
gwas<-c('freq','egs')

#summary pve
#P<-data.frame('model'=NA,'pve'=NA,'pve_025'=NA,'pve_975'=NA,'pop'=NA,'pheno'=NA,'grm'=NA,'gwas'=NA)
#for(i in 1:length(pop)){
# for(j in 1:length(pheno)){
#  for(k in 1:length(grm)){
#   for(l in 1:length(gwas)){
#	log_lmm<-readLines(paste0('results/gemma_',pop[i],'_',pheno[j],'_',grm[k],'_lmm_',gwas[l],'_cov.log.txt'))
#	lmm_pve<-as.numeric(tail(unlist(strsplit(log_lmm[grep("pve estimate",log_lmm)],' ')),n=1))
#	lmm_se_pve<-as.numeric(tail(unlist(strsplit(log_lmm[grep("se",log_lmm)],' ')),n=1))
#	log_bslmm<-read.table(paste0('results/gemma_',pop[i],'_',pheno[j],'_',grm[k],'_',gwas[l],'_bslmm.hyp.txt'),header=T)
#	bslmm_pve<-c("PVE", mean(log_bslmm$pve),quantile(log_bslmm$pve, probs=c(0.025,0.975)))
#	params<-data.frame('model'=c('lmm','bslmm'),'pve'=c(lmm_pve,bslmm_pve[2]),'pve_025'=c(lmm_pve-1.96*lmm_se_pve,bslmm_pve[3]),'pve_975'=c(lmm_pve+1.96*lmm_se_pve,bslmm_pve[4]))
#	params$pop<-pop[i]
#	params$pheno<-pheno[j]
#	params$grm<-grm[k]
#	params$gwas<-gwas[l]
#	row.names(params)<-NULL
#	P<-rbind(P,params)
#   }
#  }		
# }
#}
#P<-subset(P,!is.na(P$pop))
#P$gwas<-factor(P$gwas,levels=c('egs','freq'))
#P$model<-factor(P$model,levels=c('lmm','bslmm'))
#P$Name<-paste0(P$pop,'_',P$model)
#P$Name<-gsub('Ligustica_Carnica','L',P$Name)
#P$Name<-gsub('Mellifera_ncorse','M',P$Name)
#P$pheno<-factor(P$pheno,levels=c('smr_brut','eb_smr','logit_varroainfestation','varroaphoretic','varroadepthmitoracine4','chargevarroaracine4','logit_taux_reop','logit_taux_reop_inf'))
#write.table(P,'results/summary_pve.txt',col.names=T,row.names=F,quote=F)

#P<-fread('results/summary_pve.txt')
#P$pheno<-factor(P$pheno,levels=c('smr_brut','eb_smr','logit_varroainfestation','varroaphoretic','varroadepthmitoracine4','chargevarroaracine4','logit_taux_reop','logit_taux_reop_inf'))
#P$Name<-factor(P$Name,levels=c('L_lmm','L_bslmm','M_lmm','M_bslmm'))
#png('results/pve.png',width=500,height=500)
#ggplot(P)+
#	geom_point(aes(x=as.numeric(pve),y=Name,col=pop))+
#	geom_segment(aes(x=as.numeric(pve_025),xend=as.numeric(pve_975),y=Name,yend=Name,col=pop))+
#	xlim(-0.2,1.2)+
#	geom_vline(aes(xintercept=0),col='red',lty=3)+
#	geom_vline(aes(xintercept=1),col='red',lty=3)+
#	theme_bw()+
#	xlab('pve')+
#	facet_grid(pheno~gwas,scale='free_y')+
#	scale_colour_manual(values=c('goldenrod','grey34'))+
#	theme(strip.text.y=element_text(angle=0))
#dev.off()

#summary gwas
#D<-data.frame('chr'=NA,'ps'=NA,'p_wald'=NA,'qvalue'=NA,'svalue'=NA,'gamma'=NA,'pop'=NA,'pheno'=NA,'grm'=NA,'gwas'=NA)
#for(i in 1:length(pop)){
# for(j in 1:length(pheno)){
#  for(k in 1:length(grm)){
#   for(l in 1:length(gwas)){
#	gwas_lmm<-fread(paste0('results/gemma_',pop[i],'_',pheno[j],'_',grm[k],'_lmm_',gwas[l],'_cov.assoc.txt'))
#	fit.ash=ash(gwas_lmm$beta,gwas_lmm$se,mixcompdist='uniform')
#	gwas_lmm=as_tibble(cbind(gwas_lmm,fit.ash$result))
#	a<-colsplit(gwas_lmm$rs,':',c('chr','ps'))
#	gwas_lmm$chr<-a$chr
#	gwas_lmm$ps<-a$ps
#	gwas_lmm<-gwas_lmm[,c('chr','ps','p_wald','qvalue','svalue')]
#	gwas_bslmm<-fread(paste0('results/gemma_',pop[i],'_',pheno[j],'_',grm[k],'_',gwas[l],'_bslmm.param.txt'))
#	a<-colsplit(gwas_bslmm$rs,':',c('chr','ps'))
#	gwas_bslmm$chr<-a$chr
#	gwas_bslmm$ps<-a$ps
#	gwas_bslmm<-gwas_bslmm[,c('chr','ps','gamma')]
#	df<-merge(gwas_lmm,gwas_bslmm,by=c('chr','ps'),all=T)
#	df$pop<-pop[i]
#	df$pheno<-pheno[j]
#	df$grm<-grm[k]
#	df$gwas<-gwas[l]
#	D<-rbind(D,df)
#   }
#  }		
# }
#}
#D<-subset(D,!is.na(D$pop))
#write.table(D,'results/summary_gwas.txt',col.names=T,row.names=F,quote=F)

D<-fread('results/summary_gwas.txt',data.table=F)
ANNOT<-data.frame('Name'=NA,'Accession'=NA,'Start'=NA,'Stop'=NA,'Strand'=NA,'GeneID'=NA,'Locus'=NA,'Locus_tag'=NA,'Protein_product'=NA,'Length'=NA,'Protein_Name'=NA,'chr'=NA,'type'=NA,'rs'=NA,'pop'=NA,'pheno'=NA,'grm'=NA,'gwas'=NA)
for(i in 1:length(pop)){
 for(j in 1:length(pheno)){
  for(k in 1:length(grm)){
   for(l in 1:length(gwas)){
print(paste0(pop[i],' ',pheno[j],' ',grm[k],' ',gwas[l]))
d<-D[D$pop==pop[i] & D$pheno==pheno[j] & D$grm==grm[k] & D$gwas==gwas[l],]
d$rs<-paste0(d$chr,':',d$ps)
qseuil<-max(d$p_wald[d$qvalue<0.1],na.rm=T)
sseuil<-max(d$p_wald[d$svalue<0.1],na.rm=T)
max_log<-round(max(-log10(d$p_wald[!is.na(d$p_wald)])))
col_chr<-rep(c('deepskyblue','orange'),8)

pval_seg<-list()
qval_seg<-list()
sval_seg<-list()
gamma_seg<-list()
for(x in 1:length(unique(d$chr))){
size<-1000000 #1Mb
a<-1
b<-size
minps<-a
maxps<-b
gam<-NA
pval<-NA
qval<-NA
sval<-NA
while(b<max(d$ps[d$chr==unique(d$chr)[x]])){
	dg<-subset(d,d$chr==unique(d$chr)[x] & d$ps>=a &d$ps<b)
	gam<-c(gam,sum(dg$gamma))
	if(qseuil>0 | sseuil>0){pval<-c(pval,nrow(dg[dg$p_wald<qseuil | dg$p_wald<sseuil,]))
	}else{pval<-c(pval,0)}
	qval<-c(qval,nrow(subset(dg,dg$qvalue<0.1)))
	sval<-c(sval,nrow(subset(dg,dg$qvalue<0.1)))
	a<-a+size
	b<-b+size
	minps<-c(minps,a)
	maxps<-c(maxps,b)
	}
pval_seg[[x]]<-data.frame('start_ps'=minps[-length(minps)],'stop_ps'=maxps[-length(maxps)],'pval_nb'=pval[-1],'chr'=unique(d$chr)[x])
qval_seg[[x]]<-data.frame('start_ps'=minps[-length(minps)],'stop_ps'=maxps[-length(maxps)],'qval_nb'=qval[-1],'chr'=unique(d$chr)[x])
sval_seg[[x]]<-data.frame('start_ps'=minps[-length(minps)],'stop_ps'=maxps[-length(maxps)],'sval_nb'=sval[-1],'chr'=unique(d$chr)[x])
gamma_seg[[x]]<-data.frame('start_ps'=minps[-length(minps)],'stop_ps'=maxps[-length(maxps)],'gamma_sum'=gam[-1],'chr'=unique(d$chr)[x])
}
pval_seg<-do.call(rbind,pval_seg)
qval_seg<-do.call(rbind,qval_seg)
sval_seg<-do.call(rbind,sval_seg)
gamma_seg<-do.call(rbind,gamma_seg)
print('regions done')

if(nrow(pval_seg[pval_seg$pval_nb>0,])>0){
	g1<-ggplot()+
		geom_rect(aes(xmin=start_ps,xmax=stop_ps,ymin=0,ymax=max_log,alpha=pval_nb/10),col='red',fill='red',data=pval_seg[pval_seg$pval_nb>0,])+
		geom_point(aes(x=ps,y=-log10(p_wald),col=as.factor(chr)),data=d)+
		geom_hline(yintercept=-log10(qseuil),col='red',lty=3)+
		geom_hline(yintercept=-log10(sseuil),col='blue',lty=3)+
		facet_grid(.~chr,scales='free',space='free')+
		theme_bw()+
		scale_colour_manual(values=col_chr)+
		xlab('Position (bp)')+ylab('-log10(p_values)')+theme(legend.position='none',axis.text.x=element_blank())+
		ggtitle(paste0('GWAS for ',pop[i],' and ',pheno[j],' with GRM on ',grm[k],' and markers as ',gwas[l]))
}else{
	g1<-ggplot()+
		geom_point(aes(x=ps,y=-log10(p_wald),col=as.factor(chr)),data=d)+
		geom_hline(yintercept=-log10(qseuil),col='red',lty=3)+
		geom_hline(yintercept=-log10(sseuil),col='blue',lty=3)+
		facet_grid(.~chr,scales='free',space='free')+
		theme_bw()+
		scale_colour_manual(values=col_chr)+
		xlab('Position (bp)')+ylab('-log10(p_values)')+theme(legend.position='none',axis.text.x=element_blank())+
		ggtitle(paste0('GWAS for ',pop[i],' and ',pheno[j],' with GRM on ',grm[k],' and markers as ',gwas[l]))}
if(nrow(qval_seg[qval_seg$qval_nb>0,])>0){
	g2<-ggplot()+
		geom_rect(aes(xmin=start_ps,xmax=stop_ps,ymin=-1,ymax=0,alpha=qval_nb/10),col='red',fill='red',data=qval_seg[qval_seg$qval_nb>0,])+
		geom_point(aes(x=ps,y=-qvalue,col=as.factor(chr)),data=d)+
		geom_hline(yintercept=-0.1,col='red',lty=3)+
		facet_grid(.~chr,scales='free',space='free')+
		theme_bw()+ylim(-max(d$qvalue),0)+
		scale_colour_manual(values=col_chr)+
		xlab('Position (bp)')+ylab('-qvalues')+theme(legend.position='none',axis.text.x=element_blank())
}else{
	g2<-ggplot()+
		geom_point(aes(x=ps,y=-qvalue,col=as.factor(chr)),data=d)+
		geom_hline(yintercept=-0.1,col='red',lty=3)+
		facet_grid(.~chr,scales='free',space='free')+
		theme_bw()+ylim(-max(d$qvalue),0)+
		scale_colour_manual(values=col_chr)+
		xlab('Position (bp)')+ylab('-qvalues')+theme(legend.position='none',axis.text.x=element_blank())}
if(nrow(sval_seg[sval_seg$sval_nb>0,])>0){
	g3<-ggplot()+
		geom_rect(aes(xmin=start_ps,xmax=stop_ps,ymin=-1,ymax=0,alpha=sval_nb/10),col='red',fill='red',data=sval_seg[sval_seg$sval_nb>0,])+
		geom_point(aes(x=ps,y=-svalue,col=as.factor(chr)),data=d)+
		geom_hline(yintercept=-0.1,col='red',lty=3)+
		facet_grid(.~chr,scales='free',space='free')+
		theme_bw()+ylim(-max(d$svalue),0)+
		scale_colour_manual(values=col_chr)+
		xlab('Position (bp)')+ylab('-svalues')+theme(legend.position='none',axis.text.x=element_blank())
}else{
	g3<-ggplot()+
		geom_point(aes(x=ps,y=-svalue,col=as.factor(chr)),data=d)+
		geom_hline(yintercept=-0.1,col='red',lty=3)+
		facet_grid(.~chr,scales='free',space='free')+
		theme_bw()+ylim(-max(d$svalue),0)+
		scale_colour_manual(values=col_chr)+
		xlab('Position (bp)')+ylab('-svalues')+theme(legend.position='none',axis.text.x=element_blank())}
g4<-ggplot()+
		geom_point(aes(x=ps,y=gamma,col=as.factor(chr)),data=d)+
		geom_segment(aes(x=start_ps,xend=stop_ps,y=gamma_sum,yend=gamma_sum),col='red',data=gamma_seg)+
		facet_grid(.~chr,scales='free',space='free')+
		theme_bw()+
		scale_colour_manual(values=col_chr)+
		xlab('Position (bp)')+ylab('gamma')+theme(legend.position='none',axis.text.x=element_blank())
png(paste0('results/gwas_',pop[i],'_',pheno[j],'_',grm[k],'_',gwas[l],'.png'),width=2000,height=500)
grid.arrange(g1,g2,g3,g4,ncol=1)
dev.off()
print('manhattan done')

sign<-subset(d,d$p_wald<qseuil & d$p_wald<sseuil & d$qvalue<0.1 & d$svalue<0.1)
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
seg_annot$pop<-pop[i]
seg_annot$pheno<-pheno[j]
seg_annot$grm<-grm[k]
seg_annot$gwas<-gwas[l]
ANNOT<-rbind(ANNOT,seg_annot)
print('annotation done')
  }
  }		
 }
}
ANNOT<-do.call(cbind,ANNOT)
ANNOT<-data.frame(ANNOT)
ANNOT<-subset(ANNOT,!is.na(ANNOT$pop))
write.table(ANNOT,'results/summary_annot.txt',col.names=T,row.names=F,quote=F)

dup<-ANNOT[c('Locus','pop','pheno')]
dup<-subset(dup,!is.na(dup$Locus))
dup<-unique(dup)
dup<-dup[duplicated(dup$Locus),]








