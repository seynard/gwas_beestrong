#module load system/R-4.1.1_gcc-9.3.0
#R
library(data.table)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(VennDiagram)
library(genetics)
library(biomaRt)
library(cowplot)
library(ggforce)
library(igraph)
library(patchwork)

dirin<-'data'
dirout<-'results'
type<-'egs'
pheno<-'logit_taux_reop_inf'
pop<-c('Mellifera','Ligustica_nus','hybrid')
version<-2

#load annotation
annot<-fread(paste0(dirin,'/proteins_48_403979.csv'))
annot$chr[annot[,1]=='linkage group LG1']<-1
annot$chr[annot[,1]=='linkage group LG2']<-2
annot$chr[annot[,1]=='linkage group LG3']<-3
annot$chr[annot[,1]=='linkage group LG4']<-4
annot$chr[annot[,1]=='linkage group LG5']<-5
annot$chr[annot[,1]=='linkage group LG6']<-6
annot$chr[annot[,1]=='linkage group LG7']<-7
annot$chr[annot[,1]=='linkage group LG8']<-8
annot$chr[annot[,1]=='linkage group LG9']<-9
annot$chr[annot[,1]=='linkage group LG10']<-10
annot$chr[annot[,1]=='linkage group LG11']<-11
annot$chr[annot[,1]=='linkage group LG12']<-12
annot$chr[annot[,1]=='linkage group LG13']<-13
annot$chr[annot[,1]=='linkage group LG14']<-14
annot$chr[annot[,1]=='linkage group LG15']<-15
annot$chr[annot[,1]=='linkage group LG16']<-16
annot<-unique(annot)
annot$start[annot$Strand=='+']<-annot$Start[annot$Strand=='+']
annot$stop[annot$Strand=='+']<-annot$Stop[annot$Strand=='+']
annot$start[annot$Strand=='-']<-annot$Stop[annot$Strand=='-']
annot$stop[annot$Strand=='-']<-annot$Start[annot$Strand=='-']

#analyse Mantra and Mash results
df_mantra<-fread(paste0(dirout,'/run_mantra_',pheno,'_',type,'_',version,'/mantra.out'),data.table=F,header=F)
colnames(df_mantra)<-c('rs','chr','ps','all_effect','all_other','nb_study','log10BF','post_proba_hetero','n','direction') 
cor_choice<-'simple_corr'
df_mash<-fread(paste0(dirout,'/result_mash_',pheno,'_',type,'_',version,'_',cor_choice,'.txt'),data.table=F)
df<-merge(df_mantra,df_mash,by=c('rs','chr','ps'),all=T)
df$pheno<-pheno
df$type<-type
colnames(df)<-c("rs","chr","ps","all_effect","all_other","nb_study","log10BF_mantra","post_proba_hetero_mantra","n","direction","log10BF_mash","n_sign","loglik_mash","pheno","type")
dat<-df
dat$pheno[dat$pheno=='logit_taux_reop_inf']<-'recap_inf'
chrsize<-dat%>%group_by(chr)%>%summarise(m=max(ps))
thresh_mantra<-5
thresh_mash<-1

datx<-subset(dat,(dat$log10BF_mash>=thresh_mash|dat$log10BF_mantra>=thresh_mantra) & dat$chr==15 & dat$ps<3000000)
gene<-annot[annot$chr==unique(datx$chr),]
size=110000
ps_min<-min(datx$ps)-size
ps_max<-max(datx$ps)+size
annotx<-gene[gene$Start<=ps_max & gene$Stop>=ps_min,]
annotx$stop[annotx$stop>ps_max & annotx$Strand=='+']<-ps_max
annotx$start[annotx$start<ps_min & annotx$Strand=='+']<-ps_min
annotx$stop[annotx$stop<ps_min & annotx$Strand=='-']<-ps_min
annotx$start[annotx$start>ps_max & annotx$Strand=='-']<-ps_max
annotx<-unique(annotx[,c('Locus','start','stop')])

pheno_ini<-'logit_taux_reop_inf'
mel<-fread(paste0(dirout,'/summary_gwas_',pop[1],'_',pheno_ini,'_freq_egs.txt.bz2'),data.table=F)
mel<-subset(mel,mel$chr==unique(datx$chr))
mel$sp<-'Mellifera'
lig<-fread(paste0(dirout,'/summary_gwas_',pop[2],'_',pheno_ini,'_freq_egs.txt.bz2'),data.table=F)
lig<-subset(lig,lig$chr==unique(datx$chr))
lig$sp<-'Ligustica & Carnica'
hyb<-fread(paste0(dirout,'/summary_gwas_',pop[3],'_',pheno_ini,'_freq_egs.txt.bz2'),data.table=F)
hyb<-subset(hyb,hyb$chr==unique(datx$chr))
hyb$sp<-'Hybrids'
df_plot<-do.call(rbind,list(mel,lig,hyb))
df_plot$sp<-factor(df_plot$sp,levels=c('Ligustica & Carnica','Mellifera','Hybrids'))
LD1<-fread(paste0(dirin,'/ld/ld_snp',datx$rs[1],'.txt'),data.table=)
colnames(LD1)<-c('chr1','bp1','rs1','chr2','bp2','rs2','r2','sp')
LD1$rs_name<-datx$rs[1]
LD1$rs[LD1$rs1==datx$rs[1]]<-LD1$rs2[LD1$rs1==datx$rs[1]]
LD1$rs[LD1$rs2==datx$rs[1]]<-LD1$rs1[LD1$rs2==datx$rs[1]]
LD1[,c('chr1','chr2','bp1','bp2','rs1','rs2')]<-NULL
LD1$sp[LD1$sp=='hybrid']<-'Hybrids'
LD1$sp[LD1$sp=='Ligustica_Carnica']<-'Ligustica & Carnica'
LD1<-merge(df_plot,LD1,by=c('rs','sp'),all=T)
LD1$rs_name[is.na(LD1$rs_name)]<-datx$rs[1]
LD1$r2[is.na(LD1$r2)]<-0
LD2<-fread(paste0(dirin,'/ld/ld_snp',datx$rs[2],'.txt'),data.table=)
colnames(LD2)<-c('chr1','bp1','rs1','chr2','bp2','rs2','r2','sp')
LD2$rs_name<-datx$rs[2]
LD2$rs[LD2$rs1==datx$rs[2]]<-LD2$rs2[LD2$rs1==datx$rs[2]]
LD2$rs[LD2$rs2==datx$rs[2]]<-LD2$rs1[LD2$rs2==datx$rs[2]]
LD2[,c('chr1','chr2','bp1','bp2','rs1','rs2')]<-NULL
LD2$sp[LD2$sp=='hybrid']<-'Hybrids'
LD2$sp[LD2$sp=='Ligustica_Carnica']<-'Ligustica & Carnica'
LD2<-merge(df_plot,LD2,by=c('rs','sp'),all=T)
LD2$rs_name[is.na(LD2$rs_name)]<-datx$rs[2]
LD2$r2[is.na(LD2$r2)]<-0
LD3<-fread(paste0(dirin,'/ld/ld_snp',datx$rs[3],'.txt'),data.table=)
colnames(LD3)<-c('chr1','bp1','rs1','chr2','bp2','rs2','r2','sp')
LD3$rs_name<-datx$rs[3]
LD3$rs[LD3$rs1==datx$rs[3]]<-LD3$rs2[LD3$rs1==datx$rs[3]]
LD3$rs[LD3$rs2==datx$rs[3]]<-LD3$rs1[LD3$rs2==datx$rs[3]]
LD3[,c('chr1','chr2','bp1','bp2','rs1','rs2')]<-NULL
LD3$sp[LD3$sp=='hybrid']<-'Hybrids'
LD3$sp[LD3$sp=='Ligustica_Carnica']<-'Ligustica & Carnica'
LD3<-merge(df_plot,LD3,by=c('rs','sp'),all=T)
LD3$rs_name[is.na(LD3$rs_name)]<-datx$rs[3]
LD3$r2[is.na(LD3$r2)]<-0
df_plot<-do.call(rbind,list(LD1,LD2,LD3))
df_plot$r2[df_plot$rs==df_plot$rs_name]<-1

sp_names<-list(
  'Ligustica & Carnica'="Ligustica \n& Carnica",
  'Mellifera'="Mellifera",
  'Hybrids'="Hybrids")
sp_labeller<-function(variable,value){return(sp_names[value])}

size<-110000
g01<-ggplot()+
	geom_rect(aes(xmin=datx$ps[1]-size,xmax=datx$ps[1]+size,ymin=min(dat$log10BF_mantra[dat$chr==datx$chr]),ymax=max(dat$log10BF_mantra[dat$chr==datx$chr])),fill='red',alpha=0.3)+
	geom_point(aes(x=ps,y=log10BF_mantra),size=0.5,alpha=0.2,dat[dat$chr==datx$chr,])+
	geom_point(aes(x=ps,y=log10BF_mantra),size=3,col='red',datx[1,])+
	geom_point(aes(x=ps,y=log10BF_mantra),size=0.5,col='red',datx[c(2,3),])+
	xlab('position in bp')+
	ylab('log10 BF Mantra')+
	ggtitle(datx$rs[1])+
	theme_bw()+
	theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_text(size=20),axis.text.y=element_text(size=12))+
	theme(plot.title=element_text(size=20))
g02<-ggplot()+
	geom_rect(aes(xmin=datx$ps[2]-size,xmax=datx$ps[2]+size,ymin=min(dat$log10BF_mantra[dat$chr==datx$chr]),ymax=max(dat$log10BF_mantra[dat$chr==datx$chr])),fill='red',alpha=0.3)+
	geom_point(aes(x=ps,y=log10BF_mantra),size=0.5,alpha=0.2,dat[dat$chr==datx$chr,])+
	geom_point(aes(x=ps,y=log10BF_mantra),size=3,col='red',datx[2,])+
	geom_point(aes(x=ps,y=log10BF_mantra),size=0.5,col='red',datx[c(1,3),])+
	xlab('position in bp')+
	ylab('log10 BF Mantra')+
	ggtitle(datx$rs[2])+
	theme_bw()+
	theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())+
	theme(plot.title=element_text(size=20))
g03<-ggplot()+
	geom_rect(aes(xmin=datx$ps[3]-size,xmax=datx$ps[3]+size,ymin=min(dat$log10BF_mantra[dat$chr==datx$chr]),ymax=max(dat$log10BF_mantra[dat$chr==datx$chr])),fill='red',alpha=0.3)+
	geom_point(aes(x=ps,y=log10BF_mantra),size=0.5,alpha=0.2,dat[dat$chr==datx$chr,])+
	geom_point(aes(x=ps,y=log10BF_mantra),size=3,col='red',datx[3,])+
	geom_point(aes(x=ps,y=log10BF_mantra),size=0.5,col='red',datx[c(1,2),])+
	xlab('position in bp')+
	ggtitle(datx$rs[3])+
	ylab('log10 BF Mantra')+
	theme_bw()+
	theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())+
	theme(plot.title=element_text(size=20))
g11<-ggplot()+
	geom_rect(aes(xmin=datx$ps[1]-size,xmax=datx$ps[1]+size,ymin=min(dat$log10BF_mash[dat$chr==datx$chr]),ymax=max(dat$log10BF_mash[dat$chr==datx$chr])),fill='red',alpha=0.3)+
	geom_point(aes(x=ps,y=log10BF_mash),size=0.5,alpha=0.2,dat[dat$chr==datx$chr,])+
	geom_point(aes(x=ps,y=log10BF_mash),size=3,col='red',datx[1,])+
	geom_point(aes(x=ps,y=log10BF_mash),size=0.5,col='red',datx[c(2,3),])+
	xlab('position in bp')+
	ylab('log10 BF Mash')+
	theme_bw()+
	theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_text(size=20),axis.text.y=element_text(size=12))
g12<-ggplot()+
	geom_rect(aes(xmin=datx$ps[2]-size,xmax=datx$ps[2]+size,ymin=min(dat$log10BF_mash[dat$chr==datx$chr]),ymax=max(dat$log10BF_mash[dat$chr==datx$chr])),fill='red',alpha=0.3)+
	geom_point(aes(x=ps,y=log10BF_mash),size=0.5,alpha=0.2,dat[dat$chr==datx$chr,])+
	geom_point(aes(x=ps,y=log10BF_mash),size=3,col='red',datx[2,])+
	geom_point(aes(x=ps,y=log10BF_mash),size=0.5,col='red',datx[c(1,3),])+
	xlab('position in bp')+
	ylab('log10 BF Mash')+
	theme_bw()+
	theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())
g13<-ggplot()+
	geom_rect(aes(xmin=datx$ps[3]-size,xmax=datx$ps[3]+size,ymin=min(dat$log10BF_mash[dat$chr==datx$chr]),ymax=max(dat$log10BF_mash[dat$chr==datx$chr])),fill='red',alpha=0.3)+
	geom_point(aes(x=ps,y=log10BF_mash),size=0.5,alpha=0.2,dat[dat$chr==datx$chr,])+
	geom_point(aes(x=ps,y=log10BF_mash),size=3,col='red',datx[3,])+
	geom_point(aes(x=ps,y=log10BF_mash),size=0.5,col='red',datx[c(1,2),])+
	xlab('position in bp')+
	ylab('log10 BF Mash')+
	theme_bw()+
	theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())
g21<-ggplot()+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=0.5,alpha=0.2,data=df_plot[df_plot$r2==0 & df_plot$rs_name==datx$rs[1],])+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=0.5,alpha=0.2,data=df_plot[df_plot$r2!=0 & df_plot$rs_name==datx$rs[1],])+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=2,data=df_plot[df_plot$rs%in%datx$rs & df_plot$rs_name==datx$rs[1],])+
	geom_point(aes(x=ps,y=-log10(p_wald)),size=2,data=df_plot[df_plot$rs%in%datx$rs & df_plot$rs_name==datx$rs[1],],pch=1,col='black')+
	xlim(datx$ps[1]-size,datx$ps[1]+size)+
	scale_alpha(guide = 'none')+
	scale_color_gradient2(midpoint=0.5,low="blue",mid="yellow",high="red")+
	xlab('position in bp')+
	ylab('-log10(p_val)')+
	facet_grid(sp~.,labeller=sp_labeller)+
	theme_bw()+
	theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_text(size=20),axis.text.y=element_text(size=12))+
	theme(legend.position='none')+
	theme(strip.text.y=element_blank())
g22<-ggplot()+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=0.5,alpha=0.2,data=df_plot[df_plot$r2==0 & df_plot$rs_name==datx$rs[2],])+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=0.5,alpha=0.2,data=df_plot[df_plot$r2!=0 & df_plot$rs_name==datx$rs[2],])+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=2,data=df_plot[df_plot$rs%in%datx$rs & df_plot$rs_name==datx$rs[2],])+
	geom_point(aes(x=ps,y=-log10(p_wald)),size=2,data=df_plot[df_plot$rs%in%datx$rs & df_plot$rs_name==datx$rs[2],],pch=1,col='black')+
	xlim(datx$ps[2]-size,datx$ps[2]+size)+
	scale_alpha(guide = 'none')+
	scale_color_gradient2(midpoint=0.5,low="blue",mid="yellow",high="red")+
	xlab('position in bp')+
	ylab('-log10(p_val)')+
	facet_grid(sp~.,labeller=sp_labeller)+
	theme_bw()+
	theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())+
	theme(legend.position='none')+
	theme(strip.text.y=element_blank())
g23<-ggplot()+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=0.5,alpha=0.2,data=df_plot[df_plot$r2==0 & df_plot$rs_name==datx$rs[3],])+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=0.5,alpha=0.2,data=df_plot[df_plot$r2!=0 & df_plot$rs_name==datx$rs[3],])+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=2,data=df_plot[df_plot$rs%in%datx$rs & df_plot$rs_name==datx$rs[3],])+
	geom_point(aes(x=ps,y=-log10(p_wald)),size=2,data=df_plot[df_plot$rs%in%datx$rs & df_plot$rs_name==datx$rs[3],],pch=1,col='black')+
	xlim(datx$ps[3]-size,datx$ps[3]+size)+
	scale_alpha(guide = 'none')+
	scale_color_gradient2(midpoint=0.5,low="blue",mid="yellow",high="red")+
	xlab('position in bp')+
	ylab('-log10(p_val)')+
	facet_grid(sp~.,labeller=sp_labeller)+
	theme_bw()+
	theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())+
	theme(legend.title=element_text(size=15),legend.text=element_text(size=12))+
	theme(strip.text.y=element_text(size=15))
g41<-ggplot()+
	geom_segment(aes(x=start/1000000,xend=stop/1000000,y=Locus,yend=Locus),arrow=arrow(length=unit(0.2,'cm')),col='black',size=1,data=annotx)+
	geom_segment(aes(x=start/1000000,xend=stop/1000000,y=Locus,yend=Locus),arrow=arrow(length=unit(0.2,'cm')),col='red',size=2,data=annotx[annotx$start>=min(datx$ps) & annotx$stop<=max(datx$ps),])+
	xlab('position in Mbp')+
	theme_bw()+ylab('Genes')+xlim((datx$ps[1]-size)/1000000,(datx$ps[1]+size)/1000000)+
	theme(axis.title.y=element_text(size=20),axis.text.y=element_text(size=12),axis.title.x=element_text(size=20),axis.text.x=element_text(size=12))
g42<-ggplot()+
	geom_segment(aes(x=start/1000000,xend=stop/1000000,y=Locus,yend=Locus),arrow=arrow(length=unit(0.2,'cm')),col='black',size=1,data=annotx)+
	geom_segment(aes(x=start/1000000,xend=stop/1000000,y=Locus,yend=Locus),arrow=arrow(length=unit(0.2,'cm')),col='red',size=2,data=annotx[annotx$start>=min(datx$ps) & annotx$stop<=max(datx$ps),])+
	xlab('position in Mbp')+
	theme_bw()+ylab('Genes')+xlim((datx$ps[2]-size)/1000000,(datx$ps[2]+size)/1000000)+
	theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=20),axis.text.x=element_text(size=12))
g43<-ggplot()+
	geom_segment(aes(x=start/1000000,xend=stop/1000000,y=Locus,yend=Locus),arrow=arrow(length=unit(0.2,'cm')),col='black',size=1,data=annotx)+
	geom_segment(aes(x=start/1000000,xend=stop/1000000,y=Locus,yend=Locus),arrow=arrow(length=unit(0.2,'cm')),col='red',size=2,data=annotx[annotx$start>=min(datx$ps) & annotx$stop<=max(datx$ps),])+
	xlab('position in Mbp')+
	theme_bw()+ylab('Genes')+xlim((datx$ps[3]-size)/1000000,(datx$ps[3]+size)/1000000)+
	theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=20),axis.text.x=element_text(size=12))

mel<-fread(paste0('results/summary_gwas_',pop[1],'_',pheno,'_freq_egs.txt.bz2'),data.table=F)
mel<-subset(mel,mel$rs%in%datx$rs)
mel<-data.frame('rs'=c(mel$rs,mel$rs),'pop'='Mellifera',
	'beta'=c(mel$beta,mel$PosteriorMean),'se'=c(mel$se,mel$PosteriorSD),
	'method'=c(rep('gemma',3),rep('ashr',3)),'lfsr'=c(rep(NA,3),mel$lfsr),'lfdr'=c(rep(NA,3),mel$lfdr))
lig<-fread(paste0('results/summary_gwas_',pop[2],'_',pheno,'_freq_egs.txt.bz2'),data.table=F)
lig<-subset(lig,lig$rs%in%datx$rs)
lig<-data.frame('rs'=c(lig$rs,lig$rs),'pop'='Ligustica & Carnica',
	'beta'=c(lig$beta,lig$PosteriorMean),'se'=c(lig$se,lig$PosteriorSD),
	'method'=c(rep('gemma',3),rep('ashr',3)),'lfsr'=c(rep(NA,3),lig$lfsr),'lfdr'=c(rep(NA,3),lig$lfdr))
hyb<-fread(paste0('results/summary_gwas_',pop[3],'_',pheno,'_freq_egs.txt.bz2'),data.table=F)
hyb<-subset(hyb,hyb$rs%in%datx$rs)
hyb<-data.frame('rs'=c(hyb$rs,hyb$rs),'pop'='Hybrids',
	'beta'=c(hyb$beta,hyb$PosteriorMean),'se'=c(hyb$se,hyb$PosteriorSD),
	'method'=c(rep('gemma',3),rep('ashr',3)),'lfsr'=c(rep(NA,3),hyb$lfsr),'lfdr'=c(rep(NA,3),hyb$lfdr))
	
beta_mantra<-readLines(paste0('results/run_mantra_',pheno,'_egs_',version,'/mantra.bet.out'))
beta_mantra1<-beta_mantra[grep(paste0('\\b',datx$rs[1],'\\b'),beta_mantra)]
beta_mantra1<-unlist(strsplit(beta_mantra1,' '))
beta_mantra1<-beta_mantra1[beta_mantra1!='']
popi<-which(beta_mantra1%in%c('Mellifera','Ligustica_nus','hybrid'))
pop<-beta_mantra1[popi]
beta<-beta_mantra1[popi+1]
se<-beta_mantra1[popi+2]
df_mantra1<-data.frame('rs'=datx$rs[1],'pop'=pop,'beta'=beta,'se'=se,'method'='mantra','lfsr'=NA,'lfdr'=NA)
df_mantra1$pop[df_mantra1$pop=='Ligustica_nus']<-'Ligustica_Carnica'
beta_mantra2<-beta_mantra[grep(paste0('\\b',datx$rs[2],'\\b'),beta_mantra)]
beta_mantra2<-unlist(strsplit(beta_mantra2,' '))
beta_mantra2<-beta_mantra2[beta_mantra2!='']
popi<-which(beta_mantra2%in%c('Mellifera','Ligustica_nus','hybrid'))
pop<-beta_mantra2[popi]
beta<-beta_mantra2[popi+1]
se<-beta_mantra2[popi+2]
df_mantra2<-data.frame('rs'=datx$rs[2],'pop'=pop,'beta'=beta,'se'=se,'method'='mantra','lfsr'=NA,'lfdr'=NA)
df_mantra2$pop[df_mantra2$pop=='Ligustica_nus']<-'Ligustica_Carnica'
beta_mantra3<-beta_mantra[grep(paste0('\\b',datx$rs[3],'\\b'),beta_mantra)]
beta_mantra3<-unlist(strsplit(beta_mantra3,' '))
beta_mantra3<-beta_mantra3[beta_mantra3!='']
popi<-which(beta_mantra3%in%c('Mellifera','Ligustica_nus','hybrid'))
pop<-beta_mantra3[popi]
beta<-beta_mantra3[popi+1]
se<-beta_mantra3[popi+2]
df_mantra3<-data.frame('rs'=datx$rs[3],'pop'=pop,'beta'=beta,'se'=se,'method'='mantra','lfsr'=NA,'lfdr'=NA)
df_mantra3$pop[df_mantra3$pop=='Ligustica_nus']<-'Ligustica_Carnica'
df_mantra<-do.call(rbind,list(df_mantra1,df_mantra2,df_mantra3))

df_mash<-fread(paste0('results/summary_mash_',pheno,'_egs_',version,'_simple_corr.txt'),data.table=F)
df_mash<-subset(df_mash,df_mash$rs%in%datx$rs)
colnames(df_mash)[colnames(df_mash)=='post_mean']<-'beta'
colnames(df_mash)[colnames(df_mash)=='post_sd']<-'se'
df_mash$method<-'mash'
colnames(df_mash)<-gsub('post_','',colnames(df_mash))
colnames(df_mash)[colnames(df_mash)=='condition']<-'pop'
df_mash$pop[df_mash$pop=='Ligustica_nus']<-'Ligustica_Carnica'
effect<-do.call(rbind,list(mel,lig,hyb,df_mantra,df_mash))
effect$pheno<-pheno
effect$beta_min<-as.numeric(effect$beta)-1.96*as.numeric(effect$se)
effect$beta_max<-as.numeric(effect$beta)+1.96*as.numeric(effect$se)
effect$sign[effect$beta_min<=0 & effect$beta_max>=0]<-0 
effect$sign[effect$beta_min<=0 & effect$beta_max<=0]<-'-' 
effect$sign[effect$beta_min>=0 & effect$beta_max>=0]<-'+' 
effect$pop[effect$pop=='Ligustica_Carnica']<-'Ligustica & Carnica'
effect$pop[effect$pop=='Mellifera']<-'Mellifera'
effect$pop[effect$pop=='hybrid']<-'Hybrids'
effect$pop<-factor(effect$pop,levels=c('Ligustica & Carnica','Mellifera','Hybrids'))
effect$method<-factor(effect$method,levels=c('gemma','ashr','mantra','mash'))

	
g51<-ggplot(effect[effect$rs==datx$rs[1],])+
	geom_point(aes(x=pop,y=as.numeric(beta),col=pop),size=2)+
	geom_segment(aes(x=pop,xend=pop,y=beta_min,yend=beta_max,col=pop))+
	xlab('Group')+ylab('SNP effect')+
	geom_hline(yintercept=0,col='red',lty=3)+
	facet_grid(.~method)+
	scale_colour_manual(values=c('goldenrod2','grey34','deepskyblue2'),name='Group')+
	theme_bw()+
	theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_text(size=20),axis.text.y=element_text(size=12))+
	theme(legend.position='none')+
	theme(strip.text.x=element_text(size=15))
g52<-ggplot(effect[effect$rs==datx$rs[2],])+
	geom_point(aes(x=pop,y=as.numeric(beta),col=pop),size=2)+
	geom_segment(aes(x=pop,xend=pop,y=beta_min,yend=beta_max,col=pop))+
	xlab('Group')+ylab('SNP effect')+
	geom_hline(yintercept=0,col='red',lty=3)+
	facet_grid(.~method)+
	scale_colour_manual(values=c('goldenrod2','grey34','deepskyblue2'),name='Group')+
	theme_bw()+
	theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=12))+
	theme(legend.position='none')+
	theme(strip.text.x=element_text(size=15))
g53<-ggplot(effect[effect$rs==datx$rs[3],])+
	geom_point(aes(x=pop,y=as.numeric(beta),col=pop),size=2)+
	geom_segment(aes(x=pop,xend=pop,y=beta_min,yend=beta_max,col=pop))+
	xlab('Group')+ylab('SNP effect')+
	geom_hline(yintercept=0,col='red',lty=3)+
	facet_grid(.~method)+
	scale_colour_manual(values=c('goldenrod2','grey34','deepskyblue2'),name='Group')+
	theme_bw()+
	theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=12))+
	theme(legend.title=element_text(size=15),legend.text=element_text(size=12))+
	theme(strip.text.x=element_text(size=15))

png(paste0('plot_',pheno,'_rs_paper_',version,'.png'),width=5000,height=4200,res=300)
print(g01+g02+g03+g11+g12+g13+g21+g22+g23+g41+g42+g43+g51+g52+g53+plot_layout(nrow=5,ncol=3,height=c(1,1,1.5,1)))
dev.off()
