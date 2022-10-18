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

dirout<-'results'
dirin<-'data'
type<-c('egs','freq')
pheno<-c('pc1','eb_smr','logit_taux_reop_inf')
version<-2
if(version==1){pop<-c('Mellifera','Ligustica_Carnica','hybrid')}else if (version==2){pop<-c('Mellifera','Ligustica_nus','hybrid')}

#analyse Mantra and Mash results
dat<-list()
for(j in 1:length(type)){
DF<-list()
for(i in 1:length(pheno)){
	df_mantra<-fread(paste0(dirout,'/run_mantra_',pheno[i],'_',type[j],'_',version,'/mantra.out'),data.table=F,header=F)
	colnames(df_mantra)<-c('rs','chr','ps','all_effect','all_other','nb_study','log10BF','post_proba_hetero','n','direction') 
	cor_choice<-'simple_corr'
	print(cor_choice)
	df_mash<-fread(paste0(dirout,'/result_mash_',pheno[i],'_',type[j],'_',version,'_',cor_choice,'.txt'),data.table=F)
	df<-merge(df_mantra,df_mash,by=c('rs','chr','ps'),all=T)
	df$pheno<-pheno[i]
	df$type<-type[j]
	colnames(df)<-c("rs","chr","ps","all_effect","all_other","nb_study","log10BF_mantra","post_proba_hetero_mantra","n","direction","log10BF_mash","n_sign","loglik_mash","pheno","type")
	DF[[i]]<-df
	}
dat[[j]]<-do.call(rbind,DF)
}
dat<-do.call(rbind,dat)
dat$pheno[dat$pheno=='pc1']<-'pc1_varroa_inf'
dat$pheno[dat$pheno=='eb_smr']<-'mnr'
dat$pheno[dat$pheno=='logit_taux_reop_inf']<-'recap_inf'
chrsize<-dat%>%group_by(chr)%>%summarise(m=max(ps))
dat$pheno<-factor(dat$pheno,levels=c('pc1_varroa_inf','mnr','recap_inf'))
datx<-subset(dat,dat$type=='egs' & dat$pheno%in%c('pc1_varroa_inf','mnr','recap_inf') & (dat$log10BF_mash>=1|dat$log10BF_mantra>=5))

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

#i=1
for(i in 1:nrow(datx)){
x<-datx[i,]
pheno<-x$pheno
print(pheno)
if(pheno=='recap_inf'){pheno_ini<-'logit_taux_reop_inf'
}else if (pheno=='mnr'){pheno_ini<-'eb_smr'
}else if(pheno=='pc1_varroa_inf'){pheno_ini<-'pc1'}
if(file.info(paste0(dirin,'/ld/ld_snp',x$rs,'.txt'))$size>0){
ld<-fread(paste0(dirin,'/ld/ld_snp',x$rs,'.txt'),data.table=F)
colnames(ld)<-c('chr1','bp1','rs1','chr2','bp2','rs2','r2','sp')
ld<-subset(ld,ld$rs1==x$rs | ld$rs2==x$rs)
ld$rs<-NA
ld$rs[ld$rs1==x$rs]<-ld$rs2[ld$rs1==x$rs]
ld$rs[ld$rs2==x$rs]<-ld$rs1[ld$rs2==x$rs]
ld_sub<-subset(ld,ld$r2>0.5)	
}else{ld<-data.frame('chr1'=NA,'bp1'=NA,'rs1'=NA,'chr2'=NA,'bp2'=NA,'rs2'=NA,'r2'=NA,'sp'=c('Mellifera','Ligustica_Carnica','hybrids'),'rs'=x$rs)
ld_sub<-ld}
if(nrow(ld_sub)>20){ld_sub<-ld%>%group_by(sp)%>%top_n(n=10,wt=r2)}
mel<-fread(paste0(dirout,'/summary_gwas_',pop[1],'_',pheno_ini,'_freq_egs.txt.bz2'),data.table=F)
mel<-subset(mel,mel$chr==x$chr)
ld_mel<-subset(ld,ld$sp=='Mellifera_corse')
mel<-merge(mel,ld_mel,by='rs',all=T)
mel$r2[mel$rs==x$rs]<-1
mel$r2[is.na(mel$r2)]<-0
mel$sp<-'Mellifera'
lig<-fread(paste0(dirout,'/summary_gwas_',pop[2],'_',pheno_ini,'_freq_egs.txt.bz2'),data.table=F)
lig<-subset(lig,lig$chr==x$chr)
ld_lig<-subset(ld,ld$sp=='Ligustica_Carnica')
lig<-merge(lig,ld_lig,by='rs',all=T)
lig$r2[lig$rs==x$rs]<-1
lig$r2[is.na(lig$r2)]<-0
lig$sp<-'Ligustica & Carnica'
hyb<-fread(paste0(dirout,'/summary_gwas_',pop[3],'_',pheno_ini,'_freq_egs.txt.bz2'),data.table=F)
hyb<-subset(hyb,hyb$chr==x$chr)
ld_hyb<-subset(ld,ld$sp=='hybrid')
hyb<-merge(hyb,ld_hyb,by='rs',all=T)
hyb$r2[hyb$rs==x$rs]<-1
hyb$r2[is.na(hyb$r2)]<-0
hyb$sp<-'Hybrids'
df_plot<-do.call(rbind,list(mel,lig,hyb))
df_plot$sp<-factor(df_plot$sp,levels=c('Ligustica & Carnica','Mellifera','Hybrids'))
gene<-annot[annot$chr==x$chr & annot$Start<=x$ps & annot$Stop>=x$ps,]
gene$type<-x$rs
size=100000
ps_min<-unique(x$ps)-size
ps_max<-unique(x$ps)+size
if(nrow(ld_sub)>0){
gene_sub<-list()
for(j in 1:nrow(ld_sub)){
	y<-c(ld_sub$bp1[j],ld_sub$bp2[j])
	y<-y[y!=x$ps]
	gene_sub[[j]]<-annot[annot$chr==x$chr & annot$Start<=y & annot$Stop>=y,]
	gene_sub[[j]]$rs1<-x$rs
	gene_sub[[j]]$rs2<-paste0(x$chr,':',y)
	gene_sub[[j]]$r2<-ld_sub$r2[j]
	gene_sub[[j]]$dist<-abs(ld_sub$bp1[j]-ld_sub$bp2[j])
	gene_sub[[j]]$sp<-ld_sub$sp[j]
	}
gene_sub<-do.call(rbind,gene_sub)
}
annotx<-annot[annot$chr==x$chr & annot$Start<=ps_max & annot$Stop>=ps_min,]
annotx$stop[annotx$stop>ps_max & annotx$Strand=='+']<-ps_max
annotx$start[annotx$start<ps_min & annotx$Strand=='+']<-ps_min
annotx$stop[annotx$stop<ps_min & annotx$Strand=='-']<-ps_min
annotx$start[annotx$start>ps_max & annotx$Strand=='-']<-ps_max
if(nrow(annotx)==0){
ps_min<-min(Gene[[i]]$Stop,Gene[[i]]$Start)-size
ps_max<-max(Gene[[i]]$Stop,Gene[[i]]$Start)+size
annotx<-annot[annot$chr==x$chr & annot$Start<=ps_max & annot$Stop>=ps_min,]
annotx$stop[annotx$stop>ps_max & annotx$Strand=='+']<-ps_max
annotx$start[annotx$start<ps_min & annotx$Strand=='+']<-ps_min
annotx$stop[annotx$stop<ps_min & annotx$Strand=='-']<-ps_min
annotx$start[annotx$start>ps_max & annotx$Strand=='-']<-ps_max
}

g0<-ggplot()+
	geom_rect(aes(xmin=ps_min,xmax=ps_max,ymin=min(dat$log10BF_mash[dat$chr==x$chr]),ymax=max(dat$log10BF_mash[dat$chr==x$chr])),fill='red',alpha=0.3)+
	geom_point(aes(x=ps,y=log10BF_mash),size=0.5,alpha=0.2,dat[dat$chr==x$chr,])+
	geom_point(aes(x=ps,y=log10BF_mash),size=3,col='red',x)+
	xlab('position in bp')+
	ylab('log10 BF Mash')+
	theme_bw()+theme(axis.text.x=element_blank(),axis.title.x=element_blank())
g1<-ggplot()+
	geom_rect(aes(xmin=ps_min,xmax=ps_max,ymin=min(dat$log10BF_mantra[dat$chr==x$chr]),ymax=max(dat$log10BF_mantra[dat$chr==x$chr])),fill='red',alpha=0.3)+
	geom_point(aes(x=ps,y=log10BF_mantra),size=0.5,alpha=0.2,dat[dat$chr==x$chr,])+
	geom_point(aes(x=ps,y=log10BF_mantra),size=3,col='red',x)+
	xlab('position in bp')+
	ylab('log10 BF Mantra')+
	theme_bw()
g12<-ggplot()+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=0.5,alpha=0.2,data=df_plot[df_plot$sp=='Ligustica & Carnica' & df_plot$r2==0,])+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=0.5,data=df_plot[df_plot$sp=='Ligustica & Carnica' & df_plot$r2!=0,])+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=2,data=df_plot[df_plot$sp=='Ligustica & Carnica' & df_plot$rs==x$rs,])+
	xlim(ps_min,ps_max)+
	scale_alpha(guide = 'none')+
	scale_color_gradient2(midpoint=0.5,low="blue",mid="yellow",high="red")+
	xlab('position in bp')+
	ylab('-log10(p_val)')+
	theme_bw()+theme(axis.text.x=element_blank(),axis.title.x=element_blank())
g22<-ggplot()+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=0.5,alpha=0.2,data=df_plot[df_plot$sp=='Mellifera' & df_plot$r2==0,])+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=0.5,data=df_plot[df_plot$sp=='Mellifera' & df_plot$r2!=0,])+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=2,data=df_plot[df_plot$sp=='Mellifera' & df_plot$rs==x$rs,])+
	xlim(ps_min,ps_max)+
	scale_alpha(guide = 'none')+
	scale_color_gradient2(midpoint=0.5,low="blue",mid="yellow",high="red")+
	xlab('position in bp')+
	ylab('-log10(p_val)')+
	theme_bw()+theme(axis.text.x=element_blank(),axis.title.x=element_blank())
g32<-ggplot()+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=0.5,alpha=0.2,data=df_plot[df_plot$sp=='Hybrids' & df_plot$r2==0,])+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=0.5,data=df_plot[df_plot$sp=='Hybrids' & df_plot$r2!=0,])+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=2,data=df_plot[df_plot$sp=='Hybrids' & df_plot$rs==x$rs,])+
	xlim(ps_min,ps_max)+
	scale_alpha(guide = 'none')+
	scale_color_gradient2(midpoint=0.5,low="blue",mid="yellow",high="red")+
	xlab('position in bp')+
	ylab('-log10(p_val)')+
	theme_bw()
g4<-ggplot()+
	geom_segment(aes(x=start,xend=stop,y=Locus,yend=Locus),arrow=arrow(length=unit(0.2,'cm')),col='black',size=0,data=annotx)+
	geom_segment(aes(x=start,xend=stop,y=Locus,yend=Locus),arrow=arrow(length=unit(0.2,'cm')),col='red',size=2,data=annotx[annotx$start>=x$ps & annotx$stop<=x$ps,])+
	geom_segment(aes(x=start,xend=stop,y=Locus,yend=Locus),arrow=arrow(length=unit(0.2,'cm')),col='black',size=2,data=annotx[annotx$start<x$ps | annotx$stop>x$ps,])+
	theme_bw()+ylab('Genes')+xlim(ps_min,ps_max)

#png('plot.png',width=3000,height=3000,res=300)
#g0+g1+g12+g22+g32+g4+plot_layout(ncol=1,height=c(2,2,1,1,1,2))
#dev.off()
study<-colsplit(x$direction,'',c('mel','lig','hyb'))
if(study$lig!='?' & study$mel=='?' & study$hyb=='?'){
	png(paste0('plot_',pheno,'_rs_',x$rs,'_',version,'.png'),width=2300,height=3000,res=300)
	g0+g1+g12+g4+plot_layout(ncol=1,height=c(2,2,1,2))
	dev.off()
}else if (study$lig=='?' & study$mel!='?' & study$hyb=='?'){
	png(paste0('plot_',pheno,'_rs_',x$rs,'_',version,'.png'),width=2300,height=3000,res=300)
	g0+g1+g22+g4+plot_layout(ncol=1,height=c(2,2,1,2))
	dev.off()
}else if (study$lig=='?' & study$mel=='?' & study$hyb!='?'){
	png(paste0('plot_',pheno,'_rs_',x$rs,'_',version,'.png'),width=2300,height=3000,res=300)
	g0+g1+g32+g4+plot_layout(ncol=1,height=c(2,2,1,2))
	dev.off()
}else if (study$lig!='?' & study$mel!='?' & study$hyb=='?'){
	png(paste0('plot_',pheno,'_rs_',x$rs,'_',version,'.png'),width=2600,height=3000,res=300)
	g0+g1+g12+g22+g4+plot_layout(ncol=1,height=c(2,2,1,1,2))
	dev.off()
}else if (study$lig=='?' & study$mel!='?' & study$hyb!='?'){
	png(paste0('plot_',pheno,'_rs_',x$rs,'_',version,'.png'),width=2600,height=3000,res=300)
	g0+g1+g22+g32+g4+plot_layout(ncol=1,height=c(2,2,1,1,2))
	dev.off()
}else if (study$lig!='?' & study$mel=='?' & study$hyb!='?'){
	png(paste0('plot_',pheno,'_rs_',x$rs,'_',version,'.png'),width=2600,height=3000,res=300)
	g0+g1+g21+g32+g4+plot_layout(ncol=1,height=c(2,2,1,1,2))
	dev.off()
}else if (study$lig!='?' & study$mel!='?' & study$hyb!='?'){
	png(paste0('plot_',pheno,'_rs_',x$rs,'_',version,'.png'),width=3000,height=3000,res=300)
	g0+g1+g12+g22+g32+g4+plot_layout(ncol=1,height=c(2,2,1,1,1,2))
	dev.off()
}
}

