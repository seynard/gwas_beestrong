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

args<-commandArgs(TRUE)
dirin<-'data'
dirout<-'results'
type<-c('egs','freq')
pheno<-c('pc1','eb_smr','logit_taux_reop_inf')
pop<-c('Mellifera_ncorse','Ligustica_Carnica','hybrid_corse')

type<-'egs'
for(p in 1:length(pheno)){
df_mantra<-fread(paste0(dirout,'/run_mantra_',pheno[p],'_',type,'/mantra.out'),data.table=F,header=F)
colnames(df_mantra)<-c('rs','chr','ps','all_effect','all_other','nb_study','log10BF','post_proba_hetero','n','direction') 
cor_choice<-'simple_corr'
df_mash<-fread(paste0(dirout,'/result_mash_',pheno[p],'_',type,'_',cor_choice,'.txt'),data.table=F)
df<-merge(df_mantra,df_mash,by=c('rs','chr','ps'),all=T)
df$pheno<-pheno[p]
df$type<-type
colnames(df)<-c("rs","chr","ps","all_effect","all_other","nb_study","log10BF_mantra","post_proba_hetero_mantra","n","direction","log10BF_mash","n_sign","loglik_mash","pheno","type")
dat<-df
dat$pheno[dat$pheno=='pc1']<-'pc1_varroa_inf'
dat$pheno[dat$pheno=='eb_smr']<-'mnr'
dat$pheno[dat$pheno=='logit_taux_reop_inf']<-'recap_inf'
chrsize<-dat%>%group_by(chr)%>%summarise(m=max(ps))
dat$pheno<-factor(dat$pheno,levels=c('v_pho','v_mito','v_inf','v_load','pc1_varroa_inf','mnr','recap_inf'))
thresh_mantra<-5
thresh_mash<-1
chrsize<-dat%>%group_by(chr)%>%summarize(m=max(ps))

g1<-ggplot()+
	geom_point(aes(x=ps,y=log10BF_mantra),data=dat,alpha=0.2)+
	geom_hline(yintercept=5,col='red',lty=3)+
	facet_grid(.~chr,scales='free_x',space='free')+
	xlab('position in bp')+
	ylab('log10BF Mantra')+
	scale_x_continuous(breaks=seq(0,chrsize$m[1],5000000))+
	theme_bw()+
	theme(axis.text.x=element_text(angle=45,hjust=1))	
g2<-ggplot()+
	geom_point(aes(x=ps,y=log10BF_mash),data=dat,alpha=0.2)+
	geom_hline(yintercept=1,col='red',lty=3)+
	facet_grid(.~chr,scales='free_x',space='free')+
	xlab('position in bp')+
	ylab('log10BF mash')+
	scale_x_continuous(breaks=seq(0,chrsize$m[1],5000000))+
	theme_bw()+
	theme(axis.text.x=element_text(angle=45,hjust=1))	
png(paste0('meta_gwas_',pheno[p],'.png'),height=1500,width=4000,res=300)
print(plot_grid(g1,g2,ncol=1,align='v',axis='lr'))
dev.off()
}

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

type<-'egs'
for(p in 1:length(pheno)){
df_mantra<-fread(paste0(dirout,'/run_mantra_',pheno[p],'_',type,'/mantra.out'),data.table=F,header=F)
colnames(df_mantra)<-c('rs','chr','ps','all_effect','all_other','nb_study','log10BF','post_proba_hetero','n','direction') 
cor_choice<-'simple_corr'
df_mash<-fread(paste0(dirout,'/result_mash_',pheno[p],'_',type,'_',cor_choice,'.txt'),data.table=F)
df<-merge(df_mantra,df_mash,by=c('rs','chr','ps'),all=T)
df$pheno<-pheno[p]
df$type<-type
colnames(df)<-c("rs","chr","ps","all_effect","all_other","nb_study","log10BF_mantra","post_proba_hetero_mantra","n","direction","log10BF_mash","n_sign","loglik_mash","pheno","type")
dat<-df
dat$pheno[dat$pheno=='pc1']<-'pc1_varroa_inf'
dat$pheno[dat$pheno=='eb_smr']<-'mnr'
dat$pheno[dat$pheno=='logit_taux_reop_inf']<-'recap_inf'
chrsize<-dat%>%group_by(chr)%>%summarise(m=max(ps))
dat$pheno<-factor(dat$pheno,levels=c('v_pho','v_mito','v_inf','v_load','pc1_varroa_inf','mnr','recap_inf'))
thresh_mantra<-5
thresh_mash<-1
x<-subset(dat,dat$log10BF_mash>thresh_mash | dat$log10BF_mantra>thresh_mantra)
size=250000

for(i in 1:nrow(x)){
RS=x$rs[i] #11:14369154, 15:8485332
PS<-as.numeric(unlist(strsplit(RS,':'))[2])
CHR<-as.numeric(unlist(strsplit(RS,':'))[1])
ps_min<-PS-size 
ps_max<-PS+size 
datx<-subset(dat,dat$chr==CHR)

ld<-fread(paste0(dirout,'/ld_snp',RS,'_egs.txt'),data.table=F)
colnames(ld)<-c('chr1','bp1','rs1','chr2','bp2','rs2','r2','sp')
ld<-subset(ld,ld$rs1==RS | ld$rs2==RS)
ld$rs<-NA
ld$rs[ld$rs1==RS]<-ld$rs2[ld$rs1==RS]
ld$rs[ld$rs2==RS]<-ld$rs1[ld$rs2==RS]
ld_sub<-subset(ld,ld$r2>0.5)
mel<-fread(paste0(dirout,'/summary_gwas_Mellifera_ncorse_',pheno[p],'_freq_egs.txt.bz2'),data.table=F)
mel<-subset(mel,mel$chr==CHR)
ld_mel<-subset(ld,ld$sp=='Mellifera_ncorse')
mel<-merge(mel,ld_mel,by='rs',all=T)
mel$r2[mel$rs==RS]<-1
mel$r2[is.na(mel$r2)]<-0
mel$sp<-'Mellifera'
lig<-fread(paste0(dirout,'/summary_gwas_Ligustica_Carnica_',pheno[p],'_freq_egs.txt.bz2'),data.table=F)
lig<-subset(lig,lig$chr==CHR)
ld_lig<-subset(ld,ld$sp=='Ligustica_Carnica')
lig<-merge(lig,ld_lig,by='rs',all=T)
lig$r2[lig$rs==RS]<-1
lig$r2[is.na(lig$r2)]<-0
lig$sp<-'Ligustica & Carnica'
hyb<-fread(paste0(dirout,'/summary_gwas_hybrid_corse_',pheno[p],'_freq_egs.txt.bz2'),data.table=F)
hyb<-subset(hyb,hyb$chr==CHR)
ld_hyb<-subset(ld,ld$sp=='hybrid_corse')
hyb<-merge(hyb,ld_hyb,by='rs',all=T)
hyb$r2[hyb$rs==RS]<-1
hyb$r2[is.na(hyb$r2)]<-0
hyb$sp<-'Hybrids'
df_plot<-do.call(rbind,list(mel,lig,hyb))
df_plot$sp<-factor(df_plot$sp,levels=c('Ligustica & Carnica','Mellifera','Hybrids'))
annotx<-annot[annot$chr==CHR & annot$Start<=ps_max & annot$Stop>=ps_min,]
annotx$stop[annotx$stop>ps_max & annotx$Strand=='+']<-ps_max
annotx$start[annotx$start<ps_min & annotx$Strand=='+']<-ps_min
annotx$stop[annotx$stop<ps_min & annotx$Strand=='-']<-ps_min
annotx$start[annotx$start>ps_max & annotx$Strand=='-']<-ps_max

g1<-ggplot()+
	geom_rect(xmin=ps_min,xmax=ps_max,ymin=min(datx$log10BF_mash),ymax=max(datx$log10BF_mash),fill='grey',colour='black',alpha=0.01,size=0.15,data=datx)+
	geom_point(aes(x=ps,y=log10BF_mash),data=datx,alpha=0.2)+
	geom_point(aes(x=ps,y=log10BF_mash),col='red',data=datx[datx$rs==RS,])+
	geom_hline(yintercept=1,col='red',lty=3)+
	xlab('position in bp')+
	ylab('log10BF mash')+
	theme_bw()
g2<-ggplot()+
	geom_rect(xmin=ps_min,xmax=ps_max,ymin=min(datx$log10BF_mantra),ymax=max(datx$log10BF_mantra),fill='grey',colour='black',alpha=0.01,size=0.15,data=datx)+
	geom_point(aes(x=ps,y=log10BF_mantra),data=datx,alpha=0.2)+
	geom_point(aes(x=ps,y=log10BF_mantra),col='red',data=datx[datx$rs==RS,])+
	geom_hline(yintercept=5,col='red',lty=3)+
	xlab('position in bp')+
	ylab('log10BF Mantra')+
	theme_bw()
g3<-ggplot()+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=0.5,alpha=0.2,data=df_plot[df_plot$sp=='Ligustica & Carnica' & df_plot$r2==0,])+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=0.5,data=df_plot[df_plot$sp=='Ligustica & Carnica' & df_plot$r2!=0,])+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=2,data=df_plot[df_plot$sp=='Ligustica & Carnica' & df_plot$r2==1,])+
	xlim(ps_min,ps_max)+
	scale_alpha(guide = 'none')+
	scale_color_gradient2(midpoint=0.5,low="blue",mid="yellow",high="red")+
	xlab('position in bp')+
	ylab('-log10(p_value)')+
	ggtitle('Ligustica & Carnica')+
	theme_bw()
g4<-ggplot()+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=0.5,alpha=0.2,data=df_plot[df_plot$sp=='Mellifera' & df_plot$r2==0,])+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=0.5,data=df_plot[df_plot$sp=='Mellifera' & df_plot$r2!=0,])+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=2,data=df_plot[df_plot$sp=='Mellifera' & df_plot$r2==1,])+
	xlim(ps_min,ps_max)+
	scale_alpha(guide = 'none')+
	scale_color_gradient2(midpoint=0.5,low="blue",mid="yellow",high="red")+
	xlab('position in bp')+
	ylab('-log10(p_value)')+
	ggtitle('Mellifera')+
	theme_bw()
g5<-ggplot()+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=0.5,alpha=0.2,data=df_plot[df_plot$sp=='Hybrids' & df_plot$r2==0,])+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=0.5,data=df_plot[df_plot$sp=='Hybrids' & df_plot$r2!=0,])+
	geom_point(aes(x=ps,y=-log10(p_wald),col=r2),size=2,data=df_plot[df_plot$sp=='Hybrids' & df_plot$r2==1,])+
	xlim(ps_min,ps_max)+
	scale_alpha(guide = 'none')+
	scale_color_gradient2(midpoint=0.5,low="blue",mid="yellow",high="red")+
	xlab('position in bp')+
	ylab('-log10(p_value)')+
	ggtitle('Hybrids')+
	theme_bw()
g6<-ggplot()+
	geom_segment(aes(x=start,xend=stop,y=Locus,yend=Locus),arrow=arrow(length=unit(0.2,'cm')),col='red',size=2,data=annotx)+
	theme_bw()+ylab('Genes')+xlim(ps_min,ps_max)

png(paste0('gwas_',pheno[p],'_',RS,'.png'),height=4000,width=2000,res=300)
print(plot_grid(g1,g2,g3,g4,g5,g6,ncol=1,align='v',axis='lr',rel_heights=c(0.2,0.4,0.2,0.2,0.2,0.5)))
dev.off()
}
}
#bmel<-fread(paste0('results/summary_gwas_Mellifera_ncorse_',pheno[p],'_freq_egs.txt.bz2'),data.table=F)
#bmel<-subset(bmel,bmel$rs==RS)
#if(nrow(bmel)==0){bmel<-data.frame('rs'=RS,'pop'='Mellifera_ncorse','beta'=NA,'se'=NA,'method'=c('gemma','ashr'),'lfsr'=NA,'lfdr'=NA)}else{bmel<-data.frame('rs'=RS,'pop'='Mellifera_ncorse','beta'=c(bmel$beta,bmel$PosteriorMean),'se'=c(bmel$se,bmel$PosteriorSD),'method'=c('gemma','ashr'),'lfsr'=c(NA,bmel$lfsr),'lfdr'=c(NA,bmel$lfdr))}
#blig<-fread(paste0('results/summary_gwas_Ligustica_Carnica_',pheno[p],'_freq_egs.txt.bz2'),data.table=F)
#blig<-subset(blig,blig$rs==RS)
#if(nrow(blig)==0){blig<-data.frame('rs'=RS,'pop'='Ligustica_Carnica','beta'=NA,'se'=NA,'method'=c('gemma','ashr'),'lfsr'=NA,'lfdr'=NA)}else{blig<-data.frame('rs'=RS,'pop'='Ligustica_Carnica','beta'=c(blig$beta,blig$PosteriorMean),'se'=c(blig$se,blig$PosteriorSD),'method'=c('gemma','ashr'),'lfsr'=c(NA,blig$lfsr),'lfdr'=c(NA,blig$lfdr))}
#bhyb<-fread(paste0('results/summary_gwas_hybrid_corse_',pheno[p],'_freq_egs.txt.bz2'),data.table=F)
#bhyb<-subset(bhyb,bhyb$rs==RS)
#if(nrow(bhyb)==0){bhyb<-data.frame('rs'=RS,'pop'='hybrid_corse','beta'=NA,'se'=NA,'method'=c('gemma','ashr'),'lfsr'=NA,'lfdr'=NA)}else{bhyb<-data.frame('rs'=RS,'pop'='hybrid_corse','beta'=c(bhyb$beta,bhyb$PosteriorMean),'se'=c(bhyb$se,bhyb$PosteriorSD),'method'=c('gemma','ashr'),'lfsr'=c(NA,bhyb$lfsr),'lfdr'=c(NA,bhyb$lfdr))}

#beta_mantra<-readLines(paste0('results/run_mantra_',pheno,'_egs/mantra.bet.out'))
#beta_mantra<-beta_mantra[grep(paste0('\\b',RS,'\\b'),beta_mantra)]
#beta_mantra<-unlist(strsplit(beta_mantra,' '))
#beta_mantra<-beta_mantra[beta_mantra!='']
#popi<-which(beta_mantra%in%c('Mellifera_ncorse','Ligustica_Carnica','hybrid_corse'))
#pop<-beta_mantra[popi]
#beta<-beta_mantra[popi+1]
#se<-beta_mantra[popi+2]
#df_mantra<-data.frame('rs'=RS,'pop'=pop,'beta'=beta,'se'=se,'method'='mantra','lfsr'=NA,'lfdr'=NA)
#df_mash<-fread(paste0('results/summary_mash_',pheno,'_egs_simple_corr.txt'),data.table=F)
#df_mash<-subset(df_mash,df_mash$rs==RS)
#colnames(df_mash)[colnames(df_mash)=='post_mean']<-'beta'
#colnames(df_mash)[colnames(df_mash)=='post_sd']<-'se'
#df_mash$method<-'mash'
#colnames(df_mash)<-gsub('post_','',colnames(df_mash))
#colnames(df_mash)[colnames(df_mash)=='condition']<-'pop'
#effect<-do.call(rbind,list(bmel,blig,bhyb,df_mantra,df_mash))
#effect$pheno<-pheno
#effect$beta_min<-as.numeric(effect$beta)-1.96*as.numeric(effect$se)
#effect$beta_max<-as.numeric(effect$beta)+1.96*as.numeric(effect$se)
#effect$sign[effect$beta_min<=0 & effect$beta_max>=0]<-0 
#effect$sign[effect$beta_min<=0 & effect$beta_max<=0]<-'-' 
#effect$sign[effect$beta_min>=0 & effect$beta_max>=0]<-'+' 
#effect$pop[effect$pop=='Ligustica_Carnica']<-'Ligustica & Carnica'
#effect$pop[effect$pop=='Mellifera_ncorse']<-'Mellifera'
#effect$pop[effect$pop=='hybrid_corse']<-'Hybrids'
#effect$pop<-factor(effect$pop,levels=c('Ligustica & Carnica','Mellifera','Hybrids'))
#effect$method<-factor(effect$method,levels=c('gemma','ashr','mantra','mash'))

#g<-ggplot(effect)+
#	geom_point(aes(x=pop,y=as.numeric(beta),col=pop),size=2)+
#	geom_segment(aes(x=pop,xend=pop,y=beta_min,yend=beta_max,col=pop))+
#	xlab('Group')+ylab('SNP effect')+
#	geom_hline(yintercept=0,col='red',lty=3)+
#	facet_grid(.~method)+
#	scale_colour_manual(values=c('goldenrod2','grey34','deepskyblue2'))+
#	theme_bw()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
	
#print(g)
