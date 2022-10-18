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
#dirin<-args[1]
#dirout<-args[2]
#pheno<-args[3]
#pheno<-unlist(strsplit(pheno,','))
#type<-args[4]
#type<-unlist(strsplit(type,','))
dirin<-'data'
dirout<-'results'
type<-c('egs','freq')
pheno<-c('varroaphoretic','varroadepthmitoracine4','logit_varroainfestation','chargevarroaracine4','pc1','eb_smr','logit_taux_reop_inf')
version<-as.numeric(args[1])
if(version==1){lpop<-c('Mellifera','Ligustica_Carnica','hybrid')}else if (version==2){lpop<-c('Mellifera','Ligustica_nus','hybrid')}

#analyse Mantra and Mash results
dat<-list()
for(j in 1:length(type)){
DF<-list()
for(i in 1:length(pheno)){
	df_mantra<-fread(paste0(dirout,'/run_mantra_',pheno[i],'_',type[j],'_',version,'/mantra.out'),data.table=F,header=F)
	colnames(df_mantra)<-c('rs','chr','ps','all_effect','all_other','nb_study','log10BF','post_proba_hetero','n','direction') 
#	ll_no_corr<-read.table(paste0(dirout,'/result_mash_',pheno[i],'_',type[j],'_',version,'_no_corr.txt'),header=T,nrow=1)
#	ll_no_corr<-ll_no_corr[,4]
#	ll_simple_corr<-read.table(paste0(dirout,'/result_mash_',pheno[i],'_',type[j],'_',version,'_simple_corr.txt'),header=T,nrow=1)
#	ll_simple_corr<-ll_simple_corr[,4]
#	ll_em_corr<-read.table(paste0(dirout,'/result_mash_',pheno[i],'_',type[j],'_',version,'_em_corr.txt'),header=T,nrow=1)
#	ll_em_corr<-ll_em_corr[,4]
#	ll<-data.frame(ll_no_corr,ll_simple_corr,ll_em_corr)
#	cor_choice<-gsub('ll_','',(names(which.max(ll))))
	cor_choice<-'simple_corr'
#	print(cor_choice)
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
dat$pheno[dat$pheno=='varroaphoretic']<-'v_pho'
dat$pheno[dat$pheno=='varroadepthmitoracine4']<-'v_mito'
dat$pheno[dat$pheno=='logit_varroainfestation']<-'v_inf'
dat$pheno[dat$pheno=='chargevarroaracine4']<-'v_load'
dat$pheno[dat$pheno=='eb_smr']<-'mnr'
dat$pheno[dat$pheno=='logit_taux_reop_inf']<-'recap_inf'
chrsize<-dat%>%group_by(chr)%>%summarise(m=max(ps))
dat$pheno<-factor(dat$pheno,levels=c('v_pho','v_mito','v_inf','v_load','pc1_varroa_inf','mnr','recap_inf'))
thresh_mantra<-ceiling(quantile(dat$log10BF_mantra,probs=0.99999))
print(paste0('sign threshold Mantra ',thresh_mantra))
thresh_mash<-ceiling(quantile(dat$log10BF_mash,probs=0.99999))
print(paste0('sign threshold Mash ',thresh_mash))

EFFECT<-list()
datx<-subset(dat,dat$type=='egs' & dat$pheno%in%c('pc1_varroa_inf','mnr','recap_inf') & (dat$log10BF_mash>=thresh_mash|dat$log10BF_mantra>=thresh_mantra))
for(i in 1:nrow(datx)){
	pheno<-datx$pheno[i]
	rsx<-datx$rs[i]
	if(pheno=='pc1_varroa_inf'){phenox<-'pc1'}
	if(pheno=='mnr'){phenox<-'eb_smr'}
	if(pheno=='recap_inf'){phenox<-'logit_taux_reop_inf'}
	mel<-fread(paste0('results/summary_gwas_',lpop[1],'_',phenox,'_freq_egs.txt.bz2'),data.table=F)
	mel<-subset(mel,mel$rs==rsx)
	if(nrow(mel)==0){mel<-data.frame('rs'=rsx,'pop'='Mellifera','beta'=NA,'se'=NA,'method'=c('gemma','ashr'),'lfsr'=NA,'lfdr'=NA)}else{mel<-data.frame('rs'=rsx,'pop'='Mellifera','beta'=c(mel$beta,mel$PosteriorMean),'se'=c(mel$se,mel$PosteriorSD),'method'=c('gemma','ashr'),'lfsr'=c(NA,mel$lfsr),'lfdr'=c(NA,mel$lfdr))}
	lig<-fread(paste0('results/summary_gwas_',lpop[2],'_',phenox,'_freq_egs.txt.bz2'),data.table=F)
	lig<-subset(lig,lig$rs==rsx)
	if(nrow(lig)==0){lig<-data.frame('rs'=rsx,'pop'='Ligustica_Carnica','beta'=NA,'se'=NA,'method'=c('gemma','ashr'),'lfsr'=NA,'lfdr'=NA)}else{lig<-data.frame('rs'=rsx,'pop'='Ligustica_Carnica','beta'=c(lig$beta,lig$PosteriorMean),'se'=c(lig$se,lig$PosteriorSD),'method'=c('gemma','ashr'),'lfsr'=c(NA,lig$lfsr),'lfdr'=c(NA,lig$lfdr))}
	hyb<-fread(paste0('results/summary_gwas_',lpop[3],'_',phenox,'_freq_egs.txt.bz2'),data.table=F)
	hyb<-subset(hyb,hyb$rs==rsx)
	if(nrow(hyb)==0){hyb<-data.frame('rs'=rsx,'pop'='hybrid','beta'=NA,'se'=NA,'method'=c('gemma','ashr'),'lfsr'=NA,'lfdr'=NA)}else{hyb<-data.frame('rs'=rsx,'pop'='hybrid','beta'=c(hyb$beta,hyb$PosteriorMean),'se'=c(hyb$se,hyb$PosteriorSD),'method'=c('gemma','ashr'),'lfsr'=c(NA,hyb$lfsr),'lfdr'=c(NA,hyb$lfdr))}
	beta_mantra<-readLines(paste0('results/run_mantra_',phenox,'_egs_',version,'/mantra.bet.out'))
	beta_mantra<-beta_mantra[grep(paste0('\\b',rsx,'\\b'),beta_mantra)]
	beta_mantra<-unlist(strsplit(beta_mantra,' '))
	beta_mantra<-beta_mantra[beta_mantra!='']
	if(version==1){popi<-which(beta_mantra%in%c('Mellifera','Ligustica_Carnica','hybrid'))}else if(version==2){popi<-which(beta_mantra%in%c('Mellifera','Ligustica_nus','hybrid'))}
	pop<-beta_mantra[popi]
	beta<-beta_mantra[popi+1]
	se<-beta_mantra[popi+2]
	df_mantra<-data.frame('rs'=rsx,'pop'=pop,'beta'=beta,'se'=se,'method'='mantra','lfsr'=NA,'lfdr'=NA)
	if(version==2){df_mantra$pop[df_mantra$pop=='Ligustica_nus']<-'Ligustica_Carnica'}	
	df_mash<-fread(paste0('results/summary_mash_',phenox,'_egs_',version,'_simple_corr.txt'),data.table=F)
	df_mash<-subset(df_mash,df_mash$rs==rsx)
	colnames(df_mash)[colnames(df_mash)=='post_mean']<-'beta'
	colnames(df_mash)[colnames(df_mash)=='post_sd']<-'se'
	df_mash$method<-'mash'
	colnames(df_mash)<-gsub('post_','',colnames(df_mash))
	colnames(df_mash)[colnames(df_mash)=='condition']<-'pop'
	if(version==2){df_mash$pop[df_mash$pop=='Ligustica_nus']<-'Ligustica_Carnica'}	
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
	g<-ggplot(effect)+
	geom_point(aes(x=pop,y=as.numeric(beta),col=pop),size=2)+
	geom_segment(aes(x=pop,xend=pop,y=beta_min,yend=beta_max,col=pop))+
	xlab('Group')+ylab('SNP effect')+
	geom_hline(yintercept=0,col='red',lty=3)+
	facet_grid(.~method)+
	scale_colour_manual(values=c('goldenrod2','grey34','deepskyblue2'),name='Group')+
	theme_bw()+
	theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_text(size=20),axis.text.y=element_text(size=12))+
	theme(legend.title=element_text(size=20),legend.text=element_text(size=12))+
	theme(strip.text.x=element_text(size=20))
	png(paste0('plot_effects_',as.character(unique(effect$pheno)),'_',unique(effect$rs),'_',version,'.png'),width=3000,height=1000,res=300)
	print(g)
	dev.off()
EFFECT[[i]]<-effect
}
EFFECT<-do.call(rbind,EFFECT)
write.table(EFFECT,paste0('effect_',version,'.txt'),col.names=T,row.names=F,quote=F)
