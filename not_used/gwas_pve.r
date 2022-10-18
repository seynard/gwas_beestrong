#module load system/R-3.6.1
library(ggplot2)
library(tidyverse)
library(dplyr)
library(reshape2)
library(RColorBrewer)
args<-commandArgs(TRUE)
dir_out<-args[1]
log_files<-list.files(path=dir_out,pattern='.log.txt')

PVE<-data.frame('pve'=NA,'2.5%'=NA,'97.5%'=NA,'pop'=NA,'pheno'=NA,'grm'=NA,'gwas'=NA,'model'=NA)
for(i in 1:length(log_files)){
	file_name<-gsub('gemma_','',log_files[i])
	file_name<-gsub('.log.txt','',file_name)
	file_name<-gsub('_cov','',file_name)
	file_name<-unlist(strsplit(file_name,'_'))
	pop<-file_name[1]
	grm<-rev(file_name)[3]
	gwas<-rev(file_name)[1]
	model<-rev(file_name)[2]
	pheno<-file_name[-c(1,length(file_name)-2,length(file_name)-1,length(file_name))]
	pheno<-paste(pheno,collapse='_')
	if(grepl('bslmm',log_files[i])){
	logf<-read.table(paste0(dir_out,'/gemma_',pop,'_',pheno,'_',grm,'_bslmm_',gwas,'.hyp.txt'),header=T)
	pve<-c('pve'=mean(logf$pve),quantile(logf$pve, probs=c(0.025,0.975)),'pop'=pop,'pheno'=pheno,'grm'=grm,'gwas'=gwas,'model'=model)
	}else{
	logf<-readLines(paste0(dir_out,'/',log_files[i]))
	x<-as.numeric(rev(unlist(strsplit(logf[grep(' pve ',logf)],' ')))[1])
	y<-as.numeric(rev(unlist(strsplit(logf[grep('se\\(pve\\)',logf)],' ')))[1])
	pve<-c('pve'=x,'2.5%'=x-1.96*y,'97.5%'=x+1.96*y,'pop'=pop,'pheno'=pheno,'grm'=grm,'gwas'=gwas,'model'=model)
	}
PVE<-rbind(PVE,pve)
}
colnames(PVE)<-c('pve','min_pve','max_pve','pop','pheno','grm','gwas','model')
PVE<-subset(PVE,!is.na(PVE$pop))
PVE$type<-paste0(PVE$model,'_',PVE$grm,'_',PVE$gwas)
PVE$pop[grep('corse_',PVE$pheno)]<-'hybrid_corse'
PVE$pop[grep('ncorse_',PVE$pheno)]<-'Mellifera_ncorse'
PVE$pop[grep('nus_',PVE$pheno)]<-'all_nus'
PVE$pop[PVE$pop=='Ligustica']<-'Ligustica_Carnica'
PVE$pheno<-gsub('ncorse_','',PVE$pheno)
PVE$pheno<-gsub('Carnica_','',PVE$pheno)
PVE$pheno<-gsub('nus_','',PVE$pheno)
PVE$pheno<-gsub('corse_','',PVE$pheno)
PVE$pop<-factor(PVE$pop,levels=c('all','all_nus','hybrid','hybrid_corse','Ligustica_Carnica','Mellifera_ncorse','Mellifera','corse','us'))
PVE$pheno<-factor(PVE$pheno,levels=c('smr_brut','eb_smr','logit_varroainfestation','varroaphoretic','varroadepthmitoracine4','chargevarroaracine4','logit_taux_reop','logit_taux_reop_inf'))
PVE$type<-factor(PVE$type,levels=c('lmm_freq_freq','lmm_freq_egs','bslmm_freq_freq','bslmm_freq_egs'))
PVE$pve<-as.numeric(PVE$pve)
PVE$min_pve<-as.numeric(PVE$min_pve)
PVE$max_pve<-as.numeric(PVE$max_pve)

g<-ggplot(PVE)+
	geom_point(aes(x=pve,y=pop,col=pop))+
	geom_vline(xintercept=0,col='red',lty=3)+
	geom_vline(xintercept=1,col='red',lty=3)+
	geom_segment(aes(x=min_pve,xend=max_pve,y=pop,yend=pop,col=pop))+
	facet_grid(type~pheno)+
	xlim(-0.2,1.2)+theme_bw()

png(paste0(dir_out,'/summary_pve.png'),width=1500,height=600)
print(g)
dev.off()
