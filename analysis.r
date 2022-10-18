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
if(version==1){pop<-c('Mellifera','Ligustica_Carnica','hybrid')}else if (version==2){pop<-c('Mellifera','Ligustica_nus','hybrid')}
#pop<-c('Mellifera','Ligustica_Carnica','hybrid')
#pop<-c('Mellifera','Ligustica_nus','hybrid')

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
dat<-list()
for(j in 1:length(type)){
DF<-list()
for(i in 1:length(pheno)){
	df_mantra<-fread(paste0(dirout,'/run_mantra_',pheno[i],'_',type[j],'_',version,'/mantra.out'),data.table=F,header=F)
	colnames(df_mantra)<-c('rs','chr','ps','all_effect','all_other','nb_study','log10BF','post_proba_hetero','n','direction') 
	ll_no_corr<-read.table(paste0(dirout,'/result_mash_',pheno[i],'_',type[j],'_',version,'_no_corr.txt'),header=T,nrow=1)
	ll_no_corr<-ll_no_corr[,4]
	ll_simple_corr<-read.table(paste0(dirout,'/result_mash_',pheno[i],'_',type[j],'_',version,'_simple_corr.txt'),header=T,nrow=1)
	ll_simple_corr<-ll_simple_corr[,4]
	ll_em_corr<-read.table(paste0(dirout,'/result_mash_',pheno[i],'_',type[j],'_',version,'_em_corr.txt'),header=T,nrow=1)
	ll_em_corr<-ll_em_corr[,4]
	ll<-data.frame(ll_no_corr,ll_simple_corr,ll_em_corr)
	cor_choice<-gsub('ll_','',(names(which.max(ll))))
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
dat$pheno[dat$pheno=='varroaphoretic']<-'v_pho'
dat$pheno[dat$pheno=='varroadepthmitoracine4']<-'v_mito'
dat$pheno[dat$pheno=='logit_varroainfestation']<-'v_inf'
dat$pheno[dat$pheno=='chargevarroaracine4']<-'v_load'
dat$pheno[dat$pheno=='eb_smr']<-'mnr'
dat$pheno[dat$pheno=='logit_taux_reop_inf']<-'recap_inf'
chrsize<-dat%>%group_by(chr)%>%summarise(m=max(ps))
dat$pheno<-factor(dat$pheno,levels=c('v_pho','v_mito','v_inf','v_load','pc1_varroa_inf','mnr','recap_inf'))
corr_logBF<-dat%>%group_by(pheno,type)%>%summarise(cor=cor.test(log10BF_mantra,log10BF_mash,na.rm=T)$estimate)
corr_logBF
cor.test(dat$log10BF_mantra,dat$log10BF_mash,na.rm=T)

thresh_mantra<-ceiling(quantile(dat$log10BF_mantra,probs=0.99999))
print(paste0('sign threshold Mantra ',thresh_mantra))
thresh_mash<-ceiling(quantile(dat$log10BF_mash,probs=0.99999))
print(paste0('sign threshold Mash ',thresh_mash))

x_mash<-dat%>%group_by(type,pheno)%>%filter(log10BF_mash>=thresh_mash)
x_mash<-x_mash[,c('rs','pheno','type')]
x_mash$group<-paste0(x_mash$pheno,'_',x_mash$type)
x_mash[,c('pheno','type')]<-NULL
x_MASH<-list()
for(i in 1:length(unique(x_mash$group))){
	x_MASH[[i]]<-x_mash$rs[x_mash$group==unique(x_mash$group)[i]]
	}
names(x_MASH)<-paste0(unique(x_mash$group),'_MASH')
x_mantra<-dat%>%group_by(type,pheno)%>%filter(log10BF_mantra>=thresh_mantra)
x_mantra<-x_mantra[,c('rs','pheno','type')]
x_mantra$group<-paste0(x_mantra$pheno,'_',x_mantra$type)
x_mantra[,c('pheno','type')]<-NULL
x_MANTRA<-list()
for(i in 1:length(unique(x_mantra$group))){
	x_MANTRA[[i]]<-x_mantra$rs[x_mantra$group==unique(x_mantra$group)[i]]
	}
names(x_MANTRA)<-paste0(unique(x_mantra$group),'_MANTRA')
X<-c(x_MASH,x_MANTRA)

flog.threshold(ERROR)
venn_plot1<-venn.diagram(X[names(X)%in%c("mnr_egs_MASH","recap_inf_egs_MASH")],filename=NULL,category.names=c("MNR","Recap"),fill=c("red","orange"))
venn_plot2<-venn.diagram(X[names(X)%in%c("pc1_varroa_inf_egs_MANTRA","mnr_egs_MANTRA","recap_inf_egs_MANTRA")],filename=NULL,category.names=c("PC1","MNR","Recap"),fill=c("light green","dark green","green"))
#venn_plot3<-venn.diagram(X[names(X)%in%c("pc1_varroa_inf_egs_MASH","pc1_varroa_inf_egs_MANTRA")],filename=NULL,category.names=c("PC1","PC1"),fill=c("yellow","light green"))
venn_plot4<-venn.diagram(X[names(X)%in%c("recap_inf_egs_MASH","recap_inf_egs_MANTRA")],filename=NULL,category.names=c("Recap","Recap"),fill=c("orange","green"))
venn_plot5<-venn.diagram(X[names(X)%in%c("mnr_egs_MASH","mnr_egs_MANTRA")],filename=NULL,category.names=c("MNR","MNR"),fill=c("red","dark green"))
png(paste0(dirout,'/venn_',version,'.png'),width=3000,height=1500,res=400)
ggarrange(ggarrange(gTree(children=venn_plot1),gTree(children=venn_plot2),ncol=2),ggarrange(gTree(children=venn_plot4),gTree(children=venn_plot5),ncol=3),nrow=2)
dev.off()
flog.threshold(ERROR)
venn_plot1<-venn.diagram(X[names(X)%in%c("pc1_varroa_inf_freq_MANTRA","pc1_varroa_inf_egs_MANTRA")],filename=NULL,category.names=c("freq_MANTRA","egs_MANTRA"),fill=c("green","dark green"))
venn_plot2<-venn.diagram(X[names(X)%in%c("mnr_egs_MASH","mnr_freq_MANTRA","mnr_egs_MANTRA")],filename=NULL,category.names=c("egs_MAHS","freq_MANTRA","egs_MANTRA"),fill=c("orange","green","dark green"))
venn_plot3<-venn.diagram(X[names(X)%in%c("recap_inf_freq_MASH","recap_inf_egs_MASH","recap_inf_freq_MANTRA","recap_inf_egs_MANTRA")],filename=NULL,category.names=c("freq_MAHS","egs_MAHS","freq_MANTRA","egs_MANTRA"),fill=c("yellow","orange","green","dark green"))
png(paste0(dirout,'/venn_egs_freq_',version,'.png'),width=5000,height=1000,res=400)
ggarrange(gTree(children=venn_plot1),gTree(children=venn_plot2),gTree(children=venn_plot3),ncol=3)
dev.off()

dat_sign<-subset(dat, dat$log10BF_mantra>=thresh_mantra | dat$log10BF_mash>=thresh_mash)
table(dat_sign$pheno)
X<-list()
for(i in 1:nrow(dat_sign)){
	x<-subset(annot,annot$chr==dat_sign$chr[i] & annot$Start<=dat_sign$ps[i] & annot$Stop>=dat_sign$ps[i])
	if(nrow(x)==0){
	x<-subset(annot,annot$chr==dat_sign$chr[i])
	x$diff_start<-abs(dat_sign$ps[i]-x$Start)
	x$diff_stop<-abs(x$Stop-dat_sign$ps[i])
	m<-min(c(x$diff_start,x$diff_stop))
	x<-x[x$diff_start==m | x$diff_stop==m,]
	x<-unique(x[which.max(x$Length),])
	}else if (nrow(x)>=1){
	x<-unique(x[which.max(x$Length),])
	x$diff_start<-0
	x$diff_stop<-0}
	X[[i]]<-merge(dat_sign[i,],x,by='chr',all=T)
	}
X<-do.call(rbind,X)
X$gene_pos<-NA
X$gene_pos[X$diff_start==0 & X$diff_stop==0]<-'in'
X$gene_pos[X$ps<X$Start]<-'before'
X$gene_pos[X$ps>X$Stop]<-'after'
X$gene_dist<-NA
X$gene_dist[X$gene_pos=='in']<-0
X$gene_dist[X$gene_pos=='before']<-X$Start[X$gene_pos=='before']-X$ps[X$gene_pos=='before']
X$gene_dist[X$gene_pos=='after']<-X$ps[X$gene_pos=='after']-X$Stop[X$gene_pos=='after']
X[,c('diff_start','diff_stop')]<-NULL

locset<-unique(X$Locus)
go_apis<-biomartr::getGO(organism="apis mellifera",genes=locset,filters="ensembl_gene_id")
colnames(go_apis)<-c('Locus','description','go_nb')
X<-merge(X,go_apis,by='Locus',all=T)
write.table(X,paste0(dirout,'/sign_locus_',version,'.txt'),col.names=T,row.names=F,quote=F,sep=';')
