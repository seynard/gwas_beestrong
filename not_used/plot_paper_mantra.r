#module load system/R-4.1.1_gcc-9.3.0
library(data.table)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(reshape2)
args<-commandArgs(TRUE)
dirout<-'results'
dirin<-'data'

annot<-fread(paste0(dirin,'/proteins_48_403979.csv'))
colnames(annot)[1]<-'LG'
colnames(annot)[8]<-'prot'
colnames(annot)[10]<-'prot_name'
annot[,7]<-NULL
annot$chr[annot$LG=='linkage group LG1']<-1
annot$chr[annot$LG=='linkage group LG2']<-2
annot$chr[annot$LG=='linkage group LG3']<-3
annot$chr[annot$LG=='linkage group LG4']<-4
annot$chr[annot$LG=='linkage group LG5']<-5
annot$chr[annot$LG=='linkage group LG6']<-6
annot$chr[annot$LG=='linkage group LG7']<-7
annot$chr[annot$LG=='linkage group LG8']<-8
annot$chr[annot$LG=='linkage group LG9']<-9
annot$chr[annot$LG=='linkage group LG10']<-10
annot$chr[annot$LG=='linkage group LG11']<-11
annot$chr[annot$LG=='linkage group LG12']<-12
annot$chr[annot$LG=='linkage group LG13']<-13
annot$chr[annot$LG=='linkage group LG14']<-14
annot$chr[annot$LG=='linkage group LG15']<-15
annot$chr[annot$LG=='linkage group LG16']<-16
annot<-unique(annot)

pheno<-c('varroaphoretic','varroadepthmitoracine4','logit_varroainfestation','chargevarroaracine4','pc1','smr_brut','eb_smr','logit_taux_reop','logit_taux_reop_inf')
type<-c('egs','freq')
DF<-data.frame('ID'=NA,'chr'=NA,'ps'=NA,'all_effect'=NA,'all_other'=NA,'nb_study'=NA,'log10BF'=NA,'post_proba_hetero'=NA,'n'=NA,'direction'=NA,'pheno'=NA,'type'=NA)
for(j in 1:length(type)){
for(i in 1:length(pheno)){
	df<-fread(paste0(dirout,'/run_mantra_',pheno[i],'_',type[j],'/mantra.out'),data.table=F,header=F)
	colnames(df)<-c('ID','chr','ps','all_effect','all_other','nb_study','log10BF','post_proba_hetero','n','direction') 
	df$pheno<-pheno[i]
	df$type<-type[j]
	DF<-rbind(DF,df)
	}
}
DF$pheno[DF$pheno=='varroaphoretic']<-'v_pho'
DF$pheno[DF$pheno=='varroadepthmitoracine4']<-'v_mito'
DF$pheno[DF$pheno=='logit_varroainfestation']<-'v_inf'
DF$pheno[DF$pheno=='chargevarroaracine4']<-'v_load'
DF$pheno[DF$pheno=='pc1']<-'pc1_varroa_inf'
DF$pheno[DF$pheno=='smr_brut']<-'smr'
DF$pheno[DF$pheno=='eb_smr']<-'mnr'
DF$pheno[DF$pheno=='logit_taux_reop']<-'reop'
DF$pheno[DF$pheno=='logit_taux_reop_inf']<-'reop_inf'
chrsize<-DF%>%group_by(chr)%>%summarise(m=max(ps))
DF$pheno<-factor(DF$pheno,levels=c('v_pho','v_mito','v_inf','v_load','pc1_varroa_inf','smr','mnr','reop','reop_inf'))
DF<-subset(DF,!is.na(DF$pheno))

#for(j in 1:length(unique(DF$type))){
#print(j)
#for(i in 1:length(unique(DF$pheno))){
#print(i)
#png(paste0(dirout,'/',as.character(unique(DF$pheno)[i]),'_',unique(DF$type)[j],'_plot_mantra.png'),width=8000,height=1000,res=400)
#print(ggplot()+
#  geom_point(aes(x=ps,y=log10BF,col=post_proba_hetero),alpha=0.4,data=DF[DF$pheno==unique(DF$pheno)[i] & DF$type==unique(DF$type)[j],])+
#  scale_colour_gradient(name='Posterior probability \nof heterogeneity')+
#  geom_hline(aes(yintercept=5),col='black',lty=3,data=DF[DF$pheno==unique(DF$pheno)[i] & DF$type==unique(DF$type)[j],])+
#  geom_vline(aes(xintercept=1),col='red',lty=3,data=chrsize)+
#  geom_vline(aes(xintercept=m),col='red',lty=3,data=chrsize)+
#  scale_x_continuous(breaks=seq(0,chrsize$m[1],5000000))+
#  facet_grid(.~chr,scales='free',space='free')+theme_bw()+
#  xlab('position in bp')+theme(axis.text.x=element_text(angle=45)))
#dev.off()
#}
#}

#sign<-DF[DF$log10BF>=5,]
#X<-list()
#for(i in 1:nrow(sign)){
#	x<-sign[i,]
#	annot_x<-annot[annot$chr==x$chr,]
#	annot_x$diff_start<-abs(annot_x$Start-x$ps)
#	annot_x$diff_stop<-abs(annot_x$Stop-x$ps)
#	diff<-data.frame(row=c(seq(1:nrow(annot_x)),seq(1:nrow(annot_x))),diff=c(annot_x$diff_start,annot_x$diff_stop),sense=c(rep('start',nrow(annot_x)),rep('stop',nrow(annot_x))))
#	n<-diff[which.min(diff$diff),]
#	best<-annot_x[n$row,]
#	best$diff_start<-NULL
#	best$diff_stop<-NULL
#	X[[i]]<-cbind(x,best)
#}

#sign_best<-do.call(rbind,X)

#sign_best$gene_location<-NA
#sign_best$dist<-NA
#for(i in 1:nrow(sign_best)){
#	if(sign_best[i,'ps']>=sign_best[i,'Start'] & sign_best[i,'ps']<=sign_best[i,'Stop']){
#		sign_best[i,'gene_location']<-'in'
#		sign_best[i,'dist']<-0
#	}else if(sign_best[i,'ps']>sign_best[i,'Stop']){
#			sign_best[i,'gene_location']<-'after'
#			sign_best[i,'dist']<-sign_best[i,'ps']-sign_best[i,'Stop']
#	}else if(sign_best[i,'ps']<sign_best[i,'Start']){
#			sign_best[i,'gene_location']<-'before'
#			sign_best[i,'dist']<-sign_best[i,'Start']-sign_best[i,'ps']}
#}
#a<-colsplit(sign_best$direction,'',c('Mellifera_ncorse','Ligustica_Carnica','hybrid_corse'))
#sign_best$dir_mel<-a$Mellifera_ncorse
#sign_best$dir_lig<-a$Ligustica_Carnica
#sign_best$dir_hyb<-a$hybrid_corse
#write.table(sign_best,'summary_sign_prot.txt',col.names=T,row.names=F,quote=F,sep=';')

#plot 6 pheno EGS
df_plot<-DF[DF$type=='egs' & DF$pheno%in%c('v_pho','v_mito','v_inf','v_load','mnr','reop_inf'),]
png('plot_paper_mantra6pheno_egs.png',width=6000,height=4000,res=400)
ggplot()+
  geom_point(aes(x=ps,y=log10BF,col=post_proba_hetero),alpha=0.4,data=df_plot)+
  scale_colour_gradient(name='Posterior probability \nof heterogeneity')+
  geom_hline(aes(yintercept=5),col='black',lty=3,data=df_plot)+
  geom_vline(aes(xintercept=1),col='red',lty=3,data=chrsize)+
  geom_vline(aes(xintercept=m),col='red',lty=3,data=chrsize)+
  scale_x_continuous(breaks=seq(0,chrsize$m[1],5000000))+
  facet_grid(pheno~chr,scales='free',space='free')+theme_bw()+
  xlab('position in bp')+theme(axis.text.x=element_text(angle=45))
dev.off()
#plot 3 pheno EGS
df_plot<-DF[DF$type=='egs' & DF$pheno%in%c('pc1_varroa_inf','mnr','reop_inf'),]
png('plot_paper_mantra3pheno_egs.png',width=6000,height=2000,res=400)
ggplot()+
  geom_point(aes(x=ps,y=log10BF,col=post_proba_hetero),alpha=0.4,data=df_plot)+
  scale_colour_gradient(name='Posterior probability \nof heterogeneity')+
  geom_hline(aes(yintercept=5),col='black',lty=3,data=df_plot)+
  geom_vline(aes(xintercept=1),col='red',lty=3,data=chrsize)+
  geom_vline(aes(xintercept=m),col='red',lty=3,data=chrsize)+
  scale_x_continuous(breaks=seq(0,chrsize$m[1],5000000))+
  facet_grid(pheno~chr,scales='free',space='free')+theme_bw()+
  xlab('position in bp')+theme(axis.text.x=element_text(angle=45))
dev.off()
#plot 6 pheno freq
df_plot<-DF[DF$type=='freq' & DF$pheno%in%c('v_pho','v_mito','v_inf','v_load','mnr','reop_inf'),]
png('plot_paper_mantra6pheno_freq.png',width=6000,height=4000,res=400)
ggplot()+
  geom_point(aes(x=ps,y=log10BF,col=post_proba_hetero),alpha=0.4,data=df_plot)+
  scale_colour_gradient(name='Posterior probability \nof heterogeneity')+
  geom_hline(aes(yintercept=5),col='black',lty=3,data=df_plot)+
  geom_vline(aes(xintercept=1),col='red',lty=3,data=chrsize)+
  geom_vline(aes(xintercept=m),col='red',lty=3,data=chrsize)+
  scale_x_continuous(breaks=seq(0,chrsize$m[1],5000000))+
  facet_grid(pheno~chr,scales='free',space='free')+theme_bw()+
  xlab('position in bp')+theme(axis.text.x=element_text(angle=45))
dev.off()
#plot 3 pheno freq
df_plot<-DF[DF$type=='freq' & DF$pheno%in%c('pc1_varroa_inf','mnr','reop_inf'),]
png('plot_paper_mantra3pheno_freq.png',width=6000,height=2000,res=400)
ggplot()+
  geom_point(aes(x=ps,y=log10BF,col=post_proba_hetero),alpha=0.4,data=df_plot)+
  scale_colour_gradient(name='Posterior probability \nof heterogeneity')+
  geom_hline(aes(yintercept=5),col='black',lty=3,data=df_plot)+
  geom_vline(aes(xintercept=1),col='red',lty=3,data=chrsize)+
  geom_vline(aes(xintercept=m),col='red',lty=3,data=chrsize)+
  scale_x_continuous(breaks=seq(0,chrsize$m[1],5000000))+
  facet_grid(pheno~chr,scales='free',space='free')+theme_bw()+
  xlab('position in bp')+theme(axis.text.x=element_text(angle=45))
dev.off()
#plot 6 pheno 
df_plot<-DF[DF$pheno%in%c('v_pho','v_mito','v_inf','v_load','mnr','reop_inf'),]
df_plot$pheno<-paste0(df_plot$pheno,'_',df_plot$type)
df_plot$pheno<-factor(df_plot$pheno,levels=c('v_pho_egs','v_pho_freq','v_mito_egs','v_mito_freq','v_inf_egs','v_inf_freq','v_load_egs','v_load_freq','mnr_egs','mnr_freq','reop_inf_egs','reop_inf_freq'))
png('plot_paper_mantra6pheno.png',width=6000,height=4500,res=400)
ggplot()+
  geom_point(aes(x=ps,y=log10BF,col=post_proba_hetero),alpha=0.4,data=df_plot)+
  scale_colour_gradient(name='Posterior probability \nof heterogeneity')+
  geom_hline(aes(yintercept=5),col='black',lty=3,data=df_plot)+
  geom_vline(aes(xintercept=1),col='red',lty=3,data=chrsize)+
  geom_vline(aes(xintercept=m),col='red',lty=3,data=chrsize)+
  scale_x_continuous(breaks=seq(0,chrsize$m[1],5000000))+
  facet_grid(pheno~chr,scales='free',space='free')+theme_bw()+
  xlab('position in bp')+theme(axis.text.x=element_text(angle=45))
dev.off()
#plot 3 pheno 
df_plot<-DF[DF$pheno%in%c('pc1_varroa_inf','mnr','reop_inf'),]
df_plot$pheno<-paste0(df_plot$pheno,'_',df_plot$type)
df_plot$pheno<-factor(df_plot$pheno,levels=c('pc1_varroa_inf_egs','pc1_varroa_inf_freq','mnr_egs','mnr_freq','reop_inf_egs','reop_inf_freq'))
png('plot_paper_mantra3pheno.png',width=6000,height=3500,res=400)
ggplot()+
  geom_point(aes(x=ps,y=log10BF,col=post_proba_hetero),alpha=0.4,data=df_plot)+
  scale_colour_gradient(name='Posterior probability \nof heterogeneity')+
  geom_hline(aes(yintercept=5),col='black',lty=3,data=df_plot)+
  geom_vline(aes(xintercept=1),col='red',lty=3,data=chrsize)+
  geom_vline(aes(xintercept=m),col='red',lty=3,data=chrsize)+
  scale_x_continuous(breaks=seq(0,chrsize$m[1],5000000))+
  facet_grid(pheno~chr,scales='free',space='free')+theme_bw()+
  xlab('position in bp')+theme(axis.text.x=element_text(angle=45))
dev.off()

