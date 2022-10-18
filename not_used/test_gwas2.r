library(data.table)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)

mel<-fread('results/gemma_Mellifera_ncorse_logit_taux_reop_inf_freq_lmm_egs_cov.assoc.txt',data.table=F)
mel[,c('chr','ps','n_miss','allele1','allele0','af','logl_H1','l_remle')]<-NULL
mel$sp<-'mel'
lig<-fread('results/gemma_Ligustica_Carnica_logit_taux_reop_inf_freq_lmm_egs_cov.assoc.txt',data.table=F)
lig[,c('chr','ps','n_miss','allele1','allele0','af','logl_H1','l_remle')]<-NULL
lig$sp<-'lig'
cor<-fread('results/gemma_corse_logit_taux_reop_inf_freq_lmm_egs_cov.assoc.txt',data.table=F)
cor[,c('chr','ps','n_miss','allele1','allele0','af','logl_H1','l_remle')]<-NULL
cor$sp<-'cor'
hyb<-fread('results/gemma_hybrid_logit_taux_reop_inf_freq_lmm_egs_cov.assoc.txt',data.table=F)
hyb[,c('chr','ps','n_miss','allele1','allele0','af','logl_H1','l_remle')]<-NULL
hyb$sp<-'hyb'

x<-do.call(rbind,list(mel,lig,cor,hyb))
x$p<--log10(x$p_wald)
a<-colsplit(x$rs,':',c('CHROM','POS'))
x$CHROM<-a$CHROM
x$POS<-a$POS
x$pheno<-'logit_taux_reop_inf'

mel<-fread('results/gemma_Mellifera_ncorse_eb_smr_freq_lmm_egs_cov.assoc.txt',data.table=F)
mel[,c('chr','ps','n_miss','allele1','allele0','af','logl_H1','l_remle')]<-NULL
mel$sp<-'mel'
lig<-fread('results/gemma_Ligustica_Carnica_eb_smr_freq_lmm_egs_cov.assoc.txt',data.table=F)
lig[,c('chr','ps','n_miss','allele1','allele0','af','logl_H1','l_remle')]<-NULL
lig$sp<-'lig'
cor<-fread('results/gemma_corse_eb_smr_freq_lmm_egs_cov.assoc.txt',data.table=F)
cor[,c('chr','ps','n_miss','allele1','allele0','af','logl_H1','l_remle')]<-NULL
cor$sp<-'cor'
hyb<-fread('results/gemma_hybrid_eb_smr_freq_lmm_egs_cov.assoc.txt',data.table=F)
hyb[,c('chr','ps','n_miss','allele1','allele0','af','logl_H1','l_remle')]<-NULL
hyb$sp<-'hyb'

y<-do.call(rbind,list(mel,lig,cor,hyb))
y$p<--log10(y$p_wald)
a<-colsplit(y$rs,':',c('CHROM','POS'))
y$CHROM<-a$CHROM
y$POS<-a$POS
y$pheno<-'eb_smr'

df<-rbind(x,y)
df$t<-NA
df$t[df$p>6]<-'sign'
sub<-subset(df,df$t=='sign')
plot<-df%>%group_by(CHROM,POS)%>%summarise(p=mean(p))

ggplot()+
	geom_point(aes(x=POS,y=p,col=sp,shape=pheno),data=sub)+
	geom_point(aes(x=POS,y=p),data=plot,alpha=0.1)+
	scale_colour_manual(values=c('brown','cyan','goldenrod2','grey34'))+
	scale_shape_manual(values=c(8, 20))+	
	theme_bw()+
	facet_grid(CHROM~.)

#ggplot()+
#	geom_point(aes(x=POS,y=p,col=sp,shape=pheno),data=sub)+
#	geom_point(aes(x=POS,y=p),data=plot,alpha=0.1)+
#	scale_colour_manual(values=c('brown','cyan','goldenrod2','grey34'))+
#	scale_shape_manual(values=c(8, 20))+	
#	theme_bw()+
#	facet_grid(.~CHROM,scales='free_x',space='free_x')

