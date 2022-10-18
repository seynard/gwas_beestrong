library(data.table)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ashr)
library(reshape2)

pop<-c('Ligustica_Carnica','Mellifera_ncorse')
pheno<-c('chargevarroaracine4','eb_smr','logit_taux_reop','logit_taux_reop_inf','logit_varroainfestation','smr_brut','varroadepthmitoracine4','varroaphoretic')
grm<-'freq'
gwas<-c('freq','egs')
args<-commandArgs(TRUE)

i=as.numeric(args[1])
j=as.numeric(args[2])
k=1
l=as.numeric(args[3])

#summary gwas
gwas_lmm<-fread(paste0('results/gemma_',pop[i],'_',pheno[j],'_',grm[k],'_lmm_',gwas[l],'_cov.assoc.txt'))
fit.ash=ash(gwas_lmm$beta,gwas_lmm$se,mixcompdist='uniform')
gwas_lmm=as_tibble(cbind(gwas_lmm,fit.ash$result))
a<-colsplit(gwas_lmm$rs,':',c('chr','ps'))
gwas_lmm$chr<-a$chr
gwas_lmm$ps<-a$ps
gwas_lmm<-gwas_lmm[,c('chr','ps','p_wald','qvalue','svalue')]
gwas_bslmm<-fread(paste0('results/gemma_',pop[i],'_',pheno[j],'_',grm[k],'_',gwas[l],'_bslmm.param.txt'))
a<-colsplit(gwas_bslmm$rs,':',c('chr','ps'))
gwas_bslmm$chr<-a$chr
gwas_bslmm$ps<-a$ps
gwas_bslmm<-gwas_bslmm[,c('chr','ps','gamma')]
df<-merge(gwas_lmm,gwas_bslmm,by=c('chr','ps'),all=T)
df$pop<-pop[i]
df$pheno<-pheno[j]
df$grm<-grm[k]
df$gwas<-gwas[l]
write.table(df,paste0('results/summary_gwas_',i,'_',j,'_',k,'_',l,'.txt'),col.names=T,row.names=F,quote=F)





