library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)
library(gridExtra)
args<-commandArgs(TRUE)
dir_out<-args[1]
pheno<-args[2]
cov_file<-args[3]
grm_type<-args[4]
type<-args[5]

COV<-fread(cov_file)
PHENO<-fread(paste0(dir_out,'/pheno_',pheno,'_',type,'.txt'))
df<-cbind(PHENO,COV)
colnames(df)<-c('pheno',paste0('cov',seq(1:ncol(COV))))

fit=lm(pheno~.,df)
corrected_phenotype=data.frame(resid(fit))
colnames(corrected_phenotype)<-'corrected_phenotype'
df$nb<-rownames(df)
corrected_phenotype$nb<-rownames(corrected_phenotype)

df<-merge(df,corrected_phenotype,by='nb',all=T)

write.table(df$corrected_phenotype,paste0(dir_out,'/corrected_',pheno,'_',type,'_',grm_type,'.txt'),col.names=F,row.names=F,quote=F)
