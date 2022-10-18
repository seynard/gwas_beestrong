library(data.table)
library(ggplot2)
library(gridExtra)
library(ashr)
library(reshape2)
library(tidyverse)
args<-commandArgs(TRUE)

pheno<-args[1]

res_all<-fread(paste0('results/gemma_all_',pheno,'_freq_lmm_egs_cov.assoc.txt'))
fit.ash=ash(res_all$beta,res_all$se,mixcompdist='uniform')
res_all=as_tibble(cbind(res_all,fit.ash$result))
res_all<-res_all[,c('rs','p_wald','svalue','qvalue')]
#res_all<-res_all[,c('rs','p_wald')]
colnames(res_all)[2:ncol(res_all)]<-paste0(colnames(res_all)[2:ncol(res_all)],'_all')

res_hybrid<-fread(paste0('results/gemma_hybrid_',pheno,'_freq_lmm_egs_cov.assoc.txt'))
fit.ash=ash(res_hybrid$beta,res_hybrid$se,mixcompdist='uniform')
res_hybrid=as_tibble(cbind(res_hybrid,fit.ash$result))
res_hybrid<-res_hybrid[,c('rs','p_wald','svalue','qvalue')]
#res_hybrid<-res_hybrid[,c('rs','p_wald')]
colnames(res_hybrid)[2:ncol(res_hybrid)]<-paste0(colnames(res_hybrid)[2:ncol(res_hybrid)],'_hybrid')

res_lig<-fread(paste0('results/gemma_Ligustica_Carnica_',pheno,'_freq_lmm_egs_cov.assoc.txt'))
fit.ash=ash(res_lig$beta,res_lig$se,mixcompdist='uniform')
res_lig=as_tibble(cbind(res_lig,fit.ash$result))
res_lig<-res_lig[,c('rs','p_wald','svalue','qvalue')]
#res_lig<-res_lig[,c('rs','p_wald')]
colnames(res_lig)[2:ncol(res_lig)]<-paste0(colnames(res_lig)[2:ncol(res_lig)],'_lig')

res_mel<-fread(paste0('results/gemma_Mellifera_ncorse_',pheno,'_freq_lmm_egs_cov.assoc.txt'))
fit.ash=ash(res_mel$beta,res_mel$se,mixcompdist='uniform')
res_mel=as_tibble(cbind(res_mel,fit.ash$result))
res_mel<-res_mel[,c('rs','p_wald','svalue','qvalue')]
#res_mel<-res_mel[,c('rs','p_wald')]
colnames(res_mel)[2:ncol(res_mel)]<-paste0(colnames(res_mel)[2:ncol(res_mel)],'_mel')

res_corse<-fread(paste0('results/gemma_corse_',pheno,'_freq_lmm_egs_cov.assoc.txt'))
fit.ash=ash(res_corse$beta,res_corse$se,mixcompdist='uniform')
res_corse=as_tibble(cbind(res_corse,fit.ash$result))
res_corse<-res_corse[,c('rs','p_wald','svalue','qvalue')]
#res_corse<-res_corse[,c('rs','p_wald')]
colnames(res_corse)[2:ncol(res_corse)]<-paste0(colnames(res_corse)[2:ncol(res_corse)],'_corse')

res<-Reduce(function(x,y) merge(x,y,by='rs',all=T),list(res_all,res_hybrid,res_lig,res_mel,res_corse))
a<-colsplit(res$rs,':',c('CHROM','POS'))
res$CHROM<-a$CHROM
res$POS<-a$POS

g1<-ggplot(res)+geom_point(aes(x=POS,y=svalue_all))+ylim(1,0)+geom_hline(yintercept=0.1,col='red',lty=3)+facet_grid(.~CHROM,scales='free',space='free')+theme_bw()+ggtitle('all')
g2<-ggplot(res)+geom_point(aes(x=POS,y=svalue_hybrid))+ylim(1,0)+geom_hline(yintercept=0.1,col='red',lty=3)+facet_grid(.~CHROM,scales='free',space='free')+theme_bw()+ggtitle('hybrid')
g3<-ggplot(res)+geom_point(aes(x=POS,y=svalue_lig))+ylim(1,0)+geom_hline(yintercept=0.1,col='red',lty=3)+facet_grid(.~CHROM,scales='free',space='free')+theme_bw()+ggtitle('lig')
g4<-ggplot(res)+geom_point(aes(x=POS,y=svalue_mel))+ylim(1,0)+geom_hline(yintercept=0.1,col='red',lty=3)+facet_grid(.~CHROM,scales='free',space='free')+theme_bw()+ggtitle('mel')
g5<-ggplot(res)+geom_point(aes(x=POS,y=svalue_corse))+ylim(1,0)+geom_hline(yintercept=0.1,col='red',lty=3)+facet_grid(.~CHROM,scales='free',space='free')+theme_bw()+ggtitle('corse')

png(paste0('results/manhattan_',pheno,'.png'),height=1000,width=1200)
grid.arrange(g1,g2,g3,g4,g5,ncol=1)
dev.off()

#sbatch --mem=200G --wrap="Rscript scripts/test_gwas.r chargevarroaracine4"
#sbatch --mem=200G --wrap="Rscript scripts/test_gwas.r eb_smr"
#sbatch --mem=200G --wrap="Rscript scripts/test_gwas.r smr_brut"
#sbatch --mem=200G --wrap="Rscript scripts/test_gwas.r logit_taux_reop"
#sbatch --mem=200G --wrap="Rscript scripts/test_gwas.r logit_taux_reop_inf"
#sbatch --mem=200G --wrap="Rscript scripts/test_gwas.r varroadepthmitoracine4"
#sbatch --mem=200G --wrap="Rscript scripts/test_gwas.r logit_varroainfestation"




