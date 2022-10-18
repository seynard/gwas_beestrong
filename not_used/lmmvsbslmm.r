library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(magrittr)
library(ashr)
library(gridExtra)
library(ggpubr)
library(viridis)
library(RColorBrewer)
library(gplots)
library(reshape2)
library(qvalue)
library(cowplot)
args<-commandArgs(TRUE)
TYPE<-args[1]
PHENO<-args[2]

pp.ggplot = function(pval) {
    dat = as_tibble(pval)
    nmax=nrow(dat)
    dat %<>% mutate(stat = -log10(value)) %>%
        mutate(rank = dense_rank(value),
               expect = -log10(rank/(nmax+1)),
               hival = -log10(qbeta(0.975,rank,nmax-rank)),
               loval = -log10(qbeta(0.025,rank,nmax-rank)) )
    dat %<>% arrange(expect)
    x = dat$expect
    x = c(x,rev(x))
    y = c(dat$loval,rev(dat$hival))
    poly.dat=tibble(x,y)
    chisq <- qchisq(1 - pval, 1)
	lambda <- median(chisq) / qchisq(0.5, 1)
    p = dat %>% ggplot(aes(x=expect,y=stat)) +
        geom_polygon(data=poly.dat,aes(x=x,y=y),fill='blue', alpha=0.1) +
        geom_path(data=poly.dat,aes(x=x,y=y), alpha=0.1) +
        geom_point() +
        geom_abline(intercept=0,slope=1,colour='blue') +
        xlab("Expected -log(p)") +
        ylab("Observed -log(p)")+ 
		annotate(geom = "text",x = -Inf,y = Inf,hjust = -0.15,vjust = 1 + 0.15 * 3,label = sprintf("λ = %.2f", lambda),size = 8) +
		theme(axis.ticks = element_line(size = 0.5),panel.grid = element_blank())
    p
}

#type<-c('Ligustica_Carnica','Mellifera_ncorse')
#pheno<-c('smr_brut','eb_smr','chargevarroaracine4','logit_taux_reop','logit_taux_reop_inf','logit_varroainfestation','varroaphoretic','varroadepthmitoracine4')

#for(i in 1:length(type)){
#	TYPE<-type[i]
	print(TYPE)
#	for(j in 1:length(pheno)){
#		PHENO<-pheno[j]
		print(PHENO)
gwas_lmm<-fread(paste0('results/gemma_',TYPE,'_',PHENO,'_freq_lmm_egs_cov.assoc.txt'),data.table=F)
log_lmm<-readLines(paste0('results/gemma_',TYPE,'_',PHENO,'_freq_lmm_egs_cov.log.txt'))
pve<-tail(unlist(strsplit(log_lmm[grep("pve estimate",log_lmm)],' ')),n=1)
se_pve<-tail(unlist(strsplit(log_lmm[grep("se",log_lmm)],' ')),n=1)
vg<-tail(unlist(strsplit(log_lmm[grep("vg",log_lmm)],' ')),n=1)
ve<-tail(unlist(strsplit(log_lmm[grep("ve",log_lmm)],' ')),n=1)
lmm_params<-data.frame("pve"=pve,"se_pve"=se_pve,"vg"=vg,"ve"=ve)
a<-colsplit(gwas_lmm$rs,':',c("chr","ps"))
gwas_lmm$chr<-a$chr
gwas_lmm$ps<-a$ps
pv=as.data.frame(gwas_lmm%>%select(p_wald))
gwas.ppplot<-pp.ggplot(pv[,1])+theme_bw()
p_hist<-ggplot(pv)+
	geom_histogram(aes(x=p_wald),fill='grey')+
	xlab('p-values')+theme_bw()
df<-gwas_lmm%>%sample_n(100000)
g11<-ggplot(df,aes(x=af,y=beta))+
    geom_point(alpha=0.01)+ 
    geom_smooth()+
    xlab("Reference Allele Frequency")+ 
    ylab("SNP Effect")
g12<-ggplot(df,aes(x=af,y=beta))+
    geom_point(alpha=0.01)+ 
    geom_smooth()+
    ylim(-10,10)+
    xlab("Reference Allele Frequency")+ 
    ylab("SNP Effect")

png(paste0('results/diag_lmm_',TYPE,'_',PHENO,'.png'),width=1000,height=200)		
grid.arrange(gwas.ppplot,p_hist,g11,g12,ncol=4)
dev.off()

fit.ash=ash(gwas_lmm$beta,gwas_lmm$se,mixcompdist='uniform')
gwas_lmm=as_tibble(cbind(gwas_lmm,fit.ash$result))
qobj_pval<-qvalue(p=gwas_lmm$p_wald)
qobj_pval<-max(gwas_lmm$p_wald[qobj_pval$qvalue<0.1 & !is.na(qobj_pval$qvalue)])
qobj_qval<-max(gwas_lmm$p_wald[gwas_lmm$qvalue<0.1 & !is.na(gwas_lmm$qvalue)])
qobj_sval<-max(gwas_lmm$p_wald[gwas_lmm$svalue<0.1 & !is.na(gwas_lmm$svalue)])

gwas_bslmm<-fread(paste0('results/gemma_',TYPE,'_',PHENO,'_freq_bslmm_egs_cov.param.txt'),data.table=F)
bslmm_log<-read.table(paste0('results/gemma_',TYPE,'_',PHENO,'_freq_bslmm_egs_cov.hyp.txt'),header=T)
pve<-c("PVE", mean(bslmm_log$pve),quantile(bslmm_log$pve, probs=c(0.5,0.025,0.975)))
pge<-c("PGE",mean(bslmm_log$pge),quantile(bslmm_log$pge, probs=c(0.5,0.025,0.975)))
rho<-c("rho", mean(bslmm_log$rho),quantile(bslmm_log$rho, probs=c(0.5,0.025,0.975)))
h<-c("h", mean(bslmm_log$h),quantile(bslmm_log$h, probs=c(0.5,0.025,0.975)))
pi<-c("pi",mean(bslmm_log$pi),quantile(bslmm_log$pi, probs=c(0.5,0.025,0.975)))
n_gamma<-c("n_gamma",mean(bslmm_log$n_gamma),quantile(bslmm_log$n_gamma, probs=c(0.5,0.025,0.975)))
bslmm_params<-as.data.frame(rbind(pve,pge,rho,h,pi,n_gamma),row.names=F)
colnames(bslmm_params)<-c("hyperparam", "mean","median","2.5%", "97.5%")
a<-colsplit(gwas_bslmm$rs,':',c("chr","ps"))
gwas_bslmm$chr<-a$chr
gwas_bslmm$ps<-a$ps
gwas_bslmm$bslmm_effect<-gwas_bslmm$beta*gwas_bslmm$gamma

gwas_bslmm2<-fread(paste0('back_up/version_nov2021/gemma_',TYPE,'_',PHENO,'_freq_lmm_egs_cov_bslmm.param.txt'),data.table=F)
bslmm_log2<-read.table(paste0('back_up/version_nov2021/gemma_',TYPE,'_',PHENO,'_freq_lmm_egs_cov_bslmm.hyp.txt'),header=T)
pve2<-c("PVE", mean(bslmm_log2$pve),quantile(bslmm_log2$pve, probs=c(0.5,0.025,0.975)))
pge2<-c("PGE",mean(bslmm_log2$pge),quantile(bslmm_log2$pge, probs=c(0.5,0.025,0.975)))
rho2<-c("rho", mean(bslmm_log2$rho),quantile(bslmm_log2$rho, probs=c(0.5,0.025,0.975)))
h2<-c("h", mean(bslmm_log2$h),quantile(bslmm_log2$h, probs=c(0.5,0.025,0.975)))
pi2<-c("pi",mean(bslmm_log2$pi),quantile(bslmm_log2$pi, probs=c(0.5,0.025,0.975)))
n_gamma2<-c("n_gamma",mean(bslmm_log2$n_gamma),quantile(bslmm_log2$n_gamma, probs=c(0.5,0.025,0.975)))
bslmm_params2<-as.data.frame(rbind(pve2,pge2,rho2,h2,pi2,n_gamma2),row.names=F)
colnames(bslmm_params2)<-c("hyperparam", "mean","median","2.5%", "97.5%")
a<-colsplit(gwas_bslmm2$rs,':',c("chr","ps"))
gwas_bslmm2$chr<-a$chr
gwas_bslmm2$ps<-a$ps
gwas_bslmm2$bslmm_effect<-gwas_bslmm2$beta*gwas_bslmm2$gamma

png(paste0('results/diag_bslmm_',TYPE,'_',PHENO,'.png'),width=1000,height=2000)		
par(mfrow=c(12,3))
plot(bslmm_log$pve, type="l", ylab="PVE", main="PVE - trace")
hist(bslmm_log$pve, main="PVE - posterior distribution", xlab="PVE")
plot(density(bslmm_log$pve), main="PVE - posterior distribution", xlab="PVE")
plot(bslmm_log2$pve, type="l", ylab="PVE", main="PVE - trace")
hist(bslmm_log2$pve, main="PVE - posterior distribution", xlab="PVE")
plot(density(bslmm_log2$pve), main="PVE - posterior distribution", xlab="PVE")
plot(bslmm_log$h, type="l", ylab="h", main="h - trace")
hist(bslmm_log$h, main="h - posterior distribution", xlab="h")
plot(density(bslmm_log$h), main="h - posterior distribution", xlab="h")
plot(bslmm_log2$h, type="l", ylab="h", main="h - trace")
hist(bslmm_log2$h, main="h - posterior distribution", xlab="h")
plot(density(bslmm_log2$h), main="h - posterior distribution", xlab="h")
plot(bslmm_log$pge, type="l", ylab="PGE", main="PGE - trace")
hist(bslmm_log$pge, main="PGE - posterior distribution", xlab="PGE")
plot(density(bslmm_log$pge), main="PGE - posterior distribution", xlab="PGE")
plot(bslmm_log2$pge, type="l", ylab="PGE", main="PGE - trace")
hist(bslmm_log2$pge, main="PGE - posterior distribution", xlab="PGE")
plot(density(bslmm_log2$pge), main="PGE - posterior distribution", xlab="PGE")
plot(bslmm_log$rho, type="l", ylab="rho", main="rho - trace")
hist(bslmm_log$rho, main="rho - posterior distribution", xlab="rho")
plot(density(bslmm_log$rho), main="rho - posterior distribution", xlab="rho")
plot(bslmm_log2$rho, type="l", ylab="rho", main="rho - trace")
hist(bslmm_log2$rho, main="rho - posterior distribution", xlab="rho")
plot(density(bslmm_log2$rho), main="rho - posterior distribution", xlab="rho")
plot(bslmm_log$pi, type="l", ylab="pi", main="pi - trace")
hist(bslmm_log$pi, main="pi - posterior distribution", xlab="pi")
plot(density(bslmm_log$pi), main="pi - posterior distribution", xlab="pi")
plot(bslmm_log2$pi, type="l", ylab="pi", main="pi - trace")
hist(bslmm_log2$pi, main="pi - posterior distribution", xlab="pi")
plot(density(bslmm_log2$pi), main="pi - posterior distribution", xlab="pi")
plot(bslmm_log$n_gamma, type="l", ylab="n_gamma", main="n_gamma - trace")
hist(bslmm_log$n_gamma, main="n_gamma - posterior distribution", xlab="n_gamma")
plot(density(bslmm_log$n_gamma), main="n_gamma - posterior distribution", xlab="n_gamma")
plot(bslmm_log2$n_gamma, type="l", ylab="n_gamma", main="n_gamma - trace")
hist(bslmm_log2$n_gamma, main="n_gamma - posterior distribution", xlab="n_gamma")
plot(density(bslmm_log2$n_gamma), main="n_gamma - posterior distribution", xlab="n_gamma")
dev.off()

print(lmm_params)
print(bslmm_params)
print(bslmm_params2)

g1<-ggplot(gwas_lmm)+
	geom_point(aes(x=ps,y=-log10(p_wald),col=beta))+
	facet_grid(.~chr,scales='free', space='free')+
	geom_hline(aes(yintercept=-log10(qobj_pval)),col='red')+
	geom_hline(aes(yintercept=-log10(qobj_qval)),col='green')+
	geom_hline(aes(yintercept=-log10(qobj_sval)),col='blue')+
	theme_bw()+theme(legend.position='none')+ggtitle('lmm covariates')
g2<-ggplot(gwas_lmm)+
	geom_point(aes(x=ps,y=qvalue,col=beta))+
	geom_hline(aes(yintercept=0.1),col='red')+
	ylim(max(gwas_lmm$qvalue),0)+
	facet_grid(.~chr,scales='free',space='free')+
	theme_bw()+theme(legend.position='none')+ggtitle('qvalues')
g3<-ggplot(gwas_lmm)+
	geom_point(aes(x=ps,y=svalue,col=beta))+
	geom_hline(aes(yintercept=0.1),col='red')+
	ylim(max(gwas_lmm$svalue),0)+
	facet_grid(.~chr,scales='free',space='free')+
	theme_bw()+theme(legend.position='none')+ggtitle('svalues')
g5<-ggplot(gwas_bslmm2)+
	geom_point(aes(x=ps,y=gamma,col=bslmm_effect))+
	facet_grid(.~chr,scales='free', space='free')+
	theme_bw()+theme(legend.position='none')+ggtitle('bslmm no covariates')
g4<-ggplot(gwas_bslmm)+
	geom_point(aes(x=ps,y=gamma,col=bslmm_effect))+
	facet_grid(.~chr,scales='free',space='free')+
	theme_bw()+theme(legend.position='none')+ggtitle('bslmm covariates')
png(paste0('results/manhattan_lmm_bslmm_',TYPE,'_',PHENO,'.png'),width=1500,height=1000)
grid.arrange(g1,g2,g3,g4,g5,nrow=5)
dev.off()
#		}
#	}

#cd /work/genphyse/dynagen/seynard/GWAS/infestation
#module load system/R-3.6.1
#sbatch --mem=20G --wrap='Rscript scripts/not_used/lmmvsbslmm.r Ligustica_Carnica smr_brut'
#sbatch --mem=20G --wrap='Rscript scripts/not_used/lmmvsbslmm.r Ligustica_Carnica eb_smr'
#sbatch --mem=20G --wrap='Rscript scripts/not_used/lmmvsbslmm.r Ligustica_Carnica chargevarroaracine4'
#sbatch --mem=20G --wrap='Rscript scripts/not_used/lmmvsbslmm.r Ligustica_Carnica logit_taux_reop'
#sbatch --mem=20G --wrap='Rscript scripts/not_used/lmmvsbslmm.r Ligustica_Carnica logit_taux_reop_inf'
#sbatch --mem=20G --wrap='Rscript scripts/not_used/lmmvsbslmm.r Ligustica_Carnica logit_varroainfestation'
#sbatch --mem=20G --wrap='Rscript scripts/not_used/lmmvsbslmm.r Ligustica_Carnica varroaphoretic'
#sbatch --mem=20G --wrap='Rscript scripts/not_used/lmmvsbslmm.r Ligustica_Carnica varroadepthmitoracine4'
#sbatch --mem=20G --wrap='Rscript scripts/not_used/lmmvsbslmm.r Mellifera_ncorse smr_brut'
#sbatch --mem=20G --wrap='Rscript scripts/not_used/lmmvsbslmm.r Mellifera_ncorse eb_smr'
#sbatch --mem=20G --wrap='Rscript scripts/not_used/lmmvsbslmm.r Mellifera_ncorse chargevarroaracine4'
#sbatch --mem=20G --wrap='Rscript scripts/not_used/lmmvsbslmm.r Mellifera_ncorse logit_taux_reop'
#sbatch --mem=20G --wrap='Rscript scripts/not_used/lmmvsbslmm.r Mellifera_ncorse logit_taux_reop_inf'
#sbatch --mem=20G --wrap='Rscript scripts/not_used/lmmvsbslmm.r Mellifera_ncorse logit_varroainfestation'
#sbatch --mem=20G --wrap='Rscript scripts/not_used/lmmvsbslmm.r Mellifera_ncorse varroaphoretic'
#sbatch --mem=20G --wrap='Rscript scripts/not_used/lmmvsbslmm.r Mellifera_ncorse varroadepthmitoracine4'








