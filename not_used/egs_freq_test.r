library(data.table)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(ashr)
library(gridExtra)
library(ggpubr)
library(viridis)
library(RColorBrewer)
library(gplots)
library(reshape2)
args<-commandArgs(TRUE)
dir_out<-'results'
population='Mellifera_ncorse'

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

phenotype='eb_smr'
population='Ligustica_Carnica'
egs<-fread(paste0('results/gemma_',population,'_',phenotype,'_lmm_egs.assoc.txt'),data.table=F)
freq<-fread(paste0('results/gemma_',population,'_',phenotype,'_lmm_freq.assoc.txt'),data.table=F)

DF<-merge(egs,freq,by=c('chr','ps','rs','allele1','allele0'),all=T)
par(mfrow=c(1,1))
plot(DF$beta.x,DF$beta.y)
par(mfrow=c(1,2))
hist(DF$p_wald.x)
hist(DF$p_wald.y)
par(mfrow=c(1,1))
plot(-log10(DF$p_wald.x),-log10(DF$p_wald.y))

pv_egs=as.data.frame(egs%>%select(p_wald))
gwas.ppplot_egs=pp.ggplot(pv_egs[,1])+theme_bw()+ggtitle(paste(population, phenotype))
pv_freq=as.data.frame(freq%>%select(p_wald))
gwas.ppplot_freq=pp.ggplot(pv_freq[,1])+theme_bw()+ggtitle(paste(population, phenotype))

grid.arrange(gwas.ppplot_egs,gwas.ppplot_freq,ncol=2)

a<-colsplit(DF$rs,':',c('CHROM','POS'))
DF$chr<-a$CHROM
DF$ps<-a$POS

ggplot(DF)+
geom_point(aes(x=ps,y=-log10(p_wald.x)),col='red')+
geom_point(aes(x=ps,y=log10(p_wald.y)),col='blue')+
facet_grid(chr~.)

### PCA
list_snp<-fread('data/HAV3_1_50000.txt')
list_snp$rs<-paste0(list_snp$CHROM,':',list_snp$POS)
Q<-fread('data/queen_geno_50k.Q')

#Mellifera ncorse
g_egs<-fread('data/in_Mellifera_ncorse.egs',data.table=F)
colnames(g_egs)[1]<-'rs'
g_egs<-subset(g_egs,g_egs$rs%in%list_snp$rs)
g_egs[,c('rs','REF','ALT')]<-NULL
g_egs<-t(g_egs)
pca_egs<-prcomp(g_egs)
PoV_egs<-pca_egs$sdev^2/sum(pca_egs$sdev^2)

g_freq<-fread('data/in_Mellifera_ncorse.freq',data.table=F)
colnames(g_freq)[1]<-'rs'
g_freq<-subset(g_freq,g_freq$rs%in%list_snp$rs)
g_freq[,c('rs','REF','ALT')]<-NULL
g_freq<-t(g_freq)
g_freq[is.na(g_freq)]<--9
pca_freq<-prcomp(g_freq)
PoV_freq<-pca_freq$sdev^2/sum(pca_freq$sdev^2)

png('pca_test_Mellifera_ncorse.png',width=6000,height=2000,res=300)
par(mfrow=c(2,6))
plot(pca_egs$x[,1],pca_egs$x[,2],pch=16,xlab=round(PoV_egs[1],digit=3),ylab=round(PoV_egs[2],digit=3))
plot(pca_egs$x[,1],pca_egs$x[,3],pch=16,xlab=round(PoV_egs[1],digit=3),ylab=round(PoV_egs[3],digit=3))
plot(pca_egs$x[,1],pca_egs$x[,4],pch=16,xlab=round(PoV_egs[1],digit=3),ylab=round(PoV_egs[4],digit=3))
plot(pca_egs$x[,2],pca_egs$x[,3],pch=16,xlab=round(PoV_egs[2],digit=3),ylab=round(PoV_egs[3],digit=3))
plot(pca_egs$x[,2],pca_egs$x[,4],pch=16,xlab=round(PoV_egs[2],digit=3),ylab=round(PoV_egs[4],digit=3))
plot(pca_egs$x[,3],pca_egs$x[,4],pch=16,xlab=round(PoV_egs[3],digit=3),ylab=round(PoV_egs[4],digit=3))
plot(pca_freq$x[,1],pca_freq$x[,2],pch=16,xlab=round(PoV_freq[1],digit=3),ylab=round(PoV_freq[2],digit=3))
plot(pca_freq$x[,1],pca_freq$x[,3],pch=16,xlab=round(PoV_freq[1],digit=3),ylab=round(PoV_freq[3],digit=3))
plot(pca_freq$x[,1],pca_freq$x[,4],pch=16,xlab=round(PoV_freq[1],digit=3),ylab=round(PoV_freq[4],digit=3))
plot(pca_freq$x[,2],pca_freq$x[,3],pch=16,xlab=round(PoV_freq[2],digit=3),ylab=round(PoV_freq[3],digit=3))
plot(pca_freq$x[,2],pca_freq$x[,4],pch=16,xlab=round(PoV_freq[2],digit=3),ylab=round(PoV_freq[4],digit=3))
plot(pca_freq$x[,3],pca_freq$x[,4],pch=16,xlab=round(PoV_freq[3],digit=3),ylab=round(PoV_freq[4],digit=3))
dev.off()

q<-subset(Q,Q$num_ruche_bs%in%rownames(g_egs))
png('q_test_Mellifera_ncorse.png',width=1000,height=1000,res=300)
plot(q$Ligustica_Carnica,q$Mellifera,pch=16)
dev.off()

#Ligustica
g_egs<-fread('data/in_Ligustica_Carnica.egs',data.table=F)
colnames(g_egs)[1]<-'rs'
g_egs<-subset(g_egs,g_egs$rs%in%list_snp$rs)
g_egs[,c('rs','REF','ALT')]<-NULL
g_egs<-t(g_egs)
pca_egs<-prcomp(g_egs)
PoV_egs<-pca_egs$sdev^2/sum(pca_egs$sdev^2)

g_freq<-fread('data/in_Ligustica_Carnica.freq',data.table=F)
colnames(g_freq)[1]<-'rs'
g_freq<-subset(g_freq,g_freq$rs%in%list_snp$rs)
g_freq[,c('rs','REF','ALT')]<-NULL
g_freq<-t(g_freq)
g_freq[is.na(g_freq)]<--9
pca_freq<-prcomp(g_freq)
PoV_freq<-pca_freq$sdev^2/sum(pca_freq$sdev^2)

png('pca_test_Ligustica_Carnica.png',width=6000,height=2000,res=300)
par(mfrow=c(2,6))
plot(pca_egs$x[,1],pca_egs$x[,2],pch=16,xlab=round(PoV_egs[1],digit=3),ylab=round(PoV_egs[2],digit=3))
plot(pca_egs$x[,1],pca_egs$x[,3],pch=16,xlab=round(PoV_egs[1],digit=3),ylab=round(PoV_egs[3],digit=3))
plot(pca_egs$x[,1],pca_egs$x[,4],pch=16,xlab=round(PoV_egs[1],digit=3),ylab=round(PoV_egs[4],digit=3))
plot(pca_egs$x[,2],pca_egs$x[,3],pch=16,xlab=round(PoV_egs[2],digit=3),ylab=round(PoV_egs[3],digit=3))
plot(pca_egs$x[,2],pca_egs$x[,4],pch=16,xlab=round(PoV_egs[2],digit=3),ylab=round(PoV_egs[4],digit=3))
plot(pca_egs$x[,3],pca_egs$x[,4],pch=16,xlab=round(PoV_egs[3],digit=3),ylab=round(PoV_egs[4],digit=3))
plot(pca_freq$x[,1],pca_freq$x[,2],pch=16,xlab=round(PoV_freq[1],digit=3),ylab=round(PoV_freq[2],digit=3))
plot(pca_freq$x[,1],pca_freq$x[,3],pch=16,xlab=round(PoV_freq[1],digit=3),ylab=round(PoV_freq[3],digit=3))
plot(pca_freq$x[,1],pca_freq$x[,4],pch=16,xlab=round(PoV_freq[1],digit=3),ylab=round(PoV_freq[4],digit=3))
plot(pca_freq$x[,2],pca_freq$x[,3],pch=16,xlab=round(PoV_freq[2],digit=3),ylab=round(PoV_freq[3],digit=3))
plot(pca_freq$x[,2],pca_freq$x[,4],pch=16,xlab=round(PoV_freq[2],digit=3),ylab=round(PoV_freq[4],digit=3))
plot(pca_freq$x[,3],pca_freq$x[,4],pch=16,xlab=round(PoV_freq[3],digit=3),ylab=round(PoV_freq[4],digit=3))
dev.off()

q<-subset(Q,Q$num_ruche_bs%in%rownames(g_egs))
png('q_test_Ligustica_Carnica.png',width=1000,height=1000,res=300)
plot(q$Ligustica_Carnica,q$Mellifera,pch=16)
dev.off()
