library(data.table)
library(reshape2)

w<-fread('data/weights_Mellifera_ncorse_freq.short',data.table=F)
colnames(w)<-c('rs','w')
snp<-fread('data/in_Mellifera_ncorse_freq.freq',data.table=F,select=1)
colnames(snp)<-'rs'

W<-merge(snp,w,by='rs',all.x=T)
W$w[is.na(W$w)]<-0
a<-colsplit(W$rs,':',c('CHROM','POS'))
W$CHROM<-a$CHROM
W$POS<-a$POS
W<-W[order(W$CHROM,W$POS),]

W<-data.frame(W$rs,W$w)
colnames(W)<-c('rs','w')
write.table(W,'data/ld_score_Mellifera_ncorse_freq.txt',col.names=T,row.names=F,quote=F)


####################################################################
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
library(qvalue)
library(factoextra)
library(nFactors)
library(PCDimension)
library(PCAtools)

df<-fread('results/LDAKgmat_Ligustica_Carnica_freq.grm.raw',data.table=F)
spca <- SamplePCA(t(df))
ag.obj <- AuerGervini(spca)
agDimension(ag.obj)
f<- makeAgCpmFun("Exponential")
agfuns<-list(twice=agDimTwiceMean,specc=agDimSpectral,km=agDimKmeans, km3=agDimKmeans3,tt=agDimTtest, tt2=agDimTtest2,cpt=agDimCPT, cpm=f)
compareAgDimMethods(ag.obj, agfuns)
plot(ag.obj, agfuns)
df<-fread('results/LDAKgmat_Ligustica_Carnica_egs.grm.raw',data.table=F)
spca <- SamplePCA(t(df))
ag.obj <- AuerGervini(spca)
agDimension(ag.obj)
f<- makeAgCpmFun("Exponential")
agfuns<-list(twice=agDimTwiceMean,specc=agDimSpectral,km=agDimKmeans, km3=agDimKmeans3,tt=agDimTtest, tt2=agDimTtest2,cpt=agDimCPT, cpm=f)
compareAgDimMethods(ag.obj, agfuns)
plot(ag.obj, agfuns)

df<-fread('results/LDAKgmat_Mellifera_ncorse_freq.grm.raw',data.table=F)
spca <- SamplePCA(t(df))
ag.obj <- AuerGervini(spca)
agDimension(ag.obj)
f<- makeAgCpmFun("Exponential")
agfuns<-list(twice=agDimTwiceMean,specc=agDimSpectral,km=agDimKmeans, km3=agDimKmeans3,tt=agDimTtest, tt2=agDimTtest2,cpt=agDimCPT, cpm=f)
compareAgDimMethods(ag.obj, agfuns)
plot(ag.obj, agfuns)
df<-fread('results/LDAKgmat_Mellifera_ncorse_egs.grm.raw',data.table=F)
spca <- SamplePCA(t(df))
ag.obj <- AuerGervini(spca)
agDimension(ag.obj)
f<- makeAgCpmFun("Exponential")
agfuns<-list(twice=agDimTwiceMean,specc=agDimSpectral,km=agDimKmeans, km3=agDimKmeans3,tt=agDimTtest, tt2=agDimTtest2,cpt=agDimCPT, cpm=f)
compareAgDimMethods(ag.obj, agfuns)
plot(ag.obj, agfuns)






paran(df, iterations = 500, centile = 0, quietly = FALSE, 
    status = TRUE, all = TRUE, cfa = TRUE, graph = TRUE, color = TRUE, 
    col = c("black", "red", "blue"), lty = c(1, 2, 3), lwd = 1, legend = TRUE, 
    file = "", width = 640, height = 640, grdevice = "png", seed = 0


pclf<-fread('results/pca_Ligustica_Carnica_freq.values',data.table=F)
colnames(pclf)<-'eigenval'
pclf$pv<-(pclf$eigenval/sum(pclf$eigenval))
eiglf<-nCng(pclf$eigenval,model="components",details=TRUE)
pcle<-fread('results/pca_Ligustica_Carnica_egs.values',data.table=F)
colnames(pcle)<-'eigenval'
pcle$pv<-(pcle$eigenval/sum(pcle$eigenval))
eigle<-nCng(pcle$eigenval,model="components",details=TRUE)

pcmf<-fread('results/pca_Mellifera_ncorse_freq.values',data.table=F)
colnames(pcmf)<-'eigenval'
pcmf$pv<-(pcmf$eigenval/sum(pcmf$eigenval))
eigmf<-nCng(pcmf$eigenval,model="components",details=TRUE)
pcme<-fread('results/pca_Mellifera_ncorse_egs.values',data.table=F)
colnames(pcme)<-'eigenval'
pcme$pv<-(pcme$eigenval/sum(pcme$eigenval))
eigme<-nCng(pcme$eigenval,model="components",details=TRUE)

par(mfrow=c(2,2))
plot(seq(1:nrow(pclf)),pclf$pv,pch=16)
points(seq(1:eiglf$nFactors),pclf[1:eiglf$nFactors,'pv'],pch=16,col='red')
plot(seq(1:nrow(pcmf)),pcmf$pv,pch=16)
points(seq(1:eigmf$nFactors),pcmf[1:eigmf$nFactors,'pv'],pch=16,col='red')
plot(seq(1:nrow(pcle)),pcle$pv,pch=16)
points(seq(1:eigle$nFactors),pcle[1:eigle$nFactors,'pv'],pch=16,col='red')
plot(seq(1:nrow(pcme)),pcme$pv,pch=16)
points(seq(1:eigme$nFactors),pcme[1:eigme$nFactors,'pv'],pch=16,col='red')



     
     
     
