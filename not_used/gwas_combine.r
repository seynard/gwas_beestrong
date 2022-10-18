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
library(VennDiagram)
args<-commandArgs(TRUE)
dir_out<-args[1]
sp<-args[2]

print(sp)
p<-c('varroaphoretic','varroadepthmitoracine4','logit_varroainfestation','chargevarroaracine4','eb_smr','logit_taux_reop_inf')

df<-list()
for(i in 1:length(p)){
	print(p[i])
	x<-fread(paste0(dir_out,'/gemma_',sp,'_',p[i],'_freq_lmm_egs_cov.assoc.txt'),data.table=F)
	a<-colsplit(x$rs,':',c('CHROM','POS'))
	x$CHROM<-a$CHROM
	x$POS<-a$POS
	fit.ash=ash(x$beta,x$se,mixcompdist='uniform')
	x=as_tibble(cbind(x,fit.ash$result))	
	x[,c('chr','ps','n_miss','allele1','allele0','af','logl_H1','l_remle','betahat','sebetahat','NegativeProb','PositiveProb','lfsr','lfdr','PosteriorMean','PosteriorSD')]<-NULL
	x$sp<-sp
	x$pheno<-p[i]
	x$type<-'freq_egs'
	y<-fread(paste0(dir_out,'/gemma_',sp,'_',p[i],'_freq_lmm_freq_cov.assoc.txt'),data.table=F)
	a<-colsplit(y$rs,':',c('CHROM','POS'))
	y$CHROM<-a$CHROM
	y$POS<-a$POS
	fit.ash=ash(y$beta,y$se,mixcompdist='uniform')
	y=as_tibble(cbind(y,fit.ash$result))	
	y[,c('chr','ps','n_miss','allele1','allele0','af','logl_H1','l_remle','betahat','sebetahat','NegativeProb','PositiveProb','lfsr','lfdr','PosteriorMean','PosteriorSD')]<-NULL
	y$sp<-sp
	y$pheno<-p[i]
	y$type<-'freq_freq'
	x<-rbind(x,y)
	df[[i]]=x
}
df<-do.call(rbind,df)
write.table(df,paste0(dir_out,'/summary_gwas_',sp,'.txt'),col.names=T,row.names=F,quote=F)


#display_venn <- function(x, ...){
#  library(VennDiagram)
#  grid.newpage()
#  venn_object <- venn.diagram(x, filename = NULL, ...)
#  grid.draw(venn_object)
#}

#annot<-fread('data/proteins_48_403979.csv',data.table=F)
#annot[,1]<-NULL
#annot$CHROM<-NA
#annot$CHROM[annot$Accession=='NC_037638.1']<-1
#annot$CHROM[annot$Accession=='NC_037639.1']<-2
#annot$CHROM[annot$Accession=='NC_037640.1']<-3
#annot$CHROM[annot$Accession=='NC_037641.1']<-4
#annot$CHROM[annot$Accession=='NC_037642.1']<-5
#annot$CHROM[annot$Accession=='NC_037643.1']<-6
#annot$CHROM[annot$Accession=='NC_037644.1']<-7
#annot$CHROM[annot$Accession=='NC_037645.1']<-8
#annot$CHROM[annot$Accession=='NC_037646.1']<-9
#annot$CHROM[annot$Accession=='NC_037647.1']<-10
#annot$CHROM[annot$Accession=='NC_037648.1']<-11
#annot$CHROM[annot$Accession=='NC_037649.1']<-12
#annot$CHROM[annot$Accession=='NC_037650.1']<-13
#annot$CHROM[annot$Accession=='NC_037651.1']<-14
#annot$CHROM[annot$Accession=='NC_037652.1']<-15
#annot$CHROM[annot$Accession=='NC_037653.1']<-16
#annot<-subset(annot,!is.na(annot$CHROM))
#annot_short<-annot%>%group_by(Locus,CHROM)%>%summarise(Start=min(Start),Stop=max(Stop))

#L<-list()
#for(i in 1:length(p)){
#	print(p[i])
#	L[[i]]<-subset(df[[i]]$rs,df[[i]]$qvalue<=0.1)}
#names(L)<-p
#display_venn(L[1:4],fill=c("#999999", "#E69F00", "#56B4E9", "#009E73"))
#display_venn(L[5:6],fill=c("#999999", "#E69F00"))
#L_all<-unlist(L)
#x<-data.frame(rs=L_all)
#a<-colsplit(x$rs,':',c('CHROM','POS'))
#x$CHROM<-a$CHROM
#x$POS<-a$POS
#z<-df[[1]][,c('rs','CHROM','POS')]
#x<-rbind(x,z)
#x<-x%>%group_by(rs,CHROM,POS)%>%summarise(n=n())

#ggplot()+
#	geom_segment(aes(x=Start,xend=Stop,y=1.5,yend=1.5),data=annot_short,alpha=0.2,col='red',size=4)+
#	geom_point(aes(x=POS,y=n),data=x)+
#	facet_grid(CHROM~.,space='free_x')+
#	ylim(1,5)+
#	theme_bw()



