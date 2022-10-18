library(data.table)
library(tidyverse)
library(reshape2)
args<-commandArgs(TRUE)
dir_out<-args[1]
list_pop<-args[2]
list_pop<-unlist(strsplit(list_pop,','))
pheno<-args[3]
type<-args[4]
version<-args[5]

L<-list()
for(i in 1:length(list_pop)){
	X<-fread(paste0(dir_out,'/gemma_',list_pop[i],'_',pheno,'_freq_lmm_',type,'_cov.assoc.txt'),data.table=F)
	x<-readLines(paste0(dir_out,'/gemma_',list_pop[i],'_',pheno,'_freq_lmm_',type,'_cov.log.txt'))
	n<-grep('analyzed individuals',x)
	n<-as.numeric(gsub('## number of analyzed individuals = ','',x[n]))
	X$n<-n-X$n_miss
	X[,c('chr','ps','n_miss','logl_H1','l_remle','p_wald')]<-NULL
	X<-gather(X,variable,val,c('allele1','allele0','af','beta','se','n'))
	colnames(X)[colnames(X)=='val']<-list_pop[i]
	L[[i]]<-data.table(X)
	}
DF<-Reduce(function(x,y) merge(x,y,by=c('rs','variable'),all=T),L)
DF<-as.data.frame(DF)
a<-colsplit(DF$rs,':',c('chr','pos'))
DF$chr<-a$chr
DF$pos<-a$pos
DF<-DF[,c('rs','chr','pos','variable',list_pop)]
DF1<-DF[DF$variable%in%c('allele1','allele0'),]
DF1$value<-apply(DF1[,list_pop],1,function(x)unique(na.omit(x)))
DF1[,list_pop]<-NULL
DF1<-spread(DF1,variable,value)
DF1<-DF1[,c('rs','chr','pos','allele1','allele0')]
DF1$val<-paste0(DF1$rs,' ',DF1$chr,' ',DF1$pos,' ',DF1$allele1,' ',DF1$allele0)
DF1[,c('chr','pos','allele1','allele0')]<-NULL
DF2<-DF[DF$variable%in%c('n','af','beta','se'),]
df2<-list()
for(i in 1:length(list_pop)){
	df2i<-DF2[,c('rs','variable',list_pop[i],'chr','pos')]
	df2i<-spread(df2i,variable,list_pop[i])
	df2i<-df2i[,c('rs','n','af','beta','se')]
	df2i$presence<-1
	df2i$presence[is.na(df2i$af)]<-0
	df2i[df2i$presence==0,c('n','af','beta','se')]<-0
	df2i<-df2i[,c('rs','presence','n','af','beta','se')]
	df2i$val<-paste0(df2i$presence,' ',df2i$n,' ',df2i$af,' ',df2i$beta,' ',df2i$se)
	df2i[,c('presence','n','af','beta','se')]<-NULL
	colnames(df2i)[2]<-list_pop[i]
	df2[[i]]<-data.table(df2i)
}
X<-Reduce(function(x,y) merge(x,y,by='rs',all=T),df2)
X<-as.data.frame(X)
X<-merge(DF1,X,by='rs',all=T)
df<-gather(X,snp,string,c(val,all_of(list_pop)))
df<-df[order(df$rs) ,]
df$snp<-ordered(df$snp,levels=c('val',list_pop))
df[,c('rs','snp')]<-NULL
write.table(df,paste0(dir_out,'/run_mantra_',pheno,'_',type,'_',version,'/mantra.dat'),col.names=F,row.names=F,quote=F)
write.table(list_pop,paste0(dir_out,'/run_mantra_',pheno,'_',type,'_',version,'/mantra.in'),col.names=F,row.names=F,quote=F)

