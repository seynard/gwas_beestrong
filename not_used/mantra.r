library(data.table)
library(tidyverse)
library(reshape2)
args<-commandArgs(TRUE)
dir_res<-args[1]
dir_mantra<-args[2]
gwas<-args[4]
pheno<-args[5]
model<-args[6]
input_data<-args[7]

list_pop<-args[3]
print(list_pop)
list_pop<-unlist(strsplit(list_pop,','))
L<-list()
for(i in 1:length(list_pop)){
	if(gwas=='gemma'){
	X<-fread(paste0(dir_res,'/',gwas,'_',list_pop[i],'_',pheno,'_',model,'_',input_data,'.assoc.txt'),data.table=F)
	x<-readLines(paste0(dir_res,'/',gwas,'_',list_pop[i],'_',pheno,'_',model,'_',input_data,'.log.txt'))
	n<-grep('analyzed individuals',x)
	n<-as.numeric(gsub('## number of analyzed individuals = ','',x[n]))
	X$n<-n-X$n_miss
	X[,c('chr','ps','n_miss','logl_H1','l_remle','p_wald')]<-NULL
	X<-gather(X,variable,val,c('allele1','allele0','af','beta','se','n'))
	}else if (gwas=='ldak'){
	X<-fread(paste0(dir_res,'/',pheno,'_',list_pop[i],'_',input_data,'.assoc'),data.table=F)
	x<-fread(paste0(dir_res,'/',pheno,'_',list_pop[i],'_',input_data,'.summaries'),data.table=F)
	X<-merge(X,x,by='Predictor',all=T)
	X[,c('Chromosome','Basepair','Wald_Stat','Wald_P','Effect_Liability','SD_Liability','MAF','A1.y','A2.y','Direction','Stat')]<-NULL
	colnames(X)[colnames(X)=='A1.x']<-'allele1'
	colnames(X)[colnames(X)=='A2.x']<-'allele0'
	colnames(X)[colnames(X)=='Effect']<-'beta'
	colnames(X)[colnames(X)=='SD']<-'se'
	colnames(X)[colnames(X)=='A1_Mean']<-'af'
	colnames(X)[colnames(X)=='Predictor']<-'rs'
	X<-gather(X,variable,val,c('allele1','allele0','af','beta','se','n'))
	}
	colnames(X)[colnames(X)=='val']<-list_pop[i]
	L[[i]]<-data.table(X)
	}

DF<-Reduce(function(x,y) merge(x,y,by=c('rs','variable'),all=T),L)
DF<-as.data.frame(DF)
a<-colsplit(DF$rs,':',c('chr','pos'))
DF$chr<-a$chr
DF$pos<-a$pos
DF<-DF[,c('rs','chr','pos','variable',list_pop)]
print('DF')

DF1<-DF[DF$variable%in%c('allele1','allele0'),]
DF1$value<-apply(DF1[,list_pop],1,function(x)unique(na.omit(x)))
DF1[,list_pop]<-NULL
DF1<-spread(DF1,variable,value)
DF1<-DF1[,c('rs','chr','pos','allele1','allele0')]
print('DF1')
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
print('DF2')

X<-Reduce(function(x,y) merge(x,y,by='rs',all=T),df2)
X<-as.data.frame(X)
X<-merge(DF1,X,by='rs',all=T)
df<-gather(X,snp,string,c(val,list_pop))
df<-df[order(df$rs) ,]
df$snp<-ordered(df$snp,levels=c('val',list_pop))
df[,c('rs','snp')]<-NULL
write.table(df,paste0(dir_mantra,'/mantra.dat'),col.names=F,row.names=F,quote=F)
write.table(list_pop,paste0(dir_mantra,'/mantra.in'),col.names=F,row.names=F,quote=F)
