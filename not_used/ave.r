setwd('C:/Users/seynard/Downloads/')
library(data.table)
library(ggplot2)
library(reshape2)

llig<-fread('list_Ligustica_Carnica.txt',data.table=F,header=F)
colnames(llig)<-'col_id'
llig$sp<-'lig'
lmel<-fread('list_Mellifera.txt',data.table=F,header=F)
colnames(lmel)<-'col_id'
lmel$sp<-'mel'
lcau<-fread('list_Caucasica.txt',data.table=F,header=F)
colnames(lcau)<-'col_id'
lcau$sp<-'cau'
lhyb<-fread('list_hybrid.txt',data.table=F,header=F)
colnames(lhyb)<-'col_id'
lhyb$sp<-'hyb'
L<-do.call(rbind,list(llig,lmel,lcau,lhyb))

df_col<-fread('average_depth_col.txt',data.table=F)
head(df_col)
df_col<-merge(df_col,L,by='col_id',all=T)
summary(df_col)
ggplot(df_col)+geom_histogram(aes(x=ave_depth,fill=sp))+facet_wrap(~sp)

df_snp<-fread('average_depth_snp.txt',data.table=F)
head(df_snp)
summary(df_snp)
hist(df_snp$ave_depth)
a<-colsplit(df_snp$rs,':',c('chr','ps'))
df_snp$chr<-a$chr
df_snp$ps<-a$ps
png('depth_snp.png',width=500,height=3000)
ggplot(df_snp)+geom_point(aes(y=max_depth,x=ps),col='red')+geom_point(aes(y=min_depth,x=ps),col='blue')+facet_grid(chr~.)
dev.off()
