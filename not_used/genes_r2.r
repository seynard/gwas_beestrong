library(data.table)
library(tidyverse)
library(reshape2)

df_r2<-fread('locus_ld.csv',data.table=F)
df_r2[,c(5,7,8)]<-NULL
df_sign<-fread('results/sign_locus.txt',data.table=F)
df_sign[,c(5,6,7,9,10,11,13,14,21,22,23,25,30,31)]<-NULL

RS<-unique(df_sign$rs)

DF<-list()
for(i in 1:length(RS)){
x<-subset(df_r2,df_r2$rs1==RS[i] | df_r2$rs2==RS[i])
x[,c(1,2,3,4,6,7,8,9,13)]<-NULL
x<-unique(x)
y<-subset(df_sign,df_sign$rs==RS[i])
DF[[i]]<-x%>%group_by(sp,Locus)%>%summarise(min_ld=min(r2),max_ld=max(r2),ave_ld=mean(r2),n_snp=n())
DF[[i]]$rs<-unique(y$rs)
DF[[i]]$Locus_ini<-unique(y$Locus)
}
DF<-do.call(rbind,DF)
a<-colsplit(DF$rs,':',c('chr','ps'))
DF$chr<-a$chr
DF$ps<-a$ps

write.table(DF,'gene_ld.txt',col.names=T,row.names=F,quote=F)

df<-subset(DF,DF$Locus!=DF$Locus_ini) 
