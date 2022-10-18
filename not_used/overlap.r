library(data.table)
annot<-fread('data/proteins_48_403979.csv',data.table=F)
annot<-unique(annot[,c(1,2,3,4,6)])
annot<-subset(annot,annot[,1]!='Un')

annot$chr[annot[,1]=='linkage group LG1']<-1
annot$chr[annot[,1]=='linkage group LG2']<-2
annot$chr[annot[,1]=='linkage group LG3']<-3
annot$chr[annot[,1]=='linkage group LG4']<-4
annot$chr[annot[,1]=='linkage group LG5']<-5
annot$chr[annot[,1]=='linkage group LG6']<-6
annot$chr[annot[,1]=='linkage group LG7']<-7
annot$chr[annot[,1]=='linkage group LG8']<-8
annot$chr[annot[,1]=='linkage group LG9']<-9
annot$chr[annot[,1]=='linkage group LG10']<-10
annot$chr[annot[,1]=='linkage group LG11']<-11
annot$chr[annot[,1]=='linkage group LG12']<-12
annot$chr[annot[,1]=='linkage group LG13']<-13
annot$chr[annot[,1]=='linkage group LG14']<-14
annot$chr[annot[,1]=='linkage group LG15']<-15
annot$chr[annot[,1]=='linkage group LG16']<-16

df<-data.frame('chr'=c(1,2,3,4,4,6,8,10,11,13,14,14),
'start'=c(14780956,3846716,1065188,3115955,10289077,9826068,9057205,1464362,8903672,1823813,2783650,5325863),
'stop'=c(16624843,4937645,3074728,4735881,12165460,11855354,10299408,2526877,10027267,3144881,4475007,6711697),
'pheno'=c('recap_pc1','recap_mnr','recap_pc1','recap_mnr','recap_pc1','recap_mnr','mnr_pc1','recap_mnr','recap_mnr','recap_mnr','recap_mnr','recap_mnr'))

Ann<-list()
for(i in 1:nrow(df)){
Ann[[i]]<-data.frame('pheno'=df$pheno[i],'chr'=df$chr[i],'Locus'=unique(annot[annot$chr==df$chr[i] & annot$Start>=df$start[i] & annot$Stop<=df$stop[i],'Locus']))
}
Ann<-do.call(rbind,Ann)
write.table(Ann,'overlap_regions.txt',col.names=T,row.names=F,quote=F)

