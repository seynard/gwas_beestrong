library(data.table)
library(tidyverse)
options(warn=2) 

# individual analysis
df<-fread('sign_gwas_ind.txt',data.table=F,fill=T)
df$sp[df$sp=='Mellifera']<-'mel'
df$sp[df$sp=='Hybrids']<-'hyb'
df$sp[df$sp=='Ligustica' & df$V34=='']<-'lig'
df$sp[df$sp=='Ligustica' & df$V34=='us)']<-'lig_nus'
df[,c('V31','V32','V33','V34')]<-NULL
df$id<-paste0(df$rs,'_',df$sp)
ind<-df
vep<-fread('ind_vep.txt',data.table=F)
vep<-unique(vep)
colnames(vep)[1]<-'id'
X<-merge(ind,vep,by='id',all=T)
DF<-list()
for(i in 1:length(unique(X$rs))){
	xi<-X[X$rs==unique(X$rs)[i],]
	if('-'%in%unique(xi$DISTANCE)){xi<-xi[xi$DISTANCE=='-',]
	}else{min_d<-min(as.numeric(xi$DISTANCE));xi<-xi[as.numeric(xi$DISTANCE)==min_d,]}
	xi<-unique(xi[,c("rs","pheno","gwas","sp","Consequence","IMPACT","SYMBOL","Gene","Feature_type","BIOTYPE","DISTANCE","STRAND")])
	DF[[i]]<-xi
	}
DF_ind<-do.call(rbind,DF)
colnames(DF_ind)[3]<-'type' 

# meta GWAS
df<-fread('results/sign_locus_1.txt',data.table=F,fill=T)
df$id<-df$rs
meta1<-df
vep<-fread('meta1_vep.txt',data.table=F)
vep<-unique(vep)
colnames(vep)[1]<-'id'
X<-merge(meta1,vep,by='id',all=T)
DF<-list()
for(i in 1:length(unique(X$rs))){
	xi<-X[X$rs==unique(X$rs)[i],]
	if('-'%in%unique(xi$DISTANCE)){xi<-xi[xi$DISTANCE=='-',]
	}else{min_d<-min(as.numeric(xi$DISTANCE));xi<-xi[as.numeric(xi$DISTANCE)==min_d,]}
	xi<-unique(xi[,c("rs","pheno","type","Consequence","IMPACT","SYMBOL","Gene","Feature_type","BIOTYPE","DISTANCE","STRAND")])
	DF[[i]]<-xi
	}
DF_meta1<-do.call(rbind,DF)
DF_meta1$sp<-'meta1'

df<-fread('results/sign_locus_2.txt',data.table=F,fill=T)
df$id<-df$rs
meta2<-df
vep<-fread('meta2_vep.txt',data.table=F)
vep<-unique(vep)
colnames(vep)[1]<-'id'
X<-merge(meta2,vep,by='id',all=T)
DF<-list()
for(i in 1:length(unique(X$rs))){
	xi<-X[X$rs==unique(X$rs)[i],]
	if('-'%in%unique(xi$DISTANCE)){xi<-xi[xi$DISTANCE=='-',]
	}else{min_d<-min(as.numeric(xi$DISTANCE));xi<-xi[as.numeric(xi$DISTANCE)==min_d,]}
	xi<-unique(xi[,c("rs","pheno","type","Consequence","IMPACT","SYMBOL","Gene","Feature_type","BIOTYPE","DISTANCE","STRAND")])
	DF[[i]]<-xi
	}
DF_meta2<-do.call(rbind,DF)
DF_meta2$sp<-'meta2'

GENES<-do.call(rbind,list(DF_ind,DF_meta1,DF_meta2))
write.table(GENES,'locus_ind_meta.txt',col.names=T,row.names=F,quote=F)





















vep<-fread('ind_vep.txt',data.table=F)
vep<-vep[,c("#Uploaded_variation","Location","Allele","Consequence","IMPACT","SYMBOL","Gene","BIOTYPE","DISTANCE","STRAND")]
vep<-unique(vep)
X<-list()
for(i in 1:nrow(df)){
	df_i<-df[i,]
	vep_i<-subset(vep,vep[,1]==df_i$id)
	vep_i$DISTANCE[vep_i$DISTANCE=='-']<-0
	vep_i$DISTANCE<-as.numeric(vep_i$DISTANCE)
	if(min(vep_i$DISTANCE)==0){n<-which(vep_i$DISTANCE==0)}else{n<-which.min(vep_i$DISTANCE)}
	vep_i<-vep_i[n,]
	X[[i]]<-cbind(df_i,vep_i)
}
X<-do.call(rbind,X)
write.table(X,'ind_snp_vep.txt',col.names=T,row.names=F,quote=F)


vep<-fread('meta_vep_1.txt',data.table=F)
vep<-vep[,c("#Uploaded_variation","Location","Allele","Consequence","IMPACT","SYMBOL","Gene","BIOTYPE","DISTANCE","STRAND")]
vep<-unique(vep)

X<-list()
for(i in 1:nrow(df)){
	df_i<-df[i,]
	vep_i<-subset(vep,vep[,1]==df_i$id)
	vep_i$DISTANCE[vep_i$DISTANCE=='-']<-0
	vep_i$DISTANCE<-as.numeric(vep_i$DISTANCE)
	if(min(vep_i$DISTANCE)==0){n<-which(vep_i$DISTANCE==0)}else{n<-which.min(vep_i$DISTANCE)}
	vep_i<-vep_i[n,]
	X[[i]]<-cbind(df_i,vep_i)
}
X<-do.call(rbind,X)
write.table(X,'meta_snp_vep_1.txt',col.names=T,row.names=F,quote=F)


vep<-fread('meta_vep_2.txt',data.table=F)
vep<-vep[,c("#Uploaded_variation","Location","Allele","Consequence","IMPACT","SYMBOL","Gene","BIOTYPE","DISTANCE","STRAND")]
vep<-unique(vep)

X<-list()
for(i in 1:nrow(df)){
	df_i<-df[i,]
	vep_i<-subset(vep,vep[,1]==df_i$id)
	vep_i$DISTANCE[vep_i$DISTANCE=='-']<-0
	vep_i$DISTANCE<-as.numeric(vep_i$DISTANCE)
	if(min(vep_i$DISTANCE)==0){n<-which(vep_i$DISTANCE==0)}else{n<-which.min(vep_i$DISTANCE)}
	vep_i<-vep_i[n,]
	X[[i]]<-cbind(df_i,vep_i)
}
X<-do.call(rbind,X)
write.table(X,'meta_snp_vep_2.txt',col.names=T,row.names=F,quote=F)
