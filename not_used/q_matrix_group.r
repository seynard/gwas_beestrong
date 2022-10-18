library(tidyverse)
library(data.table)
library(ggplot2)

l_lig<-fread('data/list_Ligustica_Carnica.txt',header=F,data.table=F)
l_mel<-fread('data/list_Mellifera_ncorse.txt',header=F,data.table=F)
l_hyb<-fread('data/list_hybrid_corse.txt',header=F,data.table=F)
l_cau<-fread('data/list_Caucasica.txt',header=F,data.table=F)

pheno<-fread('data/Data_BS.csv')
pheno$group<-NA
for(i in 1:nrow(pheno)){
	if(pheno$num_ruche_bs[i]%in%tolower(l_lig[,1])){
	pheno[i,'group']<-'lig'
	}else if (pheno$num_ruche_bs[i]%in%tolower(l_mel[,1])){
	pheno[i,'group']<-'mel'
	}else if (pheno$num_ruche_bs[i]%in%tolower(l_hyb[,1])){
	pheno[i,'group']<-'hyb'
	}else if (pheno$num_ruche_bs[i]%in%tolower(l_cau[,1])){
	pheno[i,'group']<-'cau'
	}
	}

q_mat<-fread('data/queen_geno.Q',data.table=F)
q_mat$group<-NA
for(i in 1:nrow(q_mat)){
	if(q_mat$num_ruche_bs[i]%in%l_lig[,1]){
	q_mat[i,'group']<-'lig'
	}else if (q_mat$num_ruche_bs[i]%in%l_mel[,1]){
	q_mat[i,'group']<-'mel'
	}else if (q_mat$num_ruche_bs[i]%in%l_hyb[,1]){
	q_mat[i,'group']<-'hyb'
	}else if (q_mat$num_ruche_bs[i]%in%l_cau[,1]){
	q_mat[i,'group']<-'cau'
	}
	}
q_mat<-gather(q_mat,sp,val,c('Mellifera','Ligustica_Carnica','Caucasica'))
q_mat$sp<-factor(q_mat$sp,levels=c('Ligustica_Carnica','Mellifera','Caucasica'))
q_mat$group[is.na(q_mat$group)]<-'no_group'
q_mat$group<-factor(q_mat$group,levels=c('lig','mel','cau','hyb','no_group'))

ggplot(q_mat)+
	geom_bar(aes(x=num_ruche_bs,y=val,fill=sp),stat='identity',position='stack')+
	scale_fill_manual(values=c('goldenrod2','grey34','chartreuse2'))+
	facet_grid(.~group,scales='free')

q_mat50<-fread('data/queen_geno_50k.Q',data.table=F)
q_mat50$group<-NA
for(i in 1:nrow(q_mat50)){
	if(q_mat50$num_ruche_bs[i]%in%l_lig[,1]){
	q_mat50[i,'group']<-'lig'
	}else if (q_mat50$num_ruche_bs[i]%in%l_mel[,1]){
	q_mat50[i,'group']<-'mel'
	}else if (q_mat50$num_ruche_bs[i]%in%l_hyb[,1]){
	q_mat50[i,'group']<-'hyb'
	}else if (q_mat50$num_ruche_bs[i]%in%l_cau[,1]){
	q_mat50[i,'group']<-'cau'
	}
	}
q_mat50<-gather(q_mat50,sp,val,c('Mellifera','Ligustica_Carnica','Caucasica'))
q_mat50$sp<-factor(q_mat50$sp,levels=c('Ligustica_Carnica','Mellifera','Caucasica'))
q_mat50$group[is.na(q_mat50$group)]<-'no_group'
q_mat50$group<-factor(q_mat50$group,levels=c('lig','mel','cau','hyb','no_group'))

ggplot(q_mat50)+
	geom_bar(aes(x=num_ruche_bs,y=val,fill=sp),stat='identity',position='stack')+
	scale_fill_manual(values=c('goldenrod2','grey34','chartreuse2'))+
	facet_grid(.~group,scales='free')


ggplot(q_mat)+
	geom_boxplot(aes(x=group,y=val,fill=sp))+
	scale_fill_manual(values=c('goldenrod2','grey34','chartreuse2'))
	
ggplot(q_mat50)+
	geom_boxplot(aes(x=group,y=val,fill=sp))+
	scale_fill_manual(values=c('goldenrod2','grey34','chartreuse2'))






args<-commandArgs(TRUE)
dir<-args[1]
n_pop<-as.numeric(args[2])
data<-args[3]
pop_id<-args[4]
pop_id<-unlist(strsplit(pop_id,','))
choice<-args[5]
dat<-list.files(path=dir,pattern=paste0(choice,'_st_het.Q'))

Q_matrix<-matrix(ncol=n_pop+1,nrow=length(dat))
for(i in 1:length(dat)){
Q_matrix[i,1]<-gsub(paste0(choice,'_st_het.Q'),'',gsub('sim_depth_count','',dat[i],'_'))
q<-fread(paste0(dir,'/',dat[i]),data.table=F)
for(j in 1:length(pop_id)){
	Q_matrix[i,j+1]<-q[,2][q[,1]==pop_id[j]]
	}
}
Q_matrix<-as.data.frame(Q_matrix)
colnames(Q_matrix)<-c('num_ruche_bs',pop_id)
write.table(Q_matrix,paste0(dir,'/queen_geno',choice,'.Q'),col.names=T,row.names=F,quote=F)

dat<-fread(paste0(dir,'/',data),data.table=F)
dat$num_ruche_bs<-toupper(dat$num_ruche_bs)
D<-merge(dat,Q_matrix,by='num_ruche_bs')
pdf(paste0(dir,'/plot_q.pdf'),width=50,height=15)
for(i in 1:length(unique(D$structure))){
	dx<-subset(D,D$structure==unique(D$structure)[i])
	dx<-dx[,c('num_ruche_bs',all_of(pop_id))]
	dx<-dx[!duplicated(dx),]
	rownames(dx)<-dx$num_ruche_bs
	dx$num_ruche_bs<-NULL
	colnames(dx)[colnames(dx)=='Mellifera']<-'black'
	colnames(dx)[colnames(dx)=='Caucasica']<-'chartreuse4'
	colnames(dx)[colnames(dx)=='Ligustica_Carnica']<-'gold'
	par(mar=c(7,3,2,2))
	barplot(t(dx),col=colnames(dx),border="white",xlab="",las=2)
	title(main=unique(D$structure)[i])}
dev.off()
