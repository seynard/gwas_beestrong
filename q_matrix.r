library(tidyverse)
library(data.table)
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
