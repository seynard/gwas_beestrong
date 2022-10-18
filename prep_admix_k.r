#module load system/R-3.5.1
#R
args<-commandArgs(TRUE)
dir<-args[1]
#dir='/work/project/dynagen/seynard/GWAS/BeeStrongHAV3_1'
data_type=args[2]
pop_id=args[3]
npop=length(unlist(strsplit(pop_id,',')))
pop_id=unlist(strsplit(pop_id,','))

library(data.table)
library(reshape2)
library(dplyr)

fam<-fread(paste0(dir,'/',data_type,'.fam'),header=F,data.table=F)
admix_list<-intersect(list.files(path=dir,pattern=data_type),list.files(path=dir,pattern=".Q"))
admix_list<-admix_list[grep(paste0('\\b',data_type,'\\b'),admix_list)]
x<-	grep(paste0(data_type,'.',npop,'.Q'),admix_list)
q<-fread(paste0(dir,'/',admix_list[x]),header=F,data.table=F)

seuil<-as.numeric(args[4])
dat<-cbind(fam,q)
npop<-ncol(dat)-6
popn<-paste0('pop',1:npop)
colnames(dat)<-c('ind','ind_id2','father','mother','sex','pheno',popn)
j1<-max.col(dat[,7:ncol(dat)],"first")
dat$pop<-popn[j1]
for(i in 1:npop){
if(pop_id[i]=='Mellifera'){
	pop_mel<-unique(c(dat$pop[grep('OUE',dat$ind)],dat$pop[grep('UK',dat$ind)],dat$pop[grep('ESP',dat$ind)],dat$pop[grep('PacBio',dat$ind)]))
	colnames(dat)[colnames(dat)==pop_mel]<-'Mellifera'
	dat$pop[dat$pop==pop_mel]<-'Mellifera'
}else if (pop_id[i]=='Ligustica_Carnica'){
	pop_lig<-unique(dat$pop[grep('ITA',dat$ind)])
	colnames(dat)[colnames(dat)==pop_lig]<-'Ligustica_Carnica'
	dat$pop[dat$pop==pop_lig]<-'Ligustica_Carnica'	
}else if(pop_id[i]=='Caucasica'){
	pop_cau<-table(subset(dat$pop,grepl('CAU',dat$ind)))
	pop_cau<-names(which.max(pop_cau))
	colnames(dat)[colnames(dat)==pop_cau]<-'Caucasica'
	dat$pop[dat$pop==pop_cau]<-'Caucasica'
}}

ColNames<-colnames(dat[,7:(ncol(dat)-1)])
write.table(ColNames,paste0(dir,'/ColNames.txt'),col.names=F,row.names=F,quote=F)

dat_thresholdi<-list()
for(i in 1:length(pop_id)){dat_thresholdi[[i]]<-subset(dat,dat[,pop_id[i]]>seuil)}
	
dat_thresholdi<-do.call(rbind,dat_thresholdi)
write.table(dat_thresholdi[,1:6],paste0(dir,'/Unif_k.fam'),col.names=F,row.names=F,quote=F)

dat_thresholdi$ind_id2<-dat_thresholdi$pop
write.table(dat_thresholdi[,c('ind_id2','ind')],paste0(dir,'/key_k.fam'),col.names=F,row.names=F,quote=F)

write.table(dat_thresholdi[,c(2,1,3,4,5,6)],paste0(dir,'/Unif_k_pop.fam'),col.names=F,row.names=F,quote=F)
