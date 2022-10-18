library(data.table)
args<-commandArgs(TRUE)
dir<-args[1]
qg<-args[2]
pop_id2<-args[3]
pop_id2<-unlist(strsplit(pop_id2,','))
pop_id<-args[4]
pop_id<-unlist(strsplit(pop_id,','))
choice=args[7]

threshold<-as.numeric(args[5])
q_matrix<-fread(paste0(dir,'/',qg,'/queen_geno',choice,'.Q'),data.table=F)
pheno<-args[6]
pheno<-fread(paste0(dir,'/',pheno),data.table=F)
pheno$num_ruche_bs<-toupper(pheno$num_ruche_bs)
q_matrix<-merge(q_matrix,pheno,by='num_ruche_bs',all=T)
q_matrix<-q_matrix[,c('num_ruche_bs','Ligustica_Carnica','Mellifera','Caucasica','structure','pass')]
q_matrix$pass[is.na(q_matrix$pass)]<-'yes'
q_matrix<-unique(q_matrix)
q_matrix<-subset(q_matrix,!is.na(q_matrix$Ligustica_Carnica))

#for(i in 1:length(pop_id2)){
#if(pop_id2[i]%in%pop_id){
#	list<-subset(q_matrix[,1],q_matrix[,pop_id2[i]]>=threshold & q_matrix$structure!='USDA')
#	write.table(list,paste0(dir,'/list_',pop_id2[i],'.txt'),col.names=F,row.names=F,quote=F)
#}else if (pop_id2[i]=='hybrid'){
#	list<-subset(q_matrix[,1],apply(q_matrix[,pop_id],1,max)<threshold & q_matrix$structure!='USDA')
#	write.table(list,paste0(dir,'/list_',pop_id2[i],'.txt'),col.names=F,row.names=F,quote=F)
#}else if (pop_id2[i]=='all'){
#	list<-q_matrix[,1]
#	write.table(list,paste0(dir,'/list_',pop_id2[i],'.txt'),col.names=F,row.names=F,quote=F)
#}else if (pop_id2[i]=='all_nus'){
#	list<-subset(q_matrix[,1],q_matrix$structure!='USDA')
#	write.table(list,paste0(dir,'/list_',pop_id2[i],'.txt'),col.names=F,row.names=F,quote=F)
#}else if (pop_id2[i]=='us'){
#	list<-subset(q_matrix[,1],q_matrix$structure=='USDA')
#	write.table(list,paste0(dir,'/list_',pop_id2[i],'.txt'),col.names=F,row.names=F,quote=F)
#}else if (pop_id2[i]=='corse'){
#	list<-subset(q_matrix[,1],q_matrix$structure=='AOP_Corse')
#	write.table(list,paste0(dir,'/list_',pop_id2[i],'.txt'),col.names=F,row.names=F,quote=F)
#}}

for(i in 1:length(pop_id2)){
if(pop_id2[i]%in%pop_id){
	list<-subset(q_matrix[,1],q_matrix[,pop_id2[i]]>=threshold)
	write.table(list,paste0(dir,'/list_',pop_id2[i],'.txt'),col.names=F,row.names=F,quote=F)
}else if (pop_id2[i]=='hybrid'){
	list<-subset(q_matrix[,1],apply(q_matrix[,pop_id],1,max)<threshold)
	write.table(list,paste0(dir,'/list_',pop_id2[i],'.txt'),col.names=F,row.names=F,quote=F)
}else if (pop_id2[i]=='Ligustica_nus'){
	list<-subset(q_matrix[,1],q_matrix[,'Ligustica_Carnica']>=threshold & q_matrix$structure!='USDA')
	write.table(list,paste0(dir,'/list_',pop_id2[i],'.txt'),col.names=F,row.names=F,quote=F)
}
}
