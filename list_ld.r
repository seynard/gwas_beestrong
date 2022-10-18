#module load system/R-3.5.1
#R
library(data.table)
args<-commandArgs(TRUE)
dir_in<-args[1]
dir_out<-args[2]
pop_id=args[3]
pop_id=unlist(strsplit(pop_id,','))
thresh<-as.numeric(args[4])

q<-fread(paste0(dir_in,'/seqapipop_sample.3.Q'),data.table=F)
f<-fread(paste0(dir_in,'/seqapipop_sample.fam'),data.table=F)
df<-cbind(f,q)

df[,c(2,3,4,5,6)]<-NULL
colnames(df)<-c('id','pop1','pop2','pop3')
pop<-vector()
for(i in 1:3){
x<-df[df[,i+1]>0.99,]
if('Ab-PacBio'%in%x$id){pop[i]<-'Mel'}else if('CAU10'%in%x$id){pop[i]<-'Cau'}else if('ITA1A'%in%x$id){pop[i]<-'Lig'}
}
colnames(df)<-c('id',pop)

for(i in 1:length(pop_id)){
	if(pop_id[i]=='Mellifera'){
	x<-df$id[df$Mel>thresh]
	x<-cbind(x,x)
	write.table(x,paste0(dir_out,'/list_Mellifera_ld.txt'),col.names=F,row.names=F,quote=F)
	}else if(pop_id[i]=='Ligustica_Carnica'){
	x<-df$id[df$Lig>thresh]
	x<-cbind(x,x)
	write.table(x,paste0(dir_out,'/list_Ligustica_Carnica_ld.txt'),col.names=F,row.names=F,quote=F)
	}else if(pop_id[i]=='Caucasica'){
	x<-df$id[df$Cau>thresh]
	x<-cbind(x,x)
	write.table(x,paste0(dir_out,'/list_Caucasica_ld.txt'),col.names=F,row.names=F,quote=F)
	}else{
	x<-df$id[df$Mel<=thresh & df$Lig<=thresh & df$Cau<=thresh]
	x<-cbind(x,x)
	write.table(x,paste0(dir_out,'/list_hybrid_ld.txt'),col.names=F,row.names=F,quote=F)
	}
}

