#module load system/R-3.5.1
#R
#dir_out='/work/project/dynagen/seynard/GWAS/BeeStrongHAV3_1'
args<-commandArgs(TRUE)
dir_out<-args[1]
type<-args[2]

dat<-read.table(paste0(dir_out,'/',type),sep=' ',header=F,fill=T)
n<-which(duplicated(dat[,2]))
dat[,2]<-as.character(dat[,2])
if(length(n)>0){
for (i in 1:length(n)){
	dat[n[i],2]<-paste0(dat[n[i],2],'_bis')
	}
}
write.table(dat,paste0(dir_out,'/',type),sep=' ',col.names=F,row.names=F,quote=F)
