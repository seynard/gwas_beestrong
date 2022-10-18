#module load system/R-4.1.1_gcc-9.3.0
#R
library(data.table)
library(genetics)
library(reshape2)

args<-commandArgs(TRUE)
dirin<-args[1]
dirout<-args[2]
sp<-args[3]
type<-args[4]
chr<-as.numeric(args[5])
#dirin<-'data'
#dirout<-'results'
#sp<-'Mellifera_ncorse'
#type<-'egs'
#chr<-14

#R2<-data.frame('rs1'=NA,'rs2'=NA,'r'=NA,'r2'=NA,'dprime'=NA)
#write.table(R2,paste0(dirout,'/ld_',sp,'_',chr,'_',type,'.txt'),col.names=T,row.names=F,quote=F)

geno<-fread(paste0(dirin,'/',type,'2_',sp,'.txt'),data.table=F)
genoc<-geno[geno$CHROM==chr,]
genoc[,5:ncol(genoc)]<-round(genoc[,5:ncol(genoc)],digit=0)

Geno<-matrix(nrow=nrow(genoc),ncol=ncol(genoc))
#Geno[,1]<-colnames(genoc)[5:ncol(genoc)]
#Geno[,2]<-colnames(genoc)[5:ncol(genoc)]
#Geno[,3]<-0
#Geno[,4]<-0
#Geno[,5]<-0
#Geno[,6]<--9
for(i in 1:nrow(genoc)){
print(i)
x<-genoc[i,5:ncol(genoc)]
ref<-genoc[i,'REF']
alt<-genoc[i,'ALT']
x[x==0]<-paste0(ref,' ',ref)
x[x==1]<-paste0(ref,' ',alt)
x[x==2]<-paste0(alt,' ',alt)
x[is.na(x)]<-'0 0'
x<-c(genoc[i,'CHROM'],paste0(genoc[i,'CHROM'],':',genoc[i,'POS']),0,genoc[i,'POS'],unlist(x))
#Geno[,6+i]<-t(x)
Geno[i,]<-x
}
write.table(Geno,file=paste0(dirout,'/',type,'_',sp,'_',chr,'.tped'),sep=" ",row.names=F,col.names=F,quote=F)

Geno_fam<-data.frame(colnames(genoc)[5:ncol(genoc)],colnames(genoc)[5:ncol(genoc)],0,0,0,-9)
write.table(Geno_fam,file=paste0(dirout,'/',type,'_',sp,'_',chr,'.tfam'),sep=" ",row.names=F,col.names=F,quote=F)

#n<-1
#while(n<nrow(genoc)){
#	m<-2
#	while(m<nrow(genoc)){
#		if(n!=m){
#		genox<-genoc[c(n,m),]
#		rs<-paste0(genox$CHROM,':',genox$POS)
#		genox[,c('CHROM','POS','REF','ALT')]<-NULL
#		genox<-round(genox,digit=0)
#		genox<-t(genox)
#		colnames(genox)<-rs
#		rownames(genox)<-NULL
#		data<-makeGenotypes(genox,convert=colnames(genox),method=as.genotype.allele.count)
#		if(length(unique(data[,1][!is.na(data[,1])]))>1 & length(unique(data[,2][!is.na(data[,2])]))>1){
#			ld<-LD(data)
#			r<-ld["r"]
#			r<-melt(r)
#			r$L1<-NULL
#			r2<-ld["R^2"]
#			r2<-melt(r2)
#			r2$L1<-NULL
#			dprime<-ld["D'"]
#			dprime<-melt(dprime)
#			dprime$L1<-NULL
#			x<-Reduce(function(x,y) merge(x,y,by=c('Var1','Var2'),all=T),list(r,r2,dprime))
#			x<-subset(x,!is.na(x[,4]) & x[,4]>=0.1)
##			print(x)
#			write.table(x,file=paste0(dirout,'/ld_',sp,'_',chr,'_',type,'.txt'),append=TRUE,sep=" ",row.names=FALSE,col.names=F,quote=F)
#		}}
#		m<-m+1
#	}
#	n<-n+1
#}
