library(data.table)
library(dplyr)
args<-commandArgs(TRUE)
dirin=args[1]
dirout=args[2]
geno_in=args[3]
snp=args[4]
chrom=as.numeric(args[5])

#dirin='/work/project/dynagen/seynard/GWAS/BeeStrongHAV3_1'
#dirout='/work/project/dynagen/seynard/GWAS/BeeStrongHAV3_1/ldak'
#geno_in='geno_hom_Mellifera'
#snp='allele_id_final.txt'
#chrom=1

print(chrom)
snp<-fread(paste0(dirin,'/',snp),data.table=F)
colnames(snp)<-c('CHROM','POS','REF','ALT')
snp<-subset(snp,snp$CHROM==chrom)

ped<-fread(paste0(dirin,'/',geno_in,'_',chrom,'.txt'),data.table=F)
if(nrow(ped)==nrow(snp)){
ped$x0<-paste0(snp$ALT,' ',snp$ALT)
ped$x1<-paste0(snp$REF,' ',snp$ALT)
ped$x2<-paste0(snp$REF,' ',snp$REF)
}else{print('not same length')}

ped_geno<-ped[,3:ncol(ped)]
a<-ncol(ped_geno)-2
b<-ncol(ped_geno)-1
c<-ncol(ped_geno)
ped_geno<- t(apply(ped_geno, 1, function(x) {
                    x[x=="0"] <- x[a]
                    return(x)
                 }))
ped_geno<- t(apply(ped_geno, 1, function(x) {
                    x[x == "1"] <- x[b]
                    return(x)
                 }))
ped_geno<- t(apply(ped_geno, 1, function(x) {
                    x[x == "2"] <- x[c]
                    return(x)
                 }))
ped_geno<-as.data.frame(ped_geno)
ped_geno[,c('x0','x1','x2')]<-NULL
ped<-data.frame(ped[,1],ped[,2],ped_geno)
ped$rs<-paste0(ped[,1],':',ped[,2])
x<-data.frame('chr'=ped[,1],'rs'=ped$rs,'cm'=0,'pos'=ped[,2],ped_geno)
write.table(x,paste0(dirin,'/',geno_in,'_',chrom,'.tped'),col.names=F,row.names=F,quote=F)

if(chrom==max(snp$CHROM)){
	con <- file(paste0(dirin,'/',geno_in,'.txt'),"r")
	first_line <- readLines(con,n=1)
	close(con)
	first_line<-unlist(strsplit(first_line,' '))
	n<-c(grep('CHROM',first_line),grep('POS',first_line),grep('REF',first_line),grep('ALT',first_line))
	first_line<-first_line[-n]
	fam<-data.frame('v1'=1,'v2'=first_line,'v3'=0,'v4'=0,'v5'=0,'v6'=-9)
	write.table(fam,paste0(dirin,'/',geno_in,'.tfam'),col.names=F,row.names=F,quote=F)
	}
	

	
	
