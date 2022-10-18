#module load system/R-3.5.1
#R
#dir='/work/project/dynagen/seynard/GWAS/bam_BeeStrongMOSAR'
args<-commandArgs(TRUE)
dir<-args[1]

library(data.table)
library(gdata)
library(ggfortify)
library(gridExtra)
library(grid)
library(gplots)
library(reshape2)
library(ggdendro)
library(ggcorrplot)
library(factoextra)
library(FactoMineR)
library(corrplot)
library(RColorBrewer)
library(viridis)
library(cluster)
library(tidyr)
library(NbClust)
#library(lga)
library(clValid)
#library(bclust)
library(mclust)
library(cowplot)
library(ape)
library(readr)

data_type=args[2]
pop=args[3]

ind<-fread(paste0(dir,'/',data_type,'.fam'),header=F,data.table=F)
colnames(ind)[1]<-'ind'
Ind<-ind[,1]

admix_list<-intersect(list.files(path=dir,pattern=data_type),list.files(path=dir,pattern=".Q"))
admix_list<-admix_list[grep(paste0('\\b',data_type,'\\b'),admix_list)]
log_list<-intersect(list.files(path=paste0(dir,'/log/'),pattern=data_type),list.files(path=paste0(dir,'/log/'),pattern = ".out"))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colour=sample(col_vector,20)

for(i in 1:length(admix_list)){
	admix<-fread(paste0(dir,'/',admix_list[i]),header=F,data.table=F)
	colnames(admix)<-seq(1:ncol(admix))
	admix<-cbind(Ind,admix)
	admix$pop<-NA
	admix$pop[is.na(admix$pop)]<-'-'
	admix$pop[admix$pop=='-' & grepl('AOC',admix$Ind)]<-'-_AOC'
	admix$pop[admix$pop=='-' & grepl('BR',admix$Ind)]<-'-_BR'
	admix$pop[admix$pop=='-' & grepl('BER',admix$Ind)]<-'-_BER'
	admix$pop[admix$pop=='-' & grepl('Buck',admix$Ind)]<-'-_Buck'
	admix$pop[admix$pop=='-' & grepl('BS',admix$Ind)]<-'-_BS'
	admix$pop[admix$pop=='-' & grepl('CIReine',admix$Ind)]<-'-_CIreine'
	admix$pop[admix$pop=='-' & grepl('ITSAP',admix$Ind)]<-'-_ITSAP'
	admix$pop[admix$pop=='-' & grepl('JFM',admix$Ind)]<-'-_JFM'
	admix$pop[admix$pop=='-' & grepl('Sav',admix$Ind)]<-'-_Savoie'
	admix$pop[admix$pop=='-' & grepl('Sar',admix$Ind)]<-'-_YLC'
	admix$pop[admix$pop=='-' & grepl('YC',admix$Ind)]<-'-_YLC'
	admix$pop[admix$pop=='-' & grepl('JSC',admix$Ind)]<-'-_JSC'
	admix$pop[admix$pop=='-' & grepl('KF',admix$Ind)]<-'-_KF'
	admix$pop[admix$pop=='-' & grepl('PacBio',admix$Ind)]<-'-_OUE'
	admix$pop[admix$pop=='-' & grepl('OUE',admix$Ind)]<-'-_OUE'
	admix$pop[admix$pop=='-' & grepl('CAR',admix$Ind)]<-'-_CAR'
	admix$pop[admix$pop=='-' & grepl('CAU',admix$Ind)]<-'-_CAU'
	admix$pop[admix$pop=='-' & grepl('CHI',admix$Ind)]<-'-_CHI'
	admix$pop[admix$pop=='-' & grepl('ESP',admix$Ind)]<-'-_ESP'
	admix$pop[admix$pop=='-' & grepl('ITA',admix$Ind)]<-'-_ITA'
	admix$pop[admix$pop=='-' & grepl('NM',admix$Ind)]<-'-_NM'
	admix$pop[admix$pop=='-' & grepl('POL',admix$Ind)]<-'-_POL'
	admix$pop[admix$pop=='-' & grepl('POR',admix$Ind)]<-'-_POR'
	admix$pop[admix$pop=='-' & grepl('SLO',admix$Ind)]<-'-_SLO'
	admix$pop[admix$pop=='-' & grepl('SOL',admix$Ind)]<-'-_SOL'
	admix$pop[admix$pop=='-' & grepl('UK',admix$Ind)]<-'-_UK'
	admix$pop[admix$pop=='-' & grepl('FL',admix$Ind)]<-'-_FL'
	k_test<-as.numeric(unlist(strsplit(admix_list[i],"\\."))[2])
	plot_name<-paste0('/plot_',data_type,'_',k_test,'.pdf')	
	logL<-log_list[grep(paste0('_',k_test,'.out'),log_list)]
	print(plot_name)
	pdf(paste0(dir,plot_name),width=30,height=10)
	group<-unique(admix$pop)
	for(j in 1:length(group)){
	dat<-subset(admix,admix$pop==group[j])
	rownames(dat)<-dat$Ind
	dat[,c('Ind','pop')]<-NULL
	barplot(t(dat),col=colour[1:k_test],border="white",xlab="",las=3,main=group[j])
	}
	if(k_test>2){
	log<-readLines(paste0(dir,'/log/',logL))	
	fst<-as.numeric(grep('Fst',log))
	Fst<-log[(fst+1):(fst+1+k_test)]
	table_fst<-matrix(nrow=(k_test),ncol=(k_test))
	diag(table_fst)<-1
	for(j in 1:(length(Fst)-2)){
		fst_val<-unlist(strsplit(Fst[j+2],'\t'))[-1]
		for(m in 1:length(fst_val)){table_fst[j+1,m]<-as.numeric(fst_val[m])}}
	par(mfrow=c(1,2))
	dd<-as.dist(table_fst, diag = FALSE, upper = FALSE)
	clust_d<-nj(dd)
	plot(unroot(clust_d),type="unrooted",no.margin=TRUE,lab4ut="axial",edge.width=2,tip.color=colour[1:k_test],cex=5)	
	q<-admix[,c(1:k_test,'Ind')]
	clust=nj(dist(q[,1:k_test]))
	k<-vector()
	for(j in 1:nrow(q)){k[j]<-which.max(q[j,1:k_test])}
	q$k<-k
	for(j in 1:max(k)){q$k[q$k==j]<-as.character(colour[j])}
	clust$tip.label<-as.character(q$Ind)
	plot(unroot(clust),type="unrooted",no.margin=TRUE,lab4ut="axial",edge.width=2,tip.color=colour[1:k_test],cex=1)}
	dev.off()
	}



