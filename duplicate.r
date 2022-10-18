library('data.table')
library('dplyr')
args<-commandArgs(TRUE)
dir<-args[1]
#dir<-'work/project/dynagen/seynard/BeeStrongHAV3_1'
test<-args[2]
test<-unlist(strsplit(test,'_'))
print(test)

if('replicate'%in%test){
	d_replicate<-fread(paste0(dir,'/depth_replicate.txt'),sep=' ',header=T,data.table=F)
	c_replicate<-fread(paste0(dir,'/count_ref_replicate.txt'),sep=' ',header=T,data.table=F)
	F_replicate<-c_replicate[,5:ncol(c_replicate)]/d_replicate[,5:ncol(d_replicate)]
	correl_replicate<-matrix(ncol=ncol(F_replicate),nrow=ncol(F_replicate))
	for(i in 1:ncol(F_replicate)){
		for(j in 1:ncol(F_replicate)){
		correl_replicate[i,j]<-cor(F_replicate[,i][!is.na(F_replicate[,i]) & !is.na(F_replicate[,j])],F_replicate[,j][!is.na(F_replicate[,i]) & !is.na(F_replicate[,j])])
		} 
	}
	diag(correl_replicate)<-NA
	correlation_replicate=mean(correl_replicate[!is.na(correl_replicate)])
	print(paste0('correlation replicates ',correlation_replicate))
}
if('random'%in%test){
	d_random<-fread(paste0(dir,'/depth_control_correl.txt'),sep=' ',header=T,data.table=F)
	c_random<-fread(paste0(dir,'/count_ref_control_correl.txt'),sep=' ',header=T,data.table=F)
	F_random<-c_random[,5:ncol(c_random)]/d_random[,5:ncol(d_random)]
	correl_random<-matrix(ncol=ncol(F_random),nrow=ncol(F_random))
	for(i in 1:ncol(F_random)){
		for(j in 1:ncol(F_random)){
		correl_random[i,j]<-cor(F_random[,i][!is.na(F_random[,i]) & !is.na(F_random[,j])],F_random[,j][!is.na(F_random[,i]) & !is.na(F_random[,j])])
		} 
	}
	diag(correl_random)<-NA
	correlation_random1=mean(correl_random[!is.na(correl_random)])
	print(paste0('correlation random colonies accross  ',correlation_random1))
	correl_random<-matrix(ncol=ncol(d_random),nrow=10)
	for(i in 1:ncol(d_random)){
		for(x in 1:10){
			samp<-sample(seq(1:nrow(d_random)),50000,replace=F)
			D<-d_random[samp,i]
			C<-c_random[samp,i]
			dat<-data.frame(D,C)
			datna<-dat[dat$D==0,]
			datna$d1<-0
			datna$d2<-0
			datna$c1<-0
			datna$c2<-0
			datnna1<-dat[dat$D!=0 & dat$D==dat$C,]
			datnna1$d1<-datnna1$D
			datnna1$d2<-datnna1$D
			datnna1$c1<-datnna1$C
			datnna1$c2<-datnna1$C
			datnna2<-dat[dat$D!=0 & dat$D!=dat$C,]
			datnna2<-datnna2%>%rowwise()%>%mutate(d1=floor(runif(1,min=1,max=D)),d2=D-d1,c1=length(which(sort(sample(seq(1,D,by=1),C,replace=F))<=d1)),c2=length(which(sort(sample(seq(1,D,by=1),C,replace=F))>d1)))
			dat<-rbind(datna,datnna1,datnna2)
			dat$f1<-dat$c1/dat$d1
			dat$f2<-dat$c2/dat$d2
			correl_random[x,i]<-cor(dat$f1[!is.na(dat$f1) & !is.na(dat$f2)],dat$f2[!is.na(dat$f1) & !is.na(dat$f2)])
		}
	}
	correlation_random2=mean(correl_random)
	print(paste0('correlation random colonies within  ',correlation_random2))
}
if('duplicate'%in%test){
	d<-fread(paste0(dir,'/depth_dup.txt'),sep=' ',header=T,data.table=F)
	r<-fread(paste0(dir,'/count_ref_dup.txt'),sep=' ',header=T,data.table=F)
	col<-colnames(d)[which(!grepl('bis',colnames(d)))]
	correl<-list()
	diff<-list()
	F_sum<-list()
	cv<-list()
	for(i in 1:length(col)){
		n<-which(grepl(col[i],colnames(d)))
		D<-d[,n]
		R<-r[,n]
		cv1<-data.frame('num_ruche'=col[i],'depth'=D[,1],'count'=R[,1],'measure'=1)
		cv2<-data.frame('num_ruche'=col[i],'depth'=D[,2],'count'=R[,2],'measure'=2)
		cv1$freq<-cv1$count/cv1$depth
		cv1$freq[is.na(cv1$freq)]<-0
		cv2$freq<-cv2$count/cv2$depth
		cv2$freq[is.na(cv2$freq)]<-0
		F_sum[[i]]<-rowSums(R)/rowSums(D)
		correl[[i]]<-data.frame('num_ruche'=col[i],'cor_freq'=cor(cv1$freq,cv2$freq))
		diff[[i]]<-data.frame(abs(cv1$freq-cv2$freq),abs(cv1$freq-F_sum[[i]]),abs(cv2$freq-F_sum[[i]]))
		cv[[i]]<-rbind(cv1,cv2)
	}
	Cor<-do.call(rbind,correl)
	print(Cor)
	CV<-do.call(rbind,cv)
	mat<-matrix(nrow=length(col),ncol=12)
	options(warn=2)
	for(i in 1:length(col)){
		mat[i,1]<-length(which(is.na(diff[[i]])))+length(which(diff[[i]]==Inf))
		mat[i,2]<-length(which(diff[[i]]==0))
		mat[i,3]<-length(which(diff[[i]]>0 & diff[[i]]<=0.1))
		mat[i,4]<-length(which(diff[[i]]>0.1 & diff[[i]]<=0.2))
		mat[i,5]<-length(which(diff[[i]]>0.2 & diff[[i]]<=0.3))
		mat[i,6]<-length(which(diff[[i]]>0.3 & diff[[i]]<=0.4))
		mat[i,7]<-length(which(diff[[i]]>0.4 & diff[[i]]<=0.5))
		mat[i,8]<-length(which(diff[[i]]>0.5 & diff[[i]]<=0.6))
		mat[i,9]<-length(which(diff[[i]]>0.6 & diff[[i]]<=0.7))
		mat[i,10]<-length(which(diff[[i]]>0.7 & diff[[i]]<=0.8))
		mat[i,11]<-length(which(diff[[i]]>0.8 & diff[[i]]<=0.9))
		mat[i,12]<-length(which(diff[[i]]>0.9 & diff[[i]]<=1))
		}
	rownames(mat)<-col
	colnames(mat)<-c('NA','0',']0,0.1]',']0.1,0.2]',']0.2,0.3]',']0.3,0.4]',']0.4,0.5]',']0.5,0.6]',']0.6,0.7]',']0.7,0.8]',']0.8,0.9]',']0.9,1]')
	print(mat)
	if('replicate'%in%test){threshold=correlation_replicate}else{threshold=0.95}
	option_duplicate<-vector()
	for(i in 1:length(col)){if(Cor[i,'cor_freq']>=threshold){option_duplicate[i]<-paste0('recode ',col[i])}else{option_duplicate[i]<-paste0('remove ',col[i])}}
	OPT<-as.data.frame(option_duplicate)
write.table(data.frame(OPT),paste0(dir,'/opt_dup.txt'),col.names=F,row.names=F,quote=F)
}






