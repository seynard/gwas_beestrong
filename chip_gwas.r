library(data.table)
library(seqinr)

setwd('Desktop/mnt')
lsnp<-fread('genotoul_dynagen_work/GWAS/infestation4/data/allele_id_final.txt')	
colnames(lsnp)<-c('chr','ps','ref','alt')
fasta0<-read.fasta('genotoul_genphyse_save/Fasta/GCA_003254395.1_Amel_HAv3_genomic.fna',seqtype="DNA",forceDNAtolower=F)
fasta1<-read.fasta('genotoul_genphyse_save/Fasta/GCA_003254395.2_Amel_HAv3.1_genomic.fna',seqtype="DNA",forceDNAtolower=F)

SNP_id<-c("1:10080627","1:20960056","1:25184394","11:9369229","12:10734707","4:11665460","5:75369","5:9190579","6:10450971","7:11806658","7:5762037","7:5772089","7:6738985","8:2468335","8:9799408","1:16327085","1:21374478","1:24201224","1:2891204","10:10400687","10:4266011","10:5359169","10:5359173","11:9527267","12:10153855","12:136634","15:4853529","2:11333369","2:4437645","2:9347464","3:6206342","5:2008472","5:2365495","8:1150346","8:9557205","1:12552506","1:15280956","1:7448807","1:7448811","10:2026877","11:14369154","13:9483955","14:3782741","14:6686131","14:8481541","14:8481589","15:2021142","15:2081876","15:2081914","15:8485332","16:1812909","16:5024160","2:12025610","2:12025647","2:16060868","2:2729874","2:8350714","3:1059430","3:12973246","3:12973248","4:10789077","4:4327611","4:7321246","4:7321247","5:6736534","5:6761414","5:8737386","7:7028040","7:7051965","7:7078376","7:8466948","8:1319815","8:1551638","8:1748209","9:11564671")

df<-list()
for(i in 1:length(SNP_id)){
print(SNP_id[i])
CHR<-unlist(strsplit(SNP_id[i],':'))[1]
POS<-as.numeric(unlist(strsplit(SNP_id[i],':'))[2])
p_min<-POS-50
p_max<-POS+50
if(CHR==1){chr_fast1<-'CM009931.2'}else if(CHR==2){chr_fast1<-'CM009932.2'}else if(CHR==3){chr_fast1<-'CM009933.2'}else if(CHR==4){chr_fast1<-'CM009934.2'}else if(CHR==5){chr_fast1<-'CM009935.2'}else if(CHR==6){chr_fast1<-'CM009936.2'}else if(CHR==7){chr_fast1<-'CM009937.2'}else if(CHR==8){chr_fast1<-'CM009938.2'}else if(CHR==9){chr_fast1<-'CM009939.2'}else if(CHR==10){chr_fast1<-'CM009940.2'}else if(CHR==11){chr_fast1<-'CM009941.2'}else if(CHR==12){chr_fast1<-'CM009942.2'}else if(CHR==13){chr_fast1<-'CM009943.2'}else if(CHR==14){chr_fast1<-'CM009944.2'}else if(CHR==15){chr_fast1<-'CM009945.2'}else if(CHR==16){chr_fast1<-'CM009946.2'}
seq_chr1<-fasta1[chr_fast1][[1]]
seq1<-seq_chr1[p_min:p_max]
x<-data.frame('id'=seq(1:101),'nt'=seq1,'ps'=seq(p_min,p_max,by=1))
snp_in<-lsnp[lsnp$chr==CHR & lsnp$ps>=p_min & lsnp$ps<=p_max,]
x$type<-NA
x$type[x$ps%in%snp_in$ps]<-'snp'
x$type[x$ps==POS]<-'ref'
if(CHR==1){chr_fast0<-'CM009931.1'}else if(CHR==2){chr_fast0<-'CM009932.1'}else if(CHR==3){chr_fast0<-'CM009933.1'}else if(CHR==4){chr_fast0<-'CM009934.1'}else if(CHR==5){chr_fast0<-'CM009935.1'}else if(CHR==6){chr_fast0<-'CM009936.1'}else if(CHR==7){chr_fast0<-'CM009937.1'}else if(CHR==8){chr_fast0<-'CM009938.1'}else if(CHR==9){chr_fast0<-'CM009939.1'}else if(CHR==10){chr_fast0<-'CM009940.1'}else if(CHR==11){chr_fast0<-'CM009941.1'}else if(CHR==12){chr_fast0<-'CM009942.1'}else if(CHR==13){chr_fast0<-'CM009943.1'}else if(CHR==14){chr_fast0<-'CM009944.1'}else if(CHR==15){chr_fast0<-'CM009945.1'}else if(CHR==16){chr_fast0<-'CM009946.1'}
seq_chr0<-paste0(fasta0[chr_fast0][[1]],collapse='')
match<-gregexpr(pattern=toupper(paste0(seq1,collapse='')),toupper(seq_chr0))
p_start<-`attributes<-`(match[[1]],NULL) 
lens<-attr(match[[1]],'match.length')
if(lens!=101){
	print('not full match')
	match_up<-gregexpr(pattern=toupper(paste0(seq1[1:51],collapse='')),toupper(seq_chr0))
	p_up<-`attributes<-`(match_up[[1]],NULL) 
	lens_up<-attr(match_up[[1]],'match.length')
	match_down<-gregexpr(pattern=toupper(paste0(seq1[51:101],collapse='')),toupper(seq_chr0))
	p_down<-`attributes<-`(match_down[[1]],NULL) 
	lens_down<-attr(match_down[[1]],'match.length')
	if(lens_up==51){
		x$ps_0<-c(seq(p_up,(p_up+50),1),rep(NA,50))
		x$nt_probe<-toupper(x$nt)
		x$nt_probe[x$type=='snp']<-'N'
		x$nt_probe[x$type=='ref']<-paste0('[',snp_in$ref[snp_in$ps==POS],'/',snp_in$alt[snp_in$ps==POS],']')
		x$nt_probe[is.na(x$ps_0)]<-NA
		seq_up<-paste0(x$nt_probe[1:50],collapse='')
		seq_down<-paste0(x$nt_probe[52:101],collapse='')
		SNP<-x$nt_probe[which(x$type=='ref')]
		nsnp_up<-length(grep('N',unlist(strsplit(seq_up,''))))
		nsnp_down<-length(grep('N',unlist(strsplit(seq_down,''))))
		df[[i]]<-data.frame('RS_v1'=SNP_id[i],'chr'=CHR,'ps_v1'=POS,'ps_v0'=x$ps_0[which(x$type=='ref')],seq_up=seq_up,seq_down=seq_down,n_snp_up=nsnp_up,n_snp_down=nsnp_down)
		df[[i]]$RS_v0<-paste0(df[[i]]$chr,':',df[[i]]$ps_v0)
		if(toupper(seq1[51])%in%snp_in[snp_in$ps==POS,c('ref','alt')]){print('ok')}else{print('error');break}
	}else if (lens_down==51){
		x$ps_0<-c(rep(NA,50),seq(p_down,(p_down+50),1))
		x$nt_probe<-toupper(x$nt)
		x$nt_probe[x$type=='snp']<-'N'
		x$nt_probe[x$type=='ref']<-paste0('[',snp_in$ref[snp_in$ps==POS],'/',snp_in$alt[snp_in$ps==POS],']')
		x$nt_probe[is.na(x$ps_0)]<-NA
		seq_up<-paste0(x$nt_probe[1:50],collapse='')
		seq_down<-paste0(x$nt_probe[52:101],collapse='')
		SNP<-x$nt_probe[which(x$type=='ref')]
		nsnp_up<-length(grep('N',unlist(strsplit(seq_up,''))))
		nsnp_down<-length(grep('N',unlist(strsplit(seq_down,''))))
		df[[i]]<-data.frame('RS_v1'=SNP_id[i],'chr'=CHR,'ps_v1'=POS,'ps_v0'=x$ps_0[which(x$type=='ref')],seq_up=seq_up,seq_down=seq_down,n_snp_up=nsnp_up,n_snp_down=nsnp_down)
		df[[i]]$RS_v0<-paste0(df[[i]]$chr,':',df[[i]]$ps_v0)
		if(toupper(seq1[51])%in%snp_in[snp_in$ps==POS,c('ref','alt')]){print('ok')}else{print('error');break}
	}else{print('not 1/2 match')}
}else{
	x$ps_0<-seq(p_start,(p_start+100),1)
	x$nt_probe<-toupper(x$nt)
	x$nt_probe[x$type=='snp']<-'N'
	x$nt_probe[x$type=='ref']<-paste0('[',snp_in$ref[snp_in$ps==POS],'/',snp_in$alt[snp_in$ps==POS],']')
	seq_up<-paste0(x$nt_probe[1:50],collapse='')
	seq_down<-paste0(x$nt_probe[52:101],collapse='')
	SNP<-x$nt_probe[which(x$type=='ref')]
	nsnp_up<-length(grep('N',unlist(strsplit(seq_up,''))))
	nsnp_down<-length(grep('N',unlist(strsplit(seq_down,''))))
	df[[i]]<-data.frame('RS_v1'=SNP_id[i],'chr'=CHR,'ps_v1'=POS,'ps_v0'=x$ps_0[which(x$type=='ref')],seq_up=seq_up,seq_down=seq_down,n_snp_up=nsnp_up,n_snp_down=nsnp_down)
	df[[i]]$RS_v0<-paste0(df[[i]]$chr,':',df[[i]]$ps_v0)
	if(toupper(seq1[51])%in%snp_in[snp_in$ps==POS,c('ref','alt')]){print('ok')}else{print('error');break}
	}
}
df<-do.call(rbind,df)
write.table(df,'genotoul_dynagen_work/GWAS/infestation4/snp_chip.txt',col.names=T,row.names=F,quote=F)

