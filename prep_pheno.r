library(data.table)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(MASS)
library(lme4)
library(lmerTest)
library(car)
library(lattice)
library(stringr)
library(tidyverse)
library(FactoMineR)
library(missMDA)
library(factoextra)
library(RColorBrewer)
library(corrplot)
library(gridGraphics)
library(ggpubr)

args<-commandArgs(TRUE)

#dir_in<-args[1]
#dir_out<-args[2]
#df=args[3]
#pop_id=args[4]
#pop_id<-unlist(strsplit(pop_id,','))
#pheno_mito=args[5]
dir_in<-'data'
dir_out<-'results'
df='Data_BS.csv'
pop_id='Mellifera,Ligustica_Carnica,Caucasica,hybrid,Ligustica_nus'
pop_id<-unlist(strsplit(pop_id,','))
pheno_mito='compareDepthsBeeVarroa.txt'

data<-fread(paste0(dir_in,'/',df),header=T,data.table=F)
for(i in 1:ncol(data)){data[,i]<-tolower(data[,i])
	data[,i]<-trimws(data[,i], which="both")
	data[,i][data[,i]=='na']<-NA
	data[,i]<-gsub(',','.',data[,i])}
colnames(data)<-gsub('\\.','_',colnames(data))
data[,4]<-gsub(' ','_',data[,4])
data[,6]<-gsub(' ','_',data[,6])
data[,9]<-gsub(' ','_',data[,9])
for(i in 1:nrow(data)){if(is.na(data$heure[i])){data$heure[i]<-NA}
	else if(nchar(data$heure[i])<5){data$heure[i]<-paste0('0',data$heure[i],collapse='')}}
for(i in 1:ncol(data)){data[,i]<-tolower(data[,i])
	data[,i]<-trimws(data[,i], which="both")
	data[,i][data[,i]=='na']<-NA
	data[,i]<-gsub(',','.',data[,i])
	data[,i][data[,i]=='']<-NA}
data$date_pre<-as.Date(data$date_pre,format='%d/%m/%Y')
data$date<-as.Date(data$date_pre,format='%d/%m/%Y')
data$year<-year(data$date)
data$date<-format(data$date,format="%m-%d")
a<-colsplit(data$date,'-',c('month','day'))
data$month<-a$month
data$week<-strftime(as.Date(data$date_pre),format="%V")
data$day<-strftime(data$date_pre,format="%j")
data$observ_smr[is.na(data$observ_smr)]<-'unknown'
data$observ_vp[is.na(data$observ_vp)]<-'unknown'
data$observ_coleval[is.na(data$observ_coleval)]<-'unknown'
data$structure[is.na(data$structure)]<-'unknown'
data$api[is.na(data$api)]<-'unknown'
data$rucher[is.na(data$rucher)]<-'unknown'
col_num<-c("nbr_abeilles","nbr_cell_cf","nbr_cell_co","nbr_cell_male","surf_res","surf_pol","poids_abeilles","nbr_varroa_pho","nbr_varroa_100abeilles","nbr_varroa_abeille","nbr_ouvert","nbr_reop","nbr_inf_reop","nbr_inf","nbr_nymph_inf","nbr_varroa_smr","nbr_inf_une_fond","nbr_inf_nr","nbr_retard","nbr_infertilite","nbr_sans_male","charge_varroa","smr")
for(i in 1:length(col_num)){data[,col_num[i]]<-as.numeric(data[,col_num[i]])}
data[data$nbr_reop>data$nbr_ouvert & !is.na(data$nbr_reop),'nbr_reop']<-data[data$nbr_reop>data$nbr_ouvert & !is.na(data$nbr_reop),'nbr_ouvert']
data$smr_brut<-data$nbr_inf_nr/data$nbr_inf_une_fond
data$smr35<-data$smr_brut
data$smr35[data$nbr_inf_une_fond<35]<-NA
eb_smr_fit=fitdistr(x=data$smr_brut[!is.na(data$smr_brut) & data$smr_brut>0 & data$smr_brut<1],densfun="beta",start=list(shape1=1,shape2=1),method="L-BFGS-B")
aprior=eb_smr_fit$estimate[1]
bprior=eb_smr_fit$estimate[2]
data$eb_smr=(aprior+data$nbr_inf_nr)/(aprior+bprior+data$nbr_inf_une_fond)
#data$smr_brut_unif<-(data$nbr_inf_nr+1)/(data$nbr_inf_une_fond+2)
#data$logit_smr<-log((data$smr_brut_unif)/(1-data$smr_brut_unif))
data$couvain<-data$nbr_cell_cf+data$nbr_cell_co
data$taux_reop_gen<-data$nbr_reop/data$nbr_ouvert
data$taux_reop_gen_unif<-(data$nbr_reop+1)/(data$nbr_ouvert+2)
data$logit_taux_reop_gen<-log((data$taux_reop_gen_unif)/(1-data$taux_reop_gen_unif))
data$taux_reop_spe<-data$nbr_inf_reop/data$nbr_inf
data$taux_reop_spe_unif<-(data$nbr_inf_reop+1)/(data$nbr_inf+2)
data$logit_taux_reop_spe<-log((data$taux_reop_spe_unif)/(1-data$taux_reop_spe_unif))
data$taux_inf<-data$nbr_inf/data$nbr_ouvert
data$taux_inf_unif<-(data$nbr_inf+1)/(data$nbr_ouvert+2)
data$logit_taux_inf<-log((data$taux_inf_unif)/(1-data$taux_inf_unif))
data$vp<-(data$nbr_varroa_pho/data$poids_abeilles)*0.14*100
pheno_mito<-fread(paste0(dir_in,'/',pheno_mito))
pheno_mito$name<-gsub('-','_',pheno_mito$name)
pheno_mito$name<-gsub('_F','',pheno_mito$name)
a<-colsplit(pheno_mito$name,'_',c('v1','v2'))
pheno_mito$name<-paste0(a$v1,'_',str_pad(a$v2,4,pad="0"))
colnames(pheno_mito)[colnames(pheno_mito)=='name']<-'num_ruche_bs'
pheno_mito$num_ruche_bs<-tolower(pheno_mito$num_ruche_bs)
data<-merge(data,pheno_mito,by='num_ruche_bs',all=T)
x<-data.frame(table(data$num_ruche_bs))
x<-subset(x[,1],x[,2]>1)
d1<-subset(data,!(data$num_ruche_bs%in%x))
d2<-list()
for(i in 1:length(x)){
	d<-subset(data,data$num_ruche_bs==x[i])
	a<-data.frame(table(data$date_pre))
	a<-subset(a,as.character(a[,1])%in%as.character(d$date_pre))
	n<-which.max(a[,2])
	d2[[i]]<-subset(d,as.character(d$date_pre)==as.character(a[n,1]))
}
d2<-do.call(rbind,d2)
data<-rbind(d1,d2)
col_num<-c("smr_brut","smr35","eb_smr","couvain","taux_reop_gen","taux_reop_gen_unif","logit_taux_reop_gen","taux_reop_spe","taux_reop_spe_unif","logit_taux_reop_spe","taux_inf","taux_inf_unif","logit_taux_inf","vp","chrDepth","beeMitoDepth","varroaMitoDepth","beeMitoratio","varroaMitoRatio")
for(i in 1:length(col_num)){data[,col_num[i]]<-as.numeric(data[,col_num[i]])}
 
data$group<-NA
list_col<-list.files(path=dir_in,pattern='list_')
list_col<-list_col[list_col!='list_ref.txt']
for(i in 1:length(list_col)){
	x<-fread(paste0(dir_in,'/',list_col[i]),header=F,data.table=F)
	print(nrow(x))
	miss<-which(!(tolower(x[,1])%in%data$num_ruche_bs))
	if(length(miss>0)){
		d<-c(tolower(x[miss,1]),rep(NA,78))
		d<-data.frame(t(d))
		colnames(d)<-colnames(data)
		data<-rbind(data,d)}
	i_id<-gsub('list_','',list_col[i])
	i_id<-gsub('.txt','',i_id)
	print(i_id)
	data$group[data$num_ruche_bs%in%tolower(x[,1])]<-i_id}
data$group[data$group=='Ligustica_nus']<-'Ligustica_Carnica'

data[data$type_corps%in%c('1/2d','1/2n','hausse_d','nuc','ruchette_d','ruchette_l','w'),c('pass','type_corps','api')]<-'corps_prob'
data$kept<-'no'
data$kept[data$pass%in%c('apisjungels','arista','no','nz','us','yes')]<-'yes'
data$kept[is.na(data$group)]<-'no'
x<-data.frame(table(data$api[data$kept=='yes']))
x<-subset(x,x$Freq<6)
data$kept[data$api%in%x$Var1]<-'no'
for(i in 1:nrow(data)){if(data$kept[i]=='no'){data[i,c("tenue_cadre","douceur","nbr_abeilles","nbr_cell_cf","nbr_cell_co","nbr_cell_male","surf_res","surf_pol","poids_abeilles","nbr_varroa_pho","nbr_varroa_100abeilles","nbr_varroa_abeille","nbr_ouvert","nbr_reop","nbr_inf_reop","nbr_inf","nbr_nymph_inf","nbr_varroa_smr","nbr_inf_une_fond","nbr_inf_nr","nbr_retard","nbr_infertilite","nbr_sans_male","charge_varroa","smr",		"smr_brut","smr35","eb_smr","couvain","taux_reop_gen","taux_reop_gen_unif","logit_taux_reop_gen","taux_reop_spe","taux_reop_spe_unif","logit_taux_reop_spe","taux_inf","taux_inf_unif","logit_taux_inf","vp","chrDepth","beeMitoDepth","varroaMitoDepth","beeMitoratio","varroaMitoRatio")]<-NA}}
col_num<-c("vp","charge_varroa","logit_taux_inf","varroaMitoRatio","smr_brut","eb_smr","logit_taux_reop_gen","logit_taux_reop_spe")
for(i in 1:length(col_num)){data[,col_num[i]]<-as.numeric(data[,col_num[i]])}

df_pc<-data.frame(phoretic_varroa=log(data$vp+1),varroa_load=(data$charge_varroa)^(1/4),varroa_brood_infestation=data$logit_taux_inf,varroa_mitochondria_sequencing_depth=(data$varroaMitoRatio)^(1/4),MNR=data$smr_brut,EB_MNR=data$eb_smr,recapping=data$logit_taux_reop_gen,recapping_infested_cells=data$logit_taux_reop_spe)
nb=estim_ncpPCA(df_pc,ncp.max=ncol(df_pc))
res.comp=imputePCA(df_pc,ncp=nb$ncp)
res.pca=PCA(res.comp$completeObs,graph=F) 
cor.mat<-round(cor(res.comp$completeObs),2)
corrplot(cor.mat,type='upper',order="hclust",tl.col="black",tl.srt=45,cl.pos='r',col=brewer.pal(n=8,name='RdYlGn'))
	
for(i in 1:length(list_col)){
	x<-fread(paste0(dir_in,'/',list_col[i]),header=F,data.table=F)
	print(nrow(x))
	i_id<-gsub('list_','',list_col[i])
	i_id<-gsub('.txt','',i_id)
	print(i_id)
	file_dir<-intersect(list.files(path=dir_in,pattern=paste0('geno_hom_',i_id)),list.files(path=dir_in,pattern='.bgs'))
	if(paste0('geno_hom_',i_id,'.bgs') %in% file_dir){
	file_dir<-intersect(list.files(path=dir_in,pattern=paste0('geno_hom_',i_id)),list.files(path=dir_in,pattern='.bgs'))
	con_geno<-file(paste0(dir_in,'/geno_hom_',i_id,'.bgs'),"r")
	first_line_geno<-readLines(con_geno,n=1)
	first_line_geno<-unlist(strsplit(first_line_geno,','))
	first_line_geno<-first_line_geno[!first_line_geno%in%c('','CHROM','POS','REF','ALT')]
	con_freq<-file(paste0(dir_in,'/freq_',i_id,'.txt.bz2'),"r")
	first_line_freq<-readLines(con_freq,n=1)
	first_line_freq<-unlist(strsplit(first_line_freq,' '))
	first_line_freq<-first_line_freq[!first_line_freq%in%c('','CHROM','POS','REF','ALT')]
	print(table(first_line_freq==first_line_geno))}
	df<-subset(data,data$num_ruche_bs%in%tolower(x[,1]))
	if(paste0('geno_hom_',i_id,'.bgs') %in% file_dir){df<-df%>%arrange(factor(num_ruche_bs,levels=tolower(first_line_geno)))}
	print(nrow(df))
	df_pc<-data.frame(phoretic_varroa=log(df$vp+1),varroa_load=(df$charge_varroa)^(1/4),varroa_brood_infestation=df$logit_taux_inf,varroa_mitochondria_sequencing_depth=(df$varroaMitoRatio)^(1/4),EB_MNR=df$eb_smr,recapping_infested_cells=df$logit_taux_reop_spe)
	df_pc[,which(colSums(is.na(df_pc))==nrow(df_pc))]<-NULL
	nb=estim_ncpPCA(df_pc,ncp.max=ncol(df_pc))
	res.comp=imputePCA(df_pc,ncp=nb$ncp)
	res.pca=PCA(res.comp$completeObs,graph=F) 
	df$pca<-res.pca$ind$coord[,1]
	df$pca[df$kept=='no']<-NA
	write.table(df,paste0(dir_out,'/pheno_',i_id,'.txt'),col.names=T,row.names=F,quote=F)
	write.table(log(df$vp+1),paste0(dir_out,'/pheno_varroaphoretic_',i_id,'.txt'),col.names=F,row.names=F,quote=F)
	write.table((df$charge_varroa)^(1/4),paste0(dir_out,'/pheno_chargevarroaracine4_',i_id,'.txt'),col.names=F,row.names=F,quote=F)
	write.table(df$logit_taux_inf,paste0(dir_out,'/pheno_logit_varroainfestation_',i_id,'.txt'),col.names=F,row.names=F,quote=F)
	write.table((df$varroaMitoRatio)^(1/4),paste0(dir_out,'/pheno_varroadepthmitoracine4_',i_id,'.txt'),col.names=F,row.names=F,quote=F)
	write.table(df$smr_brut,paste0(dir_out,'/pheno_smr_brut_',i_id,'.txt'),col.names=F,row.names=F,quote=F)
	write.table(df$eb_smr,paste0(dir_out,'/pheno_eb_smr_',i_id,'.txt'),col.names=F,row.names=F,quote=F)
	write.table(df$logit_taux_reop_gen,paste0(dir_out,'/pheno_logit_taux_reop_',i_id,'.txt'),col.names=F,row.names=F,quote=F)
	write.table(df$logit_taux_reop_spe,paste0(dir_out,'/pheno_logit_taux_reop_inf_',i_id,'.txt'),col.names=F,row.names=F,quote=F)
	write.table(df$pca,paste0(dir_out,'/pheno_pc1_',i_id,'.txt'),col.names=F,row.names=F,quote=F)
	}

length(unique(data$api))
length(unique(data$api[!(data$region%in%c('etats_unis','luxembourg','nl','nouvelle_zelande','suede','suisse'))])) 
table(data$group)

length(unique(data$api[data$region!='etats_unis']))
length(unique(data$api[!(data$region%in%c('etats_unis','luxembourg','nl','nouvelle_zelande','suede','suisse'))])) 
table(data$group)

data_sub<-subset(data,data$kept!='no' & !is.na(data$group))
nrow(data_sub)
length(unique(data_sub$api))
x<-data.frame(table(data_sub$api))
sum(x[,2])
summary(x[,2])

# varroa mito
summary(data_sub$varroaMitoRatio)
ggplot(data_sub)+geom_boxplot(aes(varroaMitoRatio,fill=group))
# phoretic varroa
summary(data_sub$vp)
ggplot(data_sub)+geom_point(aes(x=vp,y=as.numeric(poids_abeilles),col=group))
ggplot(data_sub)+geom_point(aes(x=vp,y=varroaMitoRatio,col=group))
# varroa infestation
summary(as.numeric(data_sub$taux_inf))
summary(as.numeric(data_sub$nbr_ouvert))
ggplot(data_sub)+geom_point(aes(x=as.numeric(nbr_ouvert),y=as.numeric(taux_inf),col=group))
ggplot(data_sub)+geom_point(aes(x=as.numeric(nbr_ouvert),y=as.numeric(logit_taux_inf),col=group))
# varroa load
summary(data_sub$charge_varroa)
ggplot(data_sub)+geom_boxplot(aes(charge_varroa,fill=group))
ggplot(data_sub)+geom_point(aes(x=as.numeric(charge_varroa),y=as.numeric(couvain),col=group))
ggplot(data_sub)+geom_point(aes(x=as.numeric(charge_varroa),y=as.numeric(nbr_abeilles),col=group))
#y<-as.numeric(data_sub$taux_inf)*as.numeric(data_sub$couvain)+0.01*as.numeric(data_sub$vp)*as.numeric(data_sub$nbr_abeilles)
#plot(as.numeric(data_sub$charge_varroa),y)

par(mfrow=c(2,3))
plot(as.numeric(data_sub$vp),as.numeric(data_sub$varroaMitoRatio))
plot(as.numeric(data_sub$vp),as.numeric(data_sub$taux_inf))
plot(as.numeric(data_sub$vp),as.numeric(data_sub$charge_varroa))
plot(as.numeric(data_sub$varroaMitoRatio),as.numeric(data_sub$taux_inf))
plot(as.numeric(data_sub$varroaMitoRatio),as.numeric(data_sub$charge_varroa))
plot(as.numeric(data_sub$charge_varroa),as.numeric(data_sub$taux_inf))

cor.test(as.numeric(data_sub$vp),as.numeric(data_sub$varroaMitoRatio))
cor.test(as.numeric(data_sub$vp),as.numeric(data_sub$taux_inf))
cor.test(as.numeric(data_sub$vp),as.numeric(data_sub$charge_varroa))
cor.test(as.numeric(data_sub$varroaMitoRatio),as.numeric(data_sub$taux_inf))
cor.test(as.numeric(data_sub$varroaMitoRatio),as.numeric(data_sub$charge_varroa))
cor.test(as.numeric(data_sub$charge_varroa),as.numeric(data_sub$taux_inf))

some.eu.countries<-c("France","Switzerland","Luxembourg","Netherlands","Belgium")
some.eu.maps<-map_data("world", region = some.eu.countries)
region.lab.data<-some.eu.maps%>%group_by(region)%>%summarise(long=mean(long),lat=mean(lat))
countries_lines<-tibble(x=c(0.1278,2.3522,13.4050,12.5683),xend=c(2.3522,13.4050,12.5683,0.1278),y=c(51.5074,48.8566,52.5200,55.6761),yend=c(48.8566,52.5200,55.6761,51.5074))
p1<-ggplot(some.eu.maps,aes(x=long,y=lat))+
	geom_polygon(aes(group=group),colour='black',fill='transparent')+
	geom_count(aes(x=as.numeric(coord_long_echant),y=as.numeric(coord_lat_echant),col='col2'),data=data_sub[!data_sub$region%in%c('etats_unis','nouvelle_zelande','suede'),])+
	geom_count(aes(x=as.numeric(coord_long_echant),y=as.numeric(coord_lat_echant),col='col1'),data=data[!is.na(data$group) & data$kept=='no' & !data$region%in%c('etats_unis','nouvelle_zelande','suede'),])+
	xlab('Longitude')+ylab('Latitude')+
	theme_bw()+
	theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
	ggtitle('France and neighbour countries')+
	scale_size_area(limits=c(1,200),breaks=c(5,10,50,100,150),name='# colonies',guide=guide_legend(override.aes=list(colour='black')))+
	scale_color_manual(name='',values=c('col1'='chartreuse4','col2'='orangered3'),labels=c('sequenced','sequenced \nand phenotyped'))
	
an.eu.countries<-c("Sweden")
an.eu.maps<-map_data("world", region = an.eu.countries)
region.lab.data<-an.eu.maps%>%group_by(region)%>%summarise(long=mean(long),lat=mean(lat))
countries_lines<-tibble(x=c(0.1278,2.3522,13.4050,12.5683),xend=c(2.3522,13.4050,12.5683,0.1278),y=c(51.5074,48.8566,52.5200,55.6761),yend=c(48.8566,52.5200,55.6761,51.5074))
p2<-ggplot(an.eu.maps,aes(x=long,y=lat))+
	geom_polygon(aes(group=group),colour='black',fill='transparent')+
	geom_count(aes(x=as.numeric(coord_long_echant),y=as.numeric(coord_lat_echant)),data=data_sub[data_sub$region=='suede',],color='orangered3')+
	geom_count(aes(x=as.numeric(coord_long_echant),y=as.numeric(coord_lat_echant)),data=data[!is.na(data$group) & data$kept=='no' & data$region=='suede',],color='chartreuse4')+
	xlab('Longitude')+ylab('Latitude')+
	theme_bw()+
	theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
	ggtitle('Sweden')+
	scale_size_area(limits=c(1,200),breaks=c(5,10,50,100,150))+
	theme(legend.position='none')

mapCountry<-borders('usa',colour="black",fill="transparent") 
p3<-ggplot()+
	mapCountry+ 
	coord_fixed(1.3)+
	geom_count(aes(x=as.numeric(coord_long_echant),y=as.numeric(coord_lat_echant)),data=data_sub[data_sub$region=='etats_unis',],color='orangered3')+
	geom_count(aes(x=as.numeric(coord_long_echant),y=as.numeric(coord_lat_echant)),data=data[!is.na(data$group) & data$kept=='no' & data$region=='etats_unis',],color='chartreuse4')+
	xlab('Longitude')+ylab('Latitude')+
	theme_bw()+
	theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
	ggtitle('USA')+
	scale_size_area(limits=c(1,200),breaks=c(5,10,50,100,150))+
	theme(legend.position='none')

mapCountry<-borders('nz',colour="black",fill="transparent") 
p4<-ggplot()+
	mapCountry+ 
	coord_fixed(1.3)+
	geom_count(aes(x=as.numeric(coord_long_echant),y=as.numeric(coord_lat_echant)),data=data_sub[data_sub$region=='nouvelle_zelande',],color='orangered3')+
	geom_count(aes(x=as.numeric(coord_long_echant),y=as.numeric(coord_lat_echant)),data=data[!is.na(data$group) & data$kept=='no' & data$region=='nouvelle_zelande',],color='chartreuse4')+
	xlab('Longitude')+ylab('Latitude')+
	theme_bw()+
	theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
	ggtitle('New Zealand')+
	scale_size_area(limits=c(1,200),breaks=c(5,10,50,100,150))+
	theme(legend.position='none')

png('plot_map_sample.png',width=2500,height=3000,res=400)
ggarrange(p1,ggarrange(p2,p3,p4,ncol=3,widths=c(0.3,0.7,0.3)),nrow=2,heights=c(0.72,0.3))
dev.off()
