library(data.table)
library(ggplot2)
library(tidyverse)
#library(mdthemes)
library(FactoMineR)
library(missMDA)
library(factoextra)
library(RColorBrewer)
library(corrplot)
library(gridGraphics)
library(gridExtra)
library(ggcorrplot)
library(ggpubr)
grab_grob<-function(){
  grid.echo()
  grid.grab()
}

args<-commandArgs(TRUE)
version<-as.numeric(args[1])
if(version==1){pop<-c('Ligustica_Carnica','Mellifera','hybrid')}else if (version==2){pop<-c('Ligustica_nus','Mellifera','hybrid')}

### number api + map
df_lig<-fread(paste0('results/pheno_',pop[1],'.txt'),header=T,data.table=F)
df_mel<-fread(paste0('results/pheno_',pop[2],'.txt'),header=T,data.table=F)
df_hyb<-fread(paste0('results/pheno_',pop[3],'.txt'),header=T,data.table=F)
data<-do.call(rbind,list(df_lig,df_mel,df_hyb))
length(unique(data$api))
length(unique(data$api[data$kept=='yes']))
data_sub<-subset(data,data$kept!='no' & !is.na(data$group))
nrow(data_sub)
length(unique(data_sub$api))
x<-data.frame(table(data_sub$api))
sum(x[,2])
summary(x[,2])

if(version==2){
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
	theme(plot.title=element_text(size=20),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.title=element_text(size=12),legend.text=element_text(size=10),axis.title.x=element_text(size=12),axis.title.y=element_text(size=12))+
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
	theme(plot.title=element_text(size=20),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.title=element_text(size=12),legend.text=element_text(size=10),axis.title.x=element_text(size=12),axis.title.y=element_text(size=12))+
	ggtitle('Sweden')+
	scale_size_area(limits=c(1,200),breaks=c(5,10,50,100,150))+
	theme(legend.position='none')

mapCountry<-borders('nz',colour="black",fill="transparent") 
p3<-ggplot()+
	mapCountry+ 
	coord_fixed(1.3)+
	geom_count(aes(x=as.numeric(coord_long_echant),y=as.numeric(coord_lat_echant)),data=data_sub[data_sub$region=='nouvelle_zelande',],color='orangered3')+
	geom_count(aes(x=as.numeric(coord_long_echant),y=as.numeric(coord_lat_echant)),data=data[!is.na(data$group) & data$kept=='no' & data$region=='nouvelle_zelande',],color='chartreuse4')+
	xlab('Longitude')+ylab('Latitude')+
	theme_bw()+
	theme(plot.title=element_text(size=20),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.title=element_text(size=12),legend.text=element_text(size=10),axis.title.x=element_text(size=12),axis.title.y=element_text(size=12))+
	ggtitle('New \nZealand')+
	scale_size_area(limits=c(1,200),breaks=c(5,10,50,100,150))+
	theme(legend.position='none')

png(paste0('plot_map_sample_',version,'.png'),width=2200,height=3000,res=400)
print(ggarrange(p1,ggarrange(p2,p3,ncol=2,widths=c(0.4,0.6)),nrow=2,heights=c(0.7,0.4)))
dev.off()
}else if(version==1){
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
	theme(plot.title=element_text(size=20),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.title=element_text(size=12),legend.text=element_text(size=10),axis.title.x=element_text(size=12),axis.title.y=element_text(size=12))+
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
	theme(plot.title=element_text(size=20),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.title=element_text(size=12),legend.text=element_text(size=10),axis.title.x=element_text(size=12),axis.title.y=element_text(size=12))+
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
	theme(plot.title=element_text(size=20),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.title=element_text(size=12),legend.text=element_text(size=10),axis.title.x=element_text(size=12),axis.title.y=element_text(size=12))+
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
	theme(plot.title=element_text(size=20),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.title=element_text(size=12),legend.text=element_text(size=10),axis.title.x=element_text(size=12),axis.title.y=element_text(size=12))+
	ggtitle('New \nZealand')+
	scale_size_area(limits=c(1,200),breaks=c(5,10,50,100,150))+
	theme(legend.position='none')

png(paste0('plot_map_sample_',version,'.png'),width=2500,height=3000,res=400)
print(ggarrange(p1,ggarrange(p2,p3,p4,ncol=3,widths=c(0.3,0.6,0.4)),nrow=2,heights=c(0.72,0.3)))
dev.off()
}

### admixture colonies
l_lig<-fread(paste0('data/list_',pop[1],'.txt'),header=F,data.table=F)
l_mel<-fread(paste0('data/list_',pop[2],'.txt'),header=F,data.table=F)
l_hyb<-fread(paste0('data/list_',pop[3],'.txt'),header=F,data.table=F)
q_mat50<-fread('data/queen_geno_50k.Q',data.table=F)
q_mat50$group<-NA
for(i in 1:nrow(q_mat50)){
	if(q_mat50$num_ruche_bs[i]%in%l_lig[,1]){
	q_mat50[i,'group']<-'lig'
	}else if (q_mat50$num_ruche_bs[i]%in%l_mel[,1]){
	q_mat50[i,'group']<-'mel'
	}else if (q_mat50$num_ruche_bs[i]%in%l_hyb[,1]){
	q_mat50[i,'group']<-'hyb'
	}
	}
q_mat50<-gather(q_mat50,sp,val,c('Mellifera','Ligustica_Carnica','Caucasica'))
q_mat50$sp<-factor(q_mat50$sp,levels=c('Ligustica_Carnica','Mellifera','Caucasica'))
q_mat50$group[is.na(q_mat50$group)]<-'no_group'
q_mat50$group<-factor(q_mat50$group,levels=c('lig','mel','hyb','no_group'))

group_names<-c(
  `lig` = "Ligustica & Carnica",
  `mel` = "Mellifera",
  `hyb` = "Hybrids",
  `no_group` = "no group")

png(paste0('plot_paper_admixture_barplot_',version,'.png'),width=4000,height=500,res=400)
ggplot(q_mat50[q_mat50$group!='no_group',])+
	geom_bar(aes(x=num_ruche_bs,y=val,fill=sp),stat='identity',position='stack')+
	scale_fill_manual(name='Genetic background',values=c('goldenrod2','grey34','chartreuse3'),labels=c('Ligustica & Carnica','Mellifera','Caucasia'))+
	facet_grid(.~group,scale='free',labeller=labeller(group=as_labeller(group_names)))+
	theme_bw()+
	theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+
	theme(legend.title=element_text(size=12),legend.text=element_text(size=10))+
	theme(strip.text.x=element_text(size=12))
dev.off()

png(paste0('plot_paper_admixture_boxplot_',version,'.png'),width=4000,height=1500,res=400)
ggplot(q_mat50[q_mat50$group!='no_group',])+
	geom_boxplot(aes(x=group,y=val,fill=sp))+
	scale_fill_manual(name='Genetic background',values=c('goldenrod2','grey34','chartreuse3'),labels=c('Ligustica & Carnica','Mellifera','Caucasia'))+
	theme_bw()+	facet_grid(.~group,scale='free',labeller=labeller(group=as_labeller(group_names)))+
	theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_text(size=12),axis.text.y=element_text(size=10))+ylab('Proportion of genetic background')+
	theme(legend.title=element_text(size=12),legend.text=element_text(size=10))+
	theme(strip.text.x=element_text(size=12))
dev.off()

### phenotype correlation and pca
lig<-fread(paste0('results/pheno_',pop[1],'.txt'),data.table=F)
lig$vp<-log(lig$vp+1)
lig$varroaMitoRatio<-lig$varroaMitoRatio^(1/4)
lig$charge_varroa<-lig$charge_varroa^(1/4)
lig$sp<-'lig'
mel<-fread(paste0('results/pheno_',pop[2],'.txt'),data.table=F)
mel$vp<-log(mel$vp+1)
mel$varroaMitoRatio<-mel$varroaMitoRatio^(1/4)
mel$charge_varroa<-mel$charge_varroa^(1/4)
mel$sp<-'mel'
hyb<-fread(paste0('results/pheno_',pop[3],'.txt'),data.table=F)
hyb$vp<-log(hyb$vp+1)
hyb$varroaMitoRatio<-hyb$varroaMitoRatio^(1/4)
hyb$charge_varroa<-hyb$charge_varroa^(1/4)
hyb$sp<-'hyb'
df<-do.call(rbind,list(lig,mel,hyb))

# pc + correlations 6 pheno
df_pc<-data.frame(v_pho=df$vp,v_mito=df$varroaMitoRatio,v_inf=df$logit_taux_inf,v_load=df$charge_varroa,MNR=df$eb_smr,recap_inf=df$logit_taux_reop_spe,sp=df$sp)
nb=estim_ncpPCA(df_pc[,1:6],ncp.max=ncol(df_pc))
res=imputePCA(df_pc[,1:6],ncp=nb$ncp)
res<-res$completeObs
nb=estim_ncpPCA(df_pc[df_pc$sp=='lig',1:6],ncp.max=ncol(df_pc))
res_lig=imputePCA(df_pc[df_pc$sp=='lig',1:6],ncp=nb$ncp)
res_lig<-res_lig$completeObs
nb=estim_ncpPCA(df_pc[df_pc$sp=='mel',1:6],ncp.max=ncol(df_pc))
res_mel=imputePCA(df_pc[df_pc$sp=='mel',1:6],ncp=nb$ncp)
res_mel<-res_mel$completeObs
nb=estim_ncpPCA(df_pc[df_pc$sp=='hyb',1:6],ncp.max=ncol(df_pc))
res_hyb=imputePCA(df_pc[df_pc$sp=='hyb',1:6],ncp=nb$ncp)
res_hyb<-res_hyb$completeObs
pca=PCA(res,graph=F) 
pca_lig=PCA(res_lig,graph=F) 
pca_mel=PCA(res_mel,graph=F) 
pca_hyb=PCA(res_hyb,graph=F) 
cor<-round(cor(res),2)
cor_lig<-round(cor(res_lig),2)
cor_lig<-reshape2::melt(cor_lig)
cor_lig$sp<-'lig'
cor_mel<-round(cor(res_mel),2)
cor_mel<-reshape2::melt(cor_mel)
cor_mel$sp<-'mel'
cor_hyb<-round(cor(res_hyb),2)
cor_hyb<-reshape2::melt(cor_hyb)
cor_hyb$sp<-'hyb'
cor_comb<-do.call(rbind,list(cor_lig,cor_mel,cor_hyb))
cor_comb$Var1<-paste0(cor_comb$Var1,'_',cor_comb$sp)
cor_comb$Var2<-paste0(cor_comb$Var2,'_',cor_comb$sp)
cor_comb$sp<-NULL
cor_comb<-reshape(cor_comb,idvar='Var1',timevar='Var2',direction='wide')
rownames(cor_comb)<-cor_comb$Var1
cor_comb$Var1<-NULL
colnames(cor_comb)<-gsub('value.','',colnames(cor_comb))
cor_comb[is.na(cor_comb)]<-0
cor_comb<-as.matrix(cor_comb)

#corrplot(cor,type='upper',tl.col="black",tl.srt=45,cl.pos='r',col=brewer.pal(n=8,name='RdYlGn'),order='original',cl.cex=1.5,tl.cex=1.5)
#grid.echo()
#p1<-grid.grab()
#p1<-grab_grob()
p1<-ggcorrplot(cor,type='lower',tl.srt=45,
	show.diag=TRUE,ggtheme=ggplot2::theme_minimal,
	method="circle",hc.order=FALSE,colors=c('#D73027','#FEE08B','#1A9850'),
	lab_col="black",lab_size=10,tl.cex=20,tl.col="black")+
	scale_fill_gradient2(low="#D73027",high="#1A9850",mid="#FEE08B",limit=c(-1,1))+ 
	labs(fill="Correlation")+
	theme(legend.title=element_text(size=25),legend.text=element_text(size=15))
p2<-fviz_screeplot(pca,ncp=10)+ theme_minimal()+ labs(title = "Variances PCA",x = "Principal Components", y = "% of variances")+theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),plot.title=element_text(size=30))
p3<-fviz_pca_var(pca,col.var="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel=T,axes=c(1,2),arrowsize=1,labelsize=5)+ labs(title="",x=paste0("PC1 (",round(pca$eig[1,2],digit=2),"%)"),y=paste0("PC2 (",round(pca$eig[2,2],digit=2),"%)"))+theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),legend.title=element_text(size=20),legend.text=element_text(size=15))
p4<-fviz_pca_var(pca,col.var="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel=T,axes=c(1,3),arrowsize=1,labelsize=5)+ labs(title="",x=paste0("PC1 (",round(pca$eig[1,2],digit=2),"%)"),y=paste0("PC3 (",round(pca$eig[3,2],digit=2),"%)"))+theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),legend.title=element_text(size=20),legend.text=element_text(size=15))
p5<-fviz_pca_var(pca,col.var="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel=T,axes=c(2,3),arrowsize=1,labelsize=5)+ labs(title="",x=paste0("PC2 (",round(pca$eig[2,2],digit=2),"%)"),y=paste0("PC3 (",round(pca$eig[3,2],digit=2),"%)"))+theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),legend.title=element_text(size=20),legend.text=element_text(size=15))
png(paste0('plot_paper_pca6pheno_',version,'.png'),width=5000,height=4000,res=400)
grid.arrange(p1,p2,p3,p4,p5,ncol=2,layout_matrix=cbind(c(1,1,2),c(3,4,5)))
dev.off()

#corrplot(cor_comb,type='upper',tl.col="black",tl.srt=45,cl.pos='r',col=brewer.pal(n=8,name='RdYlGn'),order='original',cl.cex=1.5,tl.cex=1.2)%>%corrRect(c(1,7,13,18))
#grid.echo()
#p1<-grid.grab()
p1<-ggcorrplot(cor_comb,type='lower',tl.srt=45,
	show.diag=TRUE,ggtheme=ggplot2::theme_minimal,
	method="circle",hc.order=FALSE,colors=c('#D73027','#FEE08B','#1A9850'),
	lab_col="black",lab_size=10,tl.cex=20,tl.col="black")+
	scale_fill_gradient2(low="#D73027",high="#1A9850",mid="#FEE08B",limit=c(-1,1))+ 
	labs(fill="Corr")+
	theme(legend.title=element_text(size=25),legend.text=element_text(size=15))
p2<-fviz_screeplot(pca_lig,ncp=10)+ theme_minimal()+ labs(title = "Variances PCA - Ligustica \n& Carnica",x = "Principal Components", y = "% of variances")+theme(axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),plot.title=element_text(size=20))
p3<-fviz_screeplot(pca_mel,ncp=10)+ theme_minimal()+ labs(title = "Variances PCA - Mellifera",x = "Principal Components", y = "% of variances")+theme(axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),plot.title=element_text(size=20))
p4<-fviz_screeplot(pca_hyb,ncp=10)+ theme_minimal()+ labs(title = "Variances PCA - Hybrids",x = "Principal Components", y = "% of variances")+theme(axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),plot.title=element_text(size=20))
p5<-fviz_pca_var(pca_lig,col.var="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel=T,axes=c(1,2),arrowsize=1,labelsize=5)+ labs(title = "PCA - Ligustica \n& Carnica",x=paste0("PC1 (",round(pca_lig$eig[1,2],digit=2),"%)"),y=paste0("PC2 (",round(pca_lig$eig[2,2],digit=2),"%)"))+theme(axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),legend.title=element_text(size=12),legend.text=element_text(size=10),plot.title=element_text(size=20))
p6<-fviz_pca_var(pca_mel,col.var="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel=T,axes=c(1,2),arrowsize=1,labelsize=5)+ labs(title = "PCA - Mellifera",x=paste0("PC1 (",round(pca_mel$eig[1,2],digit=2),"%)"),y=paste0("PC2 (",round(pca_mel$eig[2,2],digit=2),"%)"))+theme(axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),legend.title=element_text(size=12),legend.text=element_text(size=10),plot.title=element_text(size=20))
p7<-fviz_pca_var(pca_hyb,col.var="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel=T,axes=c(1,2),arrowsize=1,labelsize=5)+ labs(title = "PCA - Hybrids",x=paste0("PC1 (",round(pca_hyb$eig[1,2],digit=2),"%)"),y=paste0("PC2 (",round(pca_hyb$eig[2,2],digit=2),"%)"))+theme(axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),legend.title=element_text(size=12),legend.text=element_text(size=10),plot.title=element_text(size=20))
g1<-ggarrange(p2,p3,p4,ncol=1)
g2<-ggarrange(p5,p6,p7,ncol=1)
png(paste0('plot_paper_pca6pheno_detail_sp_',version,'.png'),width=6000,height=4000,res=400)
grid.arrange(p1,g1,g2,ncol=3,widths=c(2,1,1))
dev.off()

# pc + correlations 3 pheno
df_pc<-data.frame(pc1_varroa_inf=df$pca,MNR=df$eb_smr,recap_inf=df$logit_taux_reop_spe,sp=df$sp)
nb=estim_ncpPCA(df_pc[,1:3],ncp.max=ncol(df_pc))
res=imputePCA(df_pc[,1:3],ncp=nb$ncp)
res<-res$completeObs
nb=estim_ncpPCA(df_pc[df_pc$sp=='lig',1:3],ncp.max=ncol(df_pc))
res_lig=imputePCA(df_pc[df_pc$sp=='lig',1:3],ncp=nb$ncp)
res_lig<-res_lig$completeObs
nb=estim_ncpPCA(df_pc[df_pc$sp=='mel',1:3],ncp.max=ncol(df_pc))
res_mel=imputePCA(df_pc[df_pc$sp=='mel',1:3],ncp=nb$ncp)
res_mel<-res_mel$completeObs
nb=estim_ncpPCA(df_pc[df_pc$sp=='hyb',1:3],ncp.max=ncol(df_pc))
res_hyb=imputePCA(df_pc[df_pc$sp=='hyb',1:3],ncp=nb$ncp)
res_hyb<-res_hyb$completeObs
pca=PCA(res,graph=F) 
pca_lig=PCA(res_lig,graph=F) 
pca_mel=PCA(res_mel,graph=F) 
pca_hyb=PCA(res_hyb,graph=F) 
cor<-round(cor(res),2)
cor_lig<-round(cor(res_lig),2)
cor_lig<-reshape2::melt(cor_lig)
cor_lig$sp<-'lig'
cor_mel<-round(cor(res_mel),2)
cor_mel<-reshape2::melt(cor_mel)
cor_mel$sp<-'mel'
cor_hyb<-round(cor(res_hyb),2)
cor_hyb<-reshape2::melt(cor_hyb)
cor_hyb$sp<-'hyb'
cor_comb<-do.call(rbind,list(cor_lig,cor_mel,cor_hyb))
cor_comb$Var1<-paste0(cor_comb$Var1,'_',cor_comb$sp)
cor_comb$Var2<-paste0(cor_comb$Var2,'_',cor_comb$sp)
cor_comb$sp<-NULL
cor_comb<-reshape(cor_comb,idvar='Var1',timevar='Var2',direction='wide')
rownames(cor_comb)<-cor_comb$Var1
cor_comb$Var1<-NULL
colnames(cor_comb)<-gsub('value.','',colnames(cor_comb))
cor_comb[is.na(cor_comb)]<-0
cor_comb<-as.matrix(cor_comb)
#corrplot(cor,type='upper',tl.col="black",tl.srt=45,cl.pos='r',col=brewer.pal(n=8,name='RdYlGn'),order='original',cl.cex=1.5,tl.cex=1.5)
#grid.echo()
#p1<-grid.grab()
p1<-ggcorrplot(cor,type='lower',tl.srt=45,
	show.diag=TRUE,ggtheme=ggplot2::theme_minimal,
	method="circle",hc.order=FALSE,colors=c('#D73027','#FEE08B','#1A9850'),
	lab_col="black",lab_size=10,tl.cex=20,tl.col="black")+
	scale_fill_gradient2(low="#D73027",high="#1A9850",mid="#FEE08B",limit=c(-1,1))+ 
	labs(fill="Correlation")+
	theme(legend.title=element_text(size=25),legend.text=element_text(size=15))
p2<-fviz_screeplot(pca,ncp=10)+ theme_minimal()+ labs(title = "Variances PCA",x = "Principal Components", y = "% of variances")+theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),plot.title=element_text(size=30))
p3<-fviz_pca_var(pca,col.var="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel=T,axes=c(1,2),arrowsize=1,labelsize=5)+ labs(title="",x=paste0("PC1 (",round(pca$eig[1,2],digit=2),"%)"),y=paste0("PC2 (",round(pca$eig[2,2],digit=2),"%)"))+theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),legend.title=element_text(size=20),legend.text=element_text(size=15))
p4<-fviz_pca_var(pca,col.var="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel=T,axes=c(1,3),arrowsize=1,labelsize=5)+ labs(title="",x=paste0("PC1 (",round(pca$eig[1,2],digit=2),"%)"),y=paste0("PC3 (",round(pca$eig[3,2],digit=2),"%)"))+theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),legend.title=element_text(size=20),legend.text=element_text(size=15))
p5<-fviz_pca_var(pca,col.var="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel=T,axes=c(2,3),arrowsize=1,labelsize=5)+ labs(title="",x=paste0("PC2 (",round(pca$eig[2,2],digit=2),"%)"),y=paste0("PC3 (",round(pca$eig[3,2],digit=2),"%)"))+theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),legend.title=element_text(size=20),legend.text=element_text(size=15))
png(paste0('plot_paper_pca3pheno_',version,'.png'),width=6000,height=4000,res=400)
grid.arrange(p1,p2,p3,p4,p5,ncol=2,layout_matrix=cbind(c(1,1,2),c(3,4,5)))
dev.off()

#corrplot(cor_comb,type='upper',tl.col="black",tl.srt=45,cl.pos='r',col=brewer.pal(n=8,name='RdYlGn'),order='original',cl.cex=1.5,tl.cex=1.2)%>%corrRect(c(1,4,7,9))
#grid.echo()
#p1<-grid.grab()
p1<-ggcorrplot(cor_comb,type='lower',tl.srt=45,
	show.diag=TRUE,ggtheme=ggplot2::theme_minimal,
	method="circle",hc.order=FALSE,colors=c('#D73027','#FEE08B','#1A9850'),
	lab_col="black",lab_size=10,tl.cex=20,tl.col="black")+
	scale_fill_gradient2(low="#D73027",high="#1A9850",mid="#FEE08B",limit=c(-1,1))+ 
	labs(fill="Corr")+
	theme(legend.title=element_text(size=25),legend.text=element_text(size=15))
p2<-fviz_screeplot(pca_lig,ncp=10)+ theme_minimal()+ labs(title = "Variances PCA - Ligustica \n& Carnica",x = "Principal Components", y = "% of variances")+theme(axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),plot.title=element_text(size=20))
p3<-fviz_screeplot(pca_mel,ncp=10)+ theme_minimal()+ labs(title = "Variances PCA - Mellifera",x = "Principal Components", y = "% of variances")+theme(axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),plot.title=element_text(size=20))
p4<-fviz_screeplot(pca_hyb,ncp=10)+ theme_minimal()+ labs(title = "Variances PCA - Hybrids",x = "Principal Components", y = "% of variances")+theme(axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),plot.title=element_text(size=20))
p5<-fviz_pca_var(pca_lig,col.var="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel=T,axes=c(1,2),arrowsize=1,labelsize=5)+ labs(title = "PCA - Ligustica \n& Carnica",x=paste0("PC1 (",round(pca_lig$eig[1,2],digit=2),"%)"),y=paste0("PC2 (",round(pca_lig$eig[2,2],digit=2),"%)"))+theme(axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),legend.title=element_text(size=12),legend.text=element_text(size=10),plot.title=element_text(size=20))
p6<-fviz_pca_var(pca_mel,col.var="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel=T,axes=c(1,2),arrowsize=1,labelsize=5)+ labs(title = "PCA - Mellifera",x=paste0("PC1 (",round(pca_mel$eig[1,2],digit=2),"%)"),y=paste0("PC2 (",round(pca_mel$eig[2,2],digit=2),"%)"))+theme(axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),legend.title=element_text(size=12),legend.text=element_text(size=10),plot.title=element_text(size=20))
p7<-fviz_pca_var(pca_hyb,col.var="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel=T,axes=c(1,2),arrowsize=1,labelsize=5)+ labs(title = "PCA - Hybrids",x=paste0("PC1 (",round(pca_hyb$eig[1,2],digit=2),"%)"),y=paste0("PC2 (",round(pca_hyb$eig[2,2],digit=2),"%)"))+theme(axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),legend.title=element_text(size=12),legend.text=element_text(size=10),plot.title=element_text(size=20))
g1<-ggarrange(p2,p3,p4,ncol=1)
g2<-ggarrange(p5,p6,p7,ncol=1)
png(paste0('plot_paper_pca3pheno_detail_sp_',version,'.png'),width=6000,height=4000,res=400)
grid.arrange(p1,g1,g2,ncol=3,widths=c(2,1,1))
dev.off()

### phenotype distribution
lig<-fread(paste0('results/pheno_',pop[1],'.txt'),data.table=F)
table(lig$region)
lig<-lig[,c('num_ruche_bs','vp','varroaMitoRatio','logit_taux_inf','charge_varroa','eb_smr','logit_taux_reop_spe','pca')]
lig$vp<-log(lig$vp+1)
lig$varroaMitoRatio<-lig$varroaMitoRatio^(1/4)
lig$charge_varroa<-lig$charge_varroa^(1/4)
lig<-gather(lig,pheno,val,c('vp','varroaMitoRatio','logit_taux_inf','charge_varroa','pca','eb_smr','logit_taux_reop_spe'))
lig$sp<-'lig'
mel<-fread(paste0('results/pheno_',pop[2],'.txt'),data.table=F)
table(mel$region)
mel<-mel[,c('num_ruche_bs','vp','varroaMitoRatio','logit_taux_inf','charge_varroa','eb_smr','logit_taux_reop_spe','pca')]
mel$vp<-log(mel$vp+1)
mel$varroaMitoRatio<-mel$varroaMitoRatio^(1/4)
mel$charge_varroa<-mel$charge_varroa^(1/4)
mel<-gather(mel,pheno,val,c('vp','varroaMitoRatio','logit_taux_inf','charge_varroa','pca','eb_smr','logit_taux_reop_spe'))
mel$sp<-'mel'
hyb<-fread(paste0('results/pheno_',pop[3],'.txt'),data.table=F)
table(hyb$region)
hyb<-hyb[,c('num_ruche_bs','vp','varroaMitoRatio','logit_taux_inf','charge_varroa','eb_smr','logit_taux_reop_spe','pca')]
hyb$vp<-log(hyb$vp+1)
hyb$varroaMitoRatio<-hyb$varroaMitoRatio^(1/4)
hyb$charge_varroa<-hyb$charge_varroa^(1/4)
hyb<-gather(hyb,pheno,val,c('vp','varroaMitoRatio','logit_taux_inf','charge_varroa','pca','eb_smr','logit_taux_reop_spe'))
hyb$sp<-'hyb'

df<-rbind(lig,mel,hyb)
df$sp<-factor(df$sp,levels=c('lig','mel','hyb'))
df$pheno<-factor(df$pheno,levels=c('vp','varroaMitoRatio','logit_taux_inf','charge_varroa','pca','eb_smr','logit_taux_reop_spe'))
pheno_names<-c(
	`vp` = "v_pho *log(x+1)*",
	`varroaMitoRatio` = "v_mito *x^0.25*",
	`logit_taux_inf` = "v_inf *logit(x)*",
	`charge_varroa` = "v_load *x^0.25*",
	`pca`="varroa_inf *pc1*",
	`eb_smr` = "MNR *EB(x)*",
	`logit_taux_reop_spe` = "recap_inf *logit(x)*")
sp_names<-c(
  `lig` = "Ligustica & Carnica",
  `mel` = "Mellifera",
  `hyb` = "Hybrids")
png(paste0('plot_paper_distrib6pheno_',version,'.png'),width=6000,height=3000,res=400)
ggplot(df[df$pheno!='pc1',])+
  geom_histogram(aes(x=val))+
  facet_grid(sp~pheno,scale='free',labeller=labeller(pheno=as_labeller(pheno_names),sp=as_labeller(sp_names)))+
  theme_bw()+
  theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.text.x=element_text(size=15,angle = 45, vjust = 1, hjust=1),axis.text.y=element_text(size=15))+
  theme(strip.text.y=element_text(size=15),strip.text.x=element_text(size=15))+
#  mdthemes::md_theme_bw()+
  xlab('phenotype value after transformation')+ylab('count')
dev.off()
png(paste0('plot_paper_distrib3pheno_',version,'.png'),width=3500,height=3000,res=400)
ggplot(df[!(df$pheno%in%c('vp','varroaMitoRatio','logit_taux_inf','charge_varroa')),])+
  geom_histogram(aes(x=val))+
  facet_grid(sp~pheno,scale='free',labeller=labeller(pheno=as_labeller(pheno_names),sp=as_labeller(sp_names)))+
  theme_bw()+
  theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.text.x=element_text(size=15,angle = 45, vjust = 1, hjust=1),axis.text.y=element_text(size=15))+
  theme(strip.text.y=element_text(size=15),strip.text.x=element_text(size=15))+
#  mdthemes::md_theme_bw()+
  xlab('phenotype value after transformation')+ylab('count')
dev.off()

### pve(se)
pheno<-c('varroaphoretic','varroadepthmitoracine4','logit_varroainfestation','chargevarroaracine4','pc1','eb_smr','logit_taux_reop_inf')
#sp<-c('Ligustica_Carnica','Mellifera','hybrid')
SP<-vector()
PHENO<-vector()
PVE_lmm<-vector()
SEPVE_lmm<-vector()
PVE_bslmm<-vector()
SEPVE_bslmm<-vector()
PGE_bslmm<-vector()
SEPGE_bslmm<-vector()
h_bslmm<-vector()
SEh_bslmm<-vector()
pi_bslmm<-vector()
SEpi_bslmm<-vector()
ngamma_bslmm<-vector()
SEngamma_bslmm<-vector()
rho_bslmm<-vector()
SErho_bslmm<-vector()
for(i in 1:length(pop)){
	for(j in 1:length(pheno)){
		x<-readLines(paste0('results/gemma_',pop[i],'_',pheno[j],'_freq_lmm_egs_cov.log.txt'))
		SP<-c(SP,pop[i])
		PHENO<-c(PHENO,pheno[j])
		n<-grep('pve',x)
		pve<-as.numeric(tail(unlist(strsplit(x[n[1]],' ')),1))
		PVE_lmm<-c(PVE_lmm,pve)
		se_pve<-as.numeric(tail(unlist(strsplit(x[n[2]],' ')),1))
		SEPVE_lmm<-c(SEPVE_lmm,se_pve)
		y<-read.table(paste0('results/gemma_',pop[i],'_',pheno[j],'_freq_bslmm_egs.hyp.txt'),header=T)
		PVE_bslmm<-c(PVE_bslmm,mean(y$pve))
		SEPVE_bslmm<-c(SEPVE_bslmm,sd(y$pve))
		PGE_bslmm<-c(PGE_bslmm,mean(y$pge))
		SEPGE_bslmm<-c(SEPGE_bslmm,sd(y$pge))
		h_bslmm<-c(h_bslmm,mean(y$h))
		SEh_bslmm<-c(SEh_bslmm,sd(y$h))
		pi_bslmm<-c(pi_bslmm,mean(y$pi))
		SEpi_bslmm<-c(SEpi_bslmm,sd(y$pi))
		ngamma_bslmm<-c(ngamma_bslmm,mean(y$n_gamma))
		SEngamma_bslmm<-c(SEngamma_bslmm,sd(y$n_gamma))
		rho_bslmm<-c(rho_bslmm,mean(y$rho))
		SErho_bslmm<-c(SErho_bslmm,sd(y$rho))
		}
	}
X<-data.frame(SP,PHENO,PVE_lmm,SEPVE_lmm,PVE_bslmm,SEPVE_bslmm,PGE_bslmm,SEPGE_bslmm,h_bslmm,SEh_bslmm,pi_bslmm,SEpi_bslmm,ngamma_bslmm,SEngamma_bslmm,rho_bslmm,SErho_bslmm)
write.table(X,paste0('summary_pve_pge_h_',version,'.txt'),col.names=T,row.names=F,quote=F)
X$PVE_lmm_min<-X$PVE_lmm-1.96*X$SEPVE_lmm
X$PVE_lmm_max<-X$PVE_lmm+1.96*X$SEPVE_lmm
X$PVE_bslmm_min<-X$PVE_bslmm-1.96*X$SEPVE_bslmm
X$PVE_bslmm_max<-X$PVE_bslmm+1.96*X$SEPVE_bslmm
X$PGE_bslmm_min<-X$PGE_bslmm-1.96*X$SEPGE_bslmm
X$PGE_bslmm_max<-X$PGE_bslmm+1.96*X$SEPGE_bslmm
X$h_bslmm_min<-X$h_bslmm-1.96*X$SEh_bslmm
X$h_bslmm_max<-X$h_bslmm+1.96*X$SEh_bslmm
X$pi_bslmm_min<-X$pi_bslmm-1.96*X$SEpi_bslmm
X$pi_bslmm_max<-X$pi_bslmm+1.96*X$SEpi_bslmm
X$ngamma_bslmm_min<-X$ngamma_bslmm-1.96*X$SEngamma_bslmm
X$ngamma_bslmm_max<-X$ngamma_bslmm+1.96*X$SEngamma_bslmm
X$rho_bslmm_min<-X$rho_bslmm-1.96*X$SErho_bslmm
X$rho_bslmm_max<-X$rho_bslmm+1.96*X$SErho_bslmm
X$PHENO[X$PHENO=='varroaphoretic']<-'v_pho'
X$PHENO[X$PHENO=='varroadepthmitoracine4']<-'v_mito'
X$PHENO[X$PHENO=='logit_varroainfestation']<-'v_inf'
X$PHENO[X$PHENO=='chargevarroaracine4']<-'v_load'
X$PHENO[X$PHENO=='pc1']<-'pc1_varroa_inf'
X$PHENO[X$PHENO=='eb_smr']<-'MNR'
X$PHENO[X$PHENO=='logit_taux_reop_inf']<-'recap_inf'
X$PHENO<-factor(X$PHENO,levels=c('v_pho','v_mito','v_inf','v_load','pc1_varroa_inf','MNR','recap_inf'))
X$SP[X$SP=='Ligustica_nus']<-'Ligustica_Carnica'
X$SP<-factor(X$SP,levels=c('Ligustica_Carnica','Mellifera','hybrid'))
sp_names<-c(
  `Ligustica_Carnica` = "Ligustica & Carnica",
  `Mellifera` = "Mellifera",
  `hybrid` = "Hybrids")
  
X_plot_lmm<-X[,c("SP","PHENO","PVE_lmm","PVE_lmm_min","PVE_lmm_max")]  
colnames(X_plot_lmm)<-c("SP","PHENO","PVE","min","max")
X_plot_lmm$type<-'lmm'
X_plot_bslmm1<-X[,c("SP","PHENO","PVE_bslmm","PVE_bslmm_min","PVE_bslmm_max")]  
colnames(X_plot_bslmm1)<-c("SP","PHENO","PVE","min","max")
X_plot_bslmm1$type<-'bslmm'
X_plot_bslmm2<-X[,c("SP","PHENO","PGE_bslmm","PGE_bslmm_min","PGE_bslmm_max")]  
colnames(X_plot_bslmm2)<-c("SP","PHENO","PVE","min","max")
X_plot_bslmm2$type<-'pge'
#X_plot_bslmm3<-X[,c("SP","PHENO","h_bslmm","h_bslmm_min","h_bslmm_max")]  
#colnames(X_plot_bslmm3)<-c("SP","PHENO","PVE","min","max")
#X_plot_bslmm3$type<-'h'
#X_plot<-do.call(rbind,list(X_plot_lmm,X_plot_bslmm1,X_plot_bslmm2,X_plot_bslmm3))
X_plot<-do.call(rbind,list(X_plot_lmm,X_plot_bslmm1,X_plot_bslmm2))
X_plot$g<-paste0(X_plot$SP,'_',X_plot$type)
#X_plot$g<-factor(X_plot$g,levels=c('hybrid_corse_h','hybrid_corse_pge','hybrid_corse_bslmm','hybrid_corse_lmm','Mellifera_ncorse_h','Mellifera_ncorse_pge','Mellifera_ncorse_bslmm','Mellifera_ncorse_lmm','Ligustica_Carnica_h','Ligustica_Carnica_pge','Ligustica_Carnica_bslmm','Ligustica_Carnica_lmm'))
#X_plot$g<-factor(X_plot$g,levels=c('hybrid_corse_pge','hybrid_corse_bslmm','hybrid_corse_lmm','Mellifera_ncorse_pge','Mellifera_ncorse_bslmm','Mellifera_ncorse_lmm','Ligustica_Carnica_pge','Ligustica_Carnica_bslmm','Ligustica_Carnica_lmm'))
X_plot$g<-factor(X_plot$g,levels=c('hybrid_pge','hybrid_bslmm','hybrid_lmm','Mellifera_pge','Mellifera_bslmm','Mellifera_lmm','Ligustica_Carnica_pge','Ligustica_Carnica_bslmm','Ligustica_Carnica_lmm'))
X_plot$type<-factor(X_plot$type,levels=c('lmm','bslmm','pge'))

#png('plot_paper_pve6pheno.png',width=7000,height=1500,res=400)
#ggplot(X_plot[X_plot$PHENO!='pc1_varroa_inf',])+
#	geom_segment(aes(x=min,xend=max,y=g,yend=g,col=SP,linetype=type),size=1)+
#	geom_point(aes(x=PVE,col=SP,y=g),size=2)+
#	facet_grid(.~PHENO,scale='free_y',labeller=labeller(sp=as_labeller(sp_names)))+
#	scale_colour_manual(values=c('goldenrod2','grey34','deepskyblue2'),labels=c('Ligustica & Carnica','Mellifera','Hybrids'),name='Groups')+
#	geom_vline(aes(xintercept=0),col='red',lty=3)+
#	geom_vline(aes(xintercept=1),col='red',lty=3)+
#	xlab('value with 95% confidence interval')+	
#	scale_y_discrete(labels=c('H','PGE','PVE_bslmm','PVE_lmm','H','PGE','PVE_bslmm','PVE_lmm','H','PGE','PVE_bslmm','PVE_lmm'))+
#	scale_linetype_manual(values=c('solid','longdash','dashed','dotted'),guide=F)+
#	theme_bw()+
#	theme(axis.title.y=element_blank())
#dev.off()
#png('plot_paper_pve3pheno.png',width=5000,height=1500,res=400)
#ggplot(X_plot[!(X_plot$PHENO%in%c('v_pho','v_mito','v_inf','v_load')),])+
#	geom_segment(aes(x=min,xend=max,y=g,yend=g,col=SP,linetype=type),size=1)+
#	geom_point(aes(x=PVE,col=SP,y=g),size=2)+
#	facet_grid(.~PHENO,scale='free_y',labeller=labeller(sp=as_labeller(sp_names)))+
#	scale_colour_manual(values=c('goldenrod2','grey34','deepskyblue2'),labels=c('Ligustica & Carnica','Mellifera','Hybrids'),name='Groups')+
#	geom_vline(aes(xintercept=0),col='red',lty=3)+
#	geom_vline(aes(xintercept=1),col='red',lty=3)+
#	xlab('value with 95% confidence interval')+	
#	scale_y_discrete(labels=c('H','PGE','PVE_bslmm','PVE_lmm','H','PGE','PVE_bslmm','PVE_lmm','H','PGE','PVE_bslmm','PVE_lmm'))+
#	scale_linetype_manual(values=c('solid','longdash','dashed','dotted'),guide=F)+
#	theme_bw()+
#	theme(axis.title.y=element_blank())
#dev.off()

png(paste0('plot_paper_pve6pheno_',version,'.png'),width=7000,height=1500,res=400)
ggplot(X_plot[X_plot$PHENO!='pc1_varroa_inf',])+
	geom_segment(aes(x=min,xend=max,y=g,yend=g,col=SP,linetype=type),size=1)+
	geom_point(aes(x=PVE,col=SP,y=g),size=2)+
	facet_grid(.~PHENO,scale='free_y',labeller=labeller(sp=as_labeller(sp_names)))+
	scale_colour_manual(values=c('goldenrod2','grey34','deepskyblue2'),labels=c('Ligustica & Carnica','Mellifera','Hybrids'),name='Groups')+
	geom_vline(aes(xintercept=0),col='red',lty=3)+
	geom_vline(aes(xintercept=1),col='red',lty=3)+
	xlab('value with 95% confidence interval')+	
	scale_y_discrete(labels=c('PGE','PVE_bslmm','PVE_lmm','PGE','PVE_bslmm','PVE_lmm','PGE','PVE_bslmm','PVE_lmm'))+
	scale_linetype_manual(values=c('solid','longdash','dashed'),guide='none')+
	theme_bw()+
	theme(axis.title.y=element_blank(),axis.title.x=element_text(size=20),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))+
	theme(strip.text.x=element_text(size=15),legend.title=element_text(size=20),legend.text=element_text(size=15))
dev.off()
png(paste0('plot_paper_pve3pheno_',version,'.png'),width=5000,height=1500,res=400)
ggplot(X_plot[!(X_plot$PHENO%in%c('v_pho','v_mito','v_inf','v_load')),])+
	geom_segment(aes(x=min,xend=max,y=g,yend=g,col=SP,linetype=type),size=1)+
	geom_point(aes(x=PVE,col=SP,y=g),size=2)+
	facet_grid(.~PHENO,scale='free_y',labeller=labeller(sp=as_labeller(sp_names)))+
	scale_colour_manual(values=c('goldenrod2','grey34','deepskyblue2'),labels=c('Ligustica & Carnica','Mellifera','Hybrids'),name='Groups')+
	geom_vline(aes(xintercept=0),col='red',lty=3)+
	geom_vline(aes(xintercept=1),col='red',lty=3)+
	xlab('value with 95% confidence interval')+	
	scale_y_discrete(labels=c('PGE','PVE_bslmm','PVE_lmm','PGE','PVE_bslmm','PVE_lmm','PGE','PVE_bslmm','PVE_lmm'))+
	scale_linetype_manual(values=c('solid','longdash','dashed'),guide='none')+
	theme_bw()+
	theme(axis.title.y=element_blank(),axis.title.x=element_text(size=20),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))+
	theme(strip.text.x=element_text(size=15),legend.title=element_text(size=20),legend.text=element_text(size=15))
dev.off()
