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
args<-commandArgs(TRUE)

pdf('plot.pdf',width=15,height=15)

q<-list.files(path='.',pattern='_st_het.Q')
Q<-list()
for(i in 1:length(q)){
	name<-gsub('sim_depth_count','',q[i])
	name<-gsub('_st_het.Q','',name)
	x<-fread(q[i])
	colnames(x)<-c('subsp','q')
	x$id<-name
	Q[[i]]<-x
}
Q<-do.call(rbind,Q)

Q<-spread(Q,subsp,q)
a<-colsplit(Q$id,'_',c('v1','v2','v3'))
Q$id<-paste0(a$v1,'_',a$v2)
Q$type<-a$v3
Q$type[Q$type=='']<-'all'

Q$k<-'hybrid'
Q$k[Q$Caucasica>=0.8]<-'C'
Q$k[Q$Ligustica_Carnica>=0.8]<-'L'
Q$k[Q$Mellifera>=0.8]<-'M'

ggplot(Q)+geom_point(aes(x=Ligustica_Carnica,y=Mellifera,col=k))+facet_wrap(~type)+theme(legend.position='none')

DFall<-Q[Q$type=='all',c('id','k')]
DF50<-Q[Q$type=='50k',c('id','k')]
DF<-merge(DFall,DF50,by='id',all=T)
DF$ok<-'ok'
DF$ok[DF$k.x!=DF$k.y]<-'not_ok'
DF[,c('k.x','k.y')]<-NULL

pheno<-fread('Data_BS.csv')
a<-colsplit(pheno$num_ruche_bs,'_',c('v1','v2'))
pheno$id<-paste0(toupper(a$v1),'_',a$v2)
mito<-fread('compareDepthsBeeVarroa.txt')
a<-colsplit(mito$name,'-',c('v1','v2','v3'))
mito$id<-paste0(toupper(a$v1),'_',a$v2)

DF<-Reduce(function(x,y) merge(x,y,by='id',all=T),list(DF,Q,pheno,mito))

df<-subset(DF,DF$type=='50k')
df$k[df$pass=='us']<-'us'
df$smr_brut<-df$nbr_inf_nr/df$nbr_inf_une_fond
eb_smr_fit=fitdistr(x=df$smr_brut[!is.na(df$smr_brut) & df$smr_brut>0 & df$smr_brut<1],densfun="beta",start=list(shape1=1,shape2=1),method="L-BFGS-B")
aprior=eb_smr_fit$estimate[1]
bprior=eb_smr_fit$estimate[2]
df$eb_smr=(aprior+df$nbr_inf_nr)/(aprior+bprior+df$nbr_inf_une_fond)
m_vp<-mean(df$nbr_varroa_100abeilles[!is.na(df$nbr_varroa_100abeilles)])
vp_ini<-df$nbr_varroa_100abeilles
p<-pbinom(vp_ini,100,m_vp/100)
df$vp<-qnorm(p)
df$vp[df$vp%in%c('-Inf','Inf')]<-NA
df$taux_inf<-df$nbr_inf/df$nbr_ouvert
df$reop<-df$nbr_reop/df$nbr_ouvert
df$reop_inf<-df$nbr_inf_reop/df$nbr_inf

pca<-prcomp(df[,c('Ligustica_Carnica','Mellifera','Caucasica')])
pca<-pca$x
pca<-as.data.frame(pca)
pca$id<-df$num_ruche_bs
pca$rucher<-df$rucher
pca$date<-df$date_pre
pca$k<-df$k
ggplot(pca)+
	geom_point(aes(x=PC1,y=PC2,col=k))

df_pca<-df[,c('nbr_abeilles','nbr_cell_cf','nbr_cell_co','nbr_cell_male','surf_res','surf_pol')]
nb=estim_ncpPCA(df_pca,ncp.max=ncol(df))
res.comp=imputePCA(df_pca,ncp=nb$ncp)
pca.res=PCA(res.comp$completeObs,graph=F) 
pca_var<-pca.res$eig
pca<-pca.res$ind$coord
pca<-as.data.frame(pca)
pca$id<-df$num_ruche_bs
pca$k<-df$k
g1<-fviz_pca_var(pca.res,col.var="contrib")+theme_minimal()
g2<-ggplot(pca)+
	geom_point(aes(x=Dim.1,y=Dim.2,col=k),alpha=0.3)+
	xlab(paste0('PC1 (',round(pca_var[1,2],1),')'))+ylab(paste0('PC2 (',round(pca_var[2,2],1),')'))
grid.arrange(g1,g2,ncol=2)	
	
#phenotype c('smr_brut','eb_smr','charge_varroa','nbr_varroa_100abeilles','vp','taux_inf','varroaMitoDepth','reop','reop_inf')
df_pca<-df[,c('smr_brut','eb_smr','charge_varroa','nbr_varroa_100abeilles','vp','taux_inf','varroaMitoDepth','reop','reop_inf')]
nb=estim_ncpPCA(df_pca,ncp.max=ncol(df))
res.comp=imputePCA(df_pca,ncp=nb$ncp)
pca.res=PCA(res.comp$completeObs,graph=T) 
pca_var<-pca.res$eig
pca<-pca.res$ind$coord
pca<-as.data.frame(pca)
pca$id<-df$num_ruche_bs
pca$k<-df$k
g1<-fviz_pca_var(pca.res,col.var="contrib")+theme_minimal()
g2<-ggplot(pca)+
	geom_point(aes(x=Dim.1,y=Dim.2,col=k),alpha=0.3)+
	xlab(paste0('PC1 (',round(pca_var[1,2],1),')'))+ylab(paste0('PC2 (',round(pca_var[2,2],1),')'))
grid.arrange(g1,g2,ncol=2)	
	
ggplot(df)+
	geom_point(aes(x=smr_brut,y=eb_smr,col=nbr_inf_une_fond))+
	xlim(0,1)+
	ylim(0,1)+
	theme_bw()+
	facet_wrap(~k)

ggplot(df)+
	geom_point(aes(x=vp,y=nbr_varroa_100abeilles))+
	theme_bw()+
	facet_wrap(~k,scales="free")

g1<-ggplot(df)+
	geom_histogram(aes(x=smr_brut))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g2<-ggplot(df)+
	geom_histogram(aes(x=eb_smr))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g3<-ggplot(df)+
	geom_histogram(aes(x=charge_varroa))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g4<-ggplot(df)+
	geom_histogram(aes(x=nbr_varroa_100abeilles))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g5<-ggplot(df)+
	geom_histogram(aes(x=taux_inf))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g6<-ggplot(df)+
	geom_histogram(aes(x=varroaMitoDepth))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g7<-ggplot(df)+
	geom_histogram(aes(x=reop))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g8<-ggplot(df)+
	geom_histogram(aes(x=reop_inf))+
	theme_bw()+
	facet_grid(.~k,scales="free")
grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8,ncol=1)

g1<-ggplot(df)+
	geom_histogram(aes(x=charge_varroa))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g2<-ggplot(df)+
	geom_histogram(aes(x=log10(charge_varroa)))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g3<-ggplot(df)+
	geom_histogram(aes(x=sqrt(charge_varroa)))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g4<-ggplot(df)+
	geom_histogram(aes(x=(charge_varroa)*(1/3)))+
	theme_bw()+
	facet_grid(.~k,scales="free")
grid.arrange(g1,g2,g3,g4,ncol=2)

g1<-ggplot(df)+
	geom_histogram(aes(x=nbr_varroa_100abeilles))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g2<-ggplot(df)+
	geom_histogram(aes(x=log10(nbr_varroa_100abeilles)))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g3<-ggplot(df)+
	geom_histogram(aes(x=sqrt(nbr_varroa_100abeilles)))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g4<-ggplot(df)+
	geom_histogram(aes(x=(nbr_varroa_100abeilles)*(1/3)))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g5<-ggplot(df)+
	geom_histogram(aes(x=vp))+
	theme_bw()+
	facet_grid(.~k,scales="free")
grid.arrange(g1,g2,g3,g4,g5,ncol=2)

g1<-ggplot(df)+
	geom_histogram(aes(x=taux_inf))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g2<-ggplot(df)+
	geom_histogram(aes(x=log10(taux_inf)))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g3<-ggplot(df)+
	geom_histogram(aes(x=sqrt(taux_inf)))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g4<-ggplot(df)+
	geom_histogram(aes(x=(taux_inf)*(1/3)))+
	theme_bw()+
	facet_grid(.~k,scales="free")
grid.arrange(g1,g2,g3,g4,ncol=2)

g1<-ggplot(df)+
	geom_histogram(aes(x=varroaMitoDepth))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g2<-ggplot(df)+
	geom_histogram(aes(x=log10(varroaMitoDepth)))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g3<-ggplot(df)+
	geom_histogram(aes(x=sqrt(varroaMitoDepth)))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g4<-ggplot(df)+
	geom_histogram(aes(x=(varroaMitoDepth)*(1/3)))+
	theme_bw()+
	facet_grid(.~k,scales="free")
grid.arrange(g1,g2,g3,g4,ncol=2)

g1<-ggplot(df)+
	geom_histogram(aes(x=reop))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g2<-ggplot(df)+
	geom_histogram(aes(x=log10(reop)))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g3<-ggplot(df)+
	geom_histogram(aes(x=sqrt(reop)))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g4<-ggplot(df)+
	geom_histogram(aes(x=(reop)*(1/3)))+
	theme_bw()+
	facet_grid(.~k,scales="free")
grid.arrange(g1,g2,g3,g4,ncol=2)

m<-df%>%group_by(k)%>%summarise(m=mean(reop_inf,na.rm=T),s=sd(reop_inf,na.rm=T))
df<-merge(df,m,by='k',all=T)
df$reop_inf2<-abs(df$reop_inf-df$m)
g1<-ggplot(df)+
	geom_histogram(aes(x=reop_inf))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g2<-ggplot(df)+
	geom_histogram(aes(x=reop_inf2))+
	theme_bw()+
	facet_grid(.~k,scales="free")
grid.arrange(g1,g2,nrow=2)

g1<-ggplot(df)+
	geom_histogram(aes(x=smr_brut))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g2<-ggplot(df)+
	geom_histogram(aes(x=eb_smr))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g3<-ggplot(df)+
	geom_histogram(aes(x=log10(charge_varroa)))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g4<-ggplot(df)+
	geom_histogram(aes(x=sqrt(nbr_varroa_100abeilles)))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g5<-ggplot(df)+
	geom_histogram(aes(x=log10(taux_inf)))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g6<-ggplot(df)+
	geom_histogram(aes(x=log10(varroaMitoDepth)))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g7<-ggplot(df)+
	geom_histogram(aes(x=log10(reop)))+
	theme_bw()+
	facet_grid(.~k,scales="free")
g8<-ggplot(df)+
	geom_histogram(aes(x=reop_inf2))+
	theme_bw()+
	facet_grid(.~k,scales="free")
grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8,ncol=1)

dev.off()

