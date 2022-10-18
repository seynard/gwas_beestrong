library(data.table)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(biomartr)
library(mashr)
library(ashr)
library(flashr)

annot<-fread('data/proteins_48_403979.csv')
colnames(annot)[1]<-'LG'
colnames(annot)[8]<-'prot'
colnames(annot)[10]<-'prot_name'
annot[,7]<-NULL
annot$chr[annot$LG=='linkage group LG1']<-1
annot$chr[annot$LG=='linkage group LG2']<-2
annot$chr[annot$LG=='linkage group LG3']<-3
annot$chr[annot$LG=='linkage group LG4']<-4
annot$chr[annot$LG=='linkage group LG5']<-5
annot$chr[annot$LG=='linkage group LG6']<-6
annot$chr[annot$LG=='linkage group LG7']<-7
annot$chr[annot$LG=='linkage group LG8']<-8
annot$chr[annot$LG=='linkage group LG9']<-9
annot$chr[annot$LG=='linkage group LG10']<-10
annot$chr[annot$LG=='linkage group LG11']<-11
annot$chr[annot$LG=='linkage group LG12']<-12
annot$chr[annot$LG=='linkage group LG13']<-13
annot$chr[annot$LG=='linkage group LG14']<-14
annot$chr[annot$LG=='linkage group LG15']<-15
annot$chr[annot$LG=='linkage group LG16']<-16
annot<-unique(annot)
chrsize<-annot%>%group_by(chr)%>%summarise(m=max(Stop))
chrsize<-subset(chrsize,!is.na(chrsize$chr))

### run mash + combine mantra
pheno<-'logit_taux_reop_inf'
type<-'egs'
mel<-fread(paste0('results/summary_gwas_Mellifera_ncorse_',pheno,'_freq_',type,'.txt.bz2'),data.table=F)
lig<-fread(paste0('results/summary_gwas_Ligustica_Carnica_',pheno,'_freq_',type,'.txt.bz2'),data.table=F)
hyb<-fread(paste0('results/summary_gwas_hybrid_corse_',pheno,'_freq_',type,'.txt.bz2'),data.table=F)

RS<-intersect(intersect(mel$rs,lig$rs),hyb$rs)
rand<-sample(RS,5000,replace=F)
mel_rand<-mel[mel$rs%in%rand,]
b_mel_rand<-mel_rand[,c('rs','beta')]
colnames(b_mel_rand)[2]<-'mel'
s_mel_rand<-mel_rand[,c('rs','se')]
colnames(s_mel_rand)[2]<-'mel'
lig_rand<-lig[lig$rs%in%rand,]
b_lig_rand<-lig_rand[,c('rs','beta')]
colnames(b_lig_rand)[2]<-'lig'
s_lig_rand<-lig_rand[,c('rs','se')]
colnames(s_lig_rand)[2]<-'lig'
hyb_rand<-hyb[hyb$rs%in%rand,]
b_hyb_rand<-hyb_rand[,c('rs','beta')]
colnames(b_hyb_rand)[2]<-'hyb'
s_hyb_rand<-hyb_rand[,c('rs','se')]
colnames(s_hyb_rand)[2]<-'hyb'
brand<-Reduce(function(x,y) merge(x,y,by='rs',all=T),list(b_lig_rand,b_mel_rand,b_hyb_rand))
srand<-Reduce(function(x,y) merge(x,y,by='rs',all=T),list(s_lig_rand,s_mel_rand,s_hyb_rand))
if(all(brand$rs==srand$rs)){
	print('order ok')
	hat_rs_rand<-brand$rs
	brand$rs<-NULL
	srand$rs<-NULL
}else{
	print('reorder')
	brand<-brand[order(brand$rs),]
	srand<-srand[order(srand$rs),]
	hat_rs_rand<-brand$rs
	brand$rs<-NULL
	srand$rs<-NULL
	}
brand<-as.matrix(brand)
srand<-as.matrix(srand)
dtmp<-mash_set_data(brand,srand)
Vhat<-estimate_null_correlation_simple(dtmp)
rm(dtmp)
drand<-mash_set_data(brand,srand,V=Vhat)

mel<-subset(mel,mel$lfsr<0.5)
b_mel<-mel[,c('rs','beta')]
colnames(b_mel)[2]<-'mel'
s_mel<-mel[,c('rs','se')]
colnames(s_mel)[2]<-'mel'
lig<-subset(lig,lig$lfsr<0.5)
b_lig<-lig[,c('rs','beta')]
colnames(b_lig)[2]<-'lig'
s_lig<-lig[,c('rs','se')]
colnames(s_lig)[2]<-'lig'
hyb<-subset(hyb,hyb$lfsr<0.5)
b_hyb<-hyb[,c('rs','beta')]
colnames(b_hyb)[2]<-'hyb'
s_hyb<-hyb[,c('rs','se')]
colnames(s_hyb)[2]<-'hyb'
bhat<-Reduce(function(x,y) merge(x,y,by='rs',all=T),list(b_lig,b_mel,b_hyb))
shat<-Reduce(function(x,y) merge(x,y,by='rs',all=T),list(s_lig,s_mel,s_hyb))
if(all(bhat$rs==shat$rs)){
	print('order ok')
	hat_rs<-bhat$rs
	bhat$rs<-NULL
	shat$rs<-NULL
}else{
	print('reorder')
	bhat<-bhat[order(bhat$rs),]
	shat<-shat[order(shat$rs),]
	hat_rs<-bhat$rs
	bhat$rs<-NULL
	shat$rs<-NULL
	}
bhat<-as.matrix(bhat)
shat<-as.matrix(shat)
dstrong<-mash_set_data(bhat,shat,V=Vhat)

Upca<-cov_pca(dstrong,3)
Ued<-cov_ed(dstrong,Upca)
Uc<-cov_canonical(drand)
Umantra<-data.frame(c(0,0.3,0.2),c(0.3,0,0.1),c(0.2,0.1,0))
Umantra<-1-Umantra
Umantra<-as.matrix(Umantra)
Umantra<-list(Umantra)
#Ueffect<-matrix(nrow=ncol(bhat),ncol=ncol(bhat))
#for(i in 1:ncol(bhat)){
#	for(j in 1:ncol(bhat)){
#		Ueffect[i,j]<-cor.test(bhat[,i],bhat[,j],na.omit=T)$estimate
#		}
#	}
#Ueffect<-list(Ueffect)	
#m<-mash(drand,Ulist=c(Ued,Uc,Umantra,Ueffect),outputlevel=1)
m<-mash(drand,Ulist=c(Ued,Uc,Umantra),outputlevel=1)
m2<-mash(dstrong,g=get_fitted_g(m),fixg=TRUE)
get_lfsr(m2)
get_significant_results(m2,thresh=0.1)
get_estimated_pi(m2)
get_log10bf(m2)
get_n_significant_conditions(m2)
get_pairwise_sharing(m2)
logBF<-data.frame('rs'=hat_rs,'log10_bf'=get_log10bf(m2),'n_sign'=get_n_significant_conditions(m2))
a<-colsplit(logBF$rs,':',c('chr','ps'))
logBF$chr<-a$chr
logBF$ps<-a$ps

ggplot(logBF)+
	geom_point(aes(x=ps,y=log10_bf,col=as.factor(n_sign)))+
	facet_grid(.~chr,scales='free',space='free')+
	theme_bw()
	
mantra<-fread(paste0('results/run_mantra_',pheno,'_',type,'/mantra.out'),data.table=F)
colnames(mantra)<-c('rs','chr','ps','allele1','allele0','nb_study','log10BF','post_proba_hetero','n','direction')

DF<-merge(mantra,logBF,by=c('rs','chr','ps'),all=T)

ggplot(DF)+
	geom_point(aes(x=log10BF,y=log10_bf,col=as.factor(nb_study)),alpha=0.3)+
	geom_hline(aes(yintercept=2),col='red',lty=3)+
	geom_vline(aes(xintercept=5),col='red',lty=3)+
	theme_bw()



### run mash + combine mantra

sign<-fread('summary_sign_prot.txt',data.table=F)
sign<-subset(sign,sign$pheno%in%c('pc1_varroa_inf','mnr','reop_inf'))

sign_sub<-subset(sign,sign$dist<5000 & sign$type=='egs' & sign$nb_study==3)
ggplot()+
	geom_point(aes(x=ps,y=log10BF,col=pheno),data=sign_sub)+
	geom_vline(aes(xintercept=0),col='red',lty=3,data=chrsize)+
	geom_vline(aes(xintercept=m),col='red',lty=3,data=chrsize)+
	facet_grid(chr~.,space='free')+
	ylim(5,8)+
	theme_bw()
x<-data.frame(table(sign_sub$Locus))

p1<-c('pc1','eb_smr','logit_taux_reop_inf')
p2<-c('pc1_varroa_inf','mnr','reop_inf')
for(x in 1:length(p1)){
print(x)
pheno1=p1[x]
pheno2=p2[x]
mel<-fread(paste0('results/summary_gwas_Mellifera_ncorse_',pheno1,'_freq_egs.txt.bz2'),data.table=F)
mel$subsp<-'mel'
mel$sign<-'no'
mel$sign[mel$rs%in%sign$ID[sign$pheno==pheno2]]<-'yes'
lig<-fread(paste0('results/summary_gwas_Ligustica_Carnica_',pheno1,'_freq_egs.txt.bz2'),data.table=F)
lig$subsp<-'lig'
lig$sign<-'no'
lig$sign[lig$rs%in%sign$ID[sign$pheno==pheno2]]<-'yes'
hyb<-fread(paste0('results/summary_gwas_hybrid_corse_',pheno1,'_freq_egs.txt.bz2'),data.table=F)
hyb$subsp<-'hyb'
hyb$sign<-'no'
hyb$sign[hyb$rs%in%sign$ID[sign$pheno==pheno2]]<-'yes'
sign_x1<-subset(sign,sign$pheno==pheno2 &sign$type=='egs')
sign_x1<-sign_x1[order(sign_x1$chr,sign_x1$ps),]
p<-list()
for(i in 1:nrow(sign_x1)){
if(sign_x1$gene_location[i]=='in'){
mx<-subset(mel,mel$chr==sign_x1$chr[i] & mel$ps>=sign_x1$Start[i] & mel$ps<=sign_x1$Stop[i]) 
lx<-subset(lig,lig$chr==sign_x1$chr[i] & lig$ps>=sign_x1$Start[i] & lig$ps<=sign_x1$Stop[i]) 
hx<-subset(hyb,hyb$chr==sign_x1$chr[i] & hyb$ps>=sign_x1$Start[i] & hyb$ps<=sign_x1$Stop[i]) 
x<-do.call(rbind,list(mx,lx,hx))
x$subsp[x$subsp=='mel']<-'Mellifera'
x$subsp[x$subsp=='lig']<-'Ligustica_Carnica'
x$subsp[x$subsp=='hyb']<-'Hybrid'
x$subsp<-factor(x$subsp,levels=c('Ligustica_Carnica','Mellifera','Hybrid'))
}else if (sign_x1$gene_location[i]=='before'){
mx<-subset(mel,mel$chr==sign_x1$chr[i] & mel$ps>=(sign_x1$Start[i]-sign_x1$dist[i]-10) & mel$ps<=sign_x1$Stop[i]) 
lx<-subset(lig,lig$chr==sign_x1$chr[i] & lig$ps>=(sign_x1$Start[i]-sign_x1$dist[i]-10) & lig$ps<=sign_x1$Stop[i]) 
hx<-subset(hyb,hyb$chr==sign_x1$chr[i] & hyb$ps>=(sign_x1$Start[i]-sign_x1$dist[i]-10) & hyb$ps<=sign_x1$Stop[i]) 
x<-do.call(rbind,list(mx,lx,hx))
x$subsp[x$subsp=='mel']<-'Mellifera'
x$subsp[x$subsp=='lig']<-'Ligustica_Carnica'
x$subsp[x$subsp=='hyb']<-'Hybrid'
x$subsp<-factor(x$subsp,levels=c('Ligustica_Carnica','Mellifera','Hybrid'))
}else if (sign_x1$gene_location[i]=='after'){
mx<-subset(mel,mel$chr==sign_x1$chr[i] & mel$ps>=sign_x1$Start[i] & mel$ps<=(sign_x1$Stop[i]+sign_x1$dist[i]+10)) 
lx<-subset(lig,lig$chr==sign_x1$chr[i] & lig$ps>=sign_x1$Start[i] & lig$ps<=(sign_x1$Stop[i]+sign_x1$dist[i]+10)) 
hx<-subset(hyb,hyb$chr==sign_x1$chr[i] & hyb$ps>=sign_x1$Start[i] & hyb$ps<=(sign_x1$Stop[i]+sign_x1$dist[i]+10)) 
x<-do.call(rbind,list(mx,lx,hx))
x$subsp[x$subsp=='mel']<-'Mellifera'
x$subsp[x$subsp=='lig']<-'Ligustica_Carnica'
x$subsp[x$subsp=='hyb']<-'Hybrid'
x$subsp<-factor(x$subsp,levels=c('Ligustica_Carnica','Mellifera','Hybrid'))
}
p[[i]]<-ggplot()+
	geom_point(aes(x=ps,y=-log10(p_wald)),alpha=0.3,data=x[x$sign=='no',])+
	geom_point(aes(x=ps,y=-log10(p_wald),col=beta),size=3,data=x[x$sign=='yes',])+
	geom_segment(aes(x=Start,xend=Stop,y=0,yend=0),col='red',size=2,data=sign_x1[i,])+
	scale_colour_gradient2(low="dark blue",high="dark green",mid="orange",midpoint=0,na.value=NA)+
	facet_grid(subsp~.)+
	theme_bw()+
	ggtitle(paste0('SNP : ',sign_x1$ID[i],' closest Locus : ',sign_x1$Locus[i]))
	}
n=ceiling(sqrt(nrow(sign_x1)))
png(paste0('plot_',pheno2,'_sign_egs.png'),width=1500*n,height=1500*n,res=300)
grid.arrange(grobs=p,nrow=n,ncol=n)
dev.off()
sign_x2<-subset(sign,sign$pheno==pheno2 &sign$type=='freq')
sign_x2<-sign_x2[order(sign_x2$chr,sign_x2$ps),]
q<-list()
for(i in 1:nrow(sign_x2)){
if(sign_x2$gene_location[i]=='in'){
mx<-subset(mel,mel$chr==sign_x2$chr[i] & mel$ps>=sign_x2$Start[i] & mel$ps<=sign_x2$Stop[i]) 
lx<-subset(lig,lig$chr==sign_x2$chr[i] & lig$ps>=sign_x2$Start[i] & lig$ps<=sign_x2$Stop[i]) 
hx<-subset(hyb,hyb$chr==sign_x2$chr[i] & hyb$ps>=sign_x2$Start[i] & hyb$ps<=sign_x2$Stop[i]) 
x<-do.call(rbind,list(mx,lx,hx))
x$subsp[x$subsp=='mel']<-'Mellifera'
x$subsp[x$subsp=='lig']<-'Ligustica_Carnica'
x$subsp[x$subsp=='hyb']<-'Hybrid'
x$subsp<-factor(x$subsp,levels=c('Ligustica_Carnica','Mellifera','Hybrid'))
}else if (sign_x2$gene_location[i]=='before'){
mx<-subset(mel,mel$chr==sign_x2$chr[i] & mel$ps>=(sign_x2$Start[i]-sign_x2$dist[i]-10) & mel$ps<=sign_x2$Stop[i]) 
lx<-subset(lig,lig$chr==sign_x2$chr[i] & lig$ps>=(sign_x2$Start[i]-sign_x2$dist[i]-10) & lig$ps<=sign_x2$Stop[i]) 
hx<-subset(hyb,hyb$chr==sign_x2$chr[i] & hyb$ps>=(sign_x2$Start[i]-sign_x2$dist[i]-10) & hyb$ps<=sign_x2$Stop[i]) 
x<-do.call(rbind,list(mx,lx,hx))
x$subsp[x$subsp=='mel']<-'Mellifera'
x$subsp[x$subsp=='lig']<-'Ligustica_Carnica'
x$subsp[x$subsp=='hyb']<-'Hybrid'
x$subsp<-factor(x$subsp,levels=c('Ligustica_Carnica','Mellifera','Hybrid'))
}else if (sign_x2$gene_location[i]=='after'){
mx<-subset(mel,mel$chr==sign_x2$chr[i] & mel$ps>=sign_x2$Start[i] & mel$ps<=(sign_x2$Stop[i]+sign_x2$dist[i]+10)) 
lx<-subset(lig,lig$chr==sign_x2$chr[i] & lig$ps>=sign_x2$Start[i] & lig$ps<=(sign_x2$Stop[i]+sign_x2$dist[i]+10)) 
hx<-subset(hyb,hyb$chr==sign_x2$chr[i] & hyb$ps>=sign_x2$Start[i] & hyb$ps<=(sign_x2$Stop[i]+sign_x2$dist[i]+10)) 
x<-do.call(rbind,list(mx,lx,hx))
x$subsp[x$subsp=='mel']<-'Mellifera'
x$subsp[x$subsp=='lig']<-'Ligustica_Carnica'
x$subsp[x$subsp=='hyb']<-'Hybrid'
x$subsp<-factor(x$subsp,levels=c('Ligustica_Carnica','Mellifera','Hybrid'))
}
q[[i]]<-ggplot()+
	geom_point(aes(x=ps,y=-log10(p_wald)),alpha=0.3,data=x[x$sign=='no',])+
	geom_point(aes(x=ps,y=-log10(p_wald),col=beta),size=3,data=x[x$sign=='yes',])+
	geom_segment(aes(x=Start,xend=Stop,y=0,yend=0),col='red',size=2,data=sign_x2[i,])+
	scale_colour_gradient2(low="dark blue",high="dark green",mid="orange",midpoint=0,na.value=NA)+
	facet_grid(subsp~.)+
	theme_bw()+
	ggtitle(paste0('SNP : ',sign_x2$ID[i],' closest Locus : ',sign_x2$Locus[i]))
	}
n=ceiling(sqrt(nrow(sign_x2)))
png(paste0('plot_',pheno2,'_sign_freq.png'),width=1500*n,height=1500*n,res=300)
grid.arrange(grobs=q,nrow=n,ncol=n)
dev.off()
}

locset<-unique(sign$Locus)
go_apis<-biomartr::getGO(organism="apis mellifera",genes=locset,filters="ensembl_gene_id")
colnames(go_apis)<-c('name_loc','description','go_nb')
for(i in 1:length(unique(go_apis$name_loc))){
	x<-go_apis[go_apis$name_loc==unique(go_apis$name_loc)[i],]
	print(x)
	}

p1<-c('pc1','eb_smr','logit_taux_reop_inf')
p2<-c('pc1_varroa_inf','mnr','reop_inf')
typ<-c('egs','freq')
for(x in 1:length(p1)){
for(y in 1:length(type)){
mantra<-fread(paste0('results/run_mantra_',p1[x],'_',type[y],'/mantra.out'),data.table=F)
colnames(mantra)<-c('rs','chr','ps','allele1','allele0','nb_study','log10BF','post_proba_hetero','n','direction')
mantra<-mantra[,c('rs','chr','ps','log10BF')]
mel<-fread(paste0('results/summary_gwas_Mellifera_ncorse_',p1[x],'_freq_',type[y],'.txt.bz2'),data.table=F)
mel<-mel[,c('rs','chr','ps','p_wald','svalue')]
lig<-fread(paste0('results/summary_gwas_Ligustica_Carnica_',p1[x],'_freq_',type[y],'.txt.bz2'),data.table=F)
lig<-lig[,c('rs','chr','ps','p_wald','svalue')]
hyb<-fread(paste0('results/summary_gwas_hybrid_corse_',p1[x],'_freq_',type[y],'.txt.bz2'),data.table=F)
hyb<-hyb[,c('rs','chr','ps','p_wald','svalue')]
df<-Reduce(function(x,y) merge(x,y,by=c("rs","chr","ps"),all=T),list(mantra,mel,lig,hyb))
colnames(df)<-gsub('.x','_mel',colnames(df))
colnames(df)<-gsub('.y','_lig',colnames(df))
colnames(df)[9:ncol(df)]<-paste0(colnames(df)[9:ncol(df)],'_hyb')
g1<-ggplot()+
	geom_point(aes(x=ps,y=-log10(p_wald_mel)),col='grey34',alpha=0.3,data=df)+
	geom_point(aes(x=ps,y=-log10(p_wald_lig)),col='goldenrod2',alpha=0.3,data=df)+
	geom_point(aes(x=ps,y=-log10(p_wald_hyb)),col='cyan',alpha=0.3,data=df)+
	geom_point(aes(x=ps,y=log10BF),col='red',data=df[df$log10BF>=5,])+
	geom_hline(aes(yintercept=5),lty=3,col='red',data=df)+
	theme_bw()+facet_grid(.~chr,scales='free',space='free')+
	ylab('-log10(p_value)')+xlab('position in bp')+
	scale_x_continuous(breaks=seq(0,max(df$ps),5000000))+
	theme(axis.text.x=element_text(angle=45))
g2<-ggplot()+
	geom_point(aes(x=ps,y=1-svalue_mel),col='grey34',alpha=0.3,data=df)+
	geom_point(aes(x=ps,y=1-svalue_lig),col='goldenrod2',alpha=0.3,data=df)+
	geom_point(aes(x=ps,y=1-svalue_hyb),col='cyan',alpha=0.3,data=df)+
	geom_point(aes(x=ps,y=log10BF/6),col='red',data=df[df$log10BF>=5,])+
	geom_hline(aes(yintercept=1-0.1),lty=3,col='blue',data=df)+
	geom_hline(aes(yintercept=5/6),lty=3,col='red',data=df)+
	theme_bw()+facet_grid(.~chr,scales='free',space='free')+
	scale_y_continuous(name='1-svalues',sec.axis=sec_axis(trans=~./6,name='log10BF/6'))+
	xlab('position in bp')+
	scale_x_continuous(breaks=seq(0,max(df$ps),5000000))+
	theme(axis.text.x=element_text(angle=45),
	axis.line.y.right=element_line(color="red"),axis.ticks.y.right=element_line(color="red"),axis.text.y.right=element_text(color="red"),axis.title.y.right=element_text(color="red"),
	axis.line.y.left=element_line(color="blue"),axis.ticks.y.left=element_line(color="blue"),axis.text.y.left=element_text(color="blue"),axis.title.y.left=element_text(color="blue"))
png(paste0('plot_gwas_',p1[x],'_',type[y],'.png'),width=3500,height=1500,res=300)
grid.arrange(g1,g2,ncol=1)
dev.off()
}}




