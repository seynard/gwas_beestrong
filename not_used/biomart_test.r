library(data.table)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(biomartr)

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






























########################################################################
p1<-c('pc1','eb_smr','logit_taux_reop_inf')
p2<-c('pc1_varroa_inf','mnr','reop_inf')
x=3
type='egs'
mantra<-fread(paste0('results/run_mantra_',p1[x],'_',type,'/mantra.out'),data.table=F)
colnames(mantra)<-c('rs','chr','ps','allele1','allele0','nb_study','log10BF','post_proba_hetero','n','direction')
l_mantra<-mantra$rs[mantra$log10BF>=5]
mel<-fread(paste0('results/summary_gwas_Mellifera_ncorse_',pheno1,'_freq_',type,'.txt.bz2'),data.table=F)
mel[,c("logl_H1","l_remle","betahat","sebetahat","NegativeProb","PositiveProb","PosteriorMean","PosteriorSD","significant")]<-NULL
l_mel<-mel$rs[-log10(mel$p_wald)>=5 | mel$svalue<=0.1]
lig<-fread(paste0('results/summary_gwas_Ligustica_Carnica_',pheno1,'_freq_',gwas,'.txt.bz2'),data.table=F)
lig[,c("logl_H1","l_remle","betahat","sebetahat","NegativeProb","PositiveProb","PosteriorMean","PosteriorSD","significant")]<-NULL
l_lig<-lig$rs[-log10(lig$p_wald)>=5 | lig$svalue<=0.1]
hyb<-fread(paste0('results/summary_gwas_hybrid_corse_',pheno1,'_freq_',gwas,'.txt.bz2'),data.table=F)
hyb[,c("logl_H1","l_remle","betahat","sebetahat","NegativeProb","PositiveProb","PosteriorMean","PosteriorSD","significant")]<-NULL
l_hyb<-hyb$rs[-log10(hyb$p_wald)>=5 | hyb$svalue<=0.1]
L<-c(l_mantra,l_mel,l_lig,l_hyb)
L<-unique(L)
mantra<-mantra[mantra$rs%in%L,]
mel<-mel[mel$rs%in%L,]
lig<-lig[lig$rs%in%L,]
hyb<-hyb[hyb$rs%in%L,]
df<-Reduce(function(x,y) merge(x,y,by=c("rs","chr","ps","allele1","allele0"),all=T),list(mantra,mel,lig,hyb))
df[,c('grm.x','gwas.x','grm.y','gwas.y')]<-NULL
colnames(df)<-gsub('.x','_mel',colnames(df))
colnames(df)<-gsub('.y','_lig',colnames(df))
colnames(df)[37:(ncol(df)-2)]<-paste0(colnames(df)[37:(ncol(df)-2)],'_hyb')
df$significance<-NA
for(i in 1:nrow(df)){
significance<-vector()
if(!is.na(df$log10BF[i]) & df$log10BF[i]>=5){significance<-paste0(significance,'mantra/')}
if(!is.na(df$svalue_mel[i]) & df$svalue_mel[i]<=0.1){significance<-paste0(significance,'mash_mel/')}
if (!is.na(df$svalue_lig[i]) & df$svalue_lig[i]<=0.1){significance<-paste0(significance,'mash_lig/')}
if (!is.na(df$svalue_hyb[i]) & df$svalue_hyb[i]<=0.1){significance<-paste0(significance,'mash_hyb/')}
if (!is.na(df$p_wald_mel[i]) & -log10(df$p_wald_mel[i])>=5){significance<-paste0(significance,'gemma_mel/')}
if (!is.na(df$p_wald_lig[i]) & -log10(df$p_wald_lig[i])>=5){significance<-paste0(significance,'gemma_lig/')}
if (!is.na(df$p_wald_hyb[i]) & -log10(df$p_wald_hyb[i])>=5){significance<-paste0(significance,'gemma_hyb/')}
df$significance[i]<-significance
}

g1<-ggplot(df)+geom_point(aes(x=log10BF,y=svalue_mel,col=significance))+theme_bw()+xlim(0,7)+ylim(0,1)+geom_vline(aes(xintercept=5),lty=3,col='red')+geom_hline(aes(yintercept=0.1),lty=3,col='red')
g2<-ggplot(df)+geom_point(aes(x=log10BF,y=svalue_lig,col=significance))+theme_bw()+xlim(0,7)+ylim(0,1)+geom_vline(aes(xintercept=5),lty=3,col='red')+geom_hline(aes(yintercept=0.1),lty=3,col='red')
g3<-ggplot(df)+geom_point(aes(x=log10BF,y=svalue_hyb,col=significance))+theme_bw()+xlim(0,7)+ylim(0,1)+geom_vline(aes(xintercept=5),lty=3,col='red')+geom_hline(aes(yintercept=0.1),lty=3,col='red')
g4<-ggplot(df)+geom_point(aes(x=log10BF,y=-log10(p_wald_mel),col=significance))+theme_bw()+xlim(0,7)+ylim(0,7)+geom_vline(aes(xintercept=5),lty=3,col='red')+geom_hline(aes(yintercept=5),lty=3,col='red')
g5<-ggplot(df)+geom_point(aes(x=log10BF,y=-log10(p_wald_lig),col=significance))+theme_bw()+xlim(0,7)+ylim(0,7)+geom_vline(aes(xintercept=5),lty=3,col='red')+geom_hline(aes(yintercept=5),lty=3,col='red')
g6<-ggplot(df)+geom_point(aes(x=log10BF,y=-log10(p_wald_hyb),col=significance))+theme_bw()+xlim(0,7)+ylim(0,7)+geom_vline(aes(xintercept=5),lty=3,col='red')+geom_hline(aes(yintercept=5),lty=3,col='red')
grid.arrange(grobs=list(g1,g2,g3,g4,g5,g6),nrow=2)


























########################################################################
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











































































#mart<-useMart('metazoa_mart',host='metazoa.ensembl.org') 
#mart<-useDataset('amellifera_eg_gene',mart)
#annot<-getBM(mart=mart,attributes=c('chromosome_name','start_position','end_position','go_id','name_1006','definition_1006','beebase','ensembl_metazoa','ensembl_gene_id','ensembl_transcript_id'),uniqueRows=TRUE,values=locset)
########################################################################

sign_sub
Var1 			Freq
LOC100576198    1
LOC100576405    1
LOC100577472    1
LOC100578400    1
LOC113219185    1
LOC409774    	2	proteoglycan Cow
LOC409919    	1
LOC410480    	1
Saelao	LOC410758    	1
LOC410796    	1
LOC412865    	2
LOC412926    	1
LOC413466    	2
LOC551649    	1
LOC551889    	1
Saelao	LOC724238    	1 development eyes, antenna
LOC726948    	1

https://metazoa.ensembl.org/index.html

2 markers in Saelao et al. (selection signal for varroa resistant honeybee)
multiple SNP in region linked to Dsx (differentiation between queen and worker)
