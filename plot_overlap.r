#module load system/R-4.1.1_gcc-9.3.0
#R
library(data.table)
library(ggplot2)
library(cowplot)
library(tidyverse)
args<-commandArgs(TRUE)
version<-as.numeric(args[1])
if(version==1){pop<-c('Mellifera','Ligustica_Carnica','hybrid')}else if (version==2){pop<-c('Mellifera','Ligustica_nus','hybrid')}

reop_mantra<-fread(paste0('results/run_mantra_logit_taux_reop_inf_egs_',version,'/mantra.out'),data.table=F)
colnames(reop_mantra)<-c('rs','chr','ps','all_effect','all_0','n_study','log10_bf_mantra','proba_het','n','dir')
reop_mash<-fread(paste0('results/result_mash_logit_taux_reop_inf_egs_',version,'_simple_corr.txt'),data.table=F)
colnames(reop_mash)<-c('rs','log10_bf_mash','n_sign','loglik','chr','ps')
df_recap<-merge(reop_mantra,reop_mash,by=c('chr','ps','rs'),all=T)
df_recap$pheno<-'recap_inf'

mnr_mantra<-fread(paste0('results/run_mantra_eb_smr_egs_',version,'/mantra.out'),data.table=F)
colnames(mnr_mantra)<-c('rs','chr','ps','all_effect','all_0','n_study','log10_bf_mantra','proba_het','n','dir')
mnr_mash<-fread(paste0('results/result_mash_eb_smr_egs_',version,'_simple_corr.txt'),data.table=F)
colnames(mnr_mash)<-c('rs','log10_bf_mash','n_sign','loglik','chr','ps')
df_mnr<-merge(mnr_mantra,mnr_mash,by=c('chr','ps','rs'),all=T)
df_mnr$pheno<-'mnr'

pc1_mantra<-fread(paste0('results/run_mantra_pc1_egs_',version,'/mantra.out'),data.table=F)
colnames(pc1_mantra)<-c('rs','chr','ps','all_effect','all_0','n_study','log10_bf_mantra','proba_het','n','dir')
pc1_mash<-fread(paste0('results/result_mash_pc1_egs_',version,'_simple_corr.txt'),data.table=F)
colnames(pc1_mash)<-c('rs','log10_bf_mash','n_sign','loglik','chr','ps')
df_pc1<-merge(pc1_mantra,pc1_mash,by=c('chr','ps','rs'),all=T)
df_pc1$pheno<-'pc1_varroa_inf'

df<-do.call(rbind,list(df_recap,df_mnr,df_pc1))
df$sign[df$n_sign==0]<-'not_sign'
df$sign[df$n_sign>0]<-'sign'
df$pheno<-factor(df$pheno,levels=c('recap_inf','mnr','pc1_varroa_inf'))
df$log10_bf_mantra[df$log10_bf_mantra<0]<-0
df$log10_bf_mash[df$log10_bf_mash<0]<-0

df_sub<-df[df$log10_bf_mantra>5 | df$log10_bf_mash>1 | df$sign=='sign',c('chr','ps','pheno')]
df_sub<-df_sub[order(df_sub$chr,df_sub$ps),]
df_sub$ps_start<-df_sub$ps-500000
df_sub$ps_stop<-df_sub$ps+500000
df_sub$pheno<-as.character(df_sub$pheno)
df_sub$p[df_sub$pheno=='recap_inf']<-'p1'
df_sub$p[df_sub$pheno=='mnr']<-'p2'
df_sub$p[df_sub$pheno=='pc1_varroa_inf']<-'p3'
#df_sub$pheno<-NULL
df_sub$pheno<-factor(df_sub$pheno,levels=c('recap_inf','mnr','pc1_varroa_inf'))
df_sub$ps_start[df_sub$ps_start<0]<-0

ymin=min(c(df$log10_bf_mantra,-df$log10_bf_mash))-0.2
ymax=max(c(df$log10_bf_mantra,-df$log10_bf_mash))+0.2
chrsize<-df%>%group_by(chr)%>%summarize(m=max(ps))

#png(paste0('plot_3pheno_overlap_',version,'.png'),width=7000,height=3000,res=300)
#ggplot()+	
#	geom_rect(aes(xmin=ps_start,xmax=ps_stop,ymin=ymin,ymax=ymax,col=pheno,fill=pheno),alpha=0.2,data=df_sub)+
#	scale_fill_manual(values=c('orange','green','purple'),labels=c('recap_inf','mnr','pc1_varroa_inf'),name='1Mb regions with \nsignificant markers')+
#	scale_color_manual(values=c('orange','green','purple'),labels=c('recap_inf','mnr','pc1_varroa_inf'),name='',guide='none')+
#	geom_point(aes(x=ps,y=log10_bf_mantra),col='black',alpha=0.4,data=df[,c('chr','ps','log10_bf_mantra','proba_het','pheno')])+
#	geom_point(aes(x=ps,y=-log10_bf_mash,shape=sign),col='black',alpha=0.4,data=df[,c('chr','ps','log10_bf_mash','sign','pheno')])+
#	scale_shape_manual(values=c(19,17),guide='none')+
#	geom_hline(yintercept=-1,col='red',linetype=3,size=1.5)+
#	geom_hline(yintercept=0,col='black',linetype=1)+
#	geom_hline(yintercept=5,col='red',linetype=3,size=1.5)+
#	facet_grid(pheno~chr,scales='free_x',space='free_x')+
#	scale_x_continuous(breaks=seq(0,chrsize$m[1],5000000))+
#	xlab('position in bp')+ylab('log10(BF): positive MANTRA, negative Mash')+
#	theme_bw()+	
#	theme(axis.text.x=element_text(size=15,angle=45,hjust=1),axis.text.y=element_text(size=15),axis.title.y=element_text(size=20),axis.title.x=element_text(size=20))+
#	theme(legend.title=element_text(size=20),legend.text=element_text(size=15))+
#	theme(strip.text.x=element_text(size=20),strip.text.y=element_text(size=20))
#dev.off()	

x_recap<-df[df$pheno=='recap_inf' & (df$log10_bf_mash>1 | df$log10_bf_mantra>5),]
x_mnr<-df[df$pheno=='mnr' & (df$log10_bf_mash>1 | df$log10_bf_mantra>5),]
x_pc1<-df[df$pheno=='pc1_varroa_inf' & (df$log10_bf_mash>1 | df$log10_bf_mantra>5),]

X_recap_mnr<-list()
for(i in 1:nrow(x_recap)){
	x<-x_recap[i,]
	y<-subset(x_mnr,x_mnr$chr==x$chr)
	y$diff<-abs(x$ps-y$ps)
	y<-subset(y,y$diff<1000000)	
	if(nrow(y)>0){X_recap_mnr[[i]]<-merge(x,y,by='chr')}
}
X_recap_mnr<-do.call(rbind,X_recap_mnr)
X_recap_mnr$pheno<-'recap_mnr'
colnames(X_recap_mnr)<-gsub('.x','.1',colnames(X_recap_mnr))
colnames(X_recap_mnr)<-gsub('.y','.2',colnames(X_recap_mnr))
X_recap_pc1<-list()
for(i in 1:nrow(x_recap)){
	x<-x_recap[i,]
	y<-subset(x_pc1,x_pc1$chr==x$chr)
	y$diff<-abs(x$ps-y$ps)
	y<-subset(y,y$diff<1000000)	
	if(nrow(y)>0){X_recap_pc1[[i]]<-merge(x,y,by='chr')}
}
X_recap_pc1<-do.call(rbind,X_recap_pc1)
X_recap_pc1$pheno<-'recap_pc1'
colnames(X_recap_pc1)<-gsub('.x','.1',colnames(X_recap_pc1))
colnames(X_recap_pc1)<-gsub('.y','.2',colnames(X_recap_pc1))
X_pc1_mnr<-list()
for(i in 1:nrow(x_pc1)){
	x<-x_pc1[i,]
	y<-subset(x_mnr,x_mnr$chr==x$chr)
	y$diff<-abs(x$ps-y$ps)
	y<-subset(y,y$diff<1000000)	
	if(nrow(y)>0){X_pc1_mnr[[i]]<-merge(x,y,by='chr')}
}
X_pc1_mnr<-do.call(rbind,X_pc1_mnr)
X_pc1_mnr$pheno<-'pc1_mnr'
colnames(X_pc1_mnr)<-gsub('.x','.1',colnames(X_pc1_mnr))
colnames(X_pc1_mnr)<-gsub('.y','.2',colnames(X_pc1_mnr))

X<-do.call(rbind,list(X_recap_mnr,X_recap_pc1,X_pc1_mnr))
X<-X[,c('chr','ps.1','rs.1','ps.2','rs.2','pheno','diff')]

annot<-fread('data/proteins_48_403979.csv',data.table=F)
annot<-unique(annot[,c(1,2,3,4,6)])
annot<-subset(annot,annot[,1]!='Un')

annot$chr[annot[,1]=='linkage group LG1']<-1
annot$chr[annot[,1]=='linkage group LG2']<-2
annot$chr[annot[,1]=='linkage group LG3']<-3
annot$chr[annot[,1]=='linkage group LG4']<-4
annot$chr[annot[,1]=='linkage group LG5']<-5
annot$chr[annot[,1]=='linkage group LG6']<-6
annot$chr[annot[,1]=='linkage group LG7']<-7
annot$chr[annot[,1]=='linkage group LG8']<-8
annot$chr[annot[,1]=='linkage group LG9']<-9
annot$chr[annot[,1]=='linkage group LG10']<-10
annot$chr[annot[,1]=='linkage group LG11']<-11
annot$chr[annot[,1]=='linkage group LG12']<-12
annot$chr[annot[,1]=='linkage group LG13']<-13
annot$chr[annot[,1]=='linkage group LG14']<-14
annot$chr[annot[,1]=='linkage group LG15']<-15
annot$chr[annot[,1]=='linkage group LG16']<-16
annot<-unique(annot)
annot$start[annot$Strand=='+']<-annot$Start[annot$Strand=='+']
annot$stop[annot$Strand=='+']<-annot$Stop[annot$Strand=='+']
annot$start[annot$Strand=='-']<-annot$Stop[annot$Strand=='-']
annot$stop[annot$Strand=='-']<-annot$Start[annot$Strand=='-']

Ann<-list()
for(i in 1:nrow(X)){
	x<-X[i,]
	ps_min<-min(x$ps.1,x$ps.2)
	ps_max<-max(x$ps.1,x$ps.2)
	size<-round((1000000-(ps_max-ps_min))/2)
	ann<-unique(annot[annot$chr==x$chr & annot$Start>=ps_min-size & annot$Stop<=ps_max+size,c('Start','Stop','Strand','Locus')])
	if(nrow(ann)>0){
	Ann[[i]]<-data.frame('pheno'=x$pheno,'chr'=x$chr,'pos'=paste0(ps_min,'_',ps_max),'ps_start'=ps_min-size,'ps_stop'=ps_max+size,'locus'=ann$Locus,'Start'=ann$Start,'Stop'=ann$Stop,'Strand'=ann$Strand)}
}
Ann<-do.call(rbind,Ann)
table(Ann$pheno)
write.table(Ann,paste0('overlap_regions_',version,'.txt'),col.names=T,row.names=F,quote=F)

for(i in 1:length(unique(Ann$pheno))){
	phe=unique(Ann$pheno)[i]
	print(phe)
	a1<-Ann[Ann$pheno==phe,]
	for(j in 1:length(unique(a1$chr))){
		c=unique(a1$chr)[j]
		print(c)
		a2<-a1[a1$chr==c,]
		for(k in 1:length(unique(a2$pos))){
			p=unique(a2$pos)[k]
			print(p)
			a3<-a2[a2$pos==p,]
			s=df[df$chr==c & df$ps>=min(a3$ps_start)  & df$ps<=max(a3$ps_stop)  & (df$log10_bf_mantra>5 | df$log10_bf_mash>1),]
			png(paste0('plot_phe_',phe,'_chrom_',c,'_pos_',p,'_',version,'.png'),width=2000,height=3500,res=300)
			g1<-ggplot()+
				geom_point(aes(x=ps,y=log10_bf_mantra),col='black',alpha=0.4,data=df[df$chr==c & df$ps>=min(a3$ps_start) & df$ps<=max(a3$ps_stop),])+
				geom_point(aes(x=ps,y=-log10_bf_mash),col='black',alpha=0.4,data=df[df$chr==c & df$ps>=min(a3$ps_start) & df$ps<=max(a3$ps_stop),])+
				geom_hline(yintercept=-1,col='red',linetype=3,size=2)+
				geom_hline(yintercept=0,col='black',linetype=1,size=0.5)+
				geom_hline(yintercept=5,col='red',linetype=3,size=2)+
				geom_vline(aes(xintercept=ps),data=s,col='orange',linetype=1,size=1)+
				facet_grid(pheno~chr,scales='free_x',space='free_x')+
				theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))+
				xlim(min(a3$ps_start),max(a3$ps_stop))+
				xlab('position in bp')+ylab('log10(BF): positive MANTRA, negative Mash')+
				theme_bw()+
				theme(axis.text.x=element_blank(),axis.text.y=element_text(size=12),axis.title.y=element_text(size=20),axis.title.x=element_blank())+
				theme(strip.text.x=element_text(size=15),strip.text.y=element_text(size=15))
				annota<-annot[annot$chr==c & annot$Start>=min(a3$ps_start) & annot$Stop<=max(a3$ps_stop),]
				annotb<-annot[annot$chr==c & annot$Start>=min(s$ps) & annot$Stop<=max(s$ps),]
				annotc<-annot[(annot$chr==c & annot$Start<=min(s$ps) & annot$Stop>=min(s$ps)) | (annot$chr==c & annot$Start<=max(s$ps) & annot$Stop>=max(s$ps)),]
			g2<-ggplot()+
				geom_segment(aes(x=start,xend=stop,y=Locus,yend=Locus),arrow=arrow(length=unit(0.2,'cm')),col='black',size=1,data=annota)+
				geom_segment(aes(x=start,xend=stop,y=Locus,yend=Locus),arrow=arrow(length=unit(0.2,'cm')),col='red',size=2,data=annotb)+
				geom_segment(aes(x=start,xend=stop,y=Locus,yend=Locus),arrow=arrow(length=unit(0.2,'cm')),col='red',size=2,data=annotc)+
				geom_vline(aes(xintercept=ps),data=s,col='orange',linetype=1,size=1)+
				ylab('Genes')+xlab('position in bp')+
				xlim(min(a3$ps_start),max(a3$ps_stop))+
				theme_bw()+
				theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),axis.title.y=element_text(size=20),axis.title.x=element_text(size=20))
			print(plot_grid(g1,g2,nrow=2,ncol=1,align='v',axis='lr',rel_heights=c(0.5,0.5)))	
			dev.off()
			}
		}
}
