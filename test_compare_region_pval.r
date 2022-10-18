library(data.table)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(nlme)
library(lme4)
library(lattice)
library(ggpubr)

RESULTS<-list()
TEST_EFFECTS<-list()
KS_TEST<-list()
dfx<-data.frame(version=c(1,1,1,2,2,2),pheno1=c('mnr','recap_inf','pc1_varroa_inf','mnr','recap_inf','pc1_varroa_inf'),pheno2=c('eb_smr','logit_taux_reop_inf','pc1','eb_smr','logit_taux_reop_inf','pc1'))
for(x in 1:nrow(dfx)){
version<-dfx$version[x]
pheno1<-dfx$pheno1[x]
pheno2<-dfx$pheno2[x]
res<-fread(paste0('results/sign_locus_',version,'.txt'),data.table=F)
res<-res[res$pheno==pheno1 & res$type=='egs',]
mantra<-fread(paste0('results/run_mantra_',pheno2,'_egs_',version,'/mantra.out'),data.table=F)
colnames(mantra)<-c('rs','chr','ps','allele1','allele0','n_study','log10BF_mantra','post_proba_assoc','n','direction')
mash<-fread(paste0('results/result_mash_',pheno2,'_egs_',version,'_simple_corr.txt'),data.table=F)
colnames(mash)<-c('rs','log10BF_mash','n_sign','loglik','chr','ps')
df<-merge(mantra,mash,by=c('rs','chr','ps'),all=T)
chr_size<-df%>%group_by(chr)%>%summarise(size=max(ps))
res<-unique(res[,c("chr","rs","ps","all_effect","all_other","log10BF_mantra","post_proba_hetero_mantra","log10BF_mash","loglik_mash","pheno","type")])
res<-res[order(res$chr,res$ps),]
RES<-list()
G<-list()
DOTPLOT<-list()
RANEF<-list()
KS<-list()
size<-500
for(j in 1:nrow(res)){
print(j)	
y<-df[df$ps>=res$ps[j]-size & df$ps<=res$ps[j]+size & df$chr==res$chr[j],]
y$type<-'reg'
df_chr<-subset(df,df$chr==res$chr[j])
Z<-list()
for(i in 1:10){
	n<-sample(seq(1:chr_size$size[chr_size$chr==res$chr[j]]),1)
	a<-(n-size)
	b<-(n+size)
	choice<-df_chr$rs[df_chr$ps>=a & df_chr$ps<=b]
	while(length(choice)<(nrow(y)-0.1*nrow(y)) | length(choice)>(nrow(y)+0.1*nrow(y)) |any(res$rs%in%choice) ){
	n<-sample(seq(1:chr_size$size[chr_size$chr==res$chr[j]]),1)
	a<-(n-size)
	b<-(n+size)
	choice<-df_chr$rs[df_chr$ps>=a & df_chr$ps<=b]
	}
	z<-df_chr[df_chr$rs%in%choice,]
	z$type<-paste0('test',i)
	Z[[i]]<-z}
z<-do.call(rbind,Z)
X<-rbind(y,z)
ks<-list()
for(i in 1:10){
ks[[i]]<-round(ks.test(X$log10BF_mantra[X$type=='reg'],X$log10BF_mantra[X$type==paste0('test',i)])$p.val,digit=5)
}
ks<-as.data.frame(do.call(rbind,ks))
colnames(ks)<-'p_val_kstest'
ks$type<-paste0('test',seq(1:10))
ks$rs<-res$rs[j]
KS[[j]]<-ks
resume<-X%>%group_by(type)%>%summarise(n_snp=n(),min_ps=min(ps),max_ps=max(ps),m_mantra=mean(log10BF_mantra),sd_mantra=sd(log10BF_mantra),m_mash=mean(log10BF_mash),sd_mash=sd(log10BF_mash))
resume$size<-resume$max_ps-resume$min_ps
resume$rs<-res$rs[j]
#print(resume)
RES[[j]]<-resume
g1<-ggplot()+geom_density(aes(x=log10BF_mantra),col='black',data=X[X$type=='reg',])+geom_density(aes(x=log10BF_mantra,col=type),lty=3,data=X[X$type!='reg',])+theme_bw()+theme(legend.position='none')+ggtitle(res$rs[j])
g2<-ggplot()+geom_density(aes(x=log10BF_mash),col='black',data=X[X$type=='reg',])+geom_density(aes(x=log10BF_mash,col=type),lty=3,data=X[X$type!='reg',])+theme_bw()+theme(legend.position='none')
G[[j]]<-ggarrange(g1,g2,nrow=2,heights=c(0.55,0.45))
fit_lmer_mantra<-lmer(log10BF_mantra~1+(1|type),data=X)
#print(summary(fit_lmer_mantra))
ran_mantra<-as.data.frame(ranef(fit_lmer_mantra, condVar=TRUE))
ran_mantra$rs<-res$rs[j]
ran_mantra[,c('grpvar','term')]<-NULL
colnames(ran_mantra)<-c('type','val','sd','rs')
ran_mantra$ci_min<-ran_mantra$val-1.96*ran_mantra$sd
ran_mantra$ci_max<-ran_mantra$val+1.96*ran_mantra$sd
ran_mantra$group<-'mantra' 
fit_lmer_mash<-lmer(log10BF_mash~1+(1|type),data=X)
#print(summary(fit_lmer_mash))
ran_mash<-as.data.frame(ranef(fit_lmer_mash, condVar=TRUE))
ran_mash$rs<-res$rs[j]
ran_mash[,c('grpvar','term')]<-NULL
colnames(ran_mash)<-c('type','val','sd','rs')
ran_mash$ci_min<-ran_mash$val-1.96*ran_mash$sd
ran_mash$ci_max<-ran_mash$val+1.96*ran_mash$sd
ran_mash$group<-'mash' 
ran<-rbind(ran_mantra,ran_mash)
ran$type<-factor(ran$type,levels=c('test10','test9','test8','test7','test6','test5','test4','test3','test2','test1','reg'))
DOTPLOT[[j]]<-ggplot(ran)+
	geom_point(aes(x=val,y=type))+
	geom_vline(aes(xintercept=0),col='red',lty=3)+
	geom_segment(aes(x=ci_min,xend=ci_max,y=type,yend=type))+
	facet_wrap(~group,scale='free_x')+theme_bw()+ggtitle(res$rs[j])
RANEF[[j]]<-ran
}
RES<-do.call(rbind,RES)
RES<-as.data.frame(RES)
RES$version<-version
RES$pheno1<-pheno1
RANEF<-do.call(rbind,RANEF)
RANEF$version<-version
RANEF$pheno1<-pheno1
KS<-do.call(rbind,KS)
KS$version<-version
KS$pheno1<-pheno1
png(paste0('plot_sign_regions_density_',pheno1,'_',version,'.png'),width=4000,height=5000,res=300)
plot_density<-do.call("grid.arrange",c(G,ncol=5))            
dev.off()       
png(paste0('plot_sign_regions_effects_',pheno1,'_',version,'.png'),width=4000,height=5000,res=300)
plot_dot<-do.call("grid.arrange",c(DOTPLOT,ncol=5))                   
dev.off()       
RESULTS[[x]]<-RES
TEST_EFFECTS[[x]]<-RANEF
KS_TEST[[x]]<-KS
}
RESULTS<-do.call(rbind,RESULTS)
TEST_EFFECTS<-do.call(rbind,TEST_EFFECTS)
KS_TEST<-do.call(rbind,KS_TEST)
write.table(RESULTS,'sign_regions_results.txt',col.names=T,row.names=F,quote=F)
write.table(TEST_EFFECTS,'sign_regions_test_effects.txt',col.names=T,row.names=F,quote=F)
write.table(KS_TEST,'sign_regions_ks.txt',col.names=T,row.names=F,quote=F)


