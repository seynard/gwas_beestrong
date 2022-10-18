library(data.table)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(mashr)
library(reshape2)
library(ashr)
library(flashr)

#plot regions overlap
all<-fread('summary_gwas_all_nus.txt.bz2')
chr_pos<-all%>%group_by(CHROM)%>%summarise(chr_size=max(POS))
all<-subset(all,all$qvalue<0.2 & all$type=='freq_freq')
lig<-fread('summary_gwas_Ligustica_Carnica.txt.bz2')
lig<-subset(lig,lig$qvalue<0.2)
mel<-fread('summary_gwas_Mellifera_ncorse.txt.bz2')
mel<-subset(mel,mel$qvalue<0.2)
hyb<-fread('summary_gwas_hybrid.txt.bz2')
hyb<-subset(hyb,hyb$qvalue<0.2)
corse<-fread('summary_gwas_corse.txt.bz2')
corse<-subset(corse,corse$qvalue<0.2)
df<-do.call(rbind,list(all,lig,mel,hyb,corse))
df$pheno[df$pheno=='varroaphoretic']<-'v_pho'
df$pheno[df$pheno=='varroadepthmitoracine4']<-'v_mito'
df$pheno[df$pheno=='logit_varroainfestation']<-'v_inf'
df$pheno[df$pheno=='logit_taux_reop_inf']<-'reop_inf'
df$pheno[df$pheno=='chargevarroaracine4']<-'v_load'
df$pheno[df$pheno=='eb_smr']<-'mnr'
df$sp[df$sp=='all_nus']<-'all'
df$sp[df$sp=='Ligustica_Carnica']<-'lig'
df$sp[df$sp=='Mellifera_ncorse']<-'mel'
df$sp[df$sp=='hybrid']<-'hyb'
df$sp[df$sp=='corse']<-'corse'
df$sp<-factor(df$sp,levels=c('all','lig','mel','hyb','corse'))
df$pheno<-factor(df$pheno,levels=c('v_pho','v_mito','v_inf','v_load','reop_inf','mnr'))
png('gwas_combine.png',width=3000,height=1000)
g1<-ggplot()+
	geom_point(aes(y=qvalue,x=POS,col=sp),data=df[df$type=='freq_freq',])+
	scale_colour_manual(values=c('black','goldenrod2','grey30','cyan','brown'))+
	geom_vline(aes(xintercept=1),col='red',lty=3,data=chr_pos)+
	geom_vline(aes(xintercept=chr_size),col='red',lty=3,data=chr_pos)+
	facet_grid(pheno~CHROM,scale='free',space='free')+
	ylim(0.2,0)+
	theme_bw()
g2<-ggplot()+
	geom_point(aes(y=qvalue,x=POS,col=sp),data=df[df$type=='freq_egs',])+
	scale_colour_manual(values=c('goldenrod2','grey30','cyan','brown'))+
	geom_vline(aes(xintercept=1),col='red',lty=3,data=chr_pos)+
	geom_vline(aes(xintercept=chr_size),col='red',lty=3,data=chr_pos)+
	facet_grid(pheno~CHROM,scale='free',space='free')+
	ylim(0.2,0)+
	theme_bw()
grid.arrange(g1,g2,nrow=2)
dev.off()

png('gwas_combine_nocorse.png',width=3000,height=1000)
g1<-ggplot()+
	geom_point(aes(y=qvalue,x=POS,col=sp),data=df[df$type=='freq_freq' & df$sp!='corse',])+
	scale_colour_manual(values=c('black','goldenrod2','grey30','cyan','brown'))+
	geom_vline(aes(xintercept=1),col='red',lty=3,data=chr_pos)+
	geom_vline(aes(xintercept=chr_size),col='red',lty=3,data=chr_pos)+
	facet_grid(pheno~CHROM,scale='free',space='free')+
	ylim(0.2,0)+
	theme_bw()
g2<-ggplot()+
	geom_point(aes(y=qvalue,x=POS,col=sp),data=df[df$type=='freq_egs' & df$sp!='corse',])+
	scale_colour_manual(values=c('goldenrod2','grey30','cyan','brown'))+
	geom_vline(aes(xintercept=1),col='red',lty=3,data=chr_pos)+
	geom_vline(aes(xintercept=chr_size),col='red',lty=3,data=chr_pos)+
	facet_grid(pheno~CHROM,scale='free',space='free')+
	ylim(0.2,0)+
	theme_bw()
grid.arrange(g1,g2,nrow=2)
dev.off()

#mash on all subspecies 6 pheno

subsp='all_nus'
phe<-c('varroaphoretic','varroadepthmitoracine4','logit_varroainfestation','chargevarroaracine4','eb_smr','logit_taux_reop_inf')
i=1
df<-fread(paste0('gemma_',subsp,'_',phe[i],'_freq_lmm_freq_cov.assoc.txt'),data.table=F)
x<-ash(df$beta,df$se,'uniform')
lfsr<-get_lfsr(x)
n<-which(lfsr<0.2)




subsp='all_nus'
phe<-c('varroaphoretic','varroadepthmitoracine4','logit_varroainfestation','chargevarroaracine4','eb_smr','logit_taux_reop_inf')
b<-list()
s<-list()
for(i in 1:length(phe)){
	df<-fread(paste0('gemma_',subsp,'_',phe[i],'_freq_lmm_freq_cov.assoc.txt'),data.table=F)
	a<-colsplit(df$rs,':',c('CHROM','POS'))
	df$CHROM<-a$CHROM
	df$POS<-a$POS
	b[[i]]<-df[,c('rs','beta')]
	colnames(b[[i]])<-c('rs',paste0(subsp,'_',phe[i]))
	s[[i]]<-df[,c('rs','se')]
	colnames(s[[i]])<-c('rs',paste0(subsp,'_',phe[i]))
	}
B<-Reduce(function(x,y) merge(x,y,by='rs',all=T),b)
S<-Reduce(function(x,y) merge(x,y,by='rs',all=T),s)
if(all(B$rs==S$rs)){
	print('order ok')
}else{
	print('reorder')
	B<-B[order(B$rs),]
	S<-S[order(S$rs),]
	}
RS<-B$rs
B$rs<-NULL
S$rs<-NULL
B<-as.matrix(B)
S<-as.matrix(S)
dat<-mash_set_data(B,S)
m1by1<-mash_1by1(mash_set_data(B,S))
strong_subset<-get_significant_results(m1by1,0.5)
random_subset<-sample(1:nrow(B),5000)
dtmp<-mash_set_data(B[random_subset,],S[random_subset,])
Vhat<-estimate_null_correlation_simple(dtmp)
rm(dtmp)
drand<-mash_set_data(B[random_subset,],S[random_subset,],V=Vhat)
dstrong<-mash_set_data(B[strong_subset,],S[strong_subset,], V=Vhat)
Upca<-cov_pca(dstrong,3)
Ued<-cov_ed(dstrong,Upca)
Uc<-cov_canonical(drand)
m<-mash(drand,Ulist=c(Ued,Uc),outputlevel=1)
m2<-mash(dstrong,g=get_fitted_g(m),fixg=TRUE)
get_lfsr(m2)
get_significant_results(m2,thresh=0.1)
get_estimated_pi(m2)
get_log10bf(m2)
get_n_significant_conditions(m2)
get_pairwise_sharing(m2)
logBF<-data.frame('rs'=RS,'log10_bf'=get_log10bf(m2))
a<-colsplit(logBF$rs,':',c('chr','ps'))
logBF$chr<-a$chr
logBF$ps<-a$ps
logBF$sign<-0
logBF$sign[logBF$log10_bf>3]<-1

ggplot(logBF)+
  geom_point(aes(x=ps,y=log10_bf,col=sign))+
  facet_grid(.~chr,space='free_x',scales='free_x')+
  theme_bw()+
  theme(legend.position='none')


#mash 1 phenotype
pheno<-'logit_taux_reop_inf'
lig<-fread(paste0('gemma_Ligustica_Carnica_',pheno,'_freq_lmm_egs_cov.assoc.txt'),data.table=F)
lig<-subset(lig,-log10(lig$p_wald)>1)
b_lig<-lig[,c('rs','beta')]
colnames(b_lig)[2]<-'lig'
s_lig<-lig[,c('rs','se')]
colnames(s_lig)[2]<-'lig'
mel<-fread(paste0('gemma_Mellifera_ncorse_',pheno,'_freq_lmm_egs_cov.assoc.txt'),data.table=F)
mel<-subset(mel,-log10(mel$p_wald)>1)
b_mel<-mel[,c('rs','beta')]
colnames(b_mel)[2]<-'mel'
s_mel<-mel[,c('rs','se')]
colnames(s_mel)[2]<-'mel'
hyb<-fread(paste0('gemma_hybrid_',pheno,'_freq_lmm_egs_cov.assoc.txt'),data.table=F)
hyb<-subset(hyb,-log10(hyb$p_wald)>1)
b_hyb<-hyb[,c('rs','beta')]
colnames(b_hyb)[2]<-'hyb'
s_hyb<-hyb[,c('rs','se')]
colnames(s_hyb)[2]<-'hyb'
corse<-fread(paste0('gemma_corse_',pheno,'_freq_lmm_egs_cov.assoc.txt'),data.table=F)
corse<-subset(corse,-log10(corse$p_wald)>1)
b_corse<-corse[,c('rs','beta')]
colnames(b_corse)[2]<-'corse'
s_corse<-corse[,c('rs','se')]
colnames(s_corse)[2]<-'corse'

bhat<-Reduce(function(x,y) merge(x,y,by='rs',all=T),list(b_lig,b_mel,b_hyb,b_corse))
shat<-Reduce(function(x,y) merge(x,y,by='rs',all=T),list(s_lig,s_mel,s_hyb,s_corse))
#n<-sample(seq(1:nrow(bhat)),100000,replace=F)
#bhat_sub<-bhat[n,]
#shat_sub<-shat[n,]
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

m1by1<-mash_1by1(mash_set_data(bhat,shat))
strong_subset<-get_significant_results(m1by1,0.1)
random_subset<-sample(1:nrow(bhat),50000)
dtmp<-mash_set_data(bhat[random_subset,],shat[random_subset,])
Vhat<-estimate_null_correlation_simple(dtmp)
rm(dtmp)
drand<-mash_set_data(bhat[random_subset,],shat[random_subset,],V=Vhat)
dstrong<-mash_set_data(bhat[strong_subset,],shat[strong_subset,], V=Vhat)
Upca<-cov_pca(dstrong,3)
Ued<-cov_ed(dstrong,Upca)
Uc<-cov_canonical(drand)
m<-mash(drand,Ulist=c(Ued,Uc),outputlevel=1)
m2<-mash(dstrong,g=get_fitted_g(m),fixg=TRUE)
get_lfsr(m2)
get_significant_results(m2,thresh=0.1)
get_estimated_pi(m2)
get_log10bf(m2)
get_n_significant_conditions(m2)
get_pairwise_sharing(m2)
logBF<-data.frame('rs'=hat_rs[strong_subset],'log10_bf'=get_log10bf(m2),'n_sign'=get_n_significant_conditions(m2))
a<-colsplit(logBF$rs,':',c('chr','ps'))
logBF$chr<-a$chr
logBF$ps<-a$ps

logBF<-logBF[order(-logBF$n_sign),]
logBF_all<-logBF[logBF$n_sign==4,]
logBF_all<-logBF_all[order(-logBF_all$log10_bf),]

data<-mash_set_data(bhat,shat,alpha=0)
Uc<-cov_canonical(data)
Up<-cov_pca(data,3)
U1<-data.frame(c(1,cor.test(bhat[,1],bhat[,2],na.omit=T)$estimate,cor.test(bhat[,1],bhat[,3],na.omit=T)$estimate,cor.test(bhat[,1],bhat[,4],na.omit=T)$estimate),
c(cor.test(bhat[,1],bhat[,2],na.omit=T)$estimate,1,cor.test(bhat[,2],bhat[,3],na.omit=T)$estimate,cor.test(bhat[,2],bhat[,4],na.omit=T)$estimate),
c(cor.test(bhat[,1],bhat[,3],na.omit=T)$estimate,cor.test(bhat[,2],bhat[,3],na.omit=T)$estimate,1,cor.test(bhat[,3],bhat[,4],na.omit=T)$estimate),
c(cor.test(bhat[,1],bhat[,4],na.omit=T)$estimate,cor.test(bhat[,2],bhat[,4],na.omit=T)$estimate,cor.test(bhat[,3],bhat[,4],na.omit=T)$estimate,1))
U1<-as.matrix(U1)
U2<-data.frame(c(1,0,0.3,0.3),
c(0,1,0.3,0.6),
c(0.3,0.3,1,0.5),
c(0.3,0.6,0.5,1))
U2<-as.matrix(U2)
Um<-list(U1,U2)
#Uf<-cov_flash(data)
#Uu<-cov_udi(data)
#U<-c(Uc,Uf,Up,Uu)
U<-c(Uc,Up,Um)
res=mash(data,U)
sign<-get_significant_results(res,thresh=0.1,conditions=NULL,sig_fn=get_lfsr)
whichsign<-hat_rs[sign]
nsign<-as.data.frame(get_n_significant_conditions(res,thresh=0.1))
nsign$id<-hat_rs
nsign<-subset(nsign,nsign[,1]>0)

#mash all phenotype all pop
#https://www.sciencedirect.com/science/article/pii/S0098847221003701
POP<-'Ligustica_Carnica'
pop<-'lig'
X<-fread(paste0('summary_gwas_',POP,'.txt.bz2'),data.table=F)
X<-X[X$type=='freq_egs',]

bhat<-X[,c('rs','beta','pheno')]
bhat<-spread(bhat,pheno,beta)
colnames(bhat)<-paste0(colnames(bhat),'_',pop)
colnames(bhat)[1]<-'rs'
shat<-X[,c('rs','se','pheno')]
shat<-spread(shat,pheno,se)
colnames(shat)<-paste0(colnames(shat),'_',pop)
colnames(shat)[1]<-'rs'
if(all(bhat$rs==shat$rs)){
	print('order ok')
	rs_hat<-bhat$rs
	bhat$rs<-NULL
	shat$rs<-NULL
}else{
	print('reorder')
	bhat<-bhat[order(bhat$rs),]
	shat<-shat[order(shat$rs),]
	rs_hat<-bhat$rs
	bhat$rs<-NULL
	shat$rs<-NULL
	}
bhat<-as.matrix(bhat)
shat<-as.matrix(shat)
sign_hat<-unique(X[X$svalue<=0.1,'rs'])
sign_hat<-data.frame(sign_hat)

strong_subset<-which(rs_hat%in%sign_hat[,1])
random_subset<-sample(1:nrow(bhat),10000,replace=F)
dtmp<-mash_set_data(bhat[random_subset,],shat[random_subset,])
Vhat<-estimate_null_correlation_simple(dtmp)
rm(dtmp)
drand<-mash_set_data(bhat[random_subset,],shat[random_subset,],V=Vhat)
dstrong<-mash_set_data(bhat[strong_subset,],shat[strong_subset,], V=Vhat)
dat<-mash_set_data(bhat,shat)
Upca<-cov_pca(dstrong,3)
Ued<-cov_ed(dstrong,Upca)
Uc<-cov_canonical(drand)
U1<-matrix(nrow=ncol(bhat),ncol=ncol(bhat))
for(i in 1:ncol(bhat)){
	for(j in 1:ncol(bhat)){
		U1[i,j]<-cor.test(bhat[,i],bhat[,j],na.omt=T)$estimate
		}
	}
U1<-list(U1)	
m<-mash(drand,Ulist=c(Ued,Uc,U1),outputlevel=1)
m2<-mash(dat,g=get_fitted_g(m),fixg=TRUE)
get_lfsr(m2)
get_significant_results(m2,thresh=0.1)
get_estimated_pi(m2)
get_log10bf(m2)
get_n_significant_conditions(m2)
get_pairwise_sharing(m2)
logBF<-data.frame('rs'=hat_rs[strong_subset],'log10_bf'=get_log10bf(m2),'n_sign'=get_n_significant_conditions(m2))
a<-colsplit(logBF$rs,':',c('chr','ps'))
logBF$chr<-a$chr
logBF$ps<-a$ps

logBF<-logBF[order(-logBF$n_sign),]
logBF_all<-logBF[logBF$n_sign==4,]
logBF_all<-logBF_all[order(-logBF_all$log10_bf),]


























mel<-fread('summary_gwas_Mellifera_ncorse.txt.bz2')
b_mel<-mel[,c('rs','beta','type','pheno','sp')]
s_mel<-mel[,c('rs','se','type','pheno','sp')]
sign_mel<-mel[mel$svalue<=0.1,c('rs','type','pheno','sp')]
hyb<-fread('summary_gwas_hybrid.txt.bz2')
b_hyb<-hyb[,c('rs','beta','type','pheno','sp')]
s_hyb<-hyb[,c('rs','se','type','pheno','sp')]
sign_hyb<-hyb[hyb$svalue<=0.1,c('rs','type','pheno','sp')]
corse<-fread('summary_gwas_corse.txt.bz2')
b_corse<-corse[,c('rs','beta','type','pheno','sp')]
s_corse<-corse[,c('rs','se','type','pheno','sp')]
sign_corse<-corse[corse$svalue<=0.1,c('rs','type','pheno','sp')]

df<-do.call(rbind,list(all,lig,mel,hyb,corse))

df$pheno[df$pheno=='varroaphoretic']<-'v_pho'
df$pheno[df$pheno=='varroadepthmitoracine4']<-'v_mito'
df$pheno[df$pheno=='logit_varroainfestation']<-'v_inf'
df$pheno[df$pheno=='logit_taux_reop_inf']<-'reop_inf'
df$pheno[df$pheno=='chargevarroaracine4']<-'v_load'
df$pheno[df$pheno=='eb_smr']<-'mnr'
df$sp[df$sp=='all_nus']<-'all'
df$sp[df$sp=='Ligustica_Carnica']<-'lig'
df$sp[df$sp=='Mellifera_ncorse']<-'mel'
df$sp[df$sp=='hybrid']<-'hyb'
df$sp[df$sp=='corse']<-'corse'






















m1by1<-mash_1by1(mash_set_data(bhat,shat))
strong_subset<-get_significant_results(m1by1,0.1)
random_subset<-sample(1:nrow(bhat),50000)
dtmp<-mash_set_data(bhat[random_subset,],shat[random_subset,])
Vhat<-estimate_null_correlation_simple(dtmp)
rm(dtmp)
drand<-mash_set_data(bhat[random_subset,],shat[random_subset,],V=Vhat)
dstrong<-mash_set_data(bhat[strong_subset,],shat[strong_subset,], V=Vhat)
dat<-mash_set_data(bhat,shat)
Upca<-cov_pca(dstrong,3)
Ued<-cov_ed(dstrong,Upca)
Uc<-cov_canonical(drand)

# add matrices: correlation between effects, correlation between phenotypes, fst between pop
U1<-data.frame(c(1,cor.test(bhat[,1],bhat[,2],na.omit=T)$estimate,cor.test(bhat[,1],bhat[,3],na.omit=T)$estimate,cor.test(bhat[,1],bhat[,4],na.omit=T)$estimate),
c(cor.test(bhat[,1],bhat[,2],na.omit=T)$estimate,1,cor.test(bhat[,2],bhat[,3],na.omit=T)$estimate,cor.test(bhat[,2],bhat[,4],na.omit=T)$estimate),
c(cor.test(bhat[,1],bhat[,3],na.omit=T)$estimate,cor.test(bhat[,2],bhat[,3],na.omit=T)$estimate,1,cor.test(bhat[,3],bhat[,4],na.omit=T)$estimate),
c(cor.test(bhat[,1],bhat[,4],na.omit=T)$estimate,cor.test(bhat[,2],bhat[,4],na.omit=T)$estimate,cor.test(bhat[,3],bhat[,4],na.omit=T)$estimate,1))
U1<-as.matrix(U1)
U2<-data.frame(c(1,0,0.3,0.3),
c(0,1,0.3,0.6),
c(0.3,0.3,1,0.5),
c(0.3,0.6,0.5,1))




m<-mash(drand,Ulist=c(Ued,Uc),outputlevel=1)
m2<-mash(data,g=get_fitted_g(m),fixg=TRUE)
get_lfsr(m2)
get_significant_results(m2,thresh=0.1)
get_estimated_pi(m2)
get_log10bf(m2)
get_n_significant_conditions(m2)
get_pairwise_sharing(m2)
logBF<-data.frame('rs'=hat_rs[strong_subset],'log10_bf'=get_log10bf(m2),'n_sign'=get_n_significant_conditions(m2))
a<-colsplit(logBF$rs,':',c('chr','ps'))
logBF$chr<-a$chr
logBF$ps<-a$ps

logBF<-logBF[order(-logBF$n_sign),]
logBF_all<-logBF[logBF$n_sign==4,]
logBF_all<-logBF_all[order(-logBF_all$log10_bf),]


#######################################################################
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


