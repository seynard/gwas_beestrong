library(data.table)
library(mashr)
library(reshape2)
library(ashr)
library(rmeta)
library(tidyverse)
args<-commandArgs(TRUE)
dirout<-args[1]
pheno<-args[2]
type<-args[3]
threshold<-as.numeric(args[4])
list_pop<-args[5]
list_pop<-unlist(strsplit(list_pop,','))
version<-args[6]
#dirout<-'results'
#Pheno<-c('smr_brut','pc1','logit_taux_reop_inf','chargevarroaracine4','eb_smr','logit_taux_reop','logit_varroainfestation','varroadepthmitoracine4','varroaphoretic')
#type<-'freq' #'egs'
#threshold<-as.numeric(0.2)
#list_pop<-c('Mellifera_ncorse','Ligustica_Carnica','hybrid_corse')

print(pheno)
X<-list()
for(i in 1:length(list_pop)){X[[i]]<-fread(paste0(dirout,'/summary_gwas_',list_pop[i],'_',pheno,'_freq_',type,'.txt.bz2'),data.table=F)}

bhat_x<-list()
shat_x<-list()
for(i in 1:length(list_pop)){
	bhat_x[[i]]<-X[[i]][,c('rs','beta')]
	colnames(bhat_x[[i]])[2]<-list_pop[i]
	shat_x[[i]]<-X[[i]][,c('rs','se')]
	colnames(shat_x[[i]])[2]<-list_pop[i]
}
bhat<-Reduce(function(x,y) merge(x,y,by='rs',all=T),bhat_x)
shat<-Reduce(function(x,y) merge(x,y,by='rs',all=T),shat_x)
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

#1 no correlation
dat<-mash_set_data(bhat,shat)
ULIST<-vector()
Upca<-cov_pca(dat,3)
Ued<-cov_ed(dat,Upca)
ULIST<-c(ULIST,Ued)
Uc<-cov_canonical(dat)
ULIST<-c(ULIST,Uc)
Umantra<-data.frame(c(0,0.3,0.2),c(0.3,0,0.1),c(0.2,0.1,0))
Umantra<-1-Umantra
Umantra<-as.matrix(Umantra)
Umantra<-list(Umantra)
names(Umantra)<-'mantra_fst'
ULIST<-c(ULIST,Umantra)
bhat_com<-bhat[complete.cases(bhat),]
if(nrow(bhat_com)>3){
Ueffect<-matrix(nrow=ncol(bhat_com),ncol=ncol(bhat_com))
for(i in 1:ncol(bhat_com)){
	for(j in 1:ncol(bhat_com)){
		Ueffect[i,j]<-cor.test(bhat_com[,i],bhat_com[,j],na.omit=T)$estimate
		}
	}
Ueffect<-list(Ueffect)
names(Ueffect)<-'correlation_effect'
ULIST<-c(ULIST,Ueffect)
}
m<-mash(dat,Ulist=ULIST)
print('pi no correlation')
print(get_estimated_pi(m))
print('sharing no correlation')
print(get_pairwise_sharing(m))
logBF<-data.frame('rs'=hat_rs,'log10_bf'=get_log10bf(m),'n_sign'=get_n_significant_conditions(m),'loglik'=get_loglik(m))
a<-colsplit(logBF$rs,':',c('chr','ps'))
logBF$chr<-a$chr
logBF$ps<-a$ps
write.table(logBF,paste0(dirout,'/result_mash_',pheno,'_',type,'_',version,'_no_corr.txt'),col.names=T,row.names=F,quote=F)
post_mean<-as.data.frame(mash_compute_posterior_matrices(m,dat)$PosteriorMean)
post_mean$rs<-hat_rs
post_mean<-gather(post_mean,condition,post_mean,colnames(post_mean)[colnames(post_mean)!='rs'])
post_sd<-as.data.frame(mash_compute_posterior_matrices(m,dat)$PosteriorSD)
post_sd$rs<-hat_rs
post_sd<-gather(post_sd,condition,post_sd,colnames(post_sd)[colnames(post_sd)!='rs'])
post_lfsr<-as.data.frame(mash_compute_posterior_matrices(m,dat)$lfsr)
post_lfsr$rs<-hat_rs
post_lfsr<-gather(post_lfsr,condition,post_lfsr,colnames(post_lfsr)[colnames(post_lfsr)!='rs'])
post_lfdr<-as.data.frame(mash_compute_posterior_matrices(m,dat)$lfdr)
post_lfdr$rs<-hat_rs
post_lfdr<-gather(post_lfdr,condition,post_lfdr,colnames(post_lfdr)[colnames(post_lfdr)!='rs'])
post<-Reduce(function(x,y) merge(x,y,by=c('rs','condition'),all=T),list(post_mean,post_sd,post_lfsr,post_lfdr))
write.table(post,paste0(dirout,'/summary_mash_',pheno,'_',type,'_',version,'_no_corr.txt'),col.names=T,row.names=F,quote=F)

#2 simple correlation
V.simple=estimate_null_correlation_simple(dat)
print('correlation simple')
print(V.simple)
dat.Vsimple=mash_update_data(dat,V=V.simple)
ULIST<-vector()
Upca<-cov_pca(dat.Vsimple,3)
Ued<-cov_ed(dat.Vsimple,Upca)
ULIST<-c(ULIST,Ued)
Uc<-cov_canonical(dat.Vsimple)
ULIST<-c(ULIST,Uc)
Umantra<-data.frame(c(0,0.3,0.2),c(0.3,0,0.1),c(0.2,0.1,0))
Umantra<-1-Umantra
Umantra<-as.matrix(Umantra)
Umantra<-list(Umantra)
ULIST<-c(ULIST,Umantra)
bhat_com<-bhat[complete.cases(bhat),]
if(nrow(bhat_com)>3){
Ueffect<-matrix(nrow=ncol(bhat_com),ncol=ncol(bhat_com))
for(i in 1:ncol(bhat_com)){
	for(j in 1:ncol(bhat_com)){
		Ueffect[i,j]<-cor.test(bhat_com[,i],bhat_com[,j],na.omit=T)$estimate
		}
	}
Ueffect<-list(Ueffect)
ULIST<-c(ULIST,Ueffect)
}
m.Vsimple<-mash(dat.Vsimple,Ulist=ULIST)
print('pi simple correlation')
print(get_estimated_pi(m.Vsimple))
print('sharing simple correlation')
print(get_pairwise_sharing(m.Vsimple))
logBF.Vsimple<-data.frame('rs'=hat_rs,'log10_bf'=get_log10bf(m.Vsimple),'n_sign'=get_n_significant_conditions(m.Vsimple),'loglik'=get_loglik(m.Vsimple))
a<-colsplit(logBF.Vsimple$rs,':',c('chr','ps'))
logBF.Vsimple$chr<-a$chr
logBF.Vsimple$ps<-a$ps
write.table(logBF.Vsimple,paste0(dirout,'/result_mash_',pheno,'_',type,'_',version,'_simple_corr.txt'),col.names=T,row.names=F,quote=F)
post_mean<-as.data.frame(mash_compute_posterior_matrices(m.Vsimple,dat)$PosteriorMean)
post_mean$rs<-hat_rs
post_mean<-gather(post_mean,condition,post_mean,colnames(post_mean)[colnames(post_mean)!='rs'])
post_sd<-as.data.frame(mash_compute_posterior_matrices(m.Vsimple,dat)$PosteriorSD)
post_sd$rs<-hat_rs
post_sd<-gather(post_sd,condition,post_sd,colnames(post_sd)[colnames(post_sd)!='rs'])
post_lfsr<-as.data.frame(mash_compute_posterior_matrices(m.Vsimple,dat)$lfsr)
post_lfsr$rs<-hat_rs
post_lfsr<-gather(post_lfsr,condition,post_lfsr,colnames(post_lfsr)[colnames(post_lfsr)!='rs'])
post_lfdr<-as.data.frame(mash_compute_posterior_matrices(m.Vsimple,dat)$lfdr)
post_lfdr$rs<-hat_rs
post_lfdr<-gather(post_lfdr,condition,post_lfdr,colnames(post_lfdr)[colnames(post_lfdr)!='rs'])
post.Vsimple<-Reduce(function(x,y) merge(x,y,by=c('rs','condition'),all=T),list(post_mean,post_sd,post_lfsr,post_lfdr))
write.table(post,paste0(dirout,'/summary_mash_',pheno,'_',type,'_',version,'_simple_corr.txt'),col.names=T,row.names=F,quote=F)

#3 EM correlation
rand<-sample(seq(1:nrow(bhat)),ceiling(nrow(bhat)/10),replace=F)
bhat_rand<-bhat[rand,]
shat_rand<-shat[rand,]
dat_rand<-mash_set_data(bhat_rand,shat_rand)
V.em=mash_estimate_corr_em(dat_rand,Ulist=ULIST,V.simple,details=F,max_iter=5,est_cor=F)
#V.em=mash_estimate_corr_em(dat,Ulist=ULIST,V.simple,details=F,max_iter=5,est_cor=F)
print('covariance em (random set)')
print(V.em)
V.em=mash_estimate_corr_em(dat_rand,Ulist=ULIST,V.simple,details=T,max_iter=5,est_cor=T)
#V.em=mash_estimate_corr_em(dat,Ulist=ULIST,V.simple,details=T,max_iter=5,est_cor=T)
print('correlation em (random set)')
print(V.em$V)
m.Vem<-V.em$m
m.Vem2=mash(dat,g=get_fitted_g(m.Vem),fixg=TRUE)
m.Vem<-m.Vem2
print('pi em correlation')
print(get_estimated_pi(m.Vem))
print('sharing em correlation')
print(get_pairwise_sharing(m.Vem))
logBF.Vem<-data.frame('rs'=hat_rs,'log10_bf'=get_log10bf(m.Vem),'n_sign'=get_n_significant_conditions(m.Vem),'loglik'=get_loglik(m.Vem))
a<-colsplit(logBF.Vem$rs,':',c('chr','ps'))
logBF.Vem$chr<-a$chr
logBF.Vem$ps<-a$ps
write.table(logBF.Vem,paste0(dirout,'/result_mash_',pheno,'_',type,'_',version,'_em_corr.txt'),col.names=T,row.names=F,quote=F)
post_mean<-as.data.frame(mash_compute_posterior_matrices(m.Vem,dat)$PosteriorMean)
post_mean$rs<-hat_rs
post_mean<-gather(post_mean,condition,post_mean,colnames(post_mean)[colnames(post_mean)!='rs'])
post_sd<-as.data.frame(mash_compute_posterior_matrices(m.Vem,dat)$PosteriorSD)
post_sd$rs<-hat_rs
post_sd<-gather(post_sd,condition,post_sd,colnames(post_sd)[colnames(post_sd)!='rs'])
post_lfsr<-as.data.frame(mash_compute_posterior_matrices(m.Vem,dat)$lfsr)
post_lfsr$rs<-hat_rs
post_lfsr<-gather(post_lfsr,condition,post_lfsr,colnames(post_lfsr)[colnames(post_lfsr)!='rs'])
post_lfdr<-as.data.frame(mash_compute_posterior_matrices(m.Vem,dat)$lfdr)
post_lfdr$rs<-hat_rs
post_lfdr<-gather(post_lfdr,condition,post_lfdr,colnames(post_lfdr)[colnames(post_lfdr)!='rs'])
post.Vem<-Reduce(function(x,y) merge(x,y,by=c('rs','condition'),all=T),list(post_mean,post_sd,post_lfsr,post_lfdr))
write.table(post,paste0(dirout,'/summary_mash_',pheno,'_',type,'_',version,'_em_corr.txt'),col.names=T,row.names=F,quote=F)

loglik=c(get_loglik(m),get_loglik(m.Vsimple),get_loglik(m.Vem2))
significant=c(length(get_significant_results(m)),length(get_significant_results(m.Vsimple)),length(get_significant_results(m.Vem2)))
tb = rbind(loglik,significant)
colnames(tb)=c('without cor','V simple','V EM')
row.names(tb)=c('log likelihood', '# significance')
print(tb)
