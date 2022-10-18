library(data.table)
library(mashr)
library(reshape2)
library(ashr)
library(rmeta)
args<-commandArgs(TRUE)
dirout<-args[1]
pheno<-args[2]
type<-args[3]
threshold<-as.numeric(args[4])
list_pop<-args[5]
list_pop<-unlist(strsplit(list_pop,','))
#dirout<-'results'
#Pheno<-c('smr_brut','pc1','logit_taux_reop_inf','chargevarroaracine4','eb_smr','logit_taux_reop','logit_varroainfestation','varroadepthmitoracine4','varroaphoretic')
#type<-'freq' #'egs'
#threshold<-as.numeric(0.2)
#list_pop<-c('Mellifera_ncorse','Ligustica_Carnica','hybrid_corse')

#for(x in 1:length(Pheno)){
#pheno<-Pheno[x]
print(pheno)
X<-list()
l_rs<-list()
l_strong<-list()
for(i in 1:length(list_pop)){
	X[[i]]<-fread(paste0(dirout,'/summary_gwas_',list_pop[i],'_',pheno,'_freq_',type,'.txt.bz2'),data.table=F)
#	l_rs[[i]]<-X[[i]]$rs
#	l_strong[[i]]<-X[[i]]$rs[X[[i]]$lfsr<threshold]
	}

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
	
#RS<-Reduce(intersect,l_rs)
#rand<-sample(RS,5000,replace=F)
#brand_x<-list()
#srand_x<-list()
#for(i in 1:length(list_pop)){
#	rand_x<-X[[i]][X[[i]]$rs%in%rand,]
#	brand_x[[i]]<-rand_x[,c('rs','beta')]
#	colnames(brand_x[[i]])[2]<-list_pop[i]
#	srand_x[[i]]<-rand_x[,c('rs','se')]
#	colnames(srand_x[[i]])[2]<-list_pop[i]
#}
#brand<-Reduce(function(x,y) merge(x,y,by='rs',all=T),brand_x)
#srand<-Reduce(function(x,y) merge(x,y,by='rs',all=T),srand_x)
#if(all(brand$rs==srand$rs)){
#	print('order ok')
#	hat_rs_rand<-brand$rs
#	brand$rs<-NULL
#	srand$rs<-NULL
#}else{
#	print('reorder')
#	brand<-brand[order(brand$rs),]
#	srand<-srand[order(srand$rs),]
#	hat_rs_rand<-brand$rs
#	brand$rs<-NULL
#	srand$rs<-NULL
#	}
#brand<-as.matrix(brand)
#srand<-as.matrix(srand)
#dtmp<-mash_set_data(brand,srand)
#Vhat<-estimate_null_correlation_simple(dtmp)
#rm(dtmp)
#drand<-mash_set_data(brand,srand,V=Vhat)

#strong<-unique(unlist(l_strong))
#bstrong<-list()
#sstrong<-list()
#for(i in 1:length(list_pop)){
#	strong_x<-X[[i]][X[[i]]$rs%in%strong,]
#	bstrong[[i]]<-strong_x[,c('rs','beta')]
#	colnames(bstrong[[i]])[2]<-list_pop[i]
#	sstrong[[i]]<-strong_x[,c('rs','se')]
#	colnames(sstrong[[i]])[2]<-list_pop[i]
#}
#bhat<-Reduce(function(x,y) merge(x,y,by='rs',all=T),bstrong)
#shat<-Reduce(function(x,y) merge(x,y,by='rs',all=T),sstrong)
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
#dstrong<-mash_set_data(bhat,shat,V=Vhat)
dat<-mash_set_data(bhat,shat)
V.simple = estimate_null_correlation_simple(dat)
dat.Vsimple = mash_update_data(dat, V=V.simple)

#if(nrow(bhat)>1){
ULIST<-vector()
Upca<-cov_pca(dat.Vsimple,3)
Ued<-cov_ed(dat.Vsimple,Upca)
ULIST<-c(ULIST,Ued)
#Upca<-cov_pca(dstrong,3)
#if(all(grepl('NA',Upca)==F)){Ued<-cov_ed(dstrong,Upca)
#ULIST<-c(ULIST,Ued)}
Uc<-cov_canonical(dat.Vsimple)
#Uc<-cov_canonical(drand)
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
m<-mash(dat,Ulist=ULIST)
#m<-mash(drand,Ulist=ULIST,outputlevel=1)
#m2<-mash(dstrong,g=get_fitted_g(m),fixg=TRUE)


print(get_estimated_pi(m2))
print(hat_rs[get_significant_results(m2,thresh=threshold)])
logBF<-data.frame('rs'=hat_rs,'log10_bf'=get_log10bf(m2),'n_sign'=get_n_significant_conditions(m2))
a<-colsplit(logBF$rs,':',c('chr','ps'))
logBF$chr<-a$chr
logBF$ps<-a$ps
write.table(logBF,paste0(dirout,'/result_mash_',pheno,'_',type,'.txt'),col.names=T,row.names=F,quote=F)

x<-mash_compute_posterior_matrices(m2,dat)
x_postb


#logBF<-logBF[order(logBF$chr,logBF$ps),]
#col_pop<-data.frame(list_pop,col=NA)
#col_pop$col[grep('Mellifera',col_pop$list_pop)]<-'grey34'
#col_pop$col[grep('Ligustica',col_pop$list_pop)]<-'goldenrod2'
#col_pop$col[grep('hybrid',col_pop$list_pop)]<-'deepskyblue2'
#pdf(paste0(dirout,'/plot_mash_',pheno,'_',type,'.pdf'),width=10,height=10)
#for(i in 1:nrow(bhat)){
#mash_plot_meta(m2,i,logBF$rs[i],colors=meta.colors(box=col_pop$col))
#}
#dev.off()
#}
#}



