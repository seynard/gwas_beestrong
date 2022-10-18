#module load system/R-3.6.1
library(data.table)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(ashr)
library(gridExtra)
library(ggpubr)
library(viridis)
library(RColorBrewer)
library(gplots)
library(reshape2)
library(qvalue)
args<-commandArgs(TRUE)
dir_out<-args[1]
population<-args[2]
#dir_out<-'/work/genphyse/dynagen/seynard/GWAS/infestation/results'
#population='hybrid'

pp.ggplot = function(pval) {
    dat = as_tibble(pval)
    nmax=nrow(dat)
    dat %<>% mutate(stat = -log10(value)) %>%
        mutate(rank = dense_rank(value),
               expect = -log10(rank/(nmax+1)),
               hival = -log10(qbeta(0.975,rank,nmax-rank)),
               loval = -log10(qbeta(0.025,rank,nmax-rank)) )
    dat %<>% arrange(expect)
    x = dat$expect
    x = c(x,rev(x))
    y = c(dat$loval,rev(dat$hival))
    poly.dat=tibble(x,y)
    chisq <- qchisq(1 - pval, 1)
	lambda <- median(chisq,na.rm=T) / qchisq(0.5, 1)
    p = dat %>% ggplot(aes(x=expect,y=stat)) +
        geom_polygon(data=poly.dat,aes(x=x,y=y),fill='blue', alpha=0.1) +
        geom_path(data=poly.dat,aes(x=x,y=y), alpha=0.1) +
        geom_point() +
        geom_abline(intercept=0,slope=1,colour='blue') +
        xlab("Expected -log(p)") +
        ylab("Observed -log(p)")+ 
		annotate(geom = "text",x = -Inf,y = Inf,hjust = -0.15,vjust = 1 + 0.15 * 3,label = sprintf("λ = %.2f", lambda),size = 8) +
		theme(axis.ticks = element_line(size = 0.5),panel.grid = element_blank())
    p
}

grm_freq<-fread(paste0(dir_out,'/LDAKgmat_',population,'_freq.grm.raw'))
grm_freq_id<-fread(paste0(dir_out,'/LDAKgmat_',population,'_freq.grm.id'),header=F,data.table=F)
colnames(grm_freq)<-grm_freq_id[,1]
grm_egs<-fread(paste0(dir_out,'/LDAKgmat_',population,'_egs.grm.raw'))
grm_egs_id<-fread(paste0(dir_out,'/LDAKgmat_',population,'_egs.grm.id'),header=F,data.table=F)
colnames(grm_egs)<-grm_egs_id[,1]

png(paste0(dir_out,'/heatmap_freq_',population,'.png'),width=1000,height=1000)
heatmap.2(as.matrix(grm_freq),scale="none",col=viridis,trace='none',dendrogram='row',labRow=grm_freq_id[,1])
dev.off()
png(paste0(dir_out,'/heatmap_egs_',population,'.png'),width=1000,height=1000)
heatmap.2(as.matrix(grm_egs),scale="none",col=viridis,trace='none',dendrogram='row',labRow=grm_egs_id[,1])
dev.off()

pheno_list<-intersect(list.files(path=paste0(dir_out,'/'),pattern='gemma_'),list.files(path=paste0(dir_out,'/'),pattern=paste0('_',population,'_')))
pheno_list<-gsub(paste0('gemma_',population,'_'),'',pheno_list)
pheno_list<-gsub('_lmm','',pheno_list)
pheno_list<-gsub('_bslmm','',pheno_list)
pheno_list<-gsub('_cov','',pheno_list)
pheno_list<-gsub('.bv.txt','',pheno_list)
pheno_list<-gsub('.gamma.txt','',pheno_list)
pheno_list<-gsub('.hyp.txt','',pheno_list)
pheno_list<-gsub('.log.txt','',pheno_list)
pheno_list<-gsub('.param.txt','',pheno_list)
pheno_list<-gsub('.assoc.txt','',pheno_list)
pheno_list<-unique(pheno_list)
pheno_list<-data.frame(pheno_list)
pheno_list$grm[grepl('_freq_egs',pheno_list[,1])]<-'freq'
pheno_list$gwas[grepl('_freq_egs',pheno_list[,1])]<-'egs'
pheno_list$grm[grepl('_freq_freq',pheno_list[,1])]<-'freq'
pheno_list$gwas[grepl('_freq_freq',pheno_list[,1])]<-'freq'
pheno_list$grm[grepl('_egs_freq',pheno_list[,1])]<-'egs'
pheno_list$gwas[grepl('_egs_freq',pheno_list[,1])]<-'freq'
pheno_list$grm[grepl('_egs_egs',pheno_list[,1])]<-'egs'
pheno_list$gwas[grepl('_egs_egs',pheno_list[,1])]<-'egs'
pheno_list$pheno_list<-gsub('_freq_freq','',pheno_list$pheno_list)
pheno_list$pheno_list<-gsub('_freq_egs','',pheno_list$pheno_list)
pheno_list$pheno_list<-gsub('_egs_freq','',pheno_list$pheno_list)
pheno_list$pheno_list<-gsub('_egs_egs','',pheno_list$pheno_list)
pheno_list$pheno_list<-gsub('_cov','',pheno_list$pheno_list)
pheno_list$name<-paste0(pheno_list$pheno_list,'_',pheno_list$grm,'_',pheno_list$gwas)
plot_done<-intersect(list.files(path=paste0(dir_out,'/'),pattern='manhattan'),list.files(path=paste0(dir_out,'/'),pattern=population))
plot_done<-gsub(paste0('manhattan_',population,'_'),'',plot_done)
plot_done<-gsub('.png','',plot_done)
pheno_list<-pheno_list[!(pheno_list$name%in%plot_done),]
for(i in 1:nrow(pheno_list)){
phenotype_id=pheno_list[i,'pheno_list']
grm_id=pheno_list[i,'grm']
gwas_id=pheno_list[i,'gwas']
prfx=paste0(population,'_',phenotype_id,'_',grm_id,'_',gwas_id)
print(prfx)

input_file_lmm=paste0(dir_out,'/gemma_',population,'_',phenotype_id,'_',grm_id,'_lmm_',gwas_id,'_cov.assoc.txt')
lmm=fread(input_file_lmm)
input_file_bslmm=paste0(dir_out,'/gemma_',population,'_',phenotype_id,'_',grm_id,'_bslmm_',gwas_id,'.param.txt')
bslmm=fread(input_file_bslmm)
gwas<-merge(lmm,bslmm,by=c('chr','ps','rs'),all=T)
a<-colsplit(gwas$rs,':',c('CHROM','POS'))
gwas$chr<-a$CHROM
gwas$ps<-a$POS
pv=as.data.frame(gwas%>%select(p_wald))
pv<-subset(pv,!is.na(pv[,1]))

png(paste0(dir_out,'/g0_',prfx,'.png'),width=1500,height=300)
gwas.ppplot<-pp.ggplot(pv[,1])+theme_bw()+ggtitle(paste0(population,'_',phenotype_id,'_',grm_id,'_',gwas_id))
p_hist<-ggplot(pv)+
	geom_histogram(aes(x=p_wald),fill='grey')+
	xlab('p-values')+theme_bw()
df<-gwas%>%sample_n(100000)
g11<-ggplot(df,aes(x=af,y=beta.x))+
    geom_point(alpha=0.01)+ 
    geom_smooth()+
    xlab("Reference Allele Frequency")+ 
    ylab("SNP Effect")+
    ggtitle(paste0(population,'_',phenotype_id,'_',grm_id,'_',gwas_id,"(100,000 random SNPs)"))+theme_bw()
g12<-ggplot(df,aes(x=af,y=beta.x))+
    geom_point(alpha=0.01)+ 
    geom_smooth()+
    ylim(-10,10)+
    xlab("Reference Allele Frequency")+ 
    ylab("SNP Effect")+
    ggtitle(paste0(population,'_',phenotype_id,'_',grm_id,'_',gwas_id,"(100,000 random SNPs)"))+theme_bw()
grid.arrange(gwas.ppplot,p_hist,g11,g12,ncol=4)
dev.off()

fit.ash=ash(gwas$beta.x,gwas$se,mixcompdist='uniform')
gwas=as_tibble(cbind(gwas,fit.ash$result))

gwas%<>%mutate(significant=svalue<0.1)

png(paste0(dir_out,'/g1_',prfx,'.png'),width=1500,height=300)
g2<-ggplot(gwas)+
    geom_density(aes(x=af,color=significant),size=2)+
    theme_bw()+ 
    xlab("Allele frequency")+
	ggtitle(paste0(population,'_',phenotype_id,'_',grm_id,'_',gwas_id))
df<-gwas%>%filter(significant)
g3<-ggplot(df,aes(x=betahat,y=PosteriorMean,colour=1/se))+
        geom_point()+ 
        scale_colour_viridis_c(guide=F)+
        geom_abline(slope=1,intercept=0)+ 
        theme_bw()
g4<-ggplot(df,aes(x=af,y=beta.x,colour=1/se))+
    geom_point(alpha=0.5)+
    ##geom_smooth()+
    scale_colour_viridis_c()+ 
    geom_vline(xintercept=c(0.1,0.9))
df<-gwas%>%filter(svalue<0.2)
g5<-ggplot(df,aes(x=ps,y=-log10(p_wald),color=significant))+ 
    geom_point()+
    facet_wrap(~chr,scale='free_x')+
	scale_x_continuous(label=function(x){ x*1e-6})+
	xlab("Position (Mb)") 
if(nrow(df)>0){grid.arrange(g2,g3,g4,g5,ncol=4)}else{grid.arrange(g2,g3,g4,ncol=4)}
dev.off()

gwas$chr<-factor(gwas$chr,levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'))
qseuil<-max(gwas$p_wald[gwas$qvalue<0.1],na.rm=T)
sseuil<-max(gwas$p_wald[gwas$svalue<0.1],na.rm=T)
max_log<-round(max(-log10(gwas$p_wald[!is.na(gwas$p_wald)])))
col_chr<-rep(c('deepskyblue','orange'),8)

g61<-ggplot()+
		geom_point(aes(x=ps,y=-log10(p_wald),col=as.factor(chr)),data=gwas)+
		geom_hline(yintercept=-log10(qseuil),col='red',lty=3)+
		geom_hline(yintercept=-log10(sseuil),col='blue',lty=3)+
		facet_grid(.~chr,scales='free',space='free')+
		theme_bw()+
		scale_colour_manual(values=col_chr)+
		xlab('Position (bp)')+ylab('-log10(p_values)')+theme(legend.position='none',axis.text.x=element_blank())+
		ggtitle(paste0('GWAS for ',population,' and ',phenotype_id,' with GRM on ',grm_id,' and markers as ',gwas_id))
g62<-ggplot()+
		geom_point(aes(x=ps,y=-qvalue,col=as.factor(chr)),data=gwas)+
		geom_hline(yintercept=-0.1,col='red',lty=3)+
		facet_grid(.~chr,scales='free',space='free')+
		theme_bw()+ylim(-max(gwas$qvalue),0)+
		scale_colour_manual(values=col_chr)+
		xlab('Position (bp)')+ylab('-qvalues')+theme(legend.position='none',axis.text.x=element_blank())
g63<-ggplot()+
		geom_point(aes(x=ps,y=-svalue,col=as.factor(chr)),data=gwas)+
		geom_hline(yintercept=-0.1,col='red',lty=3)+
		facet_grid(.~chr,scales='free',space='free')+
		theme_bw()+ylim(-max(gwas$svalue),0)+
		scale_colour_manual(values=col_chr)+
		xlab('Position (bp)')+ylab('-svalues')+theme(legend.position='none',axis.text.x=element_blank())
if(max(gwas$gamma)<0.1){mgamma<-0.12}else{mgamma<-max(gwas$gamma)+0.2}
g64<-ggplot()+
		geom_point(aes(x=ps,y=gamma,col=as.factor(chr)),data=gwas)+
		geom_hline(yintercept=0.1,col='red',lty=3)+
		facet_grid(.~chr,scales='free',space='free')+
		theme_bw()+
		ylim(0,mgamma)+
		scale_colour_manual(values=col_chr)+
		xlab('Position (bp)')+ylab('gamma')+theme(legend.position='none',axis.text.x=element_blank())
png(paste0(dir_out,'/manhattan_',prfx,'.png'),width=2000,height=500)
grid.arrange(g61,g62,g63,g64,nrow=4)
dev.off()
}    

