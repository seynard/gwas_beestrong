gwas1<-fread('back_up/version2_oct2021/results/gemma_Ligustica_Carnica_chargevarroaracine4_freq_lmm_freq.assoc.txt')
gwas2<-fread('results/gemma_Ligustica_Carnica_chargevarroaracine4_freq_lmm_freq.assoc.txt')
gwas<-merge(gwas1,gwas2,by='rs')
a<-colsplit(gwas$rs,':',c('CHROM','POS'))
gwas$chr<-a$CHROM
gwas$ps<-a$POS

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
	lambda <- median(chisq) / qchisq(0.5, 1)
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

g1<-pp.ggplot(gwas$p_wald.x)+theme_bw()
g2<-pp.ggplot(gwas$p_wald.y)+theme_bw()

grid.arrange(g1,g2,ncol=2)

ggplot(gwas)+
	geom_point(aes(x=ps,y=-log10(p_wald.x)),col='red')+
	geom_point(aes(x=ps,y=log10(p_wald.y)),col='blue')+
	facet_grid(.~chr,scales='free', space='free')+
	theme_bw()+theme(legend.position='none')










pclf<-fread('results/pca_Ligustica_Carnica_freq.values',data.table=F)
pcle<-fread('results/pca_Ligustica_Carnica_egs.values',data.table=F)
pcmf<-fread('results/pca_Mellifera_ncorse_freq.values',data.table=F)
pcme<-fread('results/pca_Mellifera_ncorse_egs.values',data.table=F)

par(mfrow=c(2,2))
x<-quantile(pclf[,1],probs=0.975)
hist(pclf[,1],breaks=100)
abline(v=x,col='red')
pclf$k<-NA
pclf$k[pcle[,1]>x]<-'YES'
x<-quantile(pcle[,1],probs=0.975)
hist(pcle[,1],breaks=100)
abline(v=x,col='red')
pcle$k<-NA
pcle$k[pcle[,1]>x]<-'YES'
x<-quantile(pcmf[,1],probs=0.975)
hist(pcmf[,1],breaks=100)
abline(v=x,col='red')
pcmf$k<-NA
pcmf$k[pcmf[,1]>x]<-'YES'
x<-quantile(pcme[,1],probs=0.975)
hist(pcme[,1],breaks=100)
abline(v=x,col='red')
pcme$k<-NA
pcme$k[pcme[,1]>x]<-'YES'

ggplot()+
	geom_point(aes(x=seq(1:nrow(pclf)),y=(V1/sum(V1))*100,col=k),data=pclf)+
	geom_point(aes(x=seq(1:nrow(pcle)),y=(V1/sum(V1))*100,col=k),data=pcle)+
	geom_point(aes(x=seq(1:nrow(pcmf)),y=(V1/sum(V1))*100,col=k),data=pcmf)+
	geom_point(aes(x=seq(1:nrow(pcme)),y=(V1/sum(V1))*100,col=k),data=pcme)

ggplot()+
	geom_point(aes(x=seq(1:nrow(pclf)),y=cumsum((V1/sum(V1))*100),col=k),data=pclf)+
	geom_point(aes(x=seq(1:nrow(pcle)),y=cumsum((V1/sum(V1))*100),col=k),data=pcle)+
	geom_point(aes(x=seq(1:nrow(pcmf)),y=cumsum((V1/sum(V1))*100),col=k),data=pcmf)+
	geom_point(aes(x=seq(1:nrow(pcme)),y=cumsum((V1/sum(V1))*100),col=k),data=pcme)


