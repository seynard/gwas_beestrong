library(data.table)
library(ggplot2)
library(gridExtra)
args<-commandArgs(TRUE)
dir_sonia<-args[1]
#dir_sonia<-'work/GWAS/bam_BeeStrongMOSAR'
d_min<-as.numeric(args[2])
d_max<-as.numeric(args[3])

dat<-fread(paste0(dir_sonia,'/depth_control.txt'),sep=' ',header=T,data.table=F)
head(dat)
nrow(dat)

png(paste0(dir_sonia,'/depth_control1.png'),width=1500,height=500)
par(mfrow=c(1,3))
plot(log10(dat$mean),log10(dat$median),pch=16,xlab='mean depth',ylab='median depth')
abline(v=log10(d_min),col='red')
abline(h=log10(d_min),col='red')
abline(v=log10(d_max),col='red')
abline(h=log10(d_max),col='red')
hist(log10(dat$mean),breaks=100,xlab='mean depth')
abline(v=log10(d_min),col='red')
abline(v=log10(d_max),col='red')
hist(log10(dat$median),breaks=100,ylab='median depth')
abline(v=log10(d_min),col='red')
abline(v=log10(d_max),col='red')
dev.off()

dat_f<-subset(dat,dat$mean<d_max & dat$median<d_max & dat$mean>d_min  & dat$median>d_min)
nrow(dat_f)
png(paste0(dir_sonia,'/depth_control2.png'),width=1000,height=1000)
plot(dat_f$mean,dat_f$median,pch=16,xlab='mean depth',ylab='median depth')
dev.off()
to_remove<-subset(dat[,c('CHROM','POS')],!(paste0(dat$CHROM,'_',dat$POS) %in% paste0(dat_f$CHROM,'_',dat_f$POS)))
write.table(to_remove,paste0(dir_sonia,'/to_remove_depth.txt'),sep=' ',col.names=F,row.names=F,quote=F)

chr<-unique(dat_f$CHROM)
for(i in 1:length(chr)){
c<-chr[i]
dc<-subset(dat_f,dat_f$CHROM==c)
di<-subset(dat,dat$CHROM==c)
sc<-seq(1,max(dc$POS),by=1000000)
snpc<-vector()
snpi<-vector()
for(j in 1:length(sc)){
	h<-sc[j]
	H<-sc[j+1]
	snpc<-c(snpc,nrow(dc[dc$POS<H & dc$POS>=h,]))
	snpi<-c(snpi,nrow(di[di$POS<H & di$POS>=h,]))
	}
snpd<-data.frame(sc,snpc,snpi)
png(paste0(dir_sonia,'/depth_control',chr[i],'.png'),width=1000,height=1800)
p1<-ggplot()+geom_step(data=snpd,aes(x=sc,y=snpi),color='black')+ylim(0,max(snpd$snpc,snpd$snpi))+xlab('# snp initial')+ylab('position (bp)')
p2<-ggplot()+geom_step(data=snpd,aes(x=sc,y=snpc),color='red')+ylim(0,max(snpd$snpc,snpd$snpi))+xlab('# snp filtered')+ylab('position (bp)')
p3<-ggplot()+geom_point(data=di,aes(x=POS,y=mean),color='black')+geom_point(data=dc,aes(x=POS,y=mean),color='red')+xlab('mean depth')+ylab('position (bp)')+geom_hline(yintercept=50,col='grey',linetype=2)+geom_hline(yintercept=10,col='grey',linetype=2)
grid.arrange(p1,p2,p3,ncol=1,nrow=3,widths=1,heights=c(1,1,1))	
dev.off()
}



