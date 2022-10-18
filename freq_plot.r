#module load system/R-3.5.1
#R
#dir='/home/seynard/work/GWAS/bam_BeeStrongMOSAR'
args<-commandArgs(TRUE)
dir<-args[1]
library(data.table)
library(ggplot2)
library(gridExtra)
freq_admix=args[2]
n_pop<-as.numeric(unlist(regmatches(freq_admix, gregexpr("[[:digit:]]+", freq_admix))))
freq<-fread(paste0(dir,'/',freq_admix),data.table=F)
pop<-colnames(freq)[3:ncol(freq)]

col_subsp<-matrix(ncol=2,nrow=length(pop))
col_subsp[,1]<-pop
col_subsp[,2][col_subsp[,1]=='Mellifera']<-'black'
col_subsp[,2][col_subsp[,1]=='Caucasica']<-'green'
col_subsp[,2][col_subsp[,1]=='Ligustica_Carnica']<-'goldenrod2'

for(c in 1:length(unique(freq$CHROM))){
fc<-subset(freq,freq$CHROM==unique(freq$CHROM)[c])
png(paste0(dir,'/plot_freq',c,'_',n_pop,'.png'),width=1500,height=1000)
p1<-ggplot(fc)+geom_point(aes(x=POS,y=fc[,pop[1]]),colour=col_subsp[1,2],shape=20)+theme_classic()+xlab('position (bp)')+ylab('')
p2<-ggplot(fc)+geom_point(aes(x=POS,y=fc[,pop[2]]),colour=col_subsp[2,2],shape=20)+theme_classic()+xlab('position (bp)')+ylab('')
p3<-ggplot(fc)+geom_point(aes(x=POS,y=fc[,pop[3]]),colour=col_subsp[3,2],shape=20)+theme_classic()+xlab('position (bp)')+ylab('')
grid.arrange(p1,p2,p3,ncol=1,nrow=length(pop))
dev.off()
}
