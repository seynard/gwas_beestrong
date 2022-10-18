library(data.table)
library(ggplot2)
library(gridExtra)

mantra<-fread('test_mantra/mantra.out')
colnames(mantra)<-c('rs','CHROM','POS','A1','A0','nb_study','log10_BF','post_proba','nb_ind','sign')
gemma_mel<-fread('results/gemma_Mellifera_eb_smr_lmm_bimbam.assoc.txt')
gemma_mel<-gemma_mel[,c('rs','p_wald')]
gemma_mel$p_wald<--log10(gemma_mel$mel$p_wald)
colnames(gemma_mel)<-c('rs','mel')
gemma_lig<-fread('results/gemma_Ligustica_Carnica_eb_smr_lmm_bimbam.assoc.txt')
gemma_lig<-gemma_lig[,c('rs','p_wald')]
gemma_lig$p_wald<--log10(gemma_lig$p_wald)
colnames(gemma_lig)<-c('rs','lig')
gemma_us<-fread('results/gemma_us_eb_smr_lmm_bimbam.assoc.txt')
gemma_us<-gemma_us[,c('rs','p_wald')]
gemma_us$p_wald<--log10(gemma_us$p_wald)
colnames(gemma_us)<-c('rs','us')
gemma_hybrid<-fread('results/gemma_hybrid_eb_smr_lmm_bimbam.assoc.txt')
gemma_hybrid<-gemma_hybrid[,c('rs','p_wald')]
gemma_hybrid$p_wald<--log10(gemma_hybrid$p_wald)
colnames(gemma_hybrid)<-c('rs','hyb')
X<-Reduce(function(x,y) merge(x,y,by='rs',all=T),list(gemma_mel,gemma_lig,gemma_us,gemma_hybrid))
X<-merge(mantra,X,by='rs',all=T)
X<-as.data.frame(X)
Y<-X[order(-X$log10_BF),]
y<-Y[1:20,]

g1<-ggplot()+
	geom_point(aes(y=log10_BF,x=POS),alpha=0.3,data=Y)+
	geom_vline(aes(xintercept=POS),data=y,col='red')+
	facet_grid(.~CHROM,scales='free_x')
g2<-ggplot()+
	geom_point(aes(y=mel,x=POS),alpha=0.3,data=Y)+
	geom_vline(aes(xintercept=POS),data=y,col='red')+
	facet_grid(.~CHROM,scales='free_x')
g3<-ggplot()+
	geom_point(aes(y=lig,x=POS),alpha=0.3,data=Y)+
	geom_vline(aes(xintercept=POS),data=y,col='red')+
	facet_grid(.~CHROM,scales='free_x')
g4<-ggplot()+
	geom_point(aes(y=us,x=POS),alpha=0.3,data=Y)+
	geom_vline(aes(xintercept=POS),data=y,col='red')+
	facet_grid(.~CHROM,scales='free_x')
g5<-ggplot()+
	geom_point(aes(y=hyb,x=POS),alpha=0.3,data=Y)+
	geom_vline(aes(xintercept=POS),data=y,col='red')+
	facet_grid(.~CHROM,scales='free_x')
grid.arrange(g1,g2,g3,g4,g5,ncol=1)




g2<-ggplot(df)+geom_point(aes(y=post_proba,x=POS))+facet_grid(CHROM~.,scales='free')

grid.arrange(g1,g2,ncol=2)
