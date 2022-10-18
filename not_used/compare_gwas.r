library(data.table)
library(ggplot2)
library(tidyverse)
library(reshape2)
args<-commandArgs(TRUE)
dirin<-args[1]
dirout<-args[2]
#dirin<-'data'
#dirout<-'results'

gemma_bed<-intersect(intersect(list.files(path=dirout,pattern='gemma'),list.files(path=dirout,pattern='assoc')),list.files(path=dirout,pattern='bed'))
gemma_bimbam<-intersect(intersect(list.files(path=dirout,pattern='gemma'),list.files(path=dirout,pattern='assoc')),list.files(path=dirout,pattern='bimbam'))
gemma_freq<-intersect(intersect(list.files(path=dirout,pattern='gemma'),list.files(path=dirout,pattern='assoc')),list.files(path=dirout,pattern='freq'))

ldak_bed<-setdiff(intersect(list.files(path=dirout,pattern='bed'),list.files(path=dirout,pattern='assoc')),list.files(path=dirout,pattern='gemma'))
ldak_bimbam<-setdiff(intersect(list.files(path=dirout,pattern='bimbam'),list.files(path=dirout,pattern='assoc')),list.files(path=dirout,pattern='gemma'))
ldak_freq<-setdiff(intersect(list.files(path=dirout,pattern='freq'),list.files(path=dirout,pattern='assoc')),list.files(path=dirout,pattern='gemma'))

bslmm_bed<-intersect(intersect(list.files(path=dirout,pattern='gemma'),list.files(path=dirout,pattern='param')),list.files(path=dirout,pattern='bed'))
bslmm_bimbam<-intersect(intersect(list.files(path=dirout,pattern='gemma'),list.files(path=dirout,pattern='param')),list.files(path=dirout,pattern='bimbam'))
bslmm_freq<-intersect(intersect(list.files(path=dirout,pattern='gemma'),list.files(path=dirout,pattern='param')),list.files(path=dirout,pattern='freq'))

#pop_id<-c('hybrid','Ligustica_Carnica','Mellifera','us','all')
pop_id<-c('hybrid','Ligustica_Carnica','Mellifera','us')
pheno<-c('chargevarroa','eb_smr','lm_eb_smr','nbvarroa100abeilles','pc1','smr_brut','varroadepthmito','varroainfestation','varroaphoretic')

for(i in 1:length(pop_id)){
	for(j in 1:length(pheno)){
		print(pop_id[i])
		print(pheno[j])
GEMMA_BED<-fread(paste0(dirout,'/',intersect(gemma_bed[grep(pop_id[i],gemma_bed)],gemma_bed[grep(pheno[j],gemma_bed)])[1]))
GEMMA_BIMBAM<-fread(paste0(dirout,'/',intersect(gemma_bimbam[grep(pop_id[i],gemma_bimbam)],gemma_bimbam[grep(pheno[j],gemma_bimbam)])[1]))
GEMMA_FREQ<-fread(paste0(dirout,'/',intersect(gemma_freq[grep(pop_id[i],gemma_freq)],gemma_freq[grep(pheno[j],gemma_freq)])[1]))
LDAK_BED<-fread(paste0(dirout,'/',intersect(ldak_bed[grep(pop_id[i],ldak_bed)],ldak_bed[grep(pheno[j],ldak_bed)])[1]))
LDAK_BIMBAM<-fread(paste0(dirout,'/',intersect(ldak_bimbam[grep(pop_id[i],ldak_bimbam)],ldak_bimbam[grep(pheno[j],ldak_bimbam)])[1]))
LDAK_FREQ<-fread(paste0(dirout,'/',intersect(ldak_freq[grep(pop_id[i],ldak_freq)],ldak_freq[grep(pheno[j],ldak_freq)])[1]))
BSLMM_BED<-fread(paste0(dirout,'/',intersect(bslmm_bed[grep(pop_id[i],bslmm_bed)],bslmm_bed[grep(pheno[j],bslmm_bed)])[1]))
BSLMM_BIMBAM<-fread(paste0(dirout,'/',intersect(bslmm_bimbam[grep(pop_id[i],bslmm_bimbam)],bslmm_bimbam[grep(pheno[j],bslmm_bimbam)])[1]))
BSLMM_FREQ<-fread(paste0(dirout,'/',intersect(bslmm_freq[grep(pop_id[i],bslmm_freq)],bslmm_freq[grep(pheno[j],bslmm_freq)])[1]))
GEMMA_BED<-GEMMA_BED[,c('rs','p_wald')]
GEMMA_BED$p_wald<--log10(GEMMA_BED$p_wald)
colnames(GEMMA_BED)<-c('rs','p_val_gemma_bed')
GEMMA_BIMBAM<-GEMMA_BIMBAM[,c('rs','p_wald')]
GEMMA_BIMBAM$p_wald<--log10(GEMMA_BIMBAM$p_wald)
colnames(GEMMA_BIMBAM)<-c('rs','p_val_gemma_bimbam')
GEMMA_FREQ<-GEMMA_FREQ[,c('rs','p_wald')]
GEMMA_FREQ$p_wald<--log10(GEMMA_FREQ$p_wald)
colnames(GEMMA_FREQ)<-c('rs','p_val_gemma_freq')
LDAK_BED<-LDAK_BED[,c('Predictor','Wald_P')]
LDAK_BED$Wald_P<--log10(LDAK_BED$Wald_P)
colnames(LDAK_BED)<-c('rs','p_val_ldak_bed')
LDAK_BIMBAM<-LDAK_BIMBAM[,c('Predictor','Wald_P')]
LDAK_BIMBAM$Wald_P<--log10(LDAK_BIMBAM$Wald_P)
colnames(LDAK_BIMBAM)<-c('rs','p_val_ldak_bimbam')
LDAK_FREQ<-LDAK_FREQ[,c('Predictor','Wald_P')]
LDAK_FREQ$Wald_P<--log10(LDAK_FREQ$Wald_P)
colnames(LDAK_FREQ)<-c('rs','p_val_ldak_freq')
BSLMM_BED<-BSLMM_BED[,c('rs','gamma')]
colnames(BSLMM_BED)<-c('rs','p_val_bslmm_bed')
BSLMM_BIMBAM<-BSLMM_BIMBAM[,c('rs','gamma')]
colnames(BSLMM_BIMBAM)<-c('rs','p_val_bslmm_bimbam')
BSLMM_FREQ<-BSLMM_FREQ[,c('rs','gamma')]
colnames(BSLMM_FREQ)<-c('rs','p_val_bslmm_freq')
X<-Reduce(function(x,y) merge(x,y,by='rs',all=T),list(GEMMA_BED,GEMMA_BIMBAM,GEMMA_FREQ,LDAK_BED,LDAK_BIMBAM,LDAK_FREQ,BSLMM_BED,BSLMM_BIMBAM,BSLMM_FREQ))
X<-gather(X,type,p_val,colnames(X)[colnames(X)!='rs'])
a<-colsplit(X$rs,':',c('CHROM','POS'))
X$CHROM<-a$CHROM
X$POS<-a$POS
X$type<-gsub('p_val_','',X$type)

png(paste0(dirout,'/gwas_',pop_id[i],'_',pheno[j],'.png'),width=3000,height=1000)
p<-ggplot(X)+
geom_point(aes(x=POS,y=p_val))+
ggtitle(paste0('GWAS ',pop_id[i],' ',pheno[j]))+
facet_grid(type~CHROM,scales='free',space="free_x")
print(p)
dev.off()
	}
}
