library(data.table)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(gridExtra)
dirin<-'data'
dirout<-'results'

ldak_bed<-fread('results/varroadepthmito_Mellifera_bed.assoc')
gemma_bed<-fread('results/gemma_Mellifera_varroadepthmito_lmm_bed.assoc.txt')
colnames(ldak_bed)[colnames(ldak_bed)=='Predictor']<-'rs'
BED<-merge(gemma_bed,ldak_bed,by='rs',all=T)
BED[,c('n_miss','A1_Mean','Effect_Liability','SD_Liability','logl_H1','l_remle','Wald_Stat','chr','ps','Chromosome','Basepair')]<-NULL
a<-colsplit(BED$rs,':',c('chr','ps'))
BED$chr<-a$chr
BED$ps<-a$ps
maf<-vector()
for(i in 1:nrow(BED)){
	if(is.na(BED$af[i])){maf[i]<-NA		
	}else if(BED$af[i]<=0.5){maf[i]<-BED$af[i]
	}else{maf[i]<-1-BED$af[i]}
	}
BED$maf<-maf
BED[,c('allele1','allele0','A1','A2')]<-NULL
BED<-gather(BED,type,value_bed,c('af','beta','se','p_wald','Wald_P','Effect','SD','MAF','maf'))

ldak_bimbam<-fread('results/varroadepthmito_Mellifera_bimbam.assoc')
gemma_bimbam<-fread('results/gemma_Mellifera_varroadepthmito_lmm_bimbam.assoc.txt')
colnames(ldak_bimbam)[colnames(ldak_bimbam)=='Predictor']<-'rs'
BIMBAM<-merge(gemma_bimbam,ldak_bimbam,by='rs',all=T)
BIMBAM[,c('n_miss','A1_Mean','Effect_Liability','SD_Liability','logl_H1','l_remle','Wald_Stat','chr','ps','Chromosome','Basepair')]<-NULL
a<-colsplit(BIMBAM$rs,':',c('chr','ps'))
BIMBAM$chr<-a$chr
BIMBAM$ps<-a$ps
maf<-vector()
for(i in 1:nrow(BIMBAM)){
	if(is.na(BIMBAM$af[i])){maf[i]<-NA		
	}else if(BIMBAM$af[i]<=0.5){maf[i]<-BIMBAM$af[i]
	}else{maf[i]<-1-BIMBAM$af[i]}
	}
BIMBAM$maf<-maf
BIMBAM[,c('allele1','allele0','A1','A2')]<-NULL
BIMBAM<-gather(BIMBAM,type,value_bimbam,c('af','beta','se','p_wald','Wald_P','Effect','SD','MAF','maf'))

ldak_freq<-fread('results/varroadepthmito_Mellifera_freq.assoc')
gemma_freq<-fread('results/gemma_Mellifera_varroadepthmito_lmm_freq.assoc.txt')
colnames(ldak_freq)[colnames(ldak_freq)=='Predictor']<-'rs'
FREQ<-merge(gemma_freq,ldak_freq,by='rs',all=T)
FREQ[,c('n_miss','A1_Mean','Effect_Liability','SD_Liability','logl_H1','l_remle','Wald_Stat','chr','ps','Chromosome','Basepair')]<-NULL
a<-colsplit(FREQ$rs,':',c('chr','ps'))
FREQ$chr<-a$chr
FREQ$ps<-a$ps
maf<-vector()
for(i in 1:nrow(FREQ)){
	if(is.na(FREQ$af[i])){maf[i]<-NA		
	}else if(FREQ$af[i]<=0.5){maf[i]<-FREQ$af[i]
	}else{maf[i]<-1-FREQ$af[i]}
	}
FREQ$maf<-maf
FREQ[,c('allele1','allele0','A1','A2')]<-NULL
FREQ<-gather(FREQ,type,value_freq,c('af','beta','se','p_wald','Wald_P','Effect','SD','MAF','maf'))

DF<-Reduce(function(x,y) merge(x,y,by=c('chr','ps','rs','type'),all=T),list(BED,BIMBAM,FREQ))


par(mfrow=c(1,3))
plot(DF$value_bed[DF$type=='beta'],DF$value_bimbam[DF$type=='beta'])
plot(DF$value_bed[DF$type=='beta'],DF$value_freq[DF$type=='beta'])
plot(DF$value_bimbam[DF$type=='beta'],DF$value_freq[DF$type=='beta'])

plot(DF$value_bed[DF$type=='Effect'],DF$value_bimbam[DF$type=='Effect'])
plot(DF$value_bed[DF$type=='Effect'],DF$value_freq[DF$type=='Effect'])
plot(DF$value_bimbam[DF$type=='Effect'],DF$value_freq[DF$type=='Effect'])





plot(BED$MAF,BED$maf)
plot(BED$Effect,BED$beta)
plot(BED$SD,BED$se)
plot(BED$Wald_P,BED$p_wald)

plot(BIMBAM$MAF,BIMBAM$maf)
plot(BIMBAM$Effect,BIMBAM$beta)
plot(BIMBAM$SD,BIMBAM$se)
plot(BIMBAM$Wald_P,BIMBAM$p_wald)

plot(FREQ$MAF,FREQ$maf)
plot(FREQ$Effect,FREQ$beta)
plot(FREQ$SD,FREQ$se)
plot(FREQ$Wald_P,FREQ$p_wald)



library(KRIS)
in_gemma_freq<-fread('data/input_gemma_us_freq.txt')
in_gemma_bed<-read.bed('results/geno_hom_us_smr_brut.bed','results/geno_hom_us_smr_brut.bim','results/geno_hom_us_smr_brut.fam')
in_gemma_bimbam<-fread('data/input_gemma_us.txt')
in_ldak_freq<-fread('data/freq_imputed_us_0.05.txt.filter')
in_ldak_bed<-read.bed('data/geno_hom_us_maf.bed','data/geno_hom_us_maf.bim','data/geno_hom_us_maf.fam')
in_ldak_bimbam<-fread('data/geno_us_maf.txt')

bimbam_gemma<-fread('data/input_gemma_us.txt',data.table=F)
bimbam_ldak<-fread('data/geno_us_maf.txt',data.table=F)
list_snp<-paste0(bimbam_ldak$CHROM,':',bimbam_ldak$POS)
bimbam_gemma$freq<-rowMeans(bimbam_gemma[,4:ncol(bimbam_gemma)],na.rm=T)
bimbam_ldak$freq<-rowMeans(bimbam_ldak[,5:ncol(bimbam_ldak)],na.rm=T)
summary(abs(bimbam_gemma$freq-bimbam_ldak$freq))

freq_gemma<-fread('data/input_gemma_us_freq.txt',data.table=F)
freq_ldak<-fread('data/freq_imputed_us_0.05.txt.filter',data.table=F)
freq_gemma$freq<-rowMeans(freq_gemma[,4:ncol(freq_gemma)],na.rm=T)
freq_ldak$freq<-rowMeans(freq_ldak[,5:ncol(freq_ldak)],na.rm=T)
summary(abs(freq_gemma$freq-freq_ldak$freq))











