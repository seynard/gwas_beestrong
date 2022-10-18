#module load system/R-4.1.1_gcc-9.3.0
#R
library(data.table)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(VennDiagram)
library(genetics)
library(biomartr)
library(cowplot)
library(ggforce)
library(igraph)

args<-commandArgs(TRUE)
#dirin<-args[1]
#dirout<-args[2]
#pheno<-args[3]
#pheno<-unlist(strsplit(pheno,','))
#type<-args[4]
#type<-unlist(strsplit(type,','))
dirin<-'data'
dirout<-'results'
type<-c('egs','freq')
pheno<-c('varroaphoretic','varroadepthmitoracine4','logit_varroainfestation','chargevarroaracine4','pc1','eb_smr','logit_taux_reop_inf')
pop<-c('Mellifera','Ligustica_Carnica','hybrid','Ligustica_nus')

df_sign<-data.frame()
for(i in 1:length(pop)){
	print(pop[i])
	for(j in 1:length(pheno)){
		print(pheno[j])
		for(k in 1:length(type)){
df<-fread(paste0('results/summary_gwas_',pop[i],'_',pheno[j],'_freq_',type[k],'.txt.bz2'),data.table=F)
if(pop[i]=='Mellifera'){POP<-'Mellifera'
}else if (pop[i]=='Ligustica_Carnica'){POP<-'Ligustica & Carnica'
}else if (pop[i]=='hybrid'){POP<-'Hybrids'
}else if (pop[i]=='Ligustica_nus'){POP<-'Ligustica & Carnica (no us)'}
if(pheno[j]=='varroaphoretic'){PHENO<-'v_pho'
}else if (pheno[j]=='varroadepthmitoracine4'){PHENO<-'v_mito'
}else if (pheno[j]=='logit_varroainfestation'){PHENO<-'v_inf'
}else if (pheno[j]=='chargevarroaracine4'){PHENO<-'v_load'
}else if (pheno[j]=='pc1'){PHENO<-'pc1_varroa_inf'
}else if (pheno[j]=='eb_smr'){PHENO<-'MNR'
}else if (pheno[j]=='logit_taux_reop_inf'){PHENO<-'recap_inf'}
chrsize<-df%>%group_by(chr)%>%summarize(m=max(ps))
g1<-ggplot()+
	geom_point(aes(x=ps,y=-log10(p_wald)),alpha=0.4,data=df)+
	geom_point(aes(x=ps,y=-log10(p_wald)),col='red',size=3,data=df[df$lfdr<=0.1 & df$lfsr<=0.1,])+
	facet_grid(.~chr,scales='free_x',space='free_x')+
	theme_bw()+
	ylab('-log(p-value)')+
	scale_x_continuous(breaks=seq(0,chrsize$m[1],5000000))+
#	theme(axis.text.x=element_blank(),axis.title.x=element_blank())+
	theme(axis.text.x=element_text(angle=45,hjust=1))+	
	ggtitle(paste0('GWAS for group ',POP,' and phenotype ',PHENO))
#g2<-ggplot(df)+
#	geom_point(aes(x=ps,y=lfdr),alpha=0.4)+
#	facet_grid(.~chr,scales='free_x',space='free_x')+
#	geom_hline(aes(yintercept=0.1),col='red',lty=3)+
#	theme_bw()+
#	xlab('position in bp')+
#	ylab('LFDR')+
#	ylim(0,max(df$lfdr))+
#	scale_x_continuous(breaks=seq(0,chrsize$m[1],5000000))+
#	theme(axis.text.x=element_text(angle=45,hjust=1))	
png(paste0('results/plot_gwas_',POP,'_',PHENO,'_',type[k],'.png'),width=5200,height=800,res=400)
print(g1)
#print(plot_grid(g1,g2,nrow=2,ncol=1,align='v',axis='lr'))
dev.off()
df_sign_i<-subset(df,df$lfdr<=0.1 & df$lfsr<=0.1)
if(nrow(df_sign_i)>0){
df_sign_i$sp<-POP
df_sign<-rbind(df_sign,df_sign_i)}
			}
		}
	}
write.table(df_sign,'sign_gwas_ind.txt',col.names=T,row.names=F,quote=F)
