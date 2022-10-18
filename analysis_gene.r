#module load system/R-4.1.1_gcc-9.3.0
#R 
library(data.table)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(mashr)
library(ashr)
library(rmeta)
library(patchwork)
library(gridExtra)
library(ggvenn)
library(gprofiler2)

args<-commandArgs(TRUE)

annot<-fread(paste0('data/proteins_48_403979.csv'))
annot$chr[annot[,1]=='linkage group LG1']<-1
annot$chr[annot[,1]=='linkage group LG2']<-2
annot$chr[annot[,1]=='linkage group LG3']<-3
annot$chr[annot[,1]=='linkage group LG4']<-4
annot$chr[annot[,1]=='linkage group LG5']<-5
annot$chr[annot[,1]=='linkage group LG6']<-6
annot$chr[annot[,1]=='linkage group LG7']<-7
annot$chr[annot[,1]=='linkage group LG8']<-8
annot$chr[annot[,1]=='linkage group LG9']<-9
annot$chr[annot[,1]=='linkage group LG10']<-10
annot$chr[annot[,1]=='linkage group LG11']<-11
annot$chr[annot[,1]=='linkage group LG12']<-12
annot$chr[annot[,1]=='linkage group LG13']<-13
annot$chr[annot[,1]=='linkage group LG14']<-14
annot$chr[annot[,1]=='linkage group LG15']<-15
annot$chr[annot[,1]=='linkage group LG16']<-16
annot<-unique(annot)
annot$start[annot$Strand=='+']<-annot$Start[annot$Strand=='+']
annot$stop[annot$Strand=='+']<-annot$Stop[annot$Strand=='+']
annot$start[annot$Strand=='-']<-annot$Stop[annot$Strand=='-']
annot$stop[annot$Strand=='-']<-annot$Start[annot$Strand=='-']

vep_ind<-fread('ind_vep.txt',data.table=F)
vep_ind$type<-'ind'
a<-colsplit(vep_ind$`#Uploaded_variation`,'_',c('rs','sp'))
vep_ind$rs<-a$rs
vep_ind$sp<-a$sp
vep_ind<-vep_ind[vep_ind$rs%in%c("1:10080627","13:9483955","15:8485332","2:12025610","2:12025647","2:2729874","3:1059430","4:11665460","4:4327611"),]
vep_meta<-fread('meta2_vep.txt',data.table=F)
vep_meta$type<-'meta'
a<-colsplit(vep_meta$`#Uploaded_variation`,'_',c('rs','sp'))
vep_meta$rs<-a$rs
vep_meta$sp<-'meta'
vep_meta<-vep_meta[vep_meta$rs%in%c("1:2891204","1:7448807","1:7448811","1:12552506","1:15280956","1:16327085","1:20960056","1:21374478","1:24201224","1:25184394","2:2729874","2:4437645","2:8350714","2:9347464","2:11333369","2:12025610","2:12025647","2:16060868","3:1059430","3:6206342","3:12973246","3:12973248","4:7321246","4:7321247","4:10789077","4:11665460","5:75369","5:2008472","5:2365495","5:6736534","5:6761414","5:8737386","5:9190579","6:10450971","7:5762037","7:5772089","7:6738985","7:7028040","7:7051965","7:7078376","7:8466948","7:11806658","8:1150346","8:1319815","8:1551638","8:1748209","8:2468335","8:9557205","8:9799408","9:11564671","10:2026877","10:4266011","10:5359169","10:5359173","10:10400687","11:9369229","11:9527267","11:14369154","12:136634","12:10153855","12:10734707","14:3782741","14:6686131","14:8481541","14:8481589","15:2021142","15:2081876","15:2081914","15:4853529","15:8485332","16:1812909","16:5024160"),]
vep<-rbind(vep_ind,vep_meta)
vep$gene_name<-NA
vep$gene_name[vep$SYMBOL!='-']<-vep$SYMBOL[vep$SYMBOL!='-']
vep$gene_name[vep$SYMBOL=='-']<-vep$Gene[vep$SYMBOL=='-']
vep$DISTANCE[vep$DISTANCE=='-']<-0
vep$DISTANCE<-as.numeric(vep$DISTANCE)
vep[,c('SYMBOL','Gene','Feature','EXON','INTRON','HGVSc','HGVSp','cDNA_position','CDS_position','Protein_position','Amino_acids','Existing_variation','FLAGS','SYMBOL_SOURCE','HGNC_ID','CANONICAL','ENSP','SWISSPROT','TREMBL','UNIPARC','UNIPROT_ISOFORM','DOMAINS','HGVS_OFFSET','CLIN_SIG','SOMATIC','PHENO')]<-NULL
vep<-unique(vep)
vep<-subset(vep,vep$sp!='lig')
VEP<-vep%>%group_by(type,rs,sp)%>%slice(which.min(DISTANCE))
VEP<-as.data.frame(VEP)
VEP$pheno[VEP$rs%in%c("4:11665460","1:10080627") & VEP$sp!='meta']<-'pc1'
VEP$pheno[VEP$rs%in%c("13:9483955","2:12025610","2:12025647","15:8485332","2:2729874","3:1059430","4:4327611") & VEP$sp!='meta']<-'reop'
VEP$pheno[VEP$rs%in%c("1:20960056","1:25184394","11:9369229","12:10734707","4:11665460","5:75369","5:9190579","6:10450971","7:11806658","7:5762037","7:5772089","7:6738985","8:2468335","8:9799408") & VEP$sp=='meta']<-'pc1'
VEP$pheno[VEP$rs%in%c("1:16327085","1:21374478","1:24201224","1:2891204","10:10400687","10:4266011","10:5359169","10:5359173","11:9527267","12:10153855","12:136634","15:4853529","2:11333369","2:4437645","2:9347464","3:6206342","5:2008472","5:2365495","8:1150346","8:9557205") & VEP$sp=='meta']<-'eb_smr'
VEP$pheno[VEP$rs%in%c("1:12552506","1:15280956","1:7448807","1:7448811","10:2026877","11:14369154","14:3782741","14:6686131","14:8481541","14:8481589","15:2021142","15:2081876","15:2081914","15:8485332","16:1812909","16:5024160","2:12025610","2:12025647","2:16060868","2:2729874","2:8350714","3:1059430","3:12973246","3:12973248","4:10789077","4:7321246","4:7321247","5:6736534","5:6761414","5:8737386","7:7028040","7:7051965","7:7078376","7:8466948","8:1319815","8:1551638","8:1748209","9:11564671") & VEP$sp=='meta']<-'reop'

pheno<-c('pc1','eb_smr','reop')
version<-2
sp<-c('hyb','lignus','mel')
MBAT1<-list()
for(p in 1:length(pheno)){
	print(pheno[p])
	MBAT<-list()
	for(j in 1:length(sp)){
	mbat<-fread(paste0('mbat/',pheno[p],'_out_',sp[j],'.gene.assoc.mbat'),data.table=F)
	mbat$sp<-sp[j]
	mbat$pheno<-pheno[p]
	mbat$version<-version
	MBAT[[j]]<-mbat	
	}
	MBAT<-do.call(rbind,MBAT)
	MBAT$sp[MBAT$sp=='hyb']<-'Hybrids'
	MBAT$sp[MBAT$sp=='mel']<-'Mellifera'
	MBAT$sp[MBAT$sp%in%c('lig','lignus')]<-'Ligustica \n& Carnica'
	MBAT$sp<-factor(MBAT$sp,level=c('Ligustica \n& Carnica','Mellifera','Hybrids'))
	MBAT1[[p]]<-MBAT
}
MBAT<-do.call(rbind,MBAT1)
MBAT$version<-NULL
MBAT<-subset(MBAT,-log10(MBAT$P_mBATcombo)>5)

MASH_MBAT<-fread(paste0('mbat/out_mash_mbat.txt'),data.table=F)
MASH_MBAT<-MASH_MBAT[MASH_MBAT$version==version,]
MASH_MBAT$version<-NULL
MASH_MBAT<-subset(MASH_MBAT,MASH_MBAT$n_sign>0 & MASH_MBAT$log10_bf>4)


gene_vep_ind_pc1<-unique(VEP$gene_name[VEP$type=='ind' & VEP$pheno=='pc1'])
gene_vep_ind_eb_smr<-unique(VEP$gene_name[VEP$type=='ind' & VEP$pheno=='eb_smr'])
gene_vep_ind_reop<-unique(VEP$gene_name[VEP$type=='ind' & VEP$pheno=='reop'])
A<-list('pc1'=gene_vep_ind_pc1,'MNR'=gene_vep_ind_eb_smr,'reop'=gene_vep_ind_reop)
ggvenn(A)
gene_vep_meta_pc1<-unique(VEP$gene_name[VEP$type=='meta' & VEP$pheno=='pc1'])
gene_vep_meta_eb_smr<-unique(VEP$gene_name[VEP$type=='meta' & VEP$pheno=='eb_smr'])
gene_vep_meta_reop<-unique(VEP$gene_name[VEP$type=='meta' & VEP$pheno=='reop'])
B<-list('pc1'=gene_vep_meta_pc1,'MNR'=gene_vep_meta_eb_smr,'reop'=gene_vep_meta_reop)
ggvenn(B)
C<-list('ind'=gene_vep_ind_pc1,'meta'=gene_vep_meta_pc1)
ggvenn(C)
D<-list('ind'=gene_vep_ind_eb_smr,'meta'=gene_vep_meta_eb_smr)
ggvenn(D)
E<-list('ind'=gene_vep_ind_reop,'meta'=gene_vep_meta_reop)
ggvenn(E)

gene_meta_ind_pc1<-unique(MBAT$Gene[MBAT$pheno=='pc1'])
gene_meta_ind_eb_smr<-unique(MBAT$Gene[MBAT$pheno=='eb_smr'])
gene_meta_ind_reop<-unique(MBAT$Gene[MBAT$pheno=='reop'])
A<-list('pc1'=gene_meta_ind_pc1,'MNR'=gene_meta_ind_eb_smr,'reop'=gene_meta_ind_reop)
ggvenn(A)
gene_meta_meta_pc1<-unique(MASH_MBAT$Gene[MASH_MBAT$pheno=='pc1'])
gene_meta_meta_eb_smr<-unique(MASH_MBAT$Gene[MASH_MBAT$pheno=='eb_smr'])
gene_meta_meta_reop<-unique(MASH_MBAT$Gene[MASH_MBAT$pheno=='reop'])
B<-list('pc1'=gene_meta_meta_pc1,'MNR'=gene_meta_meta_eb_smr,'reop'=gene_meta_meta_reop)
ggvenn(B)
C<-list('ind'=gene_meta_ind_pc1,'meta'=gene_meta_meta_pc1)
ggvenn(C)
D<-list('ind'=gene_meta_ind_eb_smr,'meta'=gene_meta_meta_eb_smr)
ggvenn(D)
E<-list('ind'=gene_meta_ind_reop,'meta'=gene_meta_meta_reop)
ggvenn(E)

F<-list('vep_ind'=gene_vep_ind_pc1,'vep_meta'=gene_vep_meta_pc1,'gene_ind'=gene_meta_ind_pc1,'gene_meta'=gene_meta_meta_pc1)
ggvenn(F)
G<-list('vep_ind'=gene_vep_ind_eb_smr,'vep_meta'=gene_vep_meta_eb_smr,'gene_ind'=gene_meta_ind_eb_smr,'gene_meta'=gene_meta_meta_eb_smr)
ggvenn(G)
H<-list('vep_ind'=gene_vep_ind_reop,'vep_meta'=gene_vep_meta_reop,'gene_ind'=gene_meta_ind_reop,'gene_meta'=gene_meta_meta_reop)
ggvenn(H)


X<-VEP
X<-X[,c('sp','gene_name','pheno')]
colnames(X)[colnames(X)=='gene_name']<-'Gene'
X$analysis<-'VEP'
Y<-MBAT
Y<-Y[,c('sp','Gene','pheno')]
Y$analysis<-'MBAT'
Z<-MASH_MBAT
Z<-Z[,c('Gene','pheno')]
Z$sp<-'meta'
Z$analysis<-'MBAT'

df<-do.call(rbind,list(X,Y,Z))
df$sp[df$sp=='lig_nus']<-'Ligustica_Carnica'
df$sp[df$sp=='Ligustica \n& Carnica']<-'Ligustica_Carnica'
df$sp[df$sp=='hyb']<-'Hybrids'
df$sp[df$sp=='mel']<-'Mellifera'
df<-unique(df)
write.table(df,file='gene_analysis.txt',quote=FALSE,row.names=FALSE,col.names=TRUE,sep=';')

DF<-list()
for(i in 1:length(unique(df$Gene))){
	dfi<-subset(df,df$Gene==unique(df$Gene)[i])
	if('VEP'%in%unique(dfi$analysis) & 'MBAT'%in%unique(dfi$analysis)){
	analysis<-'VEP_MBAT'
	pheno_vep<-paste0(unique(dfi$pheno[dfi$analysis=='VEP']),collapse='_')
	sp_vep<-paste0(unique(dfi$sp[dfi$analysis=='VEP']),collapse='_')
	pheno_mbat<-paste0(unique(dfi$pheno[dfi$analysis=='MBAT']),collapse='_')
	sp_mbat<-paste0(unique(dfi$sp[dfi$analysis=='MBAT']),collapse='_')
	}else if('VEP'%in%unique(dfi$analysis) & !'MBAT'%in%unique(dfi$analysis)){
	analysis<-'VEP'
	pheno_vep<-paste0(unique(dfi$pheno[dfi$analysis=='VEP']),collapse='_')
	sp_vep<-paste0(unique(dfi$sp[dfi$analysis=='VEP']),collapse='_')
	pheno_mbat<-NA
	sp_mbat<-NA
	}else if(!'VEP'%in%unique(dfi$analysis) & 'MBAT'%in%unique(dfi$analysis)){
	analysis<-'MBAT'
	pheno_vep<-NA
	sp_vep<-NA
	pheno_mbat<-paste0(unique(dfi$pheno[dfi$analysis=='MBAT']),collapse='_')
	sp_mbat<-paste0(unique(dfi$sp[dfi$analysis=='MBAT']),collapse='_')
	}
	dfi<-data.frame('Gene'=unique(df$Gene)[i],'analysis'=analysis,'pheno_vep'=pheno_vep,'sp_vep'=sp_vep,'pheno_mbat'=pheno_mbat,'sp_mbat'=sp_mbat)
	DF[[i]]<-dfi	
	}
DF<-do.call(rbind,DF)

DF[DF$analysis=='VEP_MBAT',]
DF[DF$analysis=='VEP',]
DF[DF$analysis=='MBAT' & DF$sp_mbat%in%c("Hybrids_meta","Ligustica_Carnica_meta","Mellifera_meta","Hybrids_Ligustica_Carnica_meta","Ligustica_Carnica_Mellifera_meta"),]
DF[DF$analysis=='MBAT' & DF$pheno_mbat%in%c("eb_smr_pc1_reop","pc1_eb_smr_reop","eb_smr_reop_pc1"),]

gores<-gost(query=DF$Gene,organism="amellifera")
GO<-gores$result
ggplot(GO)+
  geom_point(aes(x=p_value,y=term_name,size=(precision*100)))+
  theme_bw()+
  scale_size_continuous(name="prop gene \nannot to function")
  
gores$result$term_name




