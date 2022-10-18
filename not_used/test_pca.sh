#! /bin/bash
##################################################
#################### PARAMETERS #################### 
dir="/work/project/dynagen/seynard/GWAS/infestation"
script="/work/project/dynagen/seynard/GWAS/infestation/scripts"
dir_out="/work/project/dynagen/seynard/GWAS/infestation/results"
dir_in="/work/project/dynagen/seynard/GWAS/infestation/data"
dir_prior="/work/project/dynagen/seynard/GWAS/BeeStrongHAV3_1/" 
fasta=${dir_save}/Fasta/GCF_003254395.2_Amel_HAv3.1_genomic.fna
vcf_name="MetaGenotypesCalled870_raw_snps"
vcf_file=${dir_in}/${vcf_name}_allfilter.vcf 
pop_id="Mellifera,Caucasica,Ligustica_Carnica"
n_pop=$(echo $(IFS=","; set -f; set -- $pop_id; echo $#))
pop_id2="Mellifera,Mellifera_ncorse,Caucasica,Ligustica_Carnica,hybrid,us,all,corse"
pheno="Data_BS.csv"
pheno_mito='compareDepthsBeeVarroa.txt'
snp50k='HAV3_1_50000.txt'
maf_min=0.01
compo_threshold=0.8
missing_rate=0.05
ncpu=10
##################################################

################################################
#################### MODULES #################### 
module load system/R-3.6.1
module load system/pandoc-2.1.3
module load bioinfo/bcftools-1.6
module load bioinfo/gemma-v0.97
module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.15
module load bioinfo/admixture_linux-1.3.0
module load bioinfo/plink-v1.90b5.3
module load bioinfo/shapeit.v2.904
module load bioinfo/samtools-1.8
module load bioinfo/bedtools-2.27.1
module load bioinfo/blast-2.2.26
module load system/Python-3.6.3
module load system/jdk-12.0.2
module load bioinfo/interproscan-5.44-79.0
module load system/Python-3.7.4
################################################
pheno='chargevarroaracine4'

mkdir -p ${dir_out}/cov_pca

type='Mellifera_ncorse'
grm_type='freq'
grm=${dir_out}/LDAKgmat_${type}_freq.grm.raw
infile="${dir_in}/in_${type}_${grm_type}.freq"

cov_file1="${dir_out}/cov_pca/cov_${type}_${grm_type}1.cov"
awk '{print $3}' ${dir_out}/pca_${type}_${grm_type}.vect > ${cov_file1}
sbatch --mem=100G --wrap="gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -k ${grm} -c ${cov_file1} -lmm 1 -n 1 -notsnp -outdir ${dir_out}/cov_pca -o ./gemma_${type}_${pheno}_${grm_type}_lmm_freq_cov1"
cov_file2="${dir_out}/cov_pca/cov_${type}_${grm_type}2.cov"
awk '{print $3" "$4}' ${dir_out}/pca_${type}_${grm_type}.vect > ${cov_file2}
sbatch --mem=100G --wrap="gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -k ${grm} -c ${cov_file2} -lmm 1 -n 1 -notsnp -outdir ${dir_out}/cov_pca -o ./gemma_${type}_${pheno}_${grm_type}_lmm_freq_cov2"
cov_file3="${dir_out}/cov_pca/cov_${type}_${grm_type}3.cov"
awk '{print $3" "$4" "$5}' ${dir_out}/pca_${type}_${grm_type}.vect > ${cov_file3}
sbatch --mem=100G --wrap="gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -k ${grm} -c ${cov_file3} -lmm 1 -n 1 -notsnp -outdir ${dir_out}/cov_pca -o ./gemma_${type}_${pheno}_${grm_type}_lmm_freq_cov3"
cov_file5="${dir_out}/cov_pca/cov_${type}_${grm_type}5.cov"
awk '{print $3" "$4" "$5" "$6" "$7}' ${dir_out}/pca_${type}_${grm_type}.vect > ${cov_file5}
sbatch --mem=100G --wrap="gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -k ${grm} -c ${cov_file5} -lmm 1 -n 1 -notsnp -outdir ${dir_out}/cov_pca -o ./gemma_${type}_${pheno}_${grm_type}_lmm_freq_cov5"
cov_file10="${dir_out}/cov_pca/cov_${type}_${grm_type}10.cov"
awk '{print $3" "$4" "$5" "$6" "$7" "$8" "$9" "$10" "$11" "$12}' ${dir_out}/pca_${type}_${grm_type}.vect > ${cov_file10}
sbatch --mem=100G --wrap="gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -k ${grm} -c ${cov_file10} -lmm 1 -n 1 -notsnp -outdir ${dir_out}/cov_pca -o ./gemma_${type}_${pheno}_${grm_type}_lmm_freq_cov10"
cov_file20="${dir_out}/cov_pca/cov_${type}_${grm_type}20.cov"
awk '{print $3" "$4" "$5" "$6" "$7" "$8" "$9" "$10" "$11" "$12" "$13" "$14" "$15" "$16" "$17" "$18" "$19" "$20" "$21" "$22}' ${dir_out}/pca_${type}_${grm_type}.vect > ${cov_file20}
sbatch --mem=100G --wrap="gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -k ${grm} -c ${cov_file20} -lmm 1 -n 1 -notsnp -outdir ${dir_out}/cov_pca -o ./gemma_${type}_${pheno}_${grm_type}_lmm_freq_cov20"
sbatch --mem=100G --wrap="gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -k ${grm} -lmm 1 -n 1 -notsnp -outdir ${dir_out}/cov_pca -o ./gemma_${type}_${pheno}_${grm_type}_lmm_freq"
sbatch --mem=100G --wrap="
gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -k ${grm} -eigen -outdir ${dir_out}/cov_pca -o ./${type}_${pheno}_${grm_type}_pca
gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -d ${dir_out}/cov_pca/${type}_${pheno}_${grm_type}_pca.eigenD.txt -u ${dir_out}/cov_pca/${type}_${pheno}_${grm_type}_pca.eigenU.txt -lmm 1 -n 1 -notsnp -outdir ${dir_out}/cov_pca -o ./gemma_${type}_${pheno}_${grm_type}_lmm_pca
"

type='Ligustica_Carnica'
grm_type='freq'
grm=${dir_out}/LDAKgmat_${type}_freq.grm.raw
infile="${dir_in}/in_${type}_${grm_type}.freq"

cov_file1="${dir_out}/cov_pca/cov_${type}_${grm_type}1.cov"
awk '{print $3}' ${dir_out}/pca_${type}_${grm_type}.vect > ${cov_file1}
sbatch --mem=100G --wrap="gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -k ${grm} -c ${cov_file1} -lmm 1 -n 1 -notsnp -outdir ${dir_out}/cov_pca -o ./gemma_${type}_${pheno}_${grm_type}_lmm_freq_cov1"
cov_file2="${dir_out}/cov_pca/cov_${type}_${grm_type}2.cov"
awk '{print $3" "$4}' ${dir_out}/pca_${type}_${grm_type}.vect > ${cov_file2}
sbatch --mem=100G --wrap="gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -k ${grm} -c ${cov_file2} -lmm 1 -n 1 -notsnp -outdir ${dir_out}/cov_pca -o ./gemma_${type}_${pheno}_${grm_type}_lmm_freq_cov2"
cov_file3="${dir_out}/cov_pca/cov_${type}_${grm_type}3.cov"
awk '{print $3" "$4" "$5}' ${dir_out}/pca_${type}_${grm_type}.vect > ${cov_file3}
sbatch --mem=100G --wrap="gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -k ${grm} -c ${cov_file3} -lmm 1 -n 1 -notsnp -outdir ${dir_out}/cov_pca -o ./gemma_${type}_${pheno}_${grm_type}_lmm_freq_cov3"
cov_file5="${dir_out}/cov_pca/cov_${type}_${grm_type}5.cov"
awk '{print $3" "$4" "$5" "$6" "$7}' ${dir_out}/pca_${type}_${grm_type}.vect > ${cov_file5}
sbatch --mem=100G --wrap="gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -k ${grm} -c ${cov_file5} -lmm 1 -n 1 -notsnp -outdir ${dir_out}/cov_pca -o ./gemma_${type}_${pheno}_${grm_type}_lmm_freq_cov5"
cov_file10="${dir_out}/cov_pca/cov_${type}_${grm_type}10.cov"
awk '{print $3" "$4" "$5" "$6" "$7" "$8" "$9" "$10" "$11" "$12}' ${dir_out}/pca_${type}_${grm_type}.vect > ${cov_file10}
sbatch --mem=100G --wrap="gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -k ${grm} -c ${cov_file10} -lmm 1 -n 1 -notsnp -outdir ${dir_out}/cov_pca -o ./gemma_${type}_${pheno}_${grm_type}_lmm_freq_cov10"
cov_file20="${dir_out}/cov_pca/cov_${type}_${grm_type}20.cov"
awk '{print $3" "$4" "$5" "$6" "$7" "$8" "$9" "$10" "$11" "$12" "$13" "$14" "$15" "$16" "$17" "$18" "$19" "$20" "$21" "$22}' ${dir_out}/pca_${type}_${grm_type}.vect > ${cov_file20}
sbatch --mem=100G --wrap="gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -k ${grm} -c ${cov_file20} -lmm 1 -n 1 -notsnp -outdir ${dir_out}/cov_pca -o ./gemma_${type}_${pheno}_${grm_type}_lmm_freq_cov20"
sbatch --mem=100G --wrap="gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -k ${grm} -lmm 1 -n 1 -notsnp -outdir ${dir_out}/cov_pca -o ./gemma_${type}_${pheno}_${grm_type}_lmm_freq"

sbatch --mem=100G --wrap="
gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -k ${grm} -eigen -outdir ${dir_out}/cov_pca -o ./${type}_${pheno}_${grm_type}_pca
gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -d ${dir_out}/cov_pca/${type}_${pheno}_${grm_type}_pca.eigenD.txt -u ${dir_out}/cov_pca/${type}_${pheno}_${grm_type}_pca.eigenU.txt -lmm 1 -n 1 -notsnp -outdir ${dir_out}/cov_pca -o ./gemma_${type}_${pheno}_${grm_type}_lmm_pca
"











module load bioinfo/EIG-7.2.1





R 
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
gwas_m0<-fread('results/cov_pca/gemma_Mellifera_ncorse_chargevarroaracine4_freq_lmm_freq.assoc.txt')
gwas_m0<-gwas_m0[,c('rs','p_wald')]
colnames(gwas_m0)[2]<-'p0'
gwas_m1<-fread('results/cov_pca/gemma_Mellifera_ncorse_chargevarroaracine4_freq_lmm_freq_cov1.assoc.txt')
gwas_m1<-gwas_m1[,c('rs','p_wald')]
colnames(gwas_m1)[2]<-'p1'
gwas_m2<-fread('results/cov_pca/gemma_Mellifera_ncorse_chargevarroaracine4_freq_lmm_freq_cov2.assoc.txt')
gwas_m2<-gwas_m2[,c('rs','p_wald')]
colnames(gwas_m2)[2]<-'p2'
gwas_m3<-fread('results/cov_pca/gemma_Mellifera_ncorse_chargevarroaracine4_freq_lmm_freq_cov3.assoc.txt')
gwas_m3<-gwas_m3[,c('rs','p_wald')]
colnames(gwas_m3)[2]<-'p3'
gwas_m5<-fread('results/cov_pca/gemma_Mellifera_ncorse_chargevarroaracine4_freq_lmm_freq_cov5.assoc.txt')
gwas_m5<-gwas_m5[,c('rs','p_wald')]
colnames(gwas_m5)[2]<-'p5'
gwas_m10<-fread('results/cov_pca/gemma_Mellifera_ncorse_chargevarroaracine4_freq_lmm_freq_cov10.assoc.txt')
gwas_m10<-gwas_m10[,c('rs','p_wald')]
colnames(gwas_m10)[2]<-'p10'
gwas_m20<-fread('results/cov_pca/gemma_Mellifera_ncorse_chargevarroaracine4_freq_lmm_freq_cov20.assoc.txt')
gwas_m20<-gwas_m20[,c('rs','p_wald')]
colnames(gwas_m20)[2]<-'p20'

pc_m<-fread('results/pca_Mellifera_ncorse_freq.values')

library('AssocTests')
tw=tw(pc_m,nrw(pc_m), criticalpoint=0.9793)






gwas_m<-Reduce(function (x,y) merge(x,y,by='rs',all=T),list(gwas_m0,gwas_m1,gwas_m2,gwas_m3,gwas_m5,gwas_m10))
#gwas_m<-Reduce(function (x,y) merge(x,y,by='rs',all=T),list(gwas_m0,gwas_m1,gwas_m2,gwas_m3,gwas_m5,gwas_m10,gwas_m20))
a<-colsplit(gwas_m$rs,':',c('CHROM','POS'))
gwas_m$chr<-a$CHROM
gwas_m$ps<-a$POS
pv0=as.data.frame(gwas_m%>%select(p0))
pv1=as.data.frame(gwas_m%>%select(p1))
pv2=as.data.frame(gwas_m%>%select(p2))
pv3=as.data.frame(gwas_m%>%select(p3))
pv5=as.data.frame(gwas_m%>%select(p5))
pv10=as.data.frame(gwas_m%>%select(p10))
#pv20=as.data.frame(gwas%>%select(p20))

g0<-pp.ggplot(pv0[,1])+theme_bw()
g1<-pp.ggplot(pv1[,1])+theme_bw()
g2<-pp.ggplot(pv2[,1])+theme_bw()
g3<-pp.ggplot(pv3[,1])+theme_bw()
g5<-pp.ggplot(pv5[,1])+theme_bw()
g10<-pp.ggplot(pv10[,1])+theme_bw()
#g20<-pp.ggplot(pv20[,1])+theme_bw()
#grid.arrange(g0,g1,g2,g3,g5,g10,g20,ncol=3)
grid.arrange(g0,g1,g2,g3,g5,g10,ncol=3)

g1<-ggplot(gwas_m)+geom_point(aes(x=-log10(p0),y=-log10(p1)))
g2<-ggplot(gwas_m)+geom_point(aes(x=-log10(p0),y=-log10(p2)))
g3<-ggplot(gwas_m)+geom_point(aes(x=-log10(p0),y=-log10(p3)))
g4<-ggplot(gwas_m)+geom_point(aes(x=-log10(p0),y=-log10(p5)))
g5<-ggplot(gwas_m)+geom_point(aes(x=-log10(p0),y=-log10(p10)))
g6<-ggplot(gwas_m)+geom_point(aes(x=-log10(p1),y=-log10(p2)))
g7<-ggplot(gwas_m)+geom_point(aes(x=-log10(p1),y=-log10(p3)))
g8<-ggplot(gwas_m)+geom_point(aes(x=-log10(p1),y=-log10(p5)))
g9<-ggplot(gwas_m)+geom_point(aes(x=-log10(p1),y=-log10(p10)))
g10<-ggplot(gwas_m)+geom_point(aes(x=-log10(p2),y=-log10(p3)))
g11<-ggplot(gwas_m)+geom_point(aes(x=-log10(p2),y=-log10(p5)))
g12<-ggplot(gwas_m)+geom_point(aes(x=-log10(p2),y=-log10(p10)))
g13<-ggplot(gwas_m)+geom_point(aes(x=-log10(p3),y=-log10(p5)))
g14<-ggplot(gwas_m)+geom_point(aes(x=-log10(p3),y=-log10(p10)))
g15<-ggplot(gwas_m)+geom_point(aes(x=-log10(p5),y=-log10(p10)))
grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,ncol=3)

ggplot(gwas_m)+
	geom_point(aes(x=ps,y=-log10(p0)))+
	geom_point(aes(x=ps,y=log10(p10)))+
	facet_grid(.~chr,scales='free', space='free')+
	theme_bw()+theme(legend.position='none')


gwas=fread('/work/project/dynagen/seynard/GWAS/infestation/results/gemma_Mellifera_ncorse_varroaphoretic_freq_lmm_freq.assoc.txt')
a<-colsplit(gwas$rs,':',c('CHROM','POS'))
gwas$chr<-a$CHROM
gwas$ps<-a$POS
pv=as.data.frame(gwas%>%select(p_wald))
gwas.ppplot1<-pp.ggplot(pv[,1])+theme_bw()
gwas=fread('/work/project/dynagen/seynard/GWAS/infestation/results/gemma_Mellifera_ncorse_varroaphoretic_egs_lmm_freq.assoc.txt')
a<-colsplit(gwas$rs,':',c('CHROM','POS'))
gwas$chr<-a$CHROM
gwas$ps<-a$POS
pv=as.data.frame(gwas%>%select(p_wald))
gwas.ppplot2<-pp.ggplot(pv[,1])+theme_bw()
gwas=fread('/work/project/dynagen/seynard/GWAS/infestation/results/gemma_Mellifera_ncorse_varroaphoretic_freq_lmm_egs.assoc.txt')
a<-colsplit(gwas$rs,':',c('CHROM','POS'))
gwas$chr<-a$CHROM
gwas$ps<-a$POS
pv=as.data.frame(gwas%>%select(p_wald))
gwas.ppplot3<-pp.ggplot(pv[,1])+theme_bw()
gwas=fread('/work/project/dynagen/seynard/GWAS/infestation/results/gemma_Mellifera_ncorse_varroaphoretic_egs_lmm_egs.assoc.txt')
a<-colsplit(gwas$rs,':',c('CHROM','POS'))
gwas$chr<-a$CHROM
gwas$ps<-a$POS
pv=as.data.frame(gwas%>%select(p_wald))
gwas.ppplot4<-pp.ggplot(pv[,1])+theme_bw()
gwas=fread('/work/project/dynagen/seynard/GWAS/infestation/results/cov_pca/gemma_Mellifera_ncorse_varroaphoretic_freq_lmm_freq_cov.assoc.txt')
a<-colsplit(gwas$rs,':',c('CHROM','POS'))
gwas$chr<-a$CHROM
gwas$ps<-a$POS
pv=as.data.frame(gwas%>%select(p_wald))
gwas.ppplot5<-pp.ggplot(pv[,1])+theme_bw()
gwas=fread('/work/project/dynagen/seynard/GWAS/infestation/results/cov_pca/gemma_Mellifera_ncorse_varroaphoretic_egs_lmm_freq_cov.assoc.txt')
a<-colsplit(gwas$rs,':',c('CHROM','POS'))
gwas$chr<-a$CHROM
gwas$ps<-a$POS
pv=as.data.frame(gwas%>%select(p_wald))
gwas.ppplot6<-pp.ggplot(pv[,1])+theme_bw()
gwas=fread('/work/project/dynagen/seynard/GWAS/infestation/results/cov_pca/gemma_Mellifera_ncorse_varroaphoretic_freq_lmm_egs_cov.assoc.txt')
a<-colsplit(gwas$rs,':',c('CHROM','POS'))
gwas$chr<-a$CHROM
gwas$ps<-a$POS
pv=as.data.frame(gwas%>%select(p_wald))
gwas.ppplot7<-pp.ggplot(pv[,1])+theme_bw()
gwas=fread('/work/project/dynagen/seynard/GWAS/infestation/results/cov_pca/gemma_Mellifera_ncorse_varroaphoretic_egs_lmm_egs_cov.assoc.txt')
a<-colsplit(gwas$rs,':',c('CHROM','POS'))
gwas$chr<-a$CHROM
gwas$ps<-a$POS
pv=as.data.frame(gwas%>%select(p_wald))
gwas.ppplot8<-pp.ggplot(pv[,1])+theme_bw()

gwas=fread('/work/project/dynagen/seynard/GWAS/infestation/results/gemma_Ligustica_Carnica_varroaphoretic_freq_lmm_freq.assoc.txt')
a<-colsplit(gwas$rs,':',c('CHROM','POS'))
gwas$chr<-a$CHROM
gwas$ps<-a$POS
pv=as.data.frame(gwas%>%select(p_wald))
gwas.ppplot9<-pp.ggplot(pv[,1])+theme_bw()
gwas=fread('/work/project/dynagen/seynard/GWAS/infestation/results/gemma_Ligustica_Carnica_varroaphoretic_egs_lmm_freq.assoc.txt')
a<-colsplit(gwas$rs,':',c('CHROM','POS'))
gwas$chr<-a$CHROM
gwas$ps<-a$POS
pv=as.data.frame(gwas%>%select(p_wald))
gwas.ppplot10<-pp.ggplot(pv[,1])+theme_bw()
gwas=fread('/work/project/dynagen/seynard/GWAS/infestation/results/gemma_Ligustica_Carnica_varroaphoretic_freq_lmm_egs.assoc.txt')
a<-colsplit(gwas$rs,':',c('CHROM','POS'))
gwas$chr<-a$CHROM
gwas$ps<-a$POS
pv=as.data.frame(gwas%>%select(p_wald))
gwas.ppplot11<-pp.ggplot(pv[,1])+theme_bw()
gwas=fread('/work/project/dynagen/seynard/GWAS/infestation/results/gemma_Ligustica_Carnica_varroaphoretic_egs_lmm_egs.assoc.txt')
a<-colsplit(gwas$rs,':',c('CHROM','POS'))
gwas$chr<-a$CHROM
gwas$ps<-a$POS
pv=as.data.frame(gwas%>%select(p_wald))
gwas.ppplot12<-pp.ggplot(pv[,1])+theme_bw()
gwas=fread('/work/project/dynagen/seynard/GWAS/infestation/results/cov_pca/gemma_Ligustica_Carnica_varroaphoretic_freq_lmm_freq_cov.assoc.txt')
a<-colsplit(gwas$rs,':',c('CHROM','POS'))
gwas$chr<-a$CHROM
gwas$ps<-a$POS
pv=as.data.frame(gwas%>%select(p_wald))
gwas.ppplot13<-pp.ggplot(pv[,1])+theme_bw()
gwas=fread('/work/project/dynagen/seynard/GWAS/infestation/results/cov_pca/gemma_Ligustica_Carnica_varroaphoretic_egs_lmm_freq_cov.assoc.txt')
a<-colsplit(gwas$rs,':',c('CHROM','POS'))
gwas$chr<-a$CHROM
gwas$ps<-a$POS
pv=as.data.frame(gwas%>%select(p_wald))
gwas.ppplot14<-pp.ggplot(pv[,1])+theme_bw()
gwas=fread('/work/project/dynagen/seynard/GWAS/infestation/results/cov_pca/gemma_Ligustica_Carnica_varroaphoretic_freq_lmm_egs_cov.assoc.txt')
a<-colsplit(gwas$rs,':',c('CHROM','POS'))
gwas$chr<-a$CHROM
gwas$ps<-a$POS
pv=as.data.frame(gwas%>%select(p_wald))
gwas.ppplot15<-pp.ggplot(pv[,1])+theme_bw()
gwas=fread('/work/project/dynagen/seynard/GWAS/infestation/results/cov_pca/gemma_Ligustica_Carnica_varroaphoretic_egs_lmm_egs_cov.assoc.txt')
a<-colsplit(gwas$rs,':',c('CHROM','POS'))
gwas$chr<-a$CHROM
gwas$ps<-a$POS
pv=as.data.frame(gwas%>%select(p_wald))
gwas.ppplot16<-pp.ggplot(pv[,1])+theme_bw()

png('/work/project/dynagen/seynard/GWAS/infestation/plot_test.png',width=1500,height=1500)
grid.arrange(gwas.ppplot1,gwas.ppplot2,gwas.ppplot3,gwas.ppplot4,
gwas.ppplot5,gwas.ppplot6,gwas.ppplot7,gwas.ppplot8,
gwas.ppplot9,gwas.ppplot10,gwas.ppplot11,gwas.ppplot12,
gwas.ppplot13,gwas.ppplot14,gwas.ppplot15,gwas.ppplot16,nrow=4)
dev.off()


gwas1=fread('/work/project/dynagen/seynard/GWAS/infestation/results/gemma_Ligustica_Carnica_varroaphoretic_freq_lmm_egs.assoc.txt')
gwas2=fread('/work/project/dynagen/seynard/GWAS/infestation/results/cov_pca/gemma_Ligustica_Carnica_varroaphoretic_freq_lmm_egs_cov.assoc.txt')

gwas<-merge(gwas1,gwas2,by='rs')
a<-colsplit(gwas$rs,':',c('chr','ps'))
gwas$chr<-a$chr
gwas$ps<-a$ps

png('/work/project/dynagen/seynard/GWAS/infestation/plot_test2.png',width=800,height=200)
ggplot(gwas)+geom_point(aes(x=ps,y=-log10(p_wald.x)),col='red')+geom_point(aes(x=ps,y=log10(p_wald.y)),col='blue')+facet_grid(.~chr,scales='free', space='free')
dev.off()


