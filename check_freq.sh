snp=('1:10080627' '1:12552506' '1:15280956' '1:16327085' '1:20960056' '1:21374478' '1:24201224' '1:25184394' '1:2891204' '1:7448807' '1:7448811' '10:10400687' '10:2026877' '10:4266011' '10:5359169' '10:5359173' '11:14369154' '11:9369229' '11:9527267' '12:10153855' '12:10734707' '12:136634' '13:9483955' '14:3782741' '14:6686131' '14:8481541' '14:8481589' '15:2021142' '15:2081876' '15:2081914' '15:4853529' '15:8485332' '16:1812909' '16:5024160' '2:11333369' '2:12025610' '2:12025647' '2:16060868' '2:2729874' '2:4437645' '2:8350714' '2:9347464' '3:1059430' '3:12973246' '3:12973248' '3:6206342' '4:10789077' '4:11665460' '4:4327611' '4:7321246' '4:7321247' '5:2008472' '5:2365495' '5:6736534' '5:6761414' '5:75369' '5:8737386' '5:9190579' '6:10450971' '7:11806658' '7:5762037' '7:5772089' '7:6738985' '7:7028040' '7:7051965' '7:7078376' '7:8466948' '8:1150346' '8:1319815' '8:1551638' '8:1748209' '8:2468335' '8:9557205' '8:9799408' '9:11564671')
head -n1 data/in_Mellifera_freq.freq  >> snp_paper_mel.txt
head -n1 data/in_Ligustica_nus_freq.freq  >> snp_paper_lignus.txt
head -n1 data/in_hybrid_freq.freq  >> snp_paper_hyb.txt

for i in ${snp[@]}
do
echo $i
grep ${i} data/in_Mellifera_freq.egs  >> snp_paper_mel.txt
grep ${i} data/in_Ligustica_nus_freq.egs  >> snp_paper_lignus.txt
grep ${i} data/in_hybrid_freq.egs  >> snp_paper_hyb.txt
done



#module load system/R-4.1.1_gcc-9.3.0
#R
library(data.table)
library(tidyverse)
library(ggplot2)
library(reshape2)
bs<-fread('data/Data_BS.csv',data.table=F)
bs<-bs[,c('num_ruche_bs','structure','api','rucher','altitude','region','departement','annee_reine','mode_fecond')]
bs$col<-toupper(bs$num_ruche_bs)

f_mel<-fread('snp_paper_mel.txt',data.table=F)
f_mel<-f_mel%>%gather('col','f',colnames(f_mel)[4:ncol(f_mel)])
f_mel$group<-'mel'
f_lignus<-fread('snp_paper_lignus.txt',data.table=F)
f_lignus<-f_lignus%>%gather('col','f',colnames(f_lignus)[4:ncol(f_lignus)])
f_lignus$group<-'lignus'
f_hyb<-fread('snp_paper_hyb.txt',data.table=F)
f_hyb<-f_hyb%>%gather('col','f',colnames(f_hyb)[4:ncol(f_hyb)])
f_hyb$group<-'hyb'
f<-do.call(rbind,list(f_mel,f_lignus,f_hyb))
colnames(f)[1]<-'rs'
mean_f<-f%>%group_by(rs,group)%>%summarise(f_m=mean(f),f_sd=sd(f))
mean_f$ci_min<-mean_f$f_m-1.96*mean_f$f_sd
mean_f$ci_max<-mean_f$f_m+1.96*mean_f$f_sd
F<-spread(mean_f[,c('rs','f_m','group')],group,f_m)
F$sd<-apply(F[,2:ncol(F)],1,sd)

ggplot(f)+
	geom_jitter(aes(x=rs,y=f,col=group))+
	facet_grid(group~rs,scale='free_x')


x<-f[f$rs=='15:2081876',c('col','f','group')]
a<-colsplit(x$col,'_',c('year','nb'))
x$year<-a$year
x$nb<-a$nb
x<-merge(x,bs,by='col')

ggplot(x)+
	geom_boxplot(aes(x=year,y=f,fill=group))+
	facet_wrap(~structure)










