library(data.table)
library(tidyverse)
library(ggplot2)
library(gridExtra)

q1<-fread('data/queen_geno.Q',data.table=F)
q1<-gather(q1,pop,q_all,c('Mellifera','Ligustica_Carnica','Caucasica'))
q2<-fread('data/queen_geno_seqapipop.Q',data.table=F)
q2<-gather(q2,pop,q_seqapipop,c('Mellifera','Ligustica_Carnica','Caucasica'))
q3<-fread('data/queen_geno_50k.Q',data.table=F)
q3<-gather(q3,pop,q_50k,c('Mellifera','Ligustica_Carnica','Caucasica'))
q<-Reduce(function(x,y) merge(x,y,by=c('num_ruche_bs','pop')),list(q1,q2,q3))

cor.test(q$q_all,q$q_seqapipop)
cor.test(q$q_all,q$q_50k)
cor.test(q$q_seqapipop,q$q_50k)
q$pop<-factor(q$pop,levels=c('Mellifera','Ligustica_Carnica','Caucasica'))

g1<-ggplot(q)+
geom_point(aes(x=q_all,y=q_seqapipop,col=pop))+
scale_colour_manual(values=c('grey34','goldenrod2','chartreuse3'))+
theme_bw()
g2<-ggplot(q)+
geom_point(aes(x=q_all,y=q_50k,col=pop))+
scale_colour_manual(values=c('grey34','goldenrod2','chartreuse3'))+
theme_bw()
g3<-ggplot(q)+
geom_point(aes(x=q_seqapipop,y=q_50k,col=pop))+
scale_colour_manual(values=c('grey34','goldenrod2','chartreuse3'))+
theme_bw()
grid.arrange(g1,g2,g3,nrow=3)

data<-fread('data/Data_BS.csv',data.table=F)
data<-data[,c('num_ruche_bs','structure','pass')]
data$num_ruche_bs<-gsub('bs','BS',data$num_ruche_bs)
data<-unique(data)
df<-merge(data,q,by='num_ruche_bs')
df_lig_all<-unique(subset(df$num_ruche_bs,df$pop=='Ligustica_Carnica' & df$q_all>=0.75))
df_mel_all<-unique(subset(df$num_ruche_bs,df$pop=='Mellifera' & df$q_all>=0.75))
df_cau_all<-unique(subset(df$num_ruche_bs,df$pop=='Caucasica' & df$q_all>=0.75))
df$type_all<-'hyb'
df$type_all[df$num_ruche_bs%in%df_lig_all]<-'lig'
df$type_all[df$num_ruche_bs%in%df_mel_all]<-'mel'
df$type_all[df$num_ruche_bs%in%df_cau_all]<-'cau'
df_lig_seqapipop<-unique(subset(df$num_ruche_bs,df$pop=='Ligustica_Carnica' & df$q_seqapipop>=0.75))
df_mel_seqapipop<-unique(subset(df$num_ruche_bs,df$pop=='Mellifera' & df$q_seqapipop>=0.75))
df_cau_seqapipop<-unique(subset(df$num_ruche_bs,df$pop=='Caucasica' & df$q_seqapipop>=0.75))
df$type_seqapipop<-'hyb'
df$type_seqapipop[df$num_ruche_bs%in%df_lig_seqapipop]<-'lig'
df$type_seqapipop[df$num_ruche_bs%in%df_mel_seqapipop]<-'mel'
df$type_seqapipop[df$num_ruche_bs%in%df_cau_seqapipop]<-'cau'
df_lig_50k<-unique(subset(df$num_ruche_bs,df$pop=='Ligustica_Carnica' & df$q_50k>=0.75))
df_mel_50k<-unique(subset(df$num_ruche_bs,df$pop=='Mellifera' & df$q_50k>=0.75))
df_cau_50k<-unique(subset(df$num_ruche_bs,df$pop=='Caucasica' & df$q_50k>=0.75))
df$type_50k<-'hyb'
df$type_50k[df$num_ruche_bs%in%df_lig_50k]<-'lig'
df$type_50k[df$num_ruche_bs%in%df_mel_50k]<-'mel'
df$type_50k[df$num_ruche_bs%in%df_cau_50k]<-'cau'
