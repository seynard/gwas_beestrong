#module load system/R-4.1.1_gcc-9.3.0
library(data.table)
library(paran)
library(PCAtools)
args<-commandArgs(TRUE)
input<-args[1]
output<-args[2]
#input<-'results/LDAKgmat_Mellifera_ncorse_freq.grm.raw'
#output<-'results/pc_Mellifera_ncorse_freq.txt'


df<-fread(input,data.table=F)
horn<-parallelPCA(df)
write.table(horn$n,output,col.names=F,row.names=F,quote=F)



