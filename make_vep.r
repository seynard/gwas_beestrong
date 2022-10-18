library(data.table)
## prep vep
# individual analysis
df<-fread('sign_gwas_ind.txt',data.table=F,fill=T)
df$sp[df$sp=='Mellifera']<-'mel'
df$sp[df$sp=='Hybrids']<-'hyb'
df$sp[df$sp=='Ligustica' & df$V34=='']<-'lig'
df$sp[df$sp=='Ligustica' & df$V34=='us)']<-'lig_nus'
df[,c('V31','V32','V33','V34')]<-NULL
df$CHROM<-NA
df$CHROM[df$chr==1]<-'NC_037638.1'
df$CHROM[df$chr==2]<-'NC_037639.1'
df$CHROM[df$chr==3]<-'NC_037640.1'
df$CHROM[df$chr==4]<-'NC_037641.1'
df$CHROM[df$chr==5]<-'NC_037642.1'
df$CHROM[df$chr==6]<-'NC_037643.1'
df$CHROM[df$chr==7]<-'NC_037644.1'
df$CHROM[df$chr==8]<-'NC_037645.1'
df$CHROM[df$chr==9]<-'NC_037646.1'
df$CHROM[df$chr==10]<-'NC_037647.1'
df$CHROM[df$chr==11]<-'NC_037648.1'
df$CHROM[df$chr==12]<-'NC_037649.1'
df$CHROM[df$chr==13]<-'NC_037650.1'
df$CHROM[df$chr==14]<-'NC_037651.1'
df$CHROM[df$chr==15]<-'NC_037652.1'
df$CHROM[df$chr==16]<-'NC_037653.1'
df$id<-paste0(df$rs,'_',df$sp)
df$snp<-paste0(df$allele1,'/',df$allele0)
df$orient<-'+'
x<-unique(df[,c('CHROM','ps','ps','snp','orient','id')])
write.table(x,'sign_gwas_ind.vep',col.names=F,row.names=F,quote=F,sep=' ')

# meta GWAS
df<-fread('results/sign_locus_1.txt',data.table=F,fill=T)
df$CHROM<-NA
df$CHROM[df$chr==1]<-'NC_037638.1'
df$CHROM[df$chr==2]<-'NC_037639.1'
df$CHROM[df$chr==3]<-'NC_037640.1'
df$CHROM[df$chr==4]<-'NC_037641.1'
df$CHROM[df$chr==5]<-'NC_037642.1'
df$CHROM[df$chr==6]<-'NC_037643.1'
df$CHROM[df$chr==7]<-'NC_037644.1'
df$CHROM[df$chr==8]<-'NC_037645.1'
df$CHROM[df$chr==9]<-'NC_037646.1'
df$CHROM[df$chr==10]<-'NC_037647.1'
df$CHROM[df$chr==11]<-'NC_037648.1'
df$CHROM[df$chr==12]<-'NC_037649.1'
df$CHROM[df$chr==13]<-'NC_037650.1'
df$CHROM[df$chr==14]<-'NC_037651.1'
df$CHROM[df$chr==15]<-'NC_037652.1'
df$CHROM[df$chr==16]<-'NC_037653.1'
df$id<-df$rs
df$snp<-paste0(df$all_effect,'/',df$all_other)
df$orient<-'+'
x<-unique(df[,c('CHROM','ps','ps','snp','orient','id')])
write.table(x,'sign_locus_1.vep',col.names=F,row.names=F,quote=F)

df<-fread('results/sign_locus_2.txt',data.table=F,fill=T)
df$CHROM<-NA
df$CHROM[df$chr==1]<-'NC_037638.1'
df$CHROM[df$chr==2]<-'NC_037639.1'
df$CHROM[df$chr==3]<-'NC_037640.1'
df$CHROM[df$chr==4]<-'NC_037641.1'
df$CHROM[df$chr==5]<-'NC_037642.1'
df$CHROM[df$chr==6]<-'NC_037643.1'
df$CHROM[df$chr==7]<-'NC_037644.1'
df$CHROM[df$chr==8]<-'NC_037645.1'
df$CHROM[df$chr==9]<-'NC_037646.1'
df$CHROM[df$chr==10]<-'NC_037647.1'
df$CHROM[df$chr==11]<-'NC_037648.1'
df$CHROM[df$chr==12]<-'NC_037649.1'
df$CHROM[df$chr==13]<-'NC_037650.1'
df$CHROM[df$chr==14]<-'NC_037651.1'
df$CHROM[df$chr==15]<-'NC_037652.1'
df$CHROM[df$chr==16]<-'NC_037653.1'
df$id<-df$rs
df$snp<-paste0(df$all_effect,'/',df$all_other)
df$orient<-'+'
x<-unique(df[,c('CHROM','ps','ps','snp','orient','id')])
write.table(x,'sign_locus_2.vep',col.names=F,row.names=F,quote=F)


# run with online tool: https://metazoa.ensembl.org/Aedes_aegypti_lvpagwg/Tools/VEP#
#./vep --biotype --buffer_size 5000 --canonical --check_existing --distance 500000 --failed 1 --hgvs --is_multispecies 1 --numbers --protein --species apis_mellifera --symbol --transcript_version --uniprot --cache --input_file [input_data] --output_file [output_file]

