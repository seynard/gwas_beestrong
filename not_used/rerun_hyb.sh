dir="/work/genphyse/dynagen/seynard/GWAS/infestation"
script="/work/genphyse/dynagen/seynard/GWAS/infestation/scripts"
dir_out="/work/genphyse/dynagen/seynard/GWAS/infestation/results"
dir_in="/work/genphyse/dynagen/seynard/GWAS/infestation/data"
dir_prior="/work/genphyse/dynagen/seynard/GWAS/BeeStrongHAV3_1/" 
fasta=${dir_save}/Fasta/GCF_003254395.2_Amel_HAv3.1_genomic.fna
vcf_name="MetaGenotypesCalled870_raw_snps"
vcf_file=${dir_in}/${vcf_name}_allfilter.vcf 
pop_id="Mellifera,Caucasica,Ligustica_Carnica"
n_pop=$(echo $(IFS=","; set -f; set -- $pop_id; echo $#))
pop_id2="Mellifera,Mellifera_ncorse,Caucasica,Ligustica_Carnica,hybrid,us,all,all_nus,corse,hybrid_corse"
pheno="Data_BS.csv"
pheno_mito='compareDepthsBeeVarroa.txt'
snp50k='HAV3_1_50000.txt'
maf_min=0.01
compo_threshold=0.8
missing_rate=0.05
ncpu=10
module load system/R-4.1.1_gcc-9.3.0
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

type='hybrid'
Ncol=$(awk "{print NF;exit}" ${dir_in}/depth.txt)
colname=($(head -n1 ${dir_in}/depth.txt))
col_non_bs=($(echo ${colname[@]//BS*/}))
n_col_snpid=$(echo ${#col_non_bs[@]})
sbatch -J hom_freq_${type} --wrap="${script}/run_hom_freq.sh ${dir_in} ${dir} ${script} ${type} ${ncpu} ${n_col_snpid} ${maf_min}"
sbatch -J ${type} -o ${dir}/log/${type}.out -e ${dir}/log/${type}.err --mem=100G --wrap="${script}/run_type.sh ${dir} ${dir_in} ${dir_out} ${script} ${type} ${maf_min} ${missing_rate}"
sbatch -J ${type} --mem=100G -o ${dir}/log/plot${type}.out -e ${dir}/log/plot${type}.err --wrap="Rscript ${script}/check_gwas.r ${dir_out} ${type}"
sbatch --mem=100G --wrap="Rscript ${script}/gwas_combine.r ${dir_out} hybrid"

##################################################
dir="/work/genphyse/dynagen/seynard/GWAS/infestation"
script="/work/genphyse/dynagen/seynard/GWAS/infestation/scripts"
dir_out="/work/genphyse/dynagen/seynard/GWAS/infestation/results"
dir_in="/work/genphyse/dynagen/seynard/GWAS/infestation/data"
dir_prior="/work/genphyse/dynagen/seynard/GWAS/BeeStrongHAV3_1/" 
fasta=${dir_save}/Fasta/GCF_003254395.2_Amel_HAv3.1_genomic.fna
vcf_name="MetaGenotypesCalled870_raw_snps"
vcf_file=${dir_in}/${vcf_name}_allfilter.vcf 
pop_id="Mellifera,Caucasica,Ligustica_Carnica"
n_pop=$(echo $(IFS=","; set -f; set -- $pop_id; echo $#))
pop_id2="Mellifera,Mellifera_ncorse,Caucasica,Ligustica_Carnica,hybrid,us,all,all_nus,corse,hybrid_corse"
pheno="Data_BS.csv"
pheno_mito='compareDepthsBeeVarroa.txt'
snp50k='HAV3_1_50000.txt'
maf_min=0.01
compo_threshold=0.8
missing_rate=0.05
ncpu=10
module load system/R-4.1.1_gcc-9.3.0
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

type='corse'
in_file="${dir_in}/freq2_${type}.txt"
infile="${dir_in}/freq_gemma_${type}.txt"
cp ${in_file} ${infile}
awk -F' ' '{$1=$1":"$2;NF=1}1' OFS='' ${infile} > ${dir_in}/tmp_freq_corse.txt
awk '{$1=$2=""; print $0}' ${infile} > ${dir_in}/tmp2_freq_corse.txt
paste --delimiter=' ' ${dir_in}/tmp_freq_corse.txt ${dir_in}/tmp2_freq_corse.txt > ${infile}
tail -n +2 ${infile} > ${dir_in}/tmp_freq_corse.txt && mv ${dir_in}/tmp_freq_corse.txt ${infile}
sed -i 's/ /,/g' ${infile}
sed -i 's/,,,,,/,/g' ${infile}
gmat="grm_gemma_${type}_freq"
sbatch --mem=20G --wrap="gemma -g ${infile} -p ${dir_out}/pheno_eb_smr_${type}.txt -gk 1 -maf ${maf_min} -notsnp -outdir ${dir_out} -o ./${gmat}"
rm ${dir_in}/tmp2_freq_corse.txt
grm="${dir_out}/grm_gemma_${type}_freq.cXX.txt"
infile="${dir_in}/in_${type}_freq.egs"
pheno='varroaphoretic'
sbatch --mem=100G --wrap="gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -k ${grm} -lmm 1 -n 1 -notsnp -outdir ${dir_out} -o ./gemma_${type}_${pheno}_grm_gemma_freq_lmm_egs"
#



