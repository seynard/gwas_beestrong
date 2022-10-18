#! /bin/bash
##################################################
#################### PARAMETERS #################### 
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
pop_id2="Mellifera,hybrid,us,all,all_nus,corse"
#pop_id2="Mellifera_ncorse,Caucasica,Ligustica_Carnica,hybrid_corse"
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
################################################

##################################################
#################### PRE-REQUISIT #################### 
# pileup and sync2count colonies to create depth, count 
# admixture on SeqApiPop with k=3, file.3.P = freq_admix
#50K snp list (HAV3_1_50000.txt)
# phenotype files
# ldak
##################################################

##############################################################
#################### IMPORT AND PREP INPUT FILES #################### 
sbatch -W --wrap="cp ${dir_prior}/temp/depth_final.txt.bz2 ${dir_in}/depth.txt.bz2; cp ${dir_prior}/temp/count_ref_final.txt.bz2 ${dir_in}/count_ref.txt.bz2; cp ${dir_prior}/freq_admix${n_pop}.txt ${dir_in}/freq_admix${n_pop}.txt"
sbatch -W --wrap="bzip2 -d ${dir_in}/depth.txt.bz2; bzip2 -d ${dir_in}/count_ref.txt.bz2"

chr=($(cut -f1 -d " " ${dir_in}/depth.txt | sort | uniq ))
chr=("${chr[@]:1}")
sbatch -W --wrap="sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dir_in}/depth.txt"
sbatch -W --wrap="sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dir_in}/count_ref.txt"
sbatch -W --wrap="sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dir_in}/freq_admix${n_pop}.txt"
chr=($(cut -f1 -d " " ${dir_in}/${snp50k} | sort | uniq ))
chr=("${chr[@]:1}")
sbatch -W --wrap="sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dir_in}/${snp50k} "
##############################################################

##########################################################################
#################### QUEEN GENETIC COMPOSITION AND GENOTYPE #################### 
# queen composition (on 50k)
Ncol=$(awk "{print NF;exit}" ${dir_in}/depth.txt)
colname=($(head -n1 ${dir_in}/depth.txt))
col_non_bs=($(echo ${colname[@]//BS*/}))
n_col_snpid=$(echo ${#col_non_bs[@]})
n_col=$((${Ncol}-${n_col_snpid}))
colname=`echo $(echo ${colname[@]}) | tr ' ' ','`
IFS=',' read -r -a array <<< "$colname"
for i in $(seq $((${n_col_snpid}+1)) 1 ${Ncol})
do
j=$((${i}-1))
echo $j
cname=${array[${j}]}
sbatch --mem=100G --job-name='het' --wrap="${script}/run_het.sh ${dir_in} ${cname} ${n_pop} ${snp50k} ${j}"
done	

job_run=($(squeue --noheader --format %i --name het))
job_n=$(echo ${#job_run[@]})
while [ ${job_n} -gt 0 ]
do
wait 10m
job_run=($(squeue --noheader --format %i --name het))
job_n=$(echo ${#job_run[@]})
done
rm ${dir_in}/tmp* ${dir_in}/*_st_het.geno.bz2
mkdir -p ${dir_in}/input
mv ${dir_in}/sim_depth_count*.txt ${dir_in}/input/
bzip2 ${dir_in}/input/*

# group of homogeneous genetic composition
choice='_50k' # or choice=''
sbatch --wrap="Rscript ${script}/q_matrix.r ${dir_in} ${n_pop} ${pheno} ${pop_id} ${choice}"
sbatch --wrap="Rscript ${script}/group_q_matrix.r ${dir_in} . ${pop_id2} ${pop_id} ${compo_threshold} ${pheno} ${choice}"

# group wise genotype
Ncol=$(awk "{print NF;exit}" ${dir_in}/depth.txt)
colname=($(head -n1 ${dir_in}/depth.txt))
col_non_bs=($(echo ${colname[@]//BS*/}))
n_col_snpid=$(echo ${#col_non_bs[@]})
pop_idi=($(echo ${pop_id2}| tr "," "\n"))
for i in ${pop_idi[@]}
do
echo ${i}
sbatch -J hom_freq_${i} --wrap="${script}/run_hom_freq.sh ${dir_in} ${dir} ${script} ${i} ${ncpu} ${n_col_snpid} ${maf_min}"
done
x=$(ls ${dir_in}/freq_imputed_*.txt|wc -l)
while [ ${x} -lt ${#pop_idi[@]} ]
do
wait 30m
x=$(ls ${dir_in}/freq_imputed_*.txt|wc -l)
done
awk '{print $1" "$2" "$3" "$4}' ${dir_in}/depth.txt > ${dir_in}/allele_id_final.txt
tail -n +2 ${dir_in}/allele_id_final.txt  > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/allele_id_final.txt
#########################################################################

#########################################################
#################### PREP PHENOTYPE FILES #################### 
sbatch -o ${dir}/log/prep_pheno.out  -e ${dir}/log/prep_pheno.err -W --wrap="Rscript ${script}/prep_pheno.r ${dir_in} ${dir_out} ${pheno} ${pop_id2} ${pheno_mito}"
#########################################################

#################### RUN GWAS #################### 
list=(`echo ${pop_id2} | sed 's/,/\n/g'`)
jobnum=()
for i in ${list[@]}
do
echo ${i}
if [[ -f ${dir_in}/geno_hom_${i}.bgs ]] 
then
echo 'pop to run' 
job=`sbatch -J ${i} --mem=100G -o ${dir}/log/${i}.out -e ${dir}/log/${i}.err --wrap="${script}/run_type.sh ${dir} ${dir_in} ${dir_out} ${script} ${i} ${maf_min} ${missing_rate}"`
job=$(echo $job | cut -d' ' -f4 )
jobnum+=(${job})
fi
done
pat=$(echo ${jobnum[@]}|tr " " "|")
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
while [ ${x} -gt 0 ]
do
sleep 30m
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
done

mkdir ${dir_in}/data_extra
mkdir ${dir_out}/results_extra
pop_idi=($(echo ${pop_id2}| tr "," "\n"))
for i in ${pop_idi[@]}
do
echo ${i}
mv ${dir_in}/list_${i}.txt ${dir_in}/depth_${i}.txt ${dir_in}/count_ref_${i}.txt ${dir_in}/geno_hom_${i}.egs ${dir_in}/geno_hom_${i}.bgs ${dir_in}/geno_hom_${i}.freqs ${dir_in}/geno_hom_${i}.prob ${dir_in}/geno_hom_${i}.txt ${dir_in}/freq_${i}.txt.bz2 ${dir_in}/freq_imputed_${i}.txt.bz2 ${dir_in}/freq2_${i}.bim ${dir_in}/freq2_${i}.fam ${dir_in}/freq2_${i}.txt ${dir_in}/egs2_${i}.bim ${dir_in}/egs2_${i}.fam ${dir_in}/egs2_${i}.txt ${dir_in}/weights_${i}_freq.short ${dir_in}/in_${i}_freq.freq ${dir_in}/data_extra
mv ${dir_out}/xx_${i}.txt ${dir_out}/results_extra
done

list=(`echo ${pop_id2} | sed 's/,/\n/g'`)
jobnum=()
for i in ${list[@]}
do
echo ${i}
if [[ -f ${dir_in}/data_extra/geno_hom_${i}.bgs ]] 
then
job=`sbatch -J ${i} --mem=100G -o ${dir}/log/plot${i}_extra.out -e ${dir}/log/plot${i}_extra.err --wrap="Rscript ${script}/check_gwas.r ${dir_out}/results_extra ${i}"`
job=$(echo $job | cut -d' ' -f4 )
jobnum+=(${job})
fi
done
pat=$(echo ${jobnum[@]}|tr " " "|")
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
while [ ${x} -gt 0 ]
do
sleep 30m
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
done

