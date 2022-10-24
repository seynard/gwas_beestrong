#! /bin/bash
##################################################
#################### PARAMETERS #################### 
dir="~/GWAS/infestation" 
script="~/GWAS/infestation/scripts"  
dir_out="~/GWAS/infestation/results" 
dir_in="~/GWAS/infestation/data" 
dir_save="~/save" 
fasta=${dir_save}/Fasta/GCF_003254395.2_Amel_HAv3.1_genomic.fna 
vcf_name="MetaGenotypesCalled870_raw_snps" 
vcf_sansfiltre=${dir_in}/${vcf_name}.vcf.gz 
vcf_file=${dir_in}/${vcf_name}_allfilter.vcf 
pop_id="Mellifera,Caucasia,Ligustica_Carnica" 
n_pop=$(echo $(IFS=","; set -f; set -- $pop_id; echo $#)) 
pop_id2="Mellifera,Caucasia,Ligustica_Carnica,hybrid" 
pheno="Data_BS.csv" 
pheno_mito='compareDepthsBeeVarroa.txt' 
snp50k='HAV3_1_50000.txt' 
maf_min=0.01 
missing_rate=0.05 
compo_threshold=0.8 
unif_threshold=0.99 
d_min=10 
d_max=50 
pool_size_def=150 
ncpu=10
nbjobs=10 
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
#################### PREP FILES #################### 
sbatch -W --mem=100G --job-name='prep_seqapipop' --wrap="${script}/prep_seqapipop.sh ${script} ${dir_in} ${vcf_file} ${n_pop} ${pop_id} ${unif_threshold} ${ncpu}"
sbatch -W --mem=200G --job-name='prep_pileup' --wrap="${script}/prep_pileup.sh ${script} ${dir_in} ${fasta} ${vcf_name} ${vcf_sansfiltre} ${vcf_file} ${snp50k} ${nbjobs} ${d_min} ${d_max} ${pool_size_def}"
##############################################################

##########################################################################
#################### QUEEN GENETIC COMPOSITION AND GENOTYPE #################### 
# queen composition (on 50K markers)
Ncol=$(awk "{print NF;exit}" ${dir_in}/depth.txt)
colname=($(head -n1 ${dir_in}/depth.txt))
col_non_bs=($(echo ${colname[@]//BS*/}))
n_col_snpid=$(echo ${#col_non_bs[@]})
n_col=$((${Ncol}-${n_col_snpid}))
colname=`echo $(echo ${colname[@]}) | tr ' ' ','`
IFS=',' read -r -a array <<< "$colname"
jobnum=()
for i in $(seq $((${n_col_snpid}+1)) 1 ${Ncol})
do
	j=$((${i}-1))
	echo $j
	cname=${array[${j}]}
	job=`sbatch --mem=100G --job-name='het' --wrap="${script}/run_het.sh ${dir_in} ${cname} ${n_pop} ${snp50k} ${j}"`
	job=$(echo $job | cut -d' ' -f4 )
	jobnum+=(${job})
done	
pat=$(echo ${jobnum[@]}|tr " " "|")
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
while [ ${x} -gt 1 ]
do
	sleep 10m
	x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
done
rm ${dir_in}/tmp* ${dir_in}/*_st_het.geno.bz2
mkdir -p ${dir_in}/input
mv ${dir_in}/sim_depth_count*.txt ${dir_in}/input/
bzip2 ${dir_in}/input/*

# group of homogeneous genetic composition
choice='_50k' 
sbatch -W --wrap="Rscript ${script}/q_matrix.r ${dir_in} ${n_pop} ${pheno} ${pop_id} ${choice}"
sbatch -W --wrap="Rscript ${script}/group_q_matrix.r ${dir_in} . ${pop_id2} ${pop_id} ${compo_threshold} ${pheno} ${choice}"

# group wise genotype
Ncol=$(awk "{print NF;exit}" ${dir_in}/depth.txt)
colname=($(head -n1 ${dir_in}/depth.txt))
col_non_bs=($(echo ${colname[@]//BS*/}))
n_col_snpid=$(echo ${#col_non_bs[@]})
pop_idi=($(echo ${pop_id2}| tr "," "\n"))
jobnum=()
for i in ${pop_idi[@]}
do
	echo ${i}
	job=`sbatch -J hom_freq_${i} --wrap="${script}/run_hom_freq.sh ${dir_in} ${dir} ${script} ${i} ${ncpu} ${n_col_snpid} ${maf_min}"`
	job=$(echo $job | cut -d' ' -f4 )
	jobnum+=(${job})
done
pat=$(echo ${jobnum[@]}|tr " " "|")
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
while [ ${x} -gt 1 ]
do
	sleep 30m
	x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
done
awk '{print $1" "$2" "$3" "$4}' ${dir_in}/depth.txt > ${dir_in}/allele_id_final.txt
tail -n +2 ${dir_in}/allele_id_final.txt  > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/allele_id_final.txt
#########################################################################

#########################################################
#################### PREP PHENOTYPE FILES #################### 
sbatch -o ${dir}/log/prep_pheno.out  -e ${dir}/log/prep_pheno.err -W --wrap="Rscript ${script}/prep_pheno.r ${dir_in} ${dir_out} ${pheno} ${pop_id2} ${pheno_mito}"
#########################################################

#########################################################
#################### ESTIMATE LD #################### 
#make list Lig, Mel, Cau, Hyb
#calculate LD for each group
mkdir ${dir_in}/ld
sbatch -J 'list_ld' --mem=20G -W --wrap="Rscript ${script}/list_ld.r ${dir_in}/seqapipop ${dir_in}/ld ${pop_id2} ${compo_threshold}"
jobnum=()
pop_idi=($(echo ${pop_id2}| tr "," "\n"))
for i in ${pop_idi[@]}
do
	echo ${i}
	sbatch --mem=200G -W --wrap="plink --bfile ${dir_in}/seqapipop/seqapipop --keep ${dir_in}/ld/list_${i}_ld.txt --keep-allele-order --make-bed --out ${dir_in}/ld/${i}"
	chr=($(cut -f1 -d " " ${dir_in}/snp_list.txt | sort | uniq ))
	unset 'chr[${#chr[@]}-1]'
	echo ${chr[@]}
	for c in ${chr[@]}
	do 
		echo $c
		job=`sbatch --mem=200G --wrap="plink --bfile ${dir_in}/ld/${i} --chr $c --keep-allele-order --make-bed --out ${dir_in}/ld/${i}_${c}; plink --bfile ${dir_in}/ld/${i}_${c} --keep-allele-order --ld-window-r2 0.1 --out ${dir_in}/ld/ld_${i}_${c} --r2 inter-chr"`
		jobnum+=(${job})
	done
done
pat=$(echo ${jobnum[@]}|tr " " "|")
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
while [ ${x} -gt 1 ]
do
	sleep 30m
	x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
done
#########################################################

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
while [ ${x} -gt 1 ]
do
	sleep 30m
	x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
done
list=(`echo ${pop_id2} | sed 's/,/\n/g'`)
jobnum=()
for i in ${list[@]}
do
	echo ${i}
	if [[ -f ${dir_in}/geno_hom_${i}.bgs ]] 
	then
		job=`sbatch -J ${i} --mem=100G -o ${dir}/log/plot${i}.out -e ${dir}/log/plot${i}.err --wrap="Rscript ${script}/check_gwas.r ${dir_out} ${i} 0.1"`
		job=$(echo $job | cut -d' ' -f4 )
		jobnum+=(${job})
	fi
done
pat=$(echo ${jobnum[@]}|tr " " "|")
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
while [ ${x} -gt 1 ]
do
	sleep 30m
	x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
done
bzip2 ${dir_out}/summary_gwas_*.txt
#########################################################

#########################################################
#################### DEFINE GROUPS #################### 
pheno_id=($(ls ${dir_out}/pheno_*.txt))
pheno_list=()
for i in ${pheno_id[@]}
do
	pheno=$(echo ${i} | cut -d '/' -f9)
	pheno=$(echo ${pheno} | awk '{gsub(/pheno_/,"")}1')
	pheno=$(echo ${pheno} | awk '{gsub(/.txt/,"")}1')
	ele=($(echo ${pop_id2} |awk -F, '{for (i=1;i<=NF;i++)print $i}'))
	for j in ${ele[@]}
	do
		pheno=$(echo ${pheno} | awk '{gsub(/_'"${j}"'/,"")}1')
		pheno=$(echo ${pheno} | awk '{gsub(/'"${j}"'/,"")}1')
	done
	pheno_list[${#pheno_list[@]}]=${pheno}
done
pheno_list=($(echo "${pheno_list[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
list_pop=($(echo ${pop_id2[@]}| tr "," "\n"))
delete=Caucasia
list_pop=("${list_pop[@]/$delete}")
#########################################################


#########################################################
#################### COMBINE GWAS WITH MANTRA #################### 
for i in ${pheno_list[@]}
do
	job=`sbatch -J 'mantra' -o ${dir}/log/mantra${i}.out -e ${dir}/log/mantra${i}.err --mem=100G --wrap="${script}/run_mantra.sh ${dir} ${script} ${dir_out} ${list_pop} ${i} egs"`
	job=$(echo $job | cut -d' ' -f4 )
	jobnum+=(${job})
done
pat=$(echo ${jobnum[@]}|tr " " "|")
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
while [ ${x} -gt 1 ]
do
	sleep 30m
	x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
done
#########################################################

#########################################################
#################### COMBINE GWAS WITH MASH #################### 
jobnum=()
for i in ${pheno_list[@]}
do
	job=`sbatch -o ${dir}/log/analysis_mash_${i}_egs1.out -e ${dir}/log/analysis_mash_${i}_egs1.err -J 'mash' --mem=200G --wrap="Rscript ${script}/run_mash.r ${dir_out} ${i} egs 0.2 ${list_pop}"`
	job=$(echo $job | cut -d' ' -f4 )
	jobnum+=(${job})
done
pat=$(echo ${jobnum[@]}|tr " " "|")
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
while [ ${x} -gt 1 ]
do
	sleep 30m
	x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
done
#########################################################

#########################################################
#################### MARKER GENE CORRESPONDANCE #################### 
var1='/'
var2='+'
awk -v var1="${var1}" -v var2="${var2}" '{print $1" "$2" "$2" "$3""var1""$4" "var2}' ${dir_in}/allele_id_final.txt > ${dir_in}/input_vep.txt	
awk '$1=="10"{$1="CM009940.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="11"{$1="CM009941.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="12"{$1="CM009942.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="13"{$1="CM009943.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="14"{$1="CM009944.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="15"{$1="CM009945.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="16"{$1="CM009946.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="1"{$1="CM009931.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="2"{$1="CM009932.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="3"{$1="CM009933.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="4"{$1="CM009934.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="5"{$1="CM009935.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="6"{$1="CM009936.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="7"{$1="CM009937.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="8"{$1="CM009938.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="9"{$1="CM009939.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
cp ${dir_in}/input_vep.txt ${dir}/ensembl-vep/
cd ${dir}/ensembl-vep/
./vep -i input_vep.txt -o output_vep.txt --species apis_mellifera --database --genomes --force_overwrite
cd
cp ${dir}/ensembl-vep/output_vep.txt ${dir_out}/output_vep.txt
#########################################################

#########################################################
#################### ANALYSIS ON GENES #################### 
awk -F' ' '$3 == "gene"' ${dir_in}/GCF_003254395.2_Amel_HAv3.1_genomic.gff > ${dir_in}/gene_list_HAV3_1.txt
awk -F' ' '{print $1" "$4" "$5" "$9}' ${dir_in}/gene_list_HAV3_1.txt > ${dir_in}/gene_list_HAV3_1.tmp && mv ${dir_in}/gene_list_HAV3_1.tmp ${dir_in}/gene_list_HAV3_1.txt
locus_name=$(cat ${dir_in}/gene_list_HAV3_1.txt | awk -F' ' '{print $4}' | cut -d';' -f1 | cut -d'-' -f2)
awk -F' ' '{print $1" "$2" "$3" "}' ${dir_in}/gene_list_HAV3_1.txt > ${dir_in}/gene_list_HAV3_1.tmp && mv ${dir_in}/gene_list_HAV3_1.tmp ${dir_in}/gene_list_HAV3_1.txt
paste -d ' ' ${dir_in}/gene_list_HAV3_1.txt <(printf "%s\n" $locus_name) > ${dir_in}/gene_list_HAV3_1.tmp && mv ${dir_in}/gene_list_HAV3_1.tmp ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NC_037638.1/1/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NC_037639.1/2/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NC_037640.1/3/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NC_037641.1/4/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NC_037642.1/5/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NC_037643.1/6/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NC_037644.1/7/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NC_037645.1/8/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NC_037646.1/9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NC_037647.1/10/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NC_037648.1/11/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NC_037649.1/12/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NC_037650.1/13/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NC_037651.1/14/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NC_037652.1/15/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NC_037653.1/16/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NC_001566.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555792.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555794.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555797.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555801.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555807.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555812.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555813.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555814.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555815.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555821.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555816.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555822.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555841.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555842.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555843.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555846.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555850.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555851.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555853.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555855.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555859.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555860.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555862.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555863.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555864.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555866.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555867.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555868.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555869.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555870.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555880.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555886.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555887.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555889.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555891.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555894.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555896.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555901.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555902.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555904.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555907.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555912.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555924.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/NW_020555938.1/-9/g' ${dir_in}/gene_list_HAV3_1.txt
sed -i -e 's/"//g' ${dir_in}/gene_list_HAV3_1.txt
awk '!seen[$0]++' ${dir_in}/gene_list_HAV3_1.txt > ${dir_in}/gene_list_HAV3_1.tmp && mv ${dir_in}/gene_list_HAV3_1.tmp ${dir_in}/gene_list_HAV3_1.txt
sed -i '/^-9/d' ${dir_in}/gene_list_HAV3_1.txt

jobnum=()
for i in ${list_pop[@]}
do
	echo ${i}
	for j in ${pheno_list[@]}
		echo ${j}
		n=$(grep 'number of analyzed individuals' ${dir_out}/gemma_${i}_${j}_freq_lmm_egs_cov.log.txt | cut -d '=' -f2)
		awk -v var="${n}" '{print $2" "$5" "$6" "$7" "$8" "$9" "$12" "var}' ${dir_out}/gemma_${i}_${j}_freq_lmm_egs_cov.assoc.txt > ${dir_out}/${j}_in_${i}.ma
		var="SNP A1 A2 freq BETA SE P N"
		sed -i "1s/.*/$var/" ${dir_out}/${j}_in_${i}.ma
		job=`sbatch --mem=8G --wrap="./gcta-1.94.1 --bfile ${i} --mBAT-combo ${dir_out}/${j}_in_${i}.ma --mBAT-gene-list ${dir_in}/gene_list_HAV3_1.txt --out ${dir_out}/${j}_in_${i} --thread-num 1 --diff-freq 1 --mBAT-wind 5 --fastBAT-ld-cutoff 0.8 --mBAT-print-all-p"`
		jobnum+=(${job})
	done
done
pat=$(echo ${jobnum[@]}|tr " " "|")
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
while [ ${x} -gt 1 ]
do
	sleep 30m
	x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
done
#########################################################


