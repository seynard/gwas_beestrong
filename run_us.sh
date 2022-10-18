#! /bin/bash
##################################################
#################### PARAMETERS #################### 
dir="/work/genphyse/dynagen/seynard/GWAS/infestation4/US"
script="/work/genphyse/dynagen/seynard/GWAS/infestation4/scripts"
dir_out="/work/genphyse/dynagen/seynard/GWAS/infestation4/US"
dir_in="/work/genphyse/dynagen/seynard/GWAS/infestation4/data"
dir_save="/genphyse/dynagen/BeeStrong"
#dir_prior="/work/genphyse/dynagen/seynard/GWAS/BeeStrongHAV3_1/" 
fasta=${dir_save}/Fasta/GCF_003254395.2_Amel_HAv3.1_genomic.fna
vcf_name="MetaGenotypesCalled870_raw_snps"
vcf_sansfiltre=${dir_in}/${vcf_name}.vcf.gz
vcf_file=${dir_in}/${vcf_name}_allfilter.vcf 
pheno="Data_BS.csv"
pheno_mito='compareDepthsBeeVarroa.txt'
maf_min=0.01
compo_threshold=0.8
missing_rate=0.05
ncpu=10
unif_threshold=0.99
nbjobs=10
d_min=10
d_max=50
pool_size_def=150
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
#################### COMBINE GWAS WITH MANTRA #################### 
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
pop_idi=($(echo ${pop_id2}| tr "," "\n"))
list_pop=
for i in ${pop_idi[@]}
do
if [[ -f ${dir_out}//pheno_${i}.txt ]] 
then
list_pop+=${i},
fi
done
list_pop=$(echo "${list_pop::-1}")
#pheno_list=(chargevarroaracine4 logit_taux_reop_inf logit_varroainfestation pc1)
#1)
list_pop='Mellifera,Ligustica_Carnica,hybrid'
jobnum=()
for i in ${pheno_list[@]}
do
job=`sbatch -J 'mantra' -o ${dir}/log/mantra${i}.out -e ${dir}/log/mantra${i}.err --mem=100G --wrap="${script}/run_mantra.sh ${dir} ${script} ${dir_out} ${list_pop} ${i} egs 1"`
job=`sbatch -J 'mantra' -o ${dir}/log/mantra${i}.out -e ${dir}/log/mantra${i}.err --mem=100G --wrap="${script}/run_mantra.sh ${dir} ${script} ${dir_out} ${list_pop} ${i} freq 1"`
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
#2)
list_pop='Mellifera,Ligustica_nus,hybrid'
jobnum=()
for i in ${pheno_list[@]}
do
job=`sbatch -J 'mantra' -o ${dir}/log/mantra${i}.out -e ${dir}/log/mantra${i}.err --mem=100G --wrap="${script}/run_mantra.sh ${dir} ${script} ${dir_out} ${list_pop} ${i} egs 2"`
job=`sbatch -J 'mantra' -o ${dir}/log/mantra${i}.out -e ${dir}/log/mantra${i}.err --mem=100G --wrap="${script}/run_mantra.sh ${dir} ${script} ${dir_out} ${list_pop} ${i} freq 2"`
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
pop_idi=($(echo ${pop_id2}| tr "," "\n"))
list_pop=
for i in ${pop_idi[@]}
do
if [[ -f ${dir_out}//pheno_${i}.txt ]] 
then
list_pop+=${i},
fi
done
list_pop=$(echo "${list_pop::-1}")
pheno_id=($(ls ${dir_out}/summary_gwas*.txt.bz2))
pheno_list=()
for i in ${pheno_id[@]}
do
pheno=$(echo ${i} | cut -d '/' -f9)
pheno=$(echo ${pheno} | awk '{gsub(/_freq_freq.txt.bz2/,"")}1')
pheno=$(echo ${pheno} | awk '{gsub(/_freq_egs.txt.bz2/,"")}1')
pheno=$(echo ${pheno} | awk '{gsub(/summary_gwas_/,"")}1')
ele=($(echo ${pop_id2} |awk -F, '{for (i=1;i<=NF;i++)print $i}'))
for j in ${ele[@]}
do
pheno=$(echo ${pheno} | awk '{gsub(/_'"${j}"'/,"")}1')
pheno=$(echo ${pheno} | awk '{gsub(/'"${j}"'_/,"")}1')
done
pheno_list[${#pheno_list[@]}]=${pheno}
done
pheno_list=($(echo "${pheno_list[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
#1)
list_pop='Mellifera,Ligustica_Carnica,hybrid'
jobnum=()
for i in ${pheno_list[@]}
do
job=`sbatch -o ${dir}/log/analysis_mash_${i}_egs1.out -e ${dir}/log/analysis_mash_${i}_egs1.err -J 'mash' --mem=200G --wrap="Rscript ${script}/run_mash.r ${dir_out} ${i} egs 0.2 ${list_pop} 1"`
job=`sbatch -o ${dir}/log/analysis_mash_${i}_freq1.out -e ${dir}/log/analysis_mash_${i}_freq1.err -J 'mash' --mem=200G --wrap="Rscript ${script}/run_mash.r ${dir_out} ${i} freq 0.2 ${list_pop} 1"`
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
#2)
list_pop='Mellifera,Ligustica_nus,hybrid'
jobnum=()
for i in ${pheno_list[@]}
do
job=`sbatch -o ${dir}/log/analysis_mash_${i}_egs2.out -e ${dir}/log/analysis_mash_${i}_egs2.err -J 'mash' --mem=200G --wrap="Rscript ${script}/run_mash.r ${dir_out} ${i} egs 0.2 ${list_pop} 2"`
job=`sbatch -o ${dir}/log/analysis_mash_${i}_freq2.out -e ${dir}/log/analysis_mash_${i}_freq2.err -J 'mash' --mem=200G --wrap="Rscript ${script}/run_mash.r ${dir_out} ${i} freq 0.2 ${list_pop} 2"`
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
job=`sbatch --mem=200G --wrap="plink --bfile ${dir_in}/ld/${i} --chr $c --make-bed --out ${dir_in}/ld/${i}_${c}; plink --bfile ${dir_in}/ld/${i}_${c} --ld-window-r2 0.1 --out ${dir_in}/ld/ld_${i}_${c} --r2 inter-chr"`
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

cd /work/genphyse/dynagen/seynard/GWAS/infestation4
plot_paper_description.r # description phenotypes and populations 
sbatch --mem=200G --wrap="Rscript scripts/plot_gwas_ind.r" # results from individual GWAS
sbatch --mem=200G --wrap="Rscript scripts/analysis.r" # results from meta-analysis
awk '{print $2}' sign_gwas_ind.txt > tmp_ind
awk -F';' '{print $3}' results/sign_locus_1.txt > tmp_meta1
awk -F';' '{print $3}' results/sign_locus_2.txt > tmp_meta2
cat tmp_ind tmp_meta1 tmp_meta2 > sign.txt
rm tmp_ind tmp_meta1 tmp_meta2
cat -n sign.txt | sort -uk2 | sort -nk1 | cut -f2- > tmp_sign && mv tmp_sign sign.txt
sed -i '/rs/d' sign.txt
N=($(awk '{print $1}' ${dir}/sign.txt | uniq))
for i in ${N[@]}
do
rs_snp=$(echo ${i} | cut -f1 -d' ')
chr=$(echo ${rs_snp} | cut -f1 -d':')
sbatch --wrap="${script}/ld_sign.sh ${dir_in}/ld ${rs_snp} ${chr}"
done
# plot_paper_description.r
sbatch --mem=200G --wrap="Rscript scripts/plot_paper_description.r 1"
sbatch --mem=200G --wrap="Rscript scripts/plot_paper_description.r 2"
# plot_overlap.r
sbatch --mem=200G --wrap="Rscript scripts/plot_overlap.r 1"
sbatch --mem=200G --wrap="Rscript scripts/plot_overlap.r 2"
# plot_LD_MantravsMash.r
sbatch --mem=200G --wrap="Rscript scripts/plot_LD_MantravsMash.r 1"
sbatch --mem=200G --wrap="Rscript scripts/plot_LD_MantravsMash.r 2"
# plot_effects_MantravsMash.r
sbatch --mem=200G --wrap="Rscript scripts/plot_effects_MantravsMash.r 1"
sbatch --mem=200G --wrap="Rscript scripts/plot_effects_MantravsMash.r 2"
