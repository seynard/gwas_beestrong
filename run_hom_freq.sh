#! /bin/bash
############################################ parameters ################################################## 
dir_in=${1}
dir=${2}
script=${3}
i=${4}
ncpu=${5}
n_col_snpid=${6}
maf_min=${7}
########################################################################################################## 
n=$(wc -l < ${dir_in}/list_${i}.txt)
echo ${i} ${n}
if [ ${n} -gt 50 ] 
then
# queen genotypes
ID=($(cat ${dir_in}/list_${i}.txt))
col_id=$(seq 1 1 ${n_col_snpid})
col_id=`echo $col_id | sed 's/ /,/g'`
for j in ${ID[@]}
do
a=$(head -n1 ${dir_in}/depth.txt | tr " " "\n"| grep -nx ${j} |  cut -d":" -f1)
if [[ ${a} == "" ]] 
then
echo $j 'not sequenced'
else
col_id+=','${a}
fi
done
cut -d" " -f${col_id} ${dir_in}/depth.txt > ${dir_in}/depth_${i}.txt
cut -d" " -f${col_id} ${dir_in}/count_ref.txt > ${dir_in}/count_ref_${i}.txt
n_col_snpid=${n_col_snpid}
ncpu=${ncpu}
sbatch -W -J hom --mem=40G -c ${ncpu} --wrap="python ${script}/beethoven/model_genoqueen_hom.py ${dir_in} depth_${i}.txt count_ref_${i}.txt ${n_col_snpid} ${ncpu} geno_hom_${i} 1000"
cut -f1 -d"," --complement ${dir_in}/geno_hom_${i}.bgs > ${dir_in}/geno_hom_${i}.txt
sed -i "s/,/ /g" ${dir_in}/geno_hom_${i}.txt
# frequency estimations
sbatch -W -J freq_${i} -o ${dir}/log/freq_${i}.out -e ${dir}/log/freq_${i}.err --wrap="python ${script}/make_freq.py ${dir_in}/depth_${i}.txt ${dir_in}/count_ref_${i}.txt ${dir_in}/freq_${i}.txt"
sbatch -W -J freq_missing_${i} -o ${dir}/log/freq_missing_${i}.out -e ${dir}/log/freq_missing_${i}.err --wrap="python ${script}/missing_freq.py ${dir_in}/freq_${i}.txt ${dir_in}/freq_imputed_${i}.txt"
sbatch -W --wrap="bzip2 -f ${dir_in}/freq_${i}.txt; bzip2 -f ${dir_in}/freq_imputed_${i}.txt"
#sbatch -W -J filter_freq_${i} -o ${dir}/log/filter_freq_${i}.out -e ${dir}/log/filter_freq_${i}.err --wrap="python ${script}/filter.py ${dir_in}/freq_imputed_${i} ${maf_min}; sed -i 's/ //g' ${dir_in}/freq_imputed_${i}_${maf_min}.txt.filter"
n=$(wc -l < ${dir_in}/geno_hom_${i}.txt)
echo $n
if [ ${n} == 0 ]
then
rm ${dir_in}/geno_hom_${i}.txt
fi
fi
