#! /bin/bash
dir=${1}
dir_in=${2}
dir_out=${3}
script=${4}
type=${5}
grm=${6}
grm_type=${7}

infile="${dir_in}/in_${type}_${grm_type}.egs"
#if [[ -f "${infile}" ]]
#then 
#echo 'infile exists'
#else
#cp ${dir_in}/egs2_${type}.txt ${infile}
#awk '{$1=$1":"$2 FS $1;}1' ${infile} > ${dir_in}/tmp_${type}_${grm_type}.egs && mv ${dir_in}/tmp_${type}_${grm_type}.egs ${infile}
#awk '{$2=$3=""; print $0}' ${infile} > ${dir_in}/tmp_${type}_${grm_type}.egs && mv ${dir_in}/tmp_${type}_${grm_type}.egs ${infile}
#sed -i 's/ /,/g' ${infile}
#sed -i 's/,,,/,/g' ${infile}
#fi
cov_file="${dir_out}/cov_${type}_${grm_type}.cov"
ncov=$(cat ${dir_out}/pc_${type}_${grm_type}.txt)
pheno_id=($(echo ${dir_out}/pheno_*_${type}.txt))
if [[ ${type} == 'corse' ]]
then
	pheno_rem=($(echo ${dir_out}/pheno_*_hybrid_${type}.txt))
	for host in ${pheno_rem[@]}; do
	pheno_id=( "${pheno_id[@]/$host}" )
	done
fi
if [[ ${ncov} -gt 0 ]]
then
N=$((ncov+2))
awk -v var=${N} '{out=""; for(i=3;i<=var;i++){out=out" "$i}; print out}' ${dir_out}/pca_${type}_${grm_type}.vect > ${cov_file}
for i in ${pheno_id[@]}
do
pheno=$(echo ${i} | cut -d '/' -f9)
pheno=$(echo ${pheno} | awk '{gsub(/pheno_/,"")}1')
pheno=$(echo ${pheno} | awk '{gsub(/.txt/,"")}1')
pheno=$(echo ${pheno} | awk '{gsub(/_'"$type"'/,"")}1')
echo ${pheno}
sbatch -o ${dir}/log/gwas_gemma_${type}_${grm_type}_egs_${pheno}.out -e ${dir}/log/gwas_gemma_${type}_${grm_type}_egs_${pheno}.err --mem=100G --wrap="gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -k ${grm} -lmm 1 -c ${cov_file} -n 1 -notsnp -outdir ${dir_out} -o ./gemma_${type}_${pheno}_${grm_type}_lmm_egs_cov"
done
else
for i in ${pheno_id[@]}
do
pheno=$(echo ${i} | cut -d '/' -f9)
pheno=$(echo ${pheno} | awk '{gsub(/pheno_/,"")}1')
pheno=$(echo ${pheno} | awk '{gsub(/.txt/,"")}1')
pheno=$(echo ${pheno} | awk '{gsub(/_'"$type"'/,"")}1')
echo ${pheno}
sbatch -o ${dir}/log/gwas_gemma_${type}_${grm_type}_egs_${pheno}.out -e ${dir}/log/gwas_gemma_${type}_${grm_type}_egs_${pheno}.err --mem=100G --wrap="gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -k ${grm} -lmm 1 -n 1 -notsnp -outdir ${dir_out} -o ./gemma_${type}_${pheno}_${grm_type}_lmm_egs_cov"
done
fi

