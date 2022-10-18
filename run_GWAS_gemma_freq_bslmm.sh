#! /bin/bash
dir=${1}
dir_in=${2}
dir_out=${3}
script=${4}
type=${5}
maf_min=${6}
missing_rate=${7}
grm=${8}
grm_type=${9}

infile="${dir_in}/in_${type}_${grm_type}.freq"
#if [[ -f "${infile}" ]]
#then 
#echo 'infile exists'
#else
#cp ${dir_in}/freq2_${type}.txt ${infile}
#awk '{$1=$1":"$2 FS $1;}1' ${infile} > ${dir_in}/tmp_${type}_${grm_type}.freq && mv ${dir_in}/tmp_${type}_${grm_type}.freq ${infile}
#awk '{$2=$3=""; print $0}' ${infile} > ${dir_in}/tmp_${type}_${grm_type}.freq && mv ${dir_in}/tmp_${type}_${grm_type}.freq ${infile}
#sed -i 's/ /,/g' ${infile}
#sed -i 's/,,,/,/g' ${infile}
#fi
pheno_id=($(echo ${dir_out}/pheno_*_${type}.txt))
if [[ ${type} == 'corse' ]]
then
	pheno_rem=($(echo ${dir_out}/pheno_*_hybrid_${type}.txt))
	for host in ${pheno_rem[@]}; do
	pheno_id=( "${pheno_id[@]/$host}" )
	done
fi
for i in ${pheno_id[@]}
do
pheno=$(echo ${i} | cut -d '/' -f9)
pheno=$(echo ${pheno} | awk '{gsub(/pheno_/,"")}1')
pheno=$(echo ${pheno} | awk '{gsub(/.txt/,"")}1')
pheno=$(echo ${pheno} | awk '{gsub(/_'"$type"'/,"")}1')
echo ${pheno}
sbatch -o ${dir}/log/gwas_gemma_${type}_${grm_type}_freq_${pheno}_bslmm.out -e ${dir}/log/gwas_gemma_${type}_${grm_type}_freq_${pheno}_bslmm.err --mem=100G --wrap="gemma -g ${infile} -p ${dir_out}/pheno_${pheno}_${type}.txt -k ${grm} -bslmm 1 -n 1 -notsnp -outdir ${dir_out} -o ./gemma_${type}_${pheno}_${grm_type}_bslmm_freq"
done
