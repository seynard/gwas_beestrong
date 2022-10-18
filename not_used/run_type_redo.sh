#! /bin/bash
dir=${1}
dir_in=${2}
dir_out=${3}
script=${4}
type=${5}
maf_min=${6}
missing_rate=${7}

# prep 
type_data='freq,egs'
list_data=(`echo ${type_data} | sed 's/,/\n/g'`)
for x in ${list_data[@]}
do
echo ${x}
if [[ ${x} == 'egs' ]]
then
awk -F',' '{$1=""; print $0}' ${dir_in}/geno_hom_${type}.${x} > ${dir_in}/in_${type}.${x}
sbatch --mem=200G -W --wrap="python ${script}/freqLDAK.py ${dir_in} in_${type}.${x} ${x}2_${type}.txt ${maf_min} ${missing_rate}"
else
sbatch --mem=200G -W --wrap="python ${script}/freqLDAK.py ${dir_in} ${x}_${type}.txt.bz2 ${x}2_${type}.txt ${maf_min} ${missing_rate}"
fi
in_type=${dir_in}/${x}2_${type}.txt
sed -i 's/nan/NA/g' ${in_type}
chr=($(awk '{print $1}' ${in_type} | sort |uniq))
chr=( "${chr[@]/CHROM/}" )
awk '{print $1" "$1":"$2" "0" "$2" "$3" "$4}' ${in_type} > ${dir_in}/${x}2_${type}.bim
tail -n +2 ${dir_in}/${x}2_${type}.bim > ${dir_in}/tmp.bim && mv ${dir_in}/tmp.bim ${dir_in}/${x}2_${type}.bim
awk -F' ' '{print $1" "$1" "0" "0" "0" "0}' ${dir_out}/pheno_${type}.txt > ${dir_in}/${x}2_${type}.fam
tail -n +2 ${dir_in}/${x}2_${type}.fam > ${dir_in}/tmp.fam && mv ${dir_in}/tmp.fam ${dir_in}/${x}2_${type}.fam
cat ${dir_in}/${x}2_${type}.fam| tr 'a-z' 'A-Z' > ${dir_in}/tmp.fam && mv ${dir_in}/tmp.fam ${dir_in}/${x}2_${type}.fam
done
