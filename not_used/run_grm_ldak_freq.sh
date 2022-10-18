#! /bin/bash
dir_in=${1}
dir_out=${2}
script=${3}
pheno=${4}
type=${5}
maf_min=${6}
missing_rate=${7}

python ${script}/freqLDAK.py ${dir_in} freq_${type}.txt.bz2 freq2_${type}.txt ${maf_min} ${missing_rate}
in_type=${dir_in}/freq2_${type}.txt
sed -i 's/nan/NA/g' ${in_type}
chr=($(awk '{print $1}' ${in_type} | sort |uniq))
chr=( "${chr[@]/CHROM/}" )
awk '{print $1" "$1":"$2" "0" "$2" "$3" "$4}' ${in_type} > ${dir_in}/freq2_${type}.bim
tail -n +2 ${dir_in}/freq2_${type}.bim > ${dir_in}/tmp.bim && mv ${dir_in}/tmp.bim ${dir_in}/freq2_${type}.bim
awk -F' ' '{print $1" "$1" "0" "0" "0" "0}' ${dir_out}/pheno_${type}.txt > ${dir_in}/freq2_${type}.fam
tail -n +2 ${dir_in}/freq2_${type}.fam > ${dir_in}/tmp.fam && mv ${dir_in}/tmp.fam ${dir_in}/freq2_${type}.fam
cat ${dir_in}/freq2_${type}.fam| tr 'a-z' 'A-Z' > ${dir_in}/tmp.fam && mv ${dir_in}/tmp.fam ${dir_in}/freq2_${type}.fam
for j in ${chr[@]}
do
./ldak5.1.linux --cut-weights ${dir_in}/sections${j}_${type}_freq --gen ${dir_in}/freq2_${type}.txt --fam ${dir_in}/freq2_${type}.fam --bim  ${dir_in}/freq2_${type}.bim --gen-skip 1 --gen-headers 4 --gen-probs 1 --chr ${j}
./ldak5.1.linux --calc-weights-all ${dir_in}/sections${j}_${type}_freq --gen ${dir_in}/freq2_${type}.txt --fam ${dir_in}/freq2_${type}.fam --bim  ${dir_in}/freq2_${type}.bim --gen-skip 1 --gen-headers 4 --gen-probs 1 --chr ${j}
done
echo -n > ${dir_in}/weights_${type}_freq.short
for j in ${chr[@]}; do
cat ${dir_in}/sections${j}_${type}_freq/weights.short >> ${dir_in}/weights_${type}_freq.short
done
./ldak5.1.linux --calc-kins-direct ${dir_in}/LDAKgmat_${type}_freq --gen ${dir_in}/freq2_${type}.txt --fam ${dir_in}/freq2_${type}.fam --bim  ${dir_in}/freq2_${type}.bim --gen-skip 1 --gen-headers 4 --gen-probs 1 --weights ${dir_in}/weights_${type}_freq.short --power -.25 --kinship-raw YES
mv ${dir_in}/LDAKgmat_${type}_freq* ${dir_out}

