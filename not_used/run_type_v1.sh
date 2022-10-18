#! /bin/bash
dir=${1}
dir_in=${2}
dir_out=${3}
script=${4}
type=${5}
maf_min=${6}
missing_rate=${7}

# GRM freq
echo 'start GRM freq'
#python ${script}/freqLDAK.py ${dir_in} freq_${type}.txt.bz2 freq2_${type}.txt ${maf_min} ${missing_rate}
sbatch --mem=200G -W --wrap="python ${script}/freqLDAK.py ${dir_in} freq_${type}.txt.bz2 freq2_${type}.txt ${maf_min} ${missing_rate}"
in_type=${dir_in}/freq2_${type}.txt
sed -i 's/nan/NA/g' ${in_type}
chr=($(awk '{print $1}' ${in_type} | sort |uniq))
chr=( "${chr[@]/CHROM/}" )
awk '{print $1" "$1":"$2" "0" "$2" "$3" "$4}' ${in_type} > ${dir_in}/freq2_${type}.bim
tail -n +2 ${dir_in}/freq2_${type}.bim > ${dir_in}/tmp.bim && mv ${dir_in}/tmp.bim ${dir_in}/freq2_${type}.bim
awk -F' ' '{print $1" "$1" "0" "0" "0" "0}' ${dir_out}/pheno_${type}.txt > ${dir_in}/freq2_${type}.fam
tail -n +2 ${dir_in}/freq2_${type}.fam > ${dir_in}/tmp.fam && mv ${dir_in}/tmp.fam ${dir_in}/freq2_${type}.fam
cat ${dir_in}/freq2_${type}.fam| tr 'a-z' 'A-Z' > ${dir_in}/tmp.fam && mv ${dir_in}/tmp.fam ${dir_in}/freq2_${type}.fam
jobnum=()
for j in ${chr[@]}
do
#./ldak5.1.linux --cut-weights ${dir_in}/sections${j}_${type}_freq --gen ${dir_in}/freq2_${type}.txt --fam ${dir_in}/freq2_${type}.fam --bim  ${dir_in}/freq2_${type}.bim --gen-skip 1 --gen-headers 4 --gen-probs 1 --chr ${j}
#./ldak5.1.linux --calc-weights-all ${dir_in}/sections${j}_${type}_freq --gen ${dir_in}/freq2_${type}.txt --fam ${dir_in}/freq2_${type}.fam --bim  ${dir_in}/freq2_${type}.bim --gen-skip 1 --gen-headers 4 --gen-probs 1 --chr ${j}
job=`sbatch -J ${type}_w --mem=20G --wrap="./ldak5.1.linux --cut-weights ${dir_in}/sections${j}_${type}_freq --gen ${dir_in}/freq2_${type}.txt --fam ${dir_in}/freq2_${type}.fam --bim  ${dir_in}/freq2_${type}.bim --gen-skip 1 --gen-headers 4 --gen-probs 1 --chr ${j};./ldak5.1.linux --calc-weights-all ${dir_in}/sections${j}_${type}_freq --gen ${dir_in}/freq2_${type}.txt --fam ${dir_in}/freq2_${type}.fam --bim  ${dir_in}/freq2_${type}.bim --gen-skip 1 --gen-headers 4 --gen-probs 1 --chr ${j}"`
job=$(echo $job | cut -d' ' -f4 )
jobnum+=(${job})
done
pat=$(echo ${jobnum[@]}|tr " " "|")
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
while [ ${x} -gt 0 ]
do
sleep 30m
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
done
echo -n > ${dir_in}/weights_${type}_freq.short
for j in ${chr[@]}; do
cat ${dir_in}/sections${j}_${type}_freq/weights.short >> ${dir_in}/weights_${type}_freq.short
done
#./ldak5.1.linux --calc-kins-direct ${dir_in}/LDAKgmat_${type}_freq --gen ${dir_in}/freq2_${type}.txt --fam ${dir_in}/freq2_${type}.fam --bim  ${dir_in}/freq2_${type}.bim --gen-skip 1 --gen-headers 4 --gen-probs 1 --weights ${dir_in}/weights_${type}_freq.short --power -.25 --kinship-raw YES
sbatch -W --wrap="./ldak5.1.linux --calc-kins-direct ${dir_in}/LDAKgmat_${type}_freq --gen ${dir_in}/freq2_${type}.txt --fam ${dir_in}/freq2_${type}.fam --bim  ${dir_in}/freq2_${type}.bim --gen-skip 1 --gen-headers 4 --gen-probs 1 --weights ${dir_in}/weights_${type}_freq.short --power -.25 --kinship-raw YES"
mv ${dir_in}/LDAKgmat_${type}_freq* ${dir_out}
n=($(wc -l ${dir_out}/pheno_${type}.txt))
n=$((n-1))
#./ldak5.1.linux --pca ${dir_in}/pca_${type}_freq --grm ${dir_out}/LDAKgmat_${type}_freq --axes ${n}
sbatch -W --wrap="./ldak5.1.linux --pca ${dir_in}/pca_${type}_freq --grm ${dir_out}/LDAKgmat_${type}_freq --axes ${n}"
mv ${dir_in}/pca_${type}_freq* ${dir_out}
# GRM egs
echo 'start GRM egs'
awk -F',' '{$1=""; print $0}' ${dir_in}/geno_hom_${type}.egs > ${dir_in}/in_${type}.egs
#python ${script}/freqLDAK.py ${dir_in} in_${type}.egs egs2_${type}.txt ${maf_min} ${missing_rate}
sbatch --mem=200G -W --wrap="python ${script}/freqLDAK.py ${dir_in} in_${type}.egs egs2_${type}.txt ${maf_min} ${missing_rate}"
in_type=${dir_in}/egs2_${type}.txt
sed -i 's/nan/NA/g' ${in_type}
chr=($(awk '{print $1}' ${in_type} | sort |uniq))
chr=( "${chr[@]/CHROM/}" )
awk '{print $1" "$1":"$2" "0" "$2" "$3" "$4}' ${in_type} > ${dir_in}/egs2_${type}.bim
tail -n +2 ${dir_in}/egs2_${type}.bim > ${dir_in}/tmp.bim && mv ${dir_in}/tmp.bim ${dir_in}/egs2_${type}.bim
awk -F' ' '{print $1" "$1" "0" "0" "0" "0}' ${dir_out}/pheno_${type}.txt > ${dir_in}/egs2_${type}.fam
tail -n +2 ${dir_in}/egs2_${type}.fam > ${dir_in}/tmp.fam && mv ${dir_in}/tmp.fam ${dir_in}/egs2_${type}.fam
cat ${dir_in}/egs2_${type}.fam| tr 'a-z' 'A-Z' > ${dir_in}/tmp.fam && mv ${dir_in}/tmp.fam ${dir_in}/egs2_${type}.fam
jobnum=()
for j in ${chr[@]}
do
#./ldak5.1.linux --cut-weights ${dir_in}/sections${j}_${type}_egs --gen ${dir_in}/egs2_${type}.txt --fam ${dir_in}/egs2_${type}.fam --bim  ${dir_in}/egs2_${type}.bim --gen-skip 1 --gen-headers 4 --gen-probs 1 --chr ${j}
#./ldak5.1.linux --calc-weights-all ${dir_in}/sections${j}_${type}_egs --gen ${dir_in}/egs2_${type}.txt --fam ${dir_in}/egs2_${type}.fam --bim  ${dir_in}/egs2_${type}.bim --gen-skip 1 --gen-headers 4 --gen-probs 1 --chr ${j}
job=`sbatch -J ${type}_w --mem=20G --wrap="./ldak5.1.linux --cut-weights ${dir_in}/sections${j}_${type}_egs --gen ${dir_in}/egs2_${type}.txt --fam ${dir_in}/egs2_${type}.fam --bim  ${dir_in}/egs2_${type}.bim --gen-skip 1 --gen-headers 4 --gen-probs 1 --chr ${j};./ldak5.1.linux --calc-weights-all ${dir_in}/sections${j}_${type}_egs --gen ${dir_in}/egs2_${type}.txt --fam ${dir_in}/egs2_${type}.fam --bim  ${dir_in}/egs2_${type}.bim --gen-skip 1 --gen-headers 4 --gen-probs 1 --chr ${j}"`
job=$(echo $job | cut -d' ' -f4 )
jobnum+=(${job})
done
pat=$(echo ${jobnum[@]}|tr " " "|")
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
while [ ${x} -gt 0 ]
do
sleep 30m
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
done
echo -n > ${dir_in}/weights_${type}_egs.short
for j in ${chr[@]}; do
cat ${dir_in}/sections${j}_${type}_egs/weights.short >> ${dir_in}/weights_${type}_egs.short
done
#./ldak5.1.linux --calc-kins-direct ${dir_in}/LDAKgmat_${type}_egs --gen ${dir_in}/egs2_${type}.txt --fam ${dir_in}/egs2_${type}.fam --bim  ${dir_in}/egs2_${type}.bim --gen-skip 1 --gen-headers 4 --gen-probs 1 --weights ${dir_in}/weights_${type}_egs.short --power -.25 --kinship-raw YES
sbatch -W --wrap="./ldak5.1.linux --calc-kins-direct ${dir_in}/LDAKgmat_${type}_egs --gen ${dir_in}/egs2_${type}.txt --fam ${dir_in}/egs2_${type}.fam --bim  ${dir_in}/egs2_${type}.bim --gen-skip 1 --gen-headers 4 --gen-probs 1 --weights ${dir_in}/weights_${type}_egs.short --power -.25 --kinship-raw YES"
mv ${dir_in}/LDAKgmat_${type}_egs* ${dir_out}
n=($(wc -l ${dir_out}/pheno_${type}.txt))
n=$((n-1))
#./ldak5.1.linux --pca ${dir_in}/pca_${type}_egs --grm ${dir_out}/LDAKgmat_${type}_egs --axes ${n}
sbatch -W --wrap="./ldak5.1.linux --pca ${dir_in}/pca_${type}_egs --grm ${dir_out}/LDAKgmat_${type}_egs --axes ${n}"
mv ${dir_in}/pca_${type}_egs* ${dir_out}
echo 'start covariate calculation'
#Rscript ${script}/choice_pc.r ${dir_out}/LDAKgmat_${type}_freq.grm.raw ${dir_out}/pc_${type}_freq.txt
sbatch -W --wrap="Rscript ${script}/choice_pc.r ${dir_out}/LDAKgmat_${type}_freq.grm.raw ${dir_out}/pc_${type}_freq.txt"
#Rscript ${script}/choice_pc.r ${dir_out}/LDAKgmat_${type}_egs.grm.raw ${dir_out}/pc_${type}_egs.txt
sbatch -W --wrap="Rscript ${script}/choice_pc.r ${dir_out}/LDAKgmat_${type}_egs.grm.raw ${dir_out}/pc_${type}_egs.txt"
echo 'start GWAS lmm'
# GRM freq GWAS freq lmm
sbatch --mem=200G --wrap="${script}/run_GWAS_gemma_freq_lmm.sh ${dir} ${dir_in} ${dir_out} ${script} ${type} ${maf_min} ${missing_rate} ${dir_out}/LDAKgmat_${type}_freq.grm.raw freq"
# GRM freq GWAS egs lmm
sbatch --mem=200G --wrap="${script}/run_GWAS_gemma_egs_lmm.sh ${dir} ${dir_in} ${dir_out} ${script} ${type} ${maf_min} ${missing_rate} ${dir_out}/LDAKgmat_${type}_freq.grm.raw freq"
## GRM egs GWAS egs lmm
#sbatch --wrap="${script}/run_GWAS_gemma_egs_lmm.sh ${dir} ${dir_in} ${dir_out} ${script} ${type} ${maf_min} ${missing_rate} ${dir_out}/LDAKgmat_${type}_egs.grm.raw egs"
## GRM egs GWAS freq lmm
#sbatch --wrap="${script}/run_GWAS_gemma_freq_lmm.sh ${dir} ${dir_in} ${dir_out} ${script} ${type} ${maf_min} ${missing_rate} ${dir_out}/LDAKgmat_${type}_egs.grm.raw egs"
sleep 30m
echo 'start GWAS bslmm'
# GRM freq GWAS freq bslmm
sbatch --mem=200G --wrap="${script}/run_GWAS_gemma_freq_bslmm.sh ${dir} ${dir_in} ${dir_out} ${script} ${type} ${maf_min} ${missing_rate} ${dir_out}/LDAKgmat_${type}_freq.grm.raw freq"
# GRM freq GWAS egs bslmm
sbatch --mem=200G --wrap="${script}/run_GWAS_gemma_egs_bslmm.sh ${dir} ${dir_in} ${dir_out} ${script} ${type} ${maf_min} ${missing_rate} ${dir_out}/LDAKgmat_${type}_freq.grm.raw freq"



