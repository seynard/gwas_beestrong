#! /bin/bash
dir=${1}
dir_in=${2}
dir_out=${3}
script=${4}
type=${5}
maf_min=${6}
missing_rate=${7}

## prep 
#type_data='freq,egs'
#list_data=(`echo ${type_data} | sed 's/,/\n/g'`)
#for x in ${list_data[@]}
#do
#echo ${x}
#if [[ ${x} == 'egs' ]]
#then
#awk -F',' '{$1=""; print $0}' ${dir_in}/geno_hom_${type}.${x} > ${dir_in}/in_${type}.${x}
#sbatch --mem=200G -W --wrap="python ${script}/freqLDAK.py ${dir_in} in_${type}.${x} ${x}2_${type}.txt ${maf_min} ${missing_rate}"
#else
#sbatch --mem=200G -W --wrap="python ${script}/freqLDAK.py ${dir_in} ${x}_${type}.txt.bz2 ${x}2_${type}.txt ${maf_min} ${missing_rate}"
#fi
#in_type=${dir_in}/${x}2_${type}.txt
#sed -i 's/nan/NA/g' ${in_type}
#chr=($(awk '{print $1}' ${in_type} | sort |uniq))
#chr=( "${chr[@]/CHROM/}" )
#awk '{print $1" "$1":"$2" "0" "$2" "$3" "$4}' ${in_type} > ${dir_in}/${x}2_${type}.bim
#tail -n +2 ${dir_in}/${x}2_${type}.bim > ${dir_in}/tmp.bim && mv ${dir_in}/tmp.bim ${dir_in}/${x}2_${type}.bim
#awk -F' ' '{print $1" "$1" "0" "0" "0" "0}' ${dir_out}/pheno_${type}.txt > ${dir_in}/${x}2_${type}.fam
#tail -n +2 ${dir_in}/${x}2_${type}.fam > ${dir_in}/tmp.fam && mv ${dir_in}/tmp.fam ${dir_in}/${x}2_${type}.fam
#cat ${dir_in}/${x}2_${type}.fam| tr 'a-z' 'A-Z' > ${dir_in}/tmp.fam && mv ${dir_in}/tmp.fam ${dir_in}/${x}2_${type}.fam
#done
## GRM freq
#echo 'start GRM freq'
#jobnum=()
#for j in ${chr[@]}
#do
#job=`sbatch -J ${type}_w --mem=20G --wrap="./ldak5.1.linux --cut-weights ${dir_in}/sections${j}_${type}_freq --gen ${dir_in}/freq2_${type}.txt --fam ${dir_in}/freq2_${type}.fam --bim  ${dir_in}/freq2_${type}.bim --gen-skip 1 --gen-headers 4 --gen-probs 1 --chr ${j};./ldak5.1.linux --calc-weights-all ${dir_in}/sections${j}_${type}_freq --gen ${dir_in}/freq2_${type}.txt --fam ${dir_in}/freq2_${type}.fam --bim  ${dir_in}/freq2_${type}.bim --gen-skip 1 --gen-headers 4 --gen-probs 1 --chr ${j}"`
#job=$(echo $job | cut -d' ' -f4 )
#jobnum+=(${job})
#done
#pat=$(echo ${jobnum[@]}|tr " " "|")
#x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
#while [ ${x} -gt 0 ]
#do
#sleep 30m
#x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
#done
#echo -n > ${dir_in}/weights_${type}_freq.short
#for j in ${chr[@]}; do
#cat ${dir_in}/sections${j}_${type}_freq/weights.short >> ${dir_in}/weights_${type}_freq.short
#done
#sbatch -W --wrap="./ldak5.1.linux --calc-kins-direct ${dir_in}/LDAKgmat_${type}_freq --gen ${dir_in}/freq2_${type}.txt --fam ${dir_in}/freq2_${type}.fam --bim  ${dir_in}/freq2_${type}.bim --gen-skip 1 --gen-headers 4 --gen-probs 1 --weights ${dir_in}/weights_${type}_freq.short --power -.25 --kinship-raw YES"
#mv ${dir_in}/LDAKgmat_${type}_freq* ${dir_out}
#n=($(wc -l ${dir_out}/pheno_${type}.txt))
#n=$((n-1))
#sbatch -W --wrap="./ldak5.1.linux --pca ${dir_in}/pca_${type}_freq --grm ${dir_out}/LDAKgmat_${type}_freq --axes ${n}"
#mv ${dir_in}/pca_${type}_freq* ${dir_out}
#echo 'start covariate calculation'
#sbatch -W --wrap="Rscript ${script}/choice_pc.r ${dir_out}/LDAKgmat_${type}_freq.grm.raw ${dir_out}/pc_${type}_freq.txt"
#echo 'start GWAS lmm'
#infile="${dir_in}/in_${type}_freq.freq"
#cp ${dir_in}/freq2_${type}.txt ${infile}
#awk '{$1=$1":"$2 FS $1;}1' ${infile} > ${dir_in}/tmp_${type}_freq.freq && mv ${dir_in}/tmp_${type}_freq.freq ${infile}
#awk '{$2=$3=""; print $0}' ${infile} > ${dir_in}/tmp_${type}_freq.freq && mv ${dir_in}/tmp_${type}_freq.freq ${infile}
#sed -i 's/ /,/g' ${infile}
#sed -i 's/,,,/,/g' ${infile}
#infile="${dir_in}/in_${type}_freq.egs"
#cp ${dir_in}/egs2_${type}.txt ${infile}
#awk '{$1=$1":"$2 FS $1;}1' ${infile} > ${dir_in}/tmp_${type}_freq.egs && mv ${dir_in}/tmp_${type}_freq.egs ${infile}
#awk '{$2=$3=""; print $0}' ${infile} > ${dir_in}/tmp_${type}_freq.egs && mv ${dir_in}/tmp_${type}_freq.egs ${infile}
#sed -i 's/ /,/g' ${infile}
#sed -i 's/,,,/,/g' ${infile}
# GRM freq GWAS freq lmm
sbatch --mem=200G --wrap="${script}/run_GWAS_gemma_freq_lmm.sh ${dir} ${dir_in} ${dir_out} ${script} ${type} ${dir_out}/LDAKgmat_${type}_freq.grm.raw freq"
# GRM freq GWAS egs lmm
sbatch --mem=200G --wrap="${script}/run_GWAS_gemma_egs_lmm.sh ${dir} ${dir_in} ${dir_out} ${script} ${type} ${dir_out}/LDAKgmat_${type}_freq.grm.raw freq"
echo 'start GWAS bslmm'
# GRM freq GWAS freq bslmm
sbatch --mem=200G --wrap="${script}/run_GWAS_gemma_freq_bslmm.sh ${dir} ${dir_in} ${dir_out} ${script} ${type} ${maf_min} ${missing_rate} ${dir_out}/LDAKgmat_${type}_freq.grm.raw freq"
# GRM freq GWAS egs bslmm
sbatch --mem=200G --wrap="${script}/run_GWAS_gemma_egs_bslmm.sh ${dir} ${dir_in} ${dir_out} ${script} ${type} ${maf_min} ${missing_rate} ${dir_out}/LDAKgmat_${type}_freq.grm.raw freq"



