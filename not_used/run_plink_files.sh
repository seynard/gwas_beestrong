#! /bin/bash
dir_in=${1}
dir_out=${2}
script=${3}
pheno=${4}
maf_min=${5}
type=${6}

chr=($(awk '{print $1}' ${dir_in}/geno_hom_${type}.txt | sort |uniq))
chr=( "${chr[@]/CHROM/}" )
rev=()
for (( i=${#chr[@]}-1; i>=0; i-- ));do rev[${#rev[@]}]=${chr[i]};done
echo  ${rev[@]}
for i in ${rev[@]}
do
awk '$1=='${i}'' ${dir_in}/geno_hom_${type}.txt >  ${dir_in}/geno_hom_${type}_${i}.txt
if [ ${i} == 1 ]
then
sbatch -W --mem=200G --wrap="Rscript ${script}/genotoped.r ${dir_in} ${dir_out} geno_hom_${type} allele_id_final.txt ${i}"
else
sbatch -W --mem=100G --wrap="Rscript ${script}/genotoped.r ${dir_in} ${dir_out} geno_hom_${type} allele_id_final.txt ${i}"
fi
done
echo -n > ${dir_in}/geno_hom_${type}.tped
for i in ${chr[@]}
do 
cat ${dir_in}/geno_hom_${type}_${i}.tped>> ${dir_in}/geno_hom_${type}.tped
done
awk '{print $1":"$2" "$3}' ${dir_in}/allele_id_final.txt > ${dir_in}/list_ref.txt
sbatch -W --wrap="plink --tfile ${dir_in}/geno_hom_${type} --reference-allele ${dir_in}/list_ref.txt --make-bed --maf ${maf_min} --out ${dir_in}/geno_hom_${type}_maf"
chr=($(awk '{print $1}' ${dir_in}/geno_hom_${type}.txt | sort |uniq))
chr=( "${chr[@]/CHROM/}" )
for i in ${chr[@]}
do
rm ${dir_in}/geno_hom_${type}_${i}.tped
rm ${dir_in}/geno_hom_${type}_${i}.txt
done

