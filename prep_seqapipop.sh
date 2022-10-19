#! /bin/bash
##################################################
#################### PARAMETERS #################### 
script=${1}
dir_in=${2}
vcf_file=${3}
n_pop=${4}
pop_id=${5}
unif_threshold=${6}
ncpu=${7}
##################################################
module load bioinfo/bcftools-1.6
module load bioinfo/plink-v1.90b5.3
module load system/R-4.1.1_gcc-9.3.0

##########################################################################################################
########################## definition of known subspecies from SeqApiPop panel ###########################
mkdir -p ${dir_in}/seqapipop
cp ${vcf_file}.gz ${dir_in}/seqapipop/vcf_seqapipop.vcf.gz; gunzip -d ${dir_in}/seqapipop/vcf_seqapipop.vcf.gz
bcftools query -f "%CHROM %POS \n" ${dir_in}/seqapipop/vcf_seqapipop.vcf > ${dir_in}/seqapipop/snp_pos.txt
chr=($(cut -f1 -d " " ${dir_in}/seqapipop/snp_pos.txt | sort | uniq ))
chr=("${chr[@]:1}")
sed -i "/NC_001566.1/d" ${dir_in}/seqapipop/vcf_seqapipop.vcf
sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dir_in}/seqapipop/vcf_seqapipop.vcf
bcftools annotate --set-id +"%CHROM:%POS" ${dir_in}/seqapipop/vcf_seqapipop.vcf > ${dir_in}/seqapipop/tmp.vcf; mv ${dir_in}/seqapipop/tmp.vcf ${dir_in}/seqapipop/vcf_seqapipop.vcf
plink --vcf ${dir_in}/seqapipop/vcf_seqapipop.vcf --keep-allele-order --make-bed --out ${dir_in}/seqapipop/seqapipop

# missing individuals
plink --bfile ${dir_in}/seqapipop/seqapipop --missing  --out ${dir_in}/seqapipop/seqapipop
nchr=${#chr[@]}
for ((i=1; i<=${nchr};i++))
do
plink --bfile ${dir_in}/seqapipop/seqapipop --make-bed --chr ${i} --out ${dir_in}/seqapipop/seqapipop_${i}
plink --bfile ${dir_in}/seqapipop/seqapipop_${i} --missing --out ${dir_in}/seqapipop/missing_seqapipop_${i}
done

# keep individuals from list seqapipop
awk "{print $1" "$1}" ${dir_in}/seqapipop/list_seqapipop.txt > ${dir_in}/seqapipop/sample_seqapipop.txt 
plink --bfile ${dir_in}/seqapipop/seqapipop --keep ${dir_in}/seqapipop/sample_seqapipop.txt --keep-allele-order --make-bed --out ${dir_in}/seqapipop/seqapipop_sample
K=${n_pop}
cd ${dir_in}/seqapipop	
admixture --cv -j${ncpu} seqapipop_sample.bed ${K}
cd
Rscript ${script}/admix_plot.r ${dir_in}/seqapipop seqapipop_sample ${pop_id}

# definition of sub-species: mellifera, caucasia, ligustica/carnica 
Rscript ${script}/prep_admix_k.r ${dir_in}/seqapipop seqapipop_sample ${pop_id} ${unif_threshold}
sample=($(awk "{print $1}" ${dir_in}/seqapipop/Unif_k.fam))
sample2=""
for i in ${sample[@]}
do
sample2+=","${i}
done
sample2="${sample2:1}"
bcftools view -s ${sample2} ${vcf_file} > ${dir_in}/seqapipop/vcf_subsp.vcf
IFS=", " read -r -a popID <<< "${pop_id}"
for i in ${popID[@]}
do
samplei=$(awk -v var=${i} -F" " "{if($1==var)print $2}" ${dir_in}/seqapipop/Unif_k_pop.fam )
samplei=$(echo $samplei | sed "y/ /,/")
bcftools view -s ${samplei} ${vcf_file} > ${dir_in}/vcf_${i}.vcf
sed -i "/NC_001566.1/d" ${dir_in}/seqapipop/vcf_${i}.vcf
cp ${dir_in}/seqapipop/vcf_${i}.vcf ${dir_in}/seqapipop/vcf_${i}_recode.vcf
bcftools annotate --set-id +"%CHROM:%POS" ${dir_in}/seqapipop/vcf_${i}_recode.vcf > ${dir_in}/seqapipop/tmp.vcf; mv ${dir_in}/seqapipop/tmp.vcf ${dir_in}/seqapipop/vcf_${i}_recode.vcf
sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dir_in}/seqapipop/vcf_${i}_recode.vcf
bgzip -c ${dir_in}/seqapipop/vcf_${i}.vcf > ${dir_in}/seqapipop/vcf_${i}.vcf.gz; bgzip -c ${dir_in}/seqapipop/vcf_${i}_recode.vcf > ${dir_in}/seqapipop/vcf_${i}_recode.vcf.gz
done
awk "{print $1}" ${dir_in}/seqapipop/Unif_k_pop.fam > ${dir_in}/seqapipop/vcf_subsp.pop
cp ${dir_in}/seqapipop/seqapipop_sample.${n_pop}.P ${dir_in}/freq_admix${n_pop}.txt
paste --delimiter=" " ${dir_in}/seqapipop/snp_pos.txt ${dir_in}/freq_admix${n_pop}.txt > ${dir_in}/temp.txt && mv ${dir_in}/temp.txt ${dir_in}/freq_admix${n_pop}.txt
col_id=$(cat ${dir_in}/ColNames.txt)
header="CHROM POS "${col_id}
header=$(echo ${header} | sed "s/\n/ /g")
awk -v v="$header" "NR==1{print v}1" ${dir_in}/freq_admix${n_pop}.txt > ${dir_in}/temp.txt && mv ${dir_in}/temp.txt ${dir_in}/freq_admix${n_pop}.txt
Rscript ${script}/freq_plot.r ${dir_in} freq_admix${n_pop}.txt
sed -i "/NC_001566.1/d" ${dir_in}/seqapipop/vcf_subsp.vcf
cp ${dir_in}/seqapipop/vcf_subsp.vcf ${dir_in}/seqapipop/vcf_subsp_recode.vcf
sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dir_in}/seqapipop/vcf_subsp_recode.vcf
sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dir_in}/freq_admix${n_pop}.txt

bgzip -c ${dir_in}/seqapipop/vcf_subsp.vcf > ${dir_in}/seqapipop/vcf_subsp.vcf.gz

