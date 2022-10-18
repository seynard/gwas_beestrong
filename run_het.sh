#! /bin/bash
############################################ parameters ################################################## 
dir_in=${1}
cname=${2}
n_pop=${3}
snp50k=${4}
j=${5}
header="CHROM POS D X"

########################################################################################################## 
echo ${j}
i=$((${j}+1))
echo ${cname}
awk '{print $1" "$2" "$"'"${i}"'"}' ${dir_in}/depth.txt > ${dir_in}/tmp_depth${cname}
awk '{print $'${i}'}' ${dir_in}/count_ref.txt > ${dir_in}/tmp_count${cname}
paste -d' ' ${dir_in}/tmp_depth${cname} ${dir_in}/tmp_count${cname} > ${dir_in}/tmp${cname}.txt && mv ${dir_in}/tmp${cname}.txt ${dir_in}/sim_depth_count${cname}.txt 	
sed -i "1s/.*/${header}/" ${dir_in}/sim_depth_count${cname}.txt
#.local/bin/qg_pool --Fmatrix ${dir_in}/freq_admix${n_pop}_recode.txt ${dir_in}/sim_depth_count${cname}.txt -o ${dir_in}
grep -wFf ${dir_in}/${snp50k} ${dir_in}/sim_depth_count${cname}.txt > ${dir_in}/sim_depth_count${cname}_50k.txt
.local/bin/qg_pool --Fmatrix ${dir_in}/freq_admix${n_pop}_recode.txt ${dir_in}/sim_depth_count${cname}_50k.txt -o ${dir_in}
