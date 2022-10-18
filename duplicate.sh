#! /bin/bash
dir_out=$1
opt=$2
l_depth=$3
depth_in=$4
count_in=$5

opt=($(echo ${opt} | tr "," "\n"))
l_depth=($(echo ${l_depth} | tr "," "\n"))

for i in ${opt[@]}
do
unset index_col
index_col=($(printf '%s\n' "${l_depth[@]}" | grep ${i}))
awk -F ' ' -v col=${index_col[0]} 'NR==1{for (i=1; i<=NF; i++) if ($i == col){c=i; break}} {print $c}' ${depth_in} > ${dir_out}/depth_0.txt
awk -F ' ' -v col=${index_col[1]} 'NR==1{for (i=1; i<=NF; i++) if ($i == col){c=i; break}} {print $c}' ${depth_in} > ${dir_out}/depth_1.txt
paste ${dir_out}/depth_0.txt ${dir_out}/depth_1.txt | awk '{ print $1 + $2; }' > ${dir_out}/depth_sum.txt
sed -i "1s/.*/${index_col[1]}/" ${dir_out}/depth_sum.txt
awk -F ' ' -v col=${index_col[0]} 'NR==1{for (i=1; i<=NF; i++) if ($i == col){c=i; break}} {print $c}' ${count_in} > ${dir_out}/count_ref_0.txt
awk -F ' ' -v col=${index_col[1]} 'NR==1{for (i=1; i<=NF; i++) if ($i == col){c=i; break}} {print $c}' ${count_in} > ${dir_out}/count_ref_1.txt
paste ${dir_out}/count_ref_0.txt ${dir_out}/count_ref_1.txt | awk '{ print $1 + $2; }' > ${dir_out}/count_ref_sum.txt;
sed -i "1s/.*/${index_col[1]}/" ${dir_out}/count_ref_sum.txt
unset dup_index
for j in ${index_col[@]}
do
for i in "${!l_depth[@]}"
do
   [[ "${l_depth[$i]}" = "${j}" ]] && break
done
dup_index+=$((i+1))','
done
dup_index=${dup_index::-1}
cut -d' ' -f${dup_index} --complement ${depth_in} > ${dir_out}/tmp && mv ${dir_out}/tmp ${depth_in}; cut -d' ' -f${dup_index} --complement ${count_in} > ${dir_out}/tmp && mv ${dir_out}/tmp ${count_in}
paste ${depth_in} ${dir_out}/depth_sum.txt > ${dir_out}/tmp && mv ${dir_out}/tmp ${depth_in}; paste ${count_in} ${dir_out}/count_ref_sum.txt > ${dir_out}/tmp && mv ${dir_out}/tmp ${count_in}
done

