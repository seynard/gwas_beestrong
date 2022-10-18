#! /bin/bash
############################################ parameters ################################################## 
i=${1}
echo $i
ID=$(awk -v col=${i} '{print $col}' data/depth.txt | head -n1)
col=($(awk -v col=${i} '{print $col}' data/depth.txt | tail -n +2))
nb_0=$(grep -o ' 0' <<< ${col[*]} | wc -l)
delete=0
col_n0=(${col[@]/$delete})
IFS='+' 
avg=$(echo "scale=3;(${col_n0[*]})/${#col_n0[@]}"|bc)
IFS=$old_IFS
min=$(printf "%d\n" "${col_n0[@]}" | sort -rn | tail -1)
max=$(printf "%d\n" "${col_n0[@]}" | sort -rn | head -1)
strdev=$(echo ${col_n0[@]} | awk 'NF {sum=0;ssq=0;for (i=1;i<=NF;i++){sum+=$i;ssq+=$i**2};print (ssq/NF-(sum/NF)**2)**0.5}')
line=$(echo ${ID} ${#col[@]} ${nb_0} ${#col_n0[@]} ${avg} ${min} ${max} ${strdev})
#echo ${line}
echo ${line} >> average_depth_col.txt
