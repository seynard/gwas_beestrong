cd /work/genphyse/dynagen/seynard/GWAS/infestation4
wc -l data/depth.txt
#6831075 depth.txt
awk '{print NF;exit}' data/depth.txt
#1514
old_IFS=$IFS

echo -n > average_depth_col.txt
line="col_id nb_snp nb_snp_depth0 nb_snp_depth_non0 ave_depth min_depth max_depth sd_depth"
echo ${line} >> average_depth_col.txt

ncol=$(awk '{print NF;exit}' data/depth.txt)
for ((i=5;i<=$ncol;i++))
do
echo $i
sbatch --mem=20G --wrap="sh scripts/seq_depth_summary_col_unit.sh ${i}"
#ID=$(awk -v col=${i} '{print $col}' data/data_extra/depth_all_nus.txt | head -n1)
#col=($(awk -v col=${i} '{print $col}' data/data_extra/depth_all_nus.txt | tail -n +2))
#nb_0=$(grep -o ' 0' <<< ${col[*]} | wc -l)
#delete=0
#col_n0=(${col[@]/$delete})
#IFS='+' 
#avg=$(echo "scale=3;(${col_n0[*]})/${#col_n0[@]}"|bc)
#IFS=$old_IFS
#min=$(printf "%d\n" "${col_n0[@]}" | sort -rn | tail -1)
#max=$(printf "%d\n" "${col_n0[@]}" | sort -rn | head -1)
#strdev=$(echo ${col_n0[@]} | awk 'NF {sum=0;ssq=0;for (i=1;i<=NF;i++){sum+=$i;ssq+=$i**2};print (ssq/NF-(sum/NF)**2)**0.5}')
#line=$(echo ${ID} ${#col[@]} ${nb_0} ${#col_n0[@]} ${avg} ${min} ${max} ${strdev})
#echo ${line}
#echo ${line} >> average_depth_col.txt
done 

#cp data/data_extra/depth_all_nus.txt data/data_extra/depth_all_nus_no_idcol.txt
#awk '{$1=$2=$3=$4=""; print $0}' data/data_extra/depth_all_nus_no_idcol.txt > data/data_extra/depth_all_nus_no_idcol.tmp && mv data/data_extra/depth_all_nus_no_idcol.tmp data/data_extra/depth_all_nus_no_idcol.txt
#sed -i 1d data/data_extra/depth_all_nus_no_idcol.txt
#echo -n > average_depth_snp.txt
#nsnp=$(wc -l < data/data_extra/depth_all_nus_no_idcol.txt)
#chr=($(awk '{print $1}' data/data_extra/depth_all_nus.txt|tail -n +2))
#pos=($(awk '{print $2}' data/data_extra/depth_all_nus.txt|tail -n +2))
#for j in $(seq 1 1 ${nsnp})
#do
#echo $j
#snp=$(echo ${chr[j]}':'${pos[j]})
#row=($(sed "${j}q;d" data/data_extra/depth_all_nus_no_idcol.txt))
#nb_0=$(grep -o ' 0' <<< ${row[*]} | wc -l)
#delete=0
#row_n0=(${row[@]/$delete})
#IFS='+' 
#avg=$(echo "scale=3;(${row_n0[*]})/${#row_n0[@]}"|bc)
#IFS=$old_IFS
#line=$(echo ${snp} ${#row[@]} ${nb_0} ${#row_n0[@]} ${avg})
#echo ${line}
#echo ${line} >> average_depth_snp.txt
#done
#rm data/data_extra/depth_all_nus_no_idcol.txt
