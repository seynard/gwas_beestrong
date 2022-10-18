dir='results'
N=($(awk -F';' '{print $3","$16}' ${dir}/sign_locus.txt | uniq))
N=("${N[@]:1}") 
i=${N[59]}
echo ${i}
rs_snp=$(echo ${i} | cut -f1 -d',')
type=$(echo ${i} | cut -f2 -d',')
chr=$(echo ${rs_snp} | cut -f1 -d':')
zgrep -w "${rs_snp}" ${dir}/${type}_Mellifera_ncorse_${chr}.ld.gz > ${dir}/ld_snp${rs_snp}_${type}_Mellifera_ncorse.txt
sed 's/$/ Mellifera_ncorse/' ${dir}/ld_snp${rs_snp}_${type}_Mellifera_ncorse.txt > ${dir}/ld_snp${rs_snp}_${type}_Mellifera_ncorse.tmp && mv ${dir}/ld_snp${rs_snp}_${type}_Mellifera_ncorse.tmp ${dir}/ld_snp${rs_snp}_${type}_Mellifera_ncorse.txt
zgrep -w "${rs_snp}" ${dir}/${type}_hybrid_corse_${chr}.ld.gz > ${dir}/ld_snp${rs_snp}_${type}_hybrid_corse.txt
sed 's/$/ hybrid_corse/' ${dir}/ld_snp${rs_snp}_${type}_hybrid_corse.txt > ${dir}/ld_snp${rs_snp}_${type}_hybrid_corse.tmp && mv ${dir}/ld_snp${rs_snp}_${type}_hybrid_corse.tmp ${dir}/ld_snp${rs_snp}_${type}_hybrid_corse.txt
cat ${dir}/ld_snp${rs_snp}_${type}_Mellifera_ncorse.txt ${dir}/ld_snp${rs_snp}_${type}_hybrid_corse.txt > ${dir}/ld_snp${rs_snp}_${type}.txt
rm ${dir}/ld_snp${rs_snp}_${type}_Mellifera_ncorse.txt ${dir}/ld_snp${rs_snp}_${type}_hybrid_corse.txt

dir='results'
N=($(awk -F';' '{print $3","$16}' ${dir}/sign_locus.txt | uniq))
N=("${N[@]:1}") 
i=${N[77]}
echo ${i}
rs_snp=$(echo ${i} | cut -f1 -d',')
type=$(echo ${i} | cut -f2 -d',')
chr=$(echo ${rs_snp} | cut -f1 -d':')
zgrep -w "${rs_snp}" ${dir}/${type}_Mellifera_ncorse_${chr}.ld.gz > ${dir}/ld_snp${rs_snp}_${type}_Mellifera_ncorse.txt
sed 's/$/ Mellifera_ncorse/' ${dir}/ld_snp${rs_snp}_${type}_Mellifera_ncorse.txt > ${dir}/ld_snp${rs_snp}_${type}_Mellifera_ncorse.tmp && mv ${dir}/ld_snp${rs_snp}_${type}_Mellifera_ncorse.tmp ${dir}/ld_snp${rs_snp}_${type}_Mellifera_ncorse.txt
zgrep -w "${rs_snp}" ${dir}/${type}_hybrid_corse_${chr}.ld.gz > ${dir}/ld_snp${rs_snp}_${type}_hybrid_corse.txt
sed 's/$/ hybrid_corse/' ${dir}/ld_snp${rs_snp}_${type}_hybrid_corse.txt > ${dir}/ld_snp${rs_snp}_${type}_hybrid_corse.tmp && mv ${dir}/ld_snp${rs_snp}_${type}_hybrid_corse.tmp ${dir}/ld_snp${rs_snp}_${type}_hybrid_corse.txt
cat ${dir}/ld_snp${rs_snp}_${type}_Mellifera_ncorse.txt ${dir}/ld_snp${rs_snp}_${type}_hybrid_corse.txt > ${dir}/ld_snp${rs_snp}_${type}.txt
rm ${dir}/ld_snp${rs_snp}_${type}_Mellifera_ncorse.txt ${dir}/ld_snp${rs_snp}_${type}_hybrid_corse.txt

dir='results'
N=($(awk -F';' '{print $3","$16}' ${dir}/sign_locus.txt | uniq))
N=("${N[@]:1}") 
i=${N[199]}
echo ${i}
rs_snp=$(echo ${i} | cut -f1 -d',')
type=$(echo ${i} | cut -f2 -d',')
chr=$(echo ${rs_snp} | cut -f1 -d':')
zgrep -w "${rs_snp}" ${dir}/${type}_Mellifera_ncorse_${chr}.ld.gz > ${dir}/ld_snp${rs_snp}_${type}_Mellifera_ncorse.txt
sed 's/$/ Mellifera_ncorse/' ${dir}/ld_snp${rs_snp}_${type}_Mellifera_ncorse.txt > ${dir}/ld_snp${rs_snp}_${type}_Mellifera_ncorse.tmp && mv ${dir}/ld_snp${rs_snp}_${type}_Mellifera_ncorse.tmp ${dir}/ld_snp${rs_snp}_${type}_Mellifera_ncorse.txt
zgrep -w "${rs_snp}" ${dir}/${type}_hybrid_corse_${chr}.ld.gz > ${dir}/ld_snp${rs_snp}_${type}_hybrid_corse.txt
sed 's/$/ hybrid_corse/' ${dir}/ld_snp${rs_snp}_${type}_hybrid_corse.txt > ${dir}/ld_snp${rs_snp}_${type}_hybrid_corse.tmp && mv ${dir}/ld_snp${rs_snp}_${type}_hybrid_corse.tmp ${dir}/ld_snp${rs_snp}_${type}_hybrid_corse.txt
cat ${dir}/ld_snp${rs_snp}_${type}_Mellifera_ncorse.txt ${dir}/ld_snp${rs_snp}_${type}_hybrid_corse.txt > ${dir}/ld_snp${rs_snp}_${type}.txt
rm ${dir}/ld_snp${rs_snp}_${type}_Mellifera_ncorse.txt ${dir}/ld_snp${rs_snp}_${type}_hybrid_corse.txt
