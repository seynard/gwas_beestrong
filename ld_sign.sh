#! /bin/bash
dir=${1}
rs_snp=${2}
chr=${3}
zgrep -w "${rs_snp}" ${dir}/ld_Mellifera_${chr}.ld > ${dir}/ld_snp${rs_snp}_Mellifera.txt
sed 's/$/ Mellifera/' ${dir}/ld_snp${rs_snp}_Mellifera.txt > ${dir}/ld_snp${rs_snp}_Mellifera.tmp && mv ${dir}/ld_snp${rs_snp}_Mellifera.tmp ${dir}/ld_snp${rs_snp}_Mellifera.txt
zgrep "${rs_snp}" ${dir}/ld_Ligustica_Carnica_${chr}.ld > ${dir}/ld_snp${rs_snp}_Ligustica_Carnica.txt
sed 's/$/ Ligustica_Carnica/' ${dir}/ld_snp${rs_snp}_Ligustica_Carnica.txt > ${dir}/ld_snp${rs_snp}_Ligustica_Carnica.tmp && mv ${dir}/ld_snp${rs_snp}_Ligustica_Carnica.tmp ${dir}/ld_snp${rs_snp}_Ligustica_Carnica.txt
zgrep -w "${rs_snp}" ${dir}/ld_hybrid_${chr}.ld > ${dir}/ld_snp${rs_snp}_hybrid.txt
sed 's/$/ hybrid/' ${dir}/ld_snp${rs_snp}_hybrid.txt > ${dir}/ld_snp${rs_snp}_hybrid.tmp && mv ${dir}/ld_snp${rs_snp}_hybrid.tmp ${dir}/ld_snp${rs_snp}_hybrid.txt
cat ${dir}/ld_snp${rs_snp}_Mellifera.txt ${dir}/ld_snp${rs_snp}_Ligustica_Carnica.txt ${dir}/ld_snp${rs_snp}_hybrid.txt > ${dir}/ld_snp${rs_snp}.txt
rm ${dir}/ld_snp${rs_snp}_Mellifera.txt ${dir}/ld_snp${rs_snp}_Ligustica_Carnica.txt ${dir}/ld_snp${rs_snp}_hybrid.txt

