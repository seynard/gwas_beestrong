#! /bin/bash
##################################################
#################### PARAMETERS #################### 
script=${1}
dir_in=${2}
dir_popoolation="/usr/local/bioinfo/src/PoPoolation2/popoolation2_1201"
fasta=${3}
vcf_name=${4}
vcf_sansfiltre=${5}
vcf_file=${6}
snp50k=${7}
nbjobs=${8}
d_min=${9}
d_max=${10}
pool_size_def=${11}
##################################################

##########################################################################################################
############################## from BAM to count, depth ##################################
# colony by colony create depth and count files
mkdir -p ${dir_in}/pileup
mkdir -p ${dir_in}/count
mkdri -p ${dir_in}/depth_count

gunzip -d ${vcf_file}.gz
if [ -e ${dir_in}/depth_count/depth.txt.bz2 ]
then
	module load bioinfo/bcftools-1.6; bcftools query -f "%CHROM %POS \n" ${vcf_file} > ${dir_in}/snp_pos.txt; bcftools query -f "%CHROM %POS %REF %ALT \n" ${vcf_file} > ${dir_in}/allele_id.txt
	bzip2 -d ${dir_in}/depth_count/depth.txt.bz2 -d ${dir_in}
	bzip2 -d ${dir_in}/depth_count/count_ref.txt.bz2 -d ${dir_in}
	bzip2 -d ${dir_in}/depth_count/count_alt.txt.bz2 -d ${dir_in}
else
	module load bioinfo/bcftools-1.6; bcftools query -f "%CHROM %POS \n" ${vcf_file} > ${dir_in}/snp_pos.txt; bcftools query -f "%CHROM %POS %REF %ALT \n" ${vcf_file} > ${dir_in}/allele_id.txt; cp ${dir_in}/allele_id.txt ${dir_in}/depth.txt
	header="CHROM POS REF ALT "
	awk -v x="${header}" "NR==1{print x} 1" ${dir_in}/depth.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/depth.txt; cp ${dir_in}/depth.txt ${dir_in}/count_ref.txt; cp ${dir_in}/depth.txt ${dir_in}/count_alt.txt
fi
sh ${script}/file_init.sh ${dir_in} ${dir_in} ${script} BamList BS
grep -v "done" ${dir_in}/BamList > ${dir_in}/Bam_not_done
bam_list2=($(awk -F" " "{print $1}" ${dir_in}/Bam_not_done))
rm ${script}/pileup_array2*
x=1
for i in ${bam_list2[@]} 
do
n=$(awk -F " " "NR==num{print $2}" num=${x} ${dir_in}/Bam_not_done)
x=$[$x + 1]
echo "sh "${script}"/make_pileup.sh "${dir_in}" "${dir_popoolation}" "${fasta}" "${vcf_file}" "${script}" "${i}" "${n}";" >> ${script}/pileup_array2.sh
done
split -l ${nbjobs} ${script}/pileup_array2.sh ${script}/pileup_array2_split_
pileup_list=(${script}/pileup_array2_split_*)
for i in ${pileup_list[@]}
do
echo -e "#! /bin/bash\n$(cat ${i})" > ${i}
sbatch ${i}
done
n_job_run=$(squeue -u seynard | grep "pileup" | wc -l)
while [ ${n_job_run} != 0 ]
do
sleep 10m 
echo ${n_job_run}
n_job_run=$(squeue -u seynard | grep "pileup" | wc -l)
done
mv ${dir_in}/*pileup.bz2 ${dir_in}/pileup

count_list=(${dir_in}/*.count)
last=$(echo ${count_list[@]:(-1)}|rev|cut -d"/" -f 1|rev | cut -d"." -f 1)
for i in ${count_list[@]} 
do
n=($(echo ${i}|rev| cut -d "/" -f 1|rev | cut -d"." -f 1))
if grep ${i} ${dir_in}/depth.txt
then
echo "already done"
else
echo "to do"
if [ $n != $last ] 
then
sbatch -J synctocount -o ${dir_in}/log/synctocount_${n}.out -e ${dir_in}/log/synctocount_${n}.err --wrap="python ${script}/synctocount.py ${i} ${dir_in}"
else
sbatch -W -J synctocount -o ${dir_in}/log/synctocount_${n}.out -e ${dir_in}/log/synctocount_${n}.err --wrap="python ${script}/synctocount.py ${i} ${dir_in}"
fi
fi
done
for i in ${count_list[@]} 
do
n=($(echo ${i}|rev| cut -d "/" -f 1|rev | cut -d"." -f 1))
bzip2 ${dir_in}/${n}.count
done
mv ${dir_in}/*.count.bz2 ${dir_in}/count/

ls -1 ${dir_in}/tmp_depth_* | split -l 100 -d - ${dir_in}/lists
cd ${dir_in}
for list in lists*; do paste $(cat $list) > merge${list##lists}; done
cd
paste ${dir_in}/merge* > ${dir_in}/depth.tmp; paste -d " " ${dir_in}/depth.txt ${dir_in}/depth.tmp > ${dir_in}/depth2.txt; mv ${dir_in}/depth2.txt ${dir_in}/depth.txt; sed -i "s/  / /g" ${dir_in}/depth.txt; sed -i "s/\t/ /g" ${dir_in}/depth.txt
rm ${dir_in}/merge* ${dir_in}/lists* ${dir_in}/depth.tmp
ls -1 ${dir_in}/tmp_count_ref_* | split -l 100 -d - ${dir_in}/lists
cd ${dir_in}
for list in lists*; do paste $(cat $list) > merge${list##lists}; done
cd
paste ${dir_in}/merge* > ${dir_in}/count_ref.tmp; paste -d " " ${dir_in}/count_ref.txt ${dir_in}/count_ref.tmp > ${dir_in}/count_ref2.txt; mv ${dir_in}/count_ref2.txt ${dir_in}/count_ref.txt; sed -i "s/  / /g" ${dir_in}/count_ref.txt; sed -i "s/\t/ /g" ${dir_in}/count_ref.txt
rm ${dir_in}/merge* ${dir_in}/lists* ${dir_in}/count_ref.tmp
ls -1 ${dir_in}/tmp_count_alt_* | split -l 100 -d - ${dir_in}/lists
cd ${dir_in}
for list in lists*; do paste $(cat $list) > merge${list##lists}; done
cd
paste ${dir_in}/merge* > ${dir_in}/count_alt.tmp; paste -d " " ${dir_in}/count_alt.txt ${dir_in}/count_alt.tmp > ${dir_in}/count_alt2.txt; mv ${dir_in}/count_alt2.txt ${dir_in}/count_alt.txt; sed -i "s/  / /g" ${dir_in}/count_alt.txt; sed -i "s/\t/ /g" ${dir_in}/count_alt.txt
rm ${dir_in}/merge* ${dir_in}/lists* ${dir_in}/count_alt.tmp ${dir_in}/tmp_depth_*.txt ${dir_in}/tmp_count_ref_*.txt ${dir_in}/tmp_count_alt_*.txt

# all colonies together
echo "number of initial SNP " $(($(wc -l ${dir_in}/depth.txt| awk "{print $1}")-1)) #number of initial SNP  7 026 976

# filter on markers
python ${script}/non_mono.py ${dir_in}/depth.txt ${dir_in}/count_ref.txt ${dir_in}/count_ref2.txt
python ${script}/non_mono.py ${dir_in}/depth.txt ${dir_in}/count_ref.txt ${dir_in}/depth2.txt
echo "number of monoallelic variants " $(wc -l ${dir_in}/log/non_mono_depth.out| awk "{print $1}") #number of monoallelic variants 24 337
echo "number of SNP kept " $(($(wc -l ${dir_in}/depth2.txt| awk "{print $1}")-1)) #number of SNP kept 6 999 639
echo "number of SNP on mitochondria " $(grep "NC_001566.1 " ${dir_in}/depth.txt | wc -l) #number of SNP on mitochondrial dna (before filters) 287
echo "number of SNP monoallelic on mitochondria" $(grep "NC_001566.1 " ${dir_in}/log/non_mono_depth.out | wc -l) #number of SNP monoallelic on mitochondrial dna 0
sed -i "/NC_001566.1/d" ${dir_in}/depth2.txt; sed -i "/NC_001566.1/d" ${dir_in}/count_ref2.txt #remove mitochondrial DNA
echo "number of SNP kept " $(($(wc -l ${dir_in}/depth2.txt| awk "{print $1}")-1)) #number of SNP kept 6 999 352

# remove real tri_allelic
grep "," ${dir_in}/depth2.txt | cut -d" " -f1,2 > ${dir_in}/alt.txt
grep -wFf ${dir_in}/alt.txt ${dir_in}/count_alt.txt > ${dir_in}/tmp_alt.txt && mv ${dir_in}/tmp_alt.txt ${dir_in}/alt.txt
sed -i "s/, /,/g" ${dir_in}/alt.txt
python ${script}/tri_allelic.py ${dir_in}/alt.txt ${dir_in}/triallelic_remove.txt remove
python ${script}/tri_allelic.py ${dir_in}/alt.txt ${dir_in}/triallelic_recode.txt recode
awk "!a[$0]++" ${dir_in}/triallelic_remove.txt > ${dir_in}/triallelic.tmp && mv ${dir_in}/triallelic.tmp ${dir_in}/triallelic_remove.txt
awk "!a[$0]++" ${dir_in}/triallelic_recode.txt > ${dir_in}/triallelic.tmp && mv ${dir_in}/triallelic.tmp ${dir_in}/triallelic_recode.txt
echo "number of tri-allelic removed variants " $(wc -l ${dir_in}/triallelic_remove.txt| awk "{print $1}") #number of tri-allelic variants removed 17 180
echo "number of tri-allelic recoded variants " $(wc -l ${dir_in}/triallelic_recode.txt| awk "{print $1}") #number of tri-allelic variants recoded 91 736
grep -vwFf ${dir_in}/triallelic_remove.txt ${dir_in}/depth2.txt > ${dir_in}/tmp.txt && mv ${dir_in}/tmp.txt ${dir_in}/depth2.txt; grep -vwFf ${dir_in}/triallelic_remove.txt ${dir_in}/count_ref2.txt > ${dir_in}/tmp.txt && mv ${dir_in}/tmp.txt ${dir_in}/count_ref2.txt
echo "number of SNP kept " $(($(wc -l ${dir_in}/depth2.txt| awk "{print $1}")-1)) #number of SNP kept 6 982 172 
cut -d" " -f1,2,3,4 ${dir_in}/depth2.txt > ${dir_in}/snp_kept.txt
Rscript ${script}/recode_triallelic.r ${dir_in}
sed -i "1s/.*/ALT/" ${dir_in}/snp_kept_recode.txt
awk "FNR==NR{a[NR]=$0;next} {sub($4,a[FNR])}1" ${dir_in}/snp_kept_recode.txt ${dir_in}/depth2.txt > ${dir_in}/tmp.txt && mv ${dir_in}/tmp.txt ${dir_in}/depth2.txt
awk "FNR==NR{a[NR]=$0;next} {sub($4,a[FNR])}1" ${dir_in}/snp_kept_recode.txt ${dir_in}/count_ref2.txt > ${dir_in}/tmp.txt && mv ${dir_in}/tmp.txt ${dir_in}/count_ref2.txt
rm ${dir_in}/snp_kept.txt ${dir_in}/snp_kept_recode.txt ${dir_in}/alt.txt
bzip2 ${dir_in}/triallelic_remove.txt; bzip2 ${dir_in}/triallelic_recode.txt

# filter on SNP depth and monomorphy
python ${script}/depth_control.py ${dir_in}/depth2.txt ${dir_in}/count_ref2.txt ${dir_in}/depth_control.txt
Rscript ${script}/depth_control_plot.r ${dir_in} ${d_min} ${d_max}
grep -vwFf ${dir_in}/to_remove_depth.txt ${dir_in}/depth2.txt > ${dir_in}/depth3.txt; grep -vwFf ${dir_in}/to_remove_depth.txt ${dir_in}/count_ref2.txt > ${dir_in}/count_ref3.txt
echo "number of SNP kept " $(($(wc -l ${dir_in}/depth3.txt| awk "{print $1}")-1)) #number of SNP kept 6 831 074
mv ${dir_in}/{depth_control*.png,to_remove_depth.txt,depth_control.txt} ${dir_in}/depth_count
awk "{print $1" "$2}" ${dir_in}/depth3.txt > ${dir_in}/snp_final.txt

# check and remove duplicated colonies
# 1) positive control on replicated colony ITSAP
sh ${script}/run_replicate.sh ${script} ${dir_in} ${dir_in} ${dir_popoolation} ${fasta} ${vcf_sansfiltre} ${vcf_name} ${pool_size_def} ${nbjobs}
grep -wFf ${dir_in}/snp_final.txt ${dir_in}/depth_replicate.txt > ${dir_in}/tmp_replicate && mv ${dir_in}/tmp_replicate ${dir_in}/depth_replicate.txt
grep -wFf ${dir_in}/snp_final.txt ${dir_in}/count_ref_replicate.txt > ${dir_in}/tmp_replicate && mv ${dir_in}/tmp_replicate ${dir_in}/count_ref_replicate.txt
grep -wFf ${dir_in}/snp_final.txt ${dir_in}/count_alt_replicate.txt > ${dir_in}/tmp_replicate && mv ${dir_in}/tmp_replicate ${dir_in}/count_alt_replicate.txt
# 2) random colonies from BeeStrong
ncol=$(awk "{print NF;exit}" ${dir_in}/depth3.txt)
N=($(shuf -i 5-${ncol} -n 10))
n=$(printf ",%s" "${N[@]}")
n=${n:1}
cut -d" " -f${n} ${dir_in}/depth3.txt > ${dir_in}/depth_control_correl.txt; cut -d" " -f${n} ${dir_in}/count_ref3.txt > ${dir_in}/count_ref_control_correl.txt
# 3) BeeStrong duplicated
l_depth=($(head -n1 ${dir_in}/depth3.txt))
dup_col=$(printf "%s\n" "${l_depth[@]}" | grep "_bis")
dup_col=($(echo "${dup_col//_bis}"))
unset dup_index
for j in ${dup_col[@]}
do
for i in "${!l_depth[@]}"
do
   [[ "${l_depth[$i]}" = "${j}" ]] && break
done
dup_index+=$((i+1))","
for i in "${!l_depth[@]}"
do
   [[ "${l_depth[$i]}" = "${j}_bis" ]] && break
done
dup_index+=$((i+1))","
done
dup_index=${dup_index::-1}
cut -d" " -f${dup_index} ${dir_in}/depth3.txt > ${dir_in}/depth_dup.txt; cut -d" " -f${dup_index} ${dir_in}/count_ref3.txt > ${dir_in}/count_ref_dup.txt
# 4) replicate (positive control), BeeStrong (negative control between colonies, positive control intra-colony), choice: duplicated sequences
Rscript ${script}/duplicate.r ${dir_in} replicate_random_duplicate
mv ${dir_in}/{depth_dup.txt,count_ref_dup.txt,depth_replicate.txt,count_ref_replicate.txt,depth_control_correl.txt,count_ref_control_correl.txt} ${dir_in}/depth_count
cp ${dir_in}/depth3.txt ${dir_in}/depth.tmp
cp ${dir_in}/count_ref3.txt ${dir_in}/count_ref.tmp
opt_dup=()
while IFS= read -r line
do
   opt_dup+=("$line")
done < ${dir_in}/opt_dup.txt
for x in "${opt_dup[@]}"
do
echo $x
n=($(echo ${x}| cut -d " " -f 2))
l_depth=($(head -n1 ${dir_in}/depth.tmp))
unset dup_index
for i in "${!l_depth[@]}"
do
   [[ "${l_depth[$i]}" = "${n}" ]] && break
done
dup_index+=$((i+1))","
for i in "${!l_depth[@]}"
do
   [[ "${l_depth[$i]}" = "${n}_bis" ]] && break
done
dup_index+=$((i+1))","
dup_index=${dup_index::-1}
if [[ ${x} == *"remove"* ]]
then
echo "to remove"
cut -d" " -f${dup_index} --complement ${dir_in}/depth.tmp > ${dir_in}/depth4.txt; cut -d" " -f${dup_index} --complement ${dir_in}/count_ref.tmp > ${dir_in}/count_ref4.txt
else 
echo "to recode"
opt=${opt_dup//recode/}
opt=($(echo ${opt} | tr " " ","))
l_depth=$(head -n1 ${dir_in}/depth.tmp)
l_depth=($(echo ${l_depth} | tr " " ","))
sh ${script}/duplicate.sh ${dir_in} ${opt} ${l_depth} ${dir_in}/depth.tmp ${dir_in}/count_ref.tmp
sed -i "s/  / /g" ${dir_in}/depth.tmp; sed -i "s/\t/ /g" ${dir_in}/depth.tmp
sed -i "s/  / /g" ${dir_in}/count_ref.tmp; sed -i "s/\t/ /g" ${dir_in}/count_ref.tmp
mv ${dir_in}/depth.tmp ${dir_in}/depth4.txt; mv ${dir_in}/count_ref.tmp ${dir_in}/count_ref4.txt
rm ${dir_in}/depth_0.txt ${dir_in}/depth_1.txt ${dir_in}/depth_sum.txt ${dir_in}/count_ref_0.txt ${dir_in}/count_ref_1.txt ${dir_in}/count_ref_sum.txt
fi
done
echo "number of final SNP " $(($(wc -l ${dir_in}/depth4.txt| awk "{print $1}")-1)) #number of final SNP  6 831 074
echo "number of colony " $(($(awk "{print NF;exit}" ${dir_in}/depth4.txt))) #number of colony 1637
rm ${dir_in}/depth.tmp ${dir_in}/count_ref.tmp
mv ${dir_in}/{depth2.txt,count_ref2.txt,depth3.txt,count_ref3.txt,depth.txt,count_ref.txt,count_alt.txt} ${dir_in}/depth_count

bzip2 ${dir_in}/depth_count/*
mv ${dir_in}/depth4.txt ${dir_in}/depth_final.txt
mv ${dir_in}/count_ref4.txt ${dir_in}/count_ref_final.txt

# rename colonies
l_depth=($(head -n1 ${dir_in}/depth_final.txt))
l_depth=(${l_depth[@]:4})
unset l_depth2
for i in ${l_depth[@]}
do
j=($(tr "-" "_" <<<"${i}"))
j=($(echo $j | sed -e "s/\(_F\)*$//g"))
bs=($(echo $j | cut -d "_" -f1))
num=($(echo $j | cut -d "_" -f2))
lnum=($(echo ${#num}))
if [ $lnum -lt 4 ]
then
toadd=$((4-$lnum))
start=1
end=$toadd
zadd=($(for ((i=$start; i<=$end; i++)); do echo -n 0; done))
num=$zadd$num
fi	 
j=$	bs"_"$num
l_depth2+=("$j"" ")
done
header="CHROM POS REF ALT "${l_depth2[@]}
sed -i "1s/.*/${header}/" ${dir_in}/depth_final.txt; sed -i "1s/.*/${header}/" ${dir_in}/count_ref_final.txt
cp ${dir_in}/depth_final.txt ${dir_in}/depth.txt; cp ${dir_in}/count_ref_final.txt ${dir_in}/count_ref.txt
bzip2 ${dir_in}/depth_final.txt; bzip2 ${dir_in}/count_ref_final.txt
mv ${dir_in}/{depth_final.txt.bz2,count_ref_final.txt.bz2} ${dir_in}/depth_count

chr=($(cut -f1 -d " " ${dir_in}/depth.txt | sort | uniq ))
chr=("${chr[@]:1}")
sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dir_in}/depth.txt
sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dir_in}/count_ref.txt
sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dir_in}/${snp50k}

