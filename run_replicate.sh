#! /bin/bash
############################################ parameters ################################################## 
script=$1
dir_out=$2
dir_in=$3 
dir_popoolation=$4
fasta=$5
vcf_sansfiltre=$6
vcf_name=$7
vcf_file=${dir_out}/${vcf_name}_allfilter.vcf 
pool_size_def=$8
nbjobs=$9 
########################################################################################################## 

############################################ modules #####################################################
module load system/R-3.5.1
module load system/Python-3.6.3
module load bioinfo/bcftools-1.6
module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.15
module load bioinfo/samtools-1.8
module load bioinfo/bedtools-2.27.1
##########################################################################################################

####################################### pre-requisites ####################################################
mkdir -p ${dir_out}/temp
mkdir -p ${dir_out}/temp/pileup_replicate
mkdir -p ${dir_out}/temp/count_replicate
mkdir -p ${dir_out}/log

#################### one BAM file per colony ####################
### script Alain ###
##########################################################################################################

################### from BAM to count, depth and frequencies ####################
### colony by colony create depth and count files, rerun when new sequences arrive ###
gzip -d ${vcf_file}.gz
if [ -e ${dir_out}/temp/depth_count/depth_replicate.txt.zip ]
then
	sbatch -W --wrap="module load bioinfo/bcftools-1.6; bcftools query -f '%CHROM %POS \n' ${vcf_file} > ${dir_out}/snp_pos.txt; bcftools query -f '%CHROM %POS %REF %ALT \n' ${vcf_file} > ${dir_out}/allele_id.txt"
	unzip ${dir_out}/temp/depth_count/depth_replicate.txt.zip -d ${dir_out}
	unzip ${dir_out}/temp/depth_count/count_ref_replicate.txt.zip -d ${dir_out}
	unzip ${dir_out}/temp/depth_count/count_alt_replicate.txt.zip -d ${dir_out}
else
	sbatch -W --wrap="module load bioinfo/bcftools-1.6; bcftools query -f '%CHROM %POS \n' ${vcf_file} > ${dir_out}/snp_pos.txt; bcftools query -f '%CHROM %POS %REF %ALT \n' ${vcf_file} > ${dir_out}/allele_id.txt; cp ${dir_out}/allele_id.txt ${dir_out}/depth_replicate.txt"
	header="CHROM POS REF ALT "
	awk -v x="${header}" 'NR==1{print x} 1' ${dir_out}/depth_replicate.txt > ${dir_out}/tmp && mv ${dir_out}/tmp ${dir_out}/depth_replicate.txt; cp ${dir_out}/depth_replicate.txt ${dir_out}/count_ref_replicate.txt; cp ${dir_out}/depth_replicate.txt ${dir_out}/count_alt_replicate.txt
fi
sbatch -W -J init -o ${dir_out}/log/init.out -e ${dir_out}/log/init.err --wrap="sh ${script}/file_init.sh ${dir_in} ${dir_out} ${script} BamList_replicate ITSAP"
sbatch -W --wrap="grep -v 'done' ${dir_out}/BamList_replicate > ${dir_out}/Bam_not_done_replicate"
bam_list2=($(awk -F' ' '{print $1}' ${dir_out}/Bam_not_done_replicate))
rm ${script}/pileup_array2_replicate*
x=1
for i in ${bam_list2[@]} 
do
n=$(awk -F ' ' 'NR==num{print $2}' num=${x} ${dir_out}/Bam_not_done_replicate)
x=$[$x + 1]
echo "sh "${script}"/make_pileup.sh "${dir_out}" "${dir_popoolation}" "${fasta}" "${vcf_file}" "${script}" "${i}" "${n}" BamList_replicate;" >> ${script}/pileup_array2_replicate.sh
done
split -l ${nbjobs} ${script}/pileup_array2_replicate.sh ${script}/pileup_array2_replicate_split_
pileup_list=(${script}/pileup_array2_replicate_split_*)
for i in ${pileup_list[@]}
do
echo -e "#! /bin/bash\n$(cat ${i})" > ${i}
sbatch ${i}
done
n_job_run=$(squeue -u seynard | grep 'pileup' | wc -l)
while [ ${n_job_run} != 0 ]
do
sleep 10m 
echo ${n_job_run}
n_job_run=$(squeue -u seynard | grep 'pileup' | wc -l)
done
mv ${dir_out}/ITSAP*.zip ${dir_out}/temp/pileup_replicate

count_list=(${dir_out}/ITSAP*.count)
last=$(echo ${count_list[@]:(-1)}|rev|cut -d'/' -f 1|rev | cut -d'.' -f 1)
for i in ${count_list[@]} 
do
n=($(echo ${i}|rev| cut -d '/' -f 1|rev | cut -d'.' -f 1))
if grep ${i} ${dir_out}/depth_replicate.txt
then
echo 'already done'
else
echo 'to do'
if [ $n != $last ] 
then
sbatch -J synctocount -o ${dir_out}/log/synctocount_replicate_${n}.out -e ${dir_out}/log/synctocount_replicate_${n}.err --wrap="python ${script}/synctocount.py ${i} ${dir_out}"
else
sbatch -W -J synctocount -o ${dir_out}/log/synctocount_replicate_${n}.out -e ${dir_out}/log/synctocount_replicate_${n}.err --wrap="python ${script}/synctocount.py ${i} ${dir_out}"
fi
fi
done
for i in ${count_list[@]} 
do
n=($(echo ${i}|rev| cut -d '/' -f 1|rev | cut -d'.' -f 1))
sbatch -W --wrap="zip -j ${dir_out}/temp/count_replicate/${n}.count.zip ${dir_out}/${n}.count"
done
rm ${dir_out}/ITSAP*.count
ls -1 ${dir_out}/tmp_depth_ITSAP* | split -l 100 -d - ${dir_out}/lists_replicate
cd ${dir_out}
for list in lists_replicate*; do paste $(cat $list) > merge${list##lists}; done
cd
sbatch -W --wrap="paste ${dir_out}/merge_replicate* > ${dir_out}/depth_replicate.tmp; paste -d ' ' ${dir_out}/depth_replicate.txt ${dir_out}/depth_replicate.tmp > ${dir_out}/depth2_replicate.txt; mv ${dir_out}/depth2_replicate.txt ${dir_out}/depth_replicate.txt; sed -i 's/  / /g' ${dir_out}/depth_replicate.txt; sed -i 's/\t/ /g' ${dir_out}/depth_replicate.txt"
rm ${dir_out}/merge* ${dir_out}/lists* ${dir_out}/depth_replicate.tmp
ls -1 ${dir_out}/tmp_count_ref_ITSAP* | split -l 100 -d - ${dir_out}/lists_replicate
cd ${dir_out}
for list in lists_replicate*; do paste $(cat $list) > merge${list##lists}; done
cd
sbatch -W --wrap="paste ${dir_out}/merge* > ${dir_out}/count_ref_replicate.tmp; paste -d ' ' ${dir_out}/count_ref_replicate.txt ${dir_out}/count_ref_replicate.tmp > ${dir_out}/count_ref2_replicate.txt; mv ${dir_out}/count_ref2_replicate.txt ${dir_out}/count_ref_replicate.txt; sed -i 's/  / /g' ${dir_out}/count_ref_replicate.txt; sed -i 's/\t/ /g' ${dir_out}/count_ref_replicate.txt"
rm ${dir_out}/merge* lists* ${dir_out}/count_ref_replicate.tmp
ls -1 ${dir_out}/tmp_count_alt_ITSAP* | split -l 100 -d - ${dir_out}/lists_replicate
cd ${dir_out}
for list in lists_replicate*; do paste $(cat $list) > ${dir_out}/merge${list##lists}; done
cd
sbatch -W --wrap="paste ${dir_out}/merge* > ${dir_out}/count_alt_replicate.tmp; paste -d ' ' ${dir_out}/count_alt_replicate.txt ${dir_out}/count_alt_replicate.tmp > ${dir_out}/count_alt2_replicate.txt; mv ${dir_out}/count_alt2_replicate.txt ${dir_out}/count_alt_replicate.txt; sed -i 's/  / /g' ${dir_out}/count_alt_replicate.txt; sed -i 's/\t/ /g' ${dir_out}/count_alt_replicate.txt"
rm ${dir_out}/merge* ${dir_out}/lists* ${dir_out}/count_alt_replicate.tmp
rm ${dir_out}/tmp_depth_ITSAP*.txt ${dir_out}/tmp_count_ref_ITSAP*.txt ${dir_out}/tmp_count_alt_ITSAP*.txt
