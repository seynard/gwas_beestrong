#! /bin/bash
###### Script to copy and make files necessary for further analyses of pool sequences ######

#dir_in='/genphyse/cytogen/BeeStrong/BeeStrongHAV3_1/Data/mapping' 
#dir_out='/work/project/dynagen/seynard/GWAS/BeeStrongHAV3_1'
dir_in=$1
dir_out=$2
script=$3
type=$4
prefix=$5

bam_list=($(find ${dir_in} -name ${prefix}*.bam))

if [ -e ${dir_out}/${type}.bz2 ]
then
    bzip2 -d ${dir_out}/${type}.bz2
	for i in ${bam_list[@]} 
	do
	if [[ ${i} != *'ToMerge'* ]]
	then
	n=($(echo ${i} | rev | cut -d '/' -f 1 | rev | cut -d'_' -f 1))
	l=($(bzgrep -R ${i} ${dir_out}/${type}))
	if [ ${#l[@]} -eq 0 ]
	then
	echo "${i} ${n}" >> ${dir_out}/${type}
	fi
	fi
	done
elif [ -e ${dir_out}/${type} ]
then
	for i in ${bam_list[@]} 
	do
	if [[ ${i} != *'ToMerge'* ]]
	then
	n=($(echo ${i} | rev | cut -d '/' -f 1 | rev | cut -d'_' -f 1))
	l=($(grep -R ${i} ${dir_sonia}/${type}))
	if [ ${#l[@]} -eq 0 ]
	then
	echo "${i} ${n}" >> ${dir_out}/${type}
	fi
	fi
	done	
else
	for i in ${bam_list[@]} 
	do
	if [[ ${i} != *'ToMerge'* ]]
	then
	n=($(echo ${i} | rev | cut -d '/' -f 1 | rev | cut -d'_' -f 1))
	printf "%s\n" "${i} ${n}" >> ${dir_out}/${type}
	fi
	done
fi
module load system/R-3.5.1
Rscript ${script}/recode_dup.r ${dir_out} ${type}
