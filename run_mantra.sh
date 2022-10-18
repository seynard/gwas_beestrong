#! /bin/bash
##################################################
#################### PARAMETERS #################### 
dir=${1}
script=${2}
dir_out=${3}
list_pop=${4}
pheno=${5}
type=${6}
version=${7}
##################################################

################################################
#################### MODULES #################### 
module load system/R-4.1.1_gcc-9.3.0
################################################

mkdir -p ${dir_out}/run_mantra_${pheno}_${type}_${version}
cp ${dir}/mantra/MANTRA_software/dmatcal ${dir_out}/run_mantra_${pheno}_${type}_${version}
cp ${dir}/mantra/MANTRA_software/mantra.v1 ${dir_out}/run_mantra_${pheno}_${type}_${version}
Rscript ${script}/prep_mantra.r ${dir_out} ${list_pop} ${pheno} ${type} ${version}

cd ${dir_out}/run_mantra_${pheno}_${type}_${version}
./dmatcal
N=$((50000*4))
split -l ${N} -d --additional-suffix=.dat mantra.dat mantra_split
arr=($(ls mantra_split*)) 
jobnum=()
for m in ${arr[@]}
do
nb=${m//mantra_split/}
nb=${nb//.dat/}
mkdir -p split_${nb}
cp mantra.in split_${nb}/
cp dmat.out split_${nb}/
printf '%s\n' 'mantra.out' 'mantra.bet.out' '156127998' > split_${nb}/fname.in
cp mantra.v1 split_${nb}/
mv ${m} split_${nb}/mantra.dat
cd split_${nb}
job=`sbatch --wrap="./mantra.v1 < fname.in"`
job=$(echo $job | cut -d' ' -f4 )
jobnum+=(${job})
cd ../
done
pat=$(echo ${jobnum[@]}|tr " " "|")
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
while [ ${x} -gt 0 ]
do
sleep 30m
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
done

cd ${dir_out}/run_mantra_${pheno}_${type}_${version}
arr=($(ls -d split*))
for m in ${arr[@]}
do
nb=${m//split_/}
if [ ${nb} == '00' ]
then
cp ${m}/mantra.out mantra.out
cp ${m}/mantra.bet.out mantra.bet.out
else
cat ${m}/mantra.out >> mantra.out
cat ${m}/mantra.bet.out >> mantra.bet.out
fi
done
rm -r ${dir_out}/run_mantra_${pheno}_${type}_${version}/split_*



