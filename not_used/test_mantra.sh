module load system/R-3.6.2
module load system/pandoc-2.1.3
module load system/Python-3.7.4
dir='/work/project/dynagen/seynard/GWAS/infestation/'
mkdir -p ${dir}/test_mantra
dir_res=${dir}'/results'
dir_mantra=${dir}'/test_mantra'
list_pop=(Ligustica_Carnica Mellifera hybrid us)
pheno='eb_smr'
model='lmm'
input_data='bimbam'
gwas='gemma'

cp ${dir}/mantra/MANTRA_software/dmatcal ${dir_mantra}
cp ${dir}/mantra/MANTRA_software/mantra.v1 ${dir_mantra}
lpop='Ligustica_Carnica,Mellifera,hybrid,us'
sbatch -W --mem=200G --wrap="Rscript ${dir}/scripts/not_used/test_mantra.r ${dir_res} ${dir_mantra} ${lpop} ${gwas} ${pheno} ${model} ${input_data} ${i}"
cd ${dir_mantra}
sbatch -W --mem=200G --wrap="./dmatcal"
N=$((50000*$(echo ${#list_pop[@]})))
split -l ${N} -d --additional-suffix=.dat mantra.dat mantra_split
arr=($(ls mantra_split*))
for i in ${arr[@]}
do
nb=${i//mantra_split/}
nb=${nb//.dat/}
mkdir -p split_${nb}
cp mantra.in split_${nb}/
cp dmat.out split_${nb}/
printf '%s\n' 'mantra.out' 'mantra.bet.out' '156127998' > split_${nb}/fname.in
cp mantra.v1 split_${nb}/
mv ${i} split_${nb}/mantra.dat
cd split_${nb}
sbatch --mem=10G --wrap="./mantra.v1 < fname.in"
cd ../
done
arr=($(ls -d split*))
for i in ${arr[@]}
do
nb=${i//split_/}
if [ ${nb} == '00' ]
then
cp ${i}/mantra.out mantra.out
cp ${i}/mantra.bet.out mantra.bet.out
else
cat ${i}/mantra.out >> mantra.out
cat ${i}/mantra.bet.out >> mantra.bet.out
fi
done
