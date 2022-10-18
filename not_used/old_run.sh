

#################### OLD #################### 
#############################################
#################### GMAT #################### 
#list=(`echo $pop_id2 | sed 's/,/\n/g'`)
#for i in ${list[@]}
#do
#echo ${i}
## ldak freq #
#sbatch -o ${dir}/log/grm_ldak_freq_${i}.out -e ${dir}/log/grm_ldak_freq_${i}.err --mem=100G --wrap="sh ${script}/run_grm_ldak_freq.sh ${dir_in} ${dir_out} ${script} ${pheno} ${i} ${maf_min} ${missing_rate}"
#done
#sbatch --wrap="Rscript ${script}/gmat.r ${dir_in} ${dir_out} ${pop_id2}"
#########################################################
##############################################
#################### GWAS #################### 
#type='Ligustica_Carnica'
#i=${type}
#sbatch -o ${dir}/log/gwas_gemma_egs_${i}.out  -e ${dir}/log/gwas_gemma_egs_${i}.err --mem=100G --wrap="${script}/run_GWAS_gemma_egs.sh ${dir_in} ${dir_out} ${i}"
#list=(`echo $pop_id2 | sed 's/,/\n/g'`)
#for i in ${list[@]}
#do
#echo ${i}
# gemma frequencies #
#sbatch -o ${dir}/log/gwas_gemma_freq_${i}.out  -e ${dir}/log/gwas_gemma_freq_${i}.err --mem=100G --wrap="${script}/run_GWAS_gemma_freq.sh ${dir_in} ${dir_out} ${i}"
# gemma bimbam #
#sbatch -o ${dir}/log/gwas_gemma_bimbam_${i}.out  -e ${dir}/log/gwas_gemma_bimbam_${i}.err --mem=100G --wrap="${script}/run_GWAS_gemma_bimbam.sh ${dir_in} ${dir_out} ${i}"
# gemma bed #
#sbatch -o ${dir}/log/gwas_gemma_bed_${i}.out  -e ${dir}/log/gwas_gemma_bed_${i}.err --mem=100G --wrap="${script}/run_GWAS_gemma_bed.sh ${dir_in} ${dir_out} ${i}"
# ldak freq #
#sbatch -o ${dir}/log/gwas_ldak_freq_${i}.out  -e ${dir}/log/gwas_ldak_freq_${i}.err --wrap="sh ${script}/run_GWAS_ldak_freq.sh ${dir_in} ${dir_out} ${script} ${pheno} ${maf_min} ${i}"
# ldak bimbam #
#sbatch -o ${dir}/log/gwas_ldak_bimbam_${i}.out  -e ${dir}/log/gwas_ldak_bimbam_${i}.err --wrap="sh ${script}/run_GWAS_ldak_bimbam.sh ${dir_in} ${dir_out} ${script} ${pheno} ${maf_min} ${i}"
# ldak bed #
#sbatch -o ${dir}/log/gwas_ldak_bed_${i}.out  -e ${dir}/log/gwas_ldak_bed_${i}.err --wrap="sh ${script}/run_GWAS_ldak_bed.sh ${dir_in} ${dir_out} ${script} ${pheno} ${maf_min} ${i}"
# gemma+ldak frequencies #
#sbatch -o ${dir}/log/gwas_gemmaldak_freq_${i}.out  -e ${dir}/log/gwas_gemmaldak_freq_${i}.err --mem=100G --wrap="${script}/run_GWAS_gemmaldak_freq.sh ${dir_in} ${dir_out} ${i}"
# gemma+ldak bimbam #
#sbatch -o ${dir}/log/gwas_gemmaldak_bimbam_${i}.out  -e ${dir}/log/gwas_gemmaldak_bimbam_${i}.err --mem=100G --wrap="${script}/run_GWAS_gemmaldak_bimbam.sh ${dir_in} ${dir_out} ${i}"
# gemma+ldak bed #
#sbatch -o ${dir}/log/gwas_gemmaldak_bed_${i}.out  -e ${dir}/log/gwas_gemmaldak_bed_${i}.err --mem=100G --wrap="${script}/run_GWAS_gemmaldak_bed.sh ${dir_in} ${dir_out} ${i}"
#done
#sbatch --mem=200G --wrap="Rscript ${script}/compare_gwas.r ${dir_in} ${dir_out}"
##############################################################
#################### Combine GWAS output (mantra) #################### 
#mkdir -p ${dir}/run_mantra
#dir_mantra=${dir}'/run_mantra'
#cp ${dir}/mantra/MANTRA_software/dmatcal ${dir_mantra}
#cp ${dir}/mantra/MANTRA_software/mantra.v1 ${dir_mantra}
#pop11=(Ligustica_Carnica Mellifera hybrid us)
#pop12='Ligustica_Carnica,Mellifera,hybrid,us'
#pop21=(Ligustica_Carnica Mellifera hybrid)
#pop22='Ligustica_Carnica,Mellifera,hybrid'
#pheno=($(head -n1 ${dir_out}/pheno_all.txt))
#pheno=("${pheno[@]:3}")
#input_data=(bimbam bed freq)
#gwas=(gemma ldak)
#for i in ${pheno[@]}
#do
#for j in ${input_data[@]}
#do
#for k in ${gwas[@]}
#do
#echo ${i} ${j} ${k} 
#i='eb_smr'
#j='bed'
#k='gemma'
## pop1
## prep input file for mantra
#sbatch -W --mem=200G --wrap="Rscript ${dir}/scripts/mantra.r ${dir_out} ${dir_mantra} ${pop12} ${k} ${i} lmm ${j}" 
#cd ${dir_mantra}
## run mantra to estimate fst between poulations
#sbatch -W --mem=200G --wrap="./dmatcal" 
## split file to run mantra
#N=$((50000*$(echo ${#pop11[@]}))) 
#split -l ${N} -d --additional-suffix=.dat mantra.dat mantra_split
## run mantra per file of 50k markers
#arr=($(ls mantra_split*)) 
#for m in ${arr[@]}
#do
#nb=${m//mantra_split/}
#nb=${nb//.dat/}
#mkdir -p split_${nb}
#cp mantra.in split_${nb}/
#cp dmat.out split_${nb}/
#printf '%s\n' 'mantra.out' 'mantra.bet.out' '156127998' > split_${nb}/fname.in
#cp mantra.v1 split_${nb}/
#mv ${m} split_${nb}/mantra.dat
#cd split_${nb}
#sbatch --mem=10G --wrap="./mantra.v1 < fname.in"
#cd ../
#done
## wait for all mantra run to be finished before combining the outputs
#arr=($(ls -d split*))
#for m in ${arr[@]}
#do
#nb=${m//split_/}
#n_ini=$(wc -l < split_${nb}/mantra.dat)
#n_pop=$((${#pop11[@]}+1))
#n=$(($n_ini/$n_pop))
#n_run=$(wc -l < split_${nb}/mantra.out)
#while [[ $n_run -lt $n ]]
#do
#sleep 10m
#n_ini=$(wc -l < split_${nb}/mantra.dat)
#n_pop=$((${#pop11[@]}+1))
#n=$(($n_ini/$n_pop))
#n_run=$(wc -l < split_${nb}/mantra.out)
#done 
#done
## combine split mantra outputs
#arr=($(ls -d split*))
#for m in ${arr[@]}
#do
#nb=${m//split_/}
#if [ ${nb} == '00' ]
#then
#cp ${m}/mantra.out mantra.out
#cp ${m}/mantra.bet.out mantra.bet.out
#else
#cat ${m}/mantra.out >> mantra.out
#cat ${m}/mantra.bet.out >> mantra.bet.out
#fi
#done
#cd
#mv ${dir_mantra}/dmat.out ${dir_out}/dmat1_${i}_${j}_${k}.out
#mv ${dir_mantra}/mantra.out ${dir_out}/mantra1_${i}_${j}_${k}.out
#mv ${dir_mantra}/mantra.bet.out ${dir_out}/mantra1_${i}_${j}_${k}.bet.out
## pop2
##sbatch -W --mem=200G --wrap="Rscript ${dir}/scripts/mantra.r ${dir_out} ${dir_mantra} ${pop22} ${k} ${i} lmm ${j}"
#cd ${dir_mantra}
#sbatch -W --mem=200G --wrap="./dmatcal"
#N=$((50000*$(echo ${#pop21[@]})))
#split -l ${N} -d --additional-suffix=.dat mantra.dat mantra_split
#arr=($(ls mantra_split*))
#for m in ${arr[@]}
#do
#nb=${m//mantra_split/}
#nb=${nb//.dat/}
#mkdir -p split_${nb}
#cp mantra.in split_${nb}/
#cp dmat.out split_${nb}/
#printf '%s\n' 'mantra.out' 'mantra.bet.out' '156129983' > split_${nb}/fname.in
#cp mantra.v1 split_${nb}/
#mv ${m} split_${nb}/mantra.dat
#cd split_${nb}
#sbatch --mem=10G --wrap="./mantra.v1 < fname.in"
#cd ../
#done
#arr=($(ls -d split*))
#for m in ${arr[@]}
#do
#nb=${m//split_/}
#n_ini=$(wc -l < split_${nb}/mantra.dat)
#n_pop=$((${#pop21[@]}+1))
#n=$(($n_ini/$n_pop))
#n_run=$(wc -l < split_${nb}/mantra.out)
#while [[ $n_run -lt $n ]]
#do
#sleep 10m
#n_ini=$(wc -l < split_${nb}/mantra.dat)
#n_pop=$((${#pop21[@]}+1))
#n=$(($n_ini/$n_pop))
#n_run=$(wc -l < split_${nb}/mantra.out)
#done 
#done
#arr=($(ls -d split*))
#for m in ${arr[@]}
#do
#nb=${m//split_/}
#if [ ${nb} == '00' ]
#then
#cp ${m}/mantra.out mantra.out
#cp ${m}/mantra.bet.out mantra.bet.out
#else
#cat ${m}/mantra.out >> mantra.out
#cat ${m}/mantra.bet.out >> mantra.bet.out
#fi
#done
#cd
#mv ${dir_mantra}/dmat.out ${dir_out}/dmat2_${i}_${j}_${k}.out
#mv ${dir_mantra}/mantra.out ${dir_out}/mantra2_${i}_${j}_${k}.out
#mv ${dir_mantra}/mantra.bet.out ${dir_out}/mantra2_${i}_${j}_${k}.bet.out
##done
##done
##done




















