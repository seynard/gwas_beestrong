# Genotype Reconstruction From Pool Sequences (Haplo/Diploid polyandric species Hymenoptera)

<!-- TOC depthFrom:2 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->
- [1. Introduction](#1-Introduction)
- [2. Data preparation](#2-Data-preparation)
	- [2.1. Initial steps](#21-Initial-steps)
	- [2.2. Making reference population](#22-Making-reference-population)
	- [2.3. Input files for simulations](#23-Input-files-for-simulations)
- [3. Run AM on simulations](#3-Run-AM-on-simulations)
	- [3.1. Simulations independent markers](#31-Simulations-independent-markers)
	- [3.2. Simulations from real data](#33-Simulations-from-real-data)
		- [3.2.1. Whole genome](#331-Whole-genome)
		- [3.2.2. 1000 markers](#331-1000-markers)
- [4. Group colonies on genetic ancestries ](#4-Group-colonies-on-genetic-ancestries)
- [5. Run HPM](#5-Run-HPM)
- [6. Post-hoc analysis](#6-Post-hoc-analysis)
- [7. Validation](#7-Validation)
	- [7.1. Real data pool and males](#71-Real-data-pool-and-males)
	- [7.2. Public data](#72-Public-data)
<!-- /TOC -->

## 1. Introduction
This repository contains the scripts to perform simulations and validation on real data that is described in the paper:
From group to individual - Genotyping by pool sequencing eusocial colonies Sonia E Eynard, Alain Vignal, Benjamin Basso, Yves Le Conte, Axel Decourtye, Lucie Genestout, Emmanuelle Labarthe, Fanny Mondet, Kamila Tabet, Bertrand Servin bioRxiv 2021.11.08.467442; doi: https://doi.org/10.1101/2021.11.08.467442\

Hymenoptera have a specific reproduction system and organisation were a single queen produces all the individuals of a colony. 
In such system the colony is a group of individuals performing different 'tasks', workers are diploid and originate genetically for 1/2 from the queen and 1/2 from a male, coming itself from a pool of males having inseminated the queen; males are haploid and thus are direct representation of 1/2 of the queen genetics, as gametes. 
Traits are in general measured at the colony level.

The pipeline combines\
i) an estimation of the genetic ancestry, colony per colony, (AM)\
ii) the grouping of colonies based on ther ancestries to homogeneous populations and\
iii) the reconstruction of honeybee queen genotypes within homogeneous groups, all using pool sequencing data as an input (HPM)

The pipeline was tested on\
Simulations: where queen genotypes, allele frequencies in the drone pool are drawn from Dirichlet distributions in which sub-species are in specific proportions. With these simulations we tested the impact of depth, population composition in sub-species and model fit.\
Simulations from real data: where males from a diversity panel (Wragg et al. 2021) were combined to create a 'queen' and a pool of drones and for which offspring were drawn from. Using real data it is possible to simulate expected linkage disequilibrium and thus test our models for this parameters.\
Validated on real data: For 34 colonies we had pool sequencing data and individual sequence for 4 male offspring of the queen.\
And compared with theroretical expectation using publicaly available data from Liu et al. 2015, for which 13 to 15 drone offspring and the queen were individually sequenced.

The genome version used is HAV3.1 (Wallberg et al. 2019, doi: 10.1186/s12864-019-5642-0), the VCF for diversity panel from Wragg et al. 2022 (doi: 10.1111/1755-0998.13665) with 870 and 628 individuals were used. 
All scripts to perform simulations and analysis are available in this GitHub repository as well as the necessary input files to perform this analysis, scripts for the statistical models are available at the GitHub repository : https://github.com/BertrandServin/beethoven. This pipeline is available in the script RUN.sh. 

## 2. Data preparation
### 2.1. Initial steps
Some pre-requisite for needed programs (and their versions)
```bash
module load bioinfo/admixture_linux-1.3.0
module load bioinfo/bcftools-1.6
module load bioinfo/bedtools-2.27.1
module load system/pandoc-2.1.3
module load bioinfo/plink-v1.90b5.3
module load system/Python-3.7.4
module load system/R-3.6.2
module load bioinfo/samtools-1.8
module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.15
```
Definition of the parameters and creation of the initial directorires
```bash
dir='~/pool2geno'
mkdir -p ${dir}/data
mkdir -p ${dir}/output
mkdir -p ${dir}/raw_data
mkdir -p ${dir}/code
dirin=${dir}/data
dirout=${dir}/output
dir_save=${dir}/raw_data
script=${dir}/code
dir_popoolation='~/PoPoolation2/popoolation2_1201'
fasta=${dir_save}/Fasta/GCF_003254395.2_Amel_HAv3.1_genomic.fna # FASTA for HAV3.1 reference genome of the honeybee
vcf_name='MetaGenotypesCalled870_raw_snps' # full initial vcf before cleaning
vcf_file=${dirin}/${vcf_name}_allfilter.vcf # name of vcf after cleaning
snp_list=${dirin}/HAV3_1_50000.txt # list of 50k SNP 
pop_id='Ligustica_Carnica,Mellifera,Caucasia' # name of genetic backgrounds/sub species
n_pop=$(echo $(IFS=','; set -f; set -- $pop_id; echo $#))
ncpu=20 # number of cpus
unif_threshold=0.99 # genetic ancestry threshold for reference individuals
pool_size_def=150 # parameter used for popoolation BAM to pileup inference
nbjobs=10 # number of jobs run in parallel
nboot=100 # number of iteration 
mkdir -p ${dirin}/log
mkdir -p ${dirout}/log
```

Part of this study relies on the use of a VCF produced by Wragg et al. 2022 from a diversity panel. 
The paper can be found doi: 10.1111/1755-0998.13665 and the scripts to produce the VCF on the github repository https://github.com/avignal5/SeqApiPop. For the purpose of this study the VCF was 'cleaned' and markers were filtered according to the pipeline describe in https://github.com/seynard/vcf_cleanup, using the version v0.\
If needed this step can be run as follows
```bash
mkdir -p ${dirin}/vcf_cleanup
sbatch -J vcf_cleanup -o ${dirin}/log/vcf_cleanup.out -e ${dirin}/log/vcf_cleanup.err -W --wrap="${script}/vcf_cleanup/run_vcfcleanup.sh seynard ${script}/vcfcleanup ${dir_save}/The870vcf ${dirin} ${vcf_name}.vcf.gz 3 -999" # script to perform cleaning
sbatch -W --wrap="mv ${dirin}/{plot_decision*,*_value.txt,venn*,*.vcf.log,geno.txt,gq.txt,info.txt,GQfiltered.txt,heterozygous.txt,list_kept.txt,missing.txt,snp_*.txt,count*.txt} ${dirin}/vcf_cleanup/" # move all files produced by vcfcleanup to specific directory
```
One additional preparation step: For the purpose of this analysis we focused on autosomes
```bash
sbatch -W --wrap="sh ${script}/1_remove_mito.sh ${dirin} ${vcf_file} NC_001566.1" # removing mitochondrial markers.
```

Final prepration of the files 
```bash
sbatch -W --wrap="bgzip -c ${vcf_file} > ${vcf_file}.gz"
sbatch -W --wrap="bcftools query -f '%CHROM %POS\n' ${dirin}/vcf.vcf > ${dirin}/snp_pos.txt; sed -i "1i\CHROM POS" ${dirin}/snp_pos.txt" # extract chromosome and position information from vcf file
chr=($(cut -f1 -d ' ' ${dirin}/snp_pos.txt | sort | uniq ))
chr=("${chr[@]:1}") # remove first value of array chr as it contains the column name 'CHROM'
sbatch -W --wrap="sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dirin}/vcf.vcf" # change chromosome name to chromosome number
```
### 2.2. Making reference population
Allele frequencies for each of the genetic background on which we focus are necessary. 
We base the definition of the reference population on the database of the diversity panel provided by Wragg et al. 2022 (10.1111/1755-0998.13665).\
Such diversity panel allows for the identification of reference individuals, 'pure' individuals representing the sub-species of interest, in our case Mellifera, Ligustica_Carnica and Caucasia
Using Admixture, in a non supervised set up with k=3 we estimate the allele frequencies in each genetic background/sub species, we can also identify 'pure' individuals for whom the genetic ancestry for one of the three genetc background is above the set threshold 'unif_threshold'
```bash
sbatch -J run_SeqApiPop -c ${ncpu} -o ${dirin}/log/run_SeqApiPop.out -e ${dirin}/log/run_SeqApiPop.err -W --wrap="${script}/2_run_SeqApiPop.sh ${script} ${dirout} ${dir_save} ${dirin} ${fasta} ${dirin}/vcf.vcf ${n_pop} ${ncpu} ${pop_id} ${unif_threshold}"
sbatch -W --wrap="bcftools query -f '%CHROM,%POS,%REF,%ALT[,%GT]\n' ${dirin}/panel_seqapipop/vcf_subsp.vcf > ${dirin}/geno_ref.txt" # create genotype file for the reference individuals
sbatch -W --wrap="sed -i -e 's:0/0:2:g' -e 's:0|0:2:g' -e 's:0/1:1:g' -e 's:0|1:1:g' -e 's:1/1:0:g' -e 's:1|1:0:g' -e 's:./.:-9:g' -e 's:\.:-9:g' ${dirin}/geno_ref.txt" # recode genotype file for the reference individuals from 0/0, 0/1, 1/1 and ./. to 0, 1, 2, -9
sbatch -W --wrap="sed -i -e 's:-91:\.1:g' ${dirin}/geno_ref.txt"
```
### 2.3. Input files for simulations
In this section we prepare the input files that will be necessary for the different simulation scenarios and estimation of the genetic composition of the queen for each colony based on the simulated genotypes and using AM 

First, in order to perform the estimation of genetic ancestry we used a subset of markers with caracteristics such as: variance across the three genetic background above 0.1, maf above 0.1 
```bash
sbatch -W="Rscript ${script}/3_marker_choice_freq.r ${dirin} ${pop_id}"
cp ${dirin}/freq_admix_mv.txt ${dirin}/sim_freq.txt 
```

For 'complete' simulations the script prepares and creates input files as well as run the simulation directly after therefore the script to prepare this data set is included into the '5_run_simul.sh'

For the simulations based on real data from haploid drones we used data from Wragg et al. (2022) (10.1111/1755-0998.13665) to create the different populations. It was first necessary to define groups based on genetic composition and to then sample individuals from the groups of interest to create our simulated colonies. For each colony we sampled 2 males to create a queen genome (2 chromosomes) and 15 other males were to create the mate group. From this we generated an offspring group for which we measured depth and frequencies.
Definition of additional parameters  that are necessary for this section of the process
```bash
prefix='simul_data'
nbmales=15 # simulated number of inseminating males
depth=30 # simulated sequencing depth
nboffspring=0 # optional to generate actual offspring from the mating of the simulated queen and the insemating drones. If set to 0 the pool will be based on an average allele frequency, if set to more than 0 we will create offspring genotypes and estimate the pool from these offspring
```
```bash
sbatch -J prep_simul_data -o ${dirin}/log/prep_simul_data.out -e ${dirin}/log/prep_simul_data.err -W --wrap="${script}/6_prep_simul_data.sh ${script} ${dirout} ${dirin} ${dirin}/vcf.vcf ${nbmales} ${n_pop}"
nbsimul=($(awk '{print $1}' $dirin/simul_data.txt | uniq))
for i in ${nbsimul[@]}
do
	n=$((${i}+1))
	simul_number=($(sed "${n}q;d" ${dirin}/choice_scenario_simul_data.txt|cut -f1))
	echo ${simul_number}
	sbatch -W --mem=20G -J prep_simul_data_${simul_number} -o ${dirin}/log/prep_simul_data_${simul_number}.out -e ${dirin}/log/prep_simul_data_${simul_number}.err --wrap="${script}/7_prep_simul_data2.sh ${script} ${dirin} ${dirout} ${nbmales} ${depth} ${nboffspring} ${simul_number} ${prefix}"
done
```

## 3. Run AM on simulations (script for AM available at https://github.com/BertrandServin/beethoven)
### 3.1. Simulations independent markers
Run simulation. Here we can adjust the number of colonies, the number of markers, the depth, the frequency mean or variance threshold for resampling, the sub-species proportions in both queen and drones pool.
Use of parameter file simul.txt columns: simulation number; dirichlet alpha parameters for queen genetic composition; dirichlet alpha parameters for inseminating drones genetic composition; number of simulated colonies; simulated sequencing depth; number of simulated markers
```bash
nbsim=$(wc -l < $dirin/simul.txt)
freq='sim_freq.txt'
nb_marker=($(sed "1q;d" ${dirin}/simul.txt|cut -f6))
sbatch -W --wrap="Rscript ${script}/4_select_marker_simul.r ${dirin} ${nb_marker} ${freq} ${n_pop} ${pop_id}"
for i in $(seq 1 1 ${nbsim})
do
	echo ${i}
	freq='sim_freq.txt'
	n_col_snpid=2
	n_col=100 # number of colonies simulated
	simul_number=($(sed "${i}q;d" ${dirin}/simul.txt|cut -f1)) # simulation number
	distrib_q=($(sed "${i}q;d" ${dirin}/simul.txt|cut -f2)) # simulated dirichlet alpha parameters for queen genetic composition
	distrib_d=($(sed "${i}q;d" ${dirin}/simul.txt|cut -f3)) # simulated dirichlet alpha parameters for inseminating drones genetic composition
	n_col_prop=($(sed "${i}q;d" ${dirin}/simul.txt|cut -f4)) # simulated population composition
	depth=($(sed "${i}q;d" ${dirin}/simul.txt|cut -f5)) # simulated sequencing depth
	nbmales=15 # simulated number of inseminating males
	sbatch -W -J simul_${simul_number} -o ${dirout}/log/simul_${simul_number}.out -e ${dirout}/log/simul_${simul_number}.err --wrap="${script}/5_run_simul.sh ${dirin} ${dirout} ${script} ${nb_marker} ${simul_number} ${freq} ${n_col} ${n_col_prop} ${depth} ${distrib_q} ${distrib_d} ${ncpu} ${n_col_snpid} ${n_pop} ${pop_id} ${nbmales}" # run simulation 
done
```

### 3.2. Simulations from real data
Definition of additional parameters necessary for this part of the pipeline. 
```bash
prefix1='simul_data'
n_col_snpid=4
nbsimul=($(awk '{print $1}' ${dirin}/simul_data.txt | uniq))
```
#### 3.2.1. Whole genome (linked markers)
```bash
freq=freq_admix${n_pop}.txt
for i in ${nbsimul[@]}
do
	echo ${i}
	prefix=$(echo ${prefix1}${i})
	echo ${prefix}
	sbatch -W --mem=50G -J ${prefix} -o ${dirout}/log/${prefix}.out -e ${dirout}/log/${prefix}.err --wrap="${script}/8_run_simul_data.sh ${script} ${dirin} ${dirout} ${freq} ${n_col_snpid} ${ncpu} ${n_pop} ${pop_id} ${prefix} ${snp_list}"
done
```
#### 3.2.2. 1000 markers (independent markers)
```bash
freq='sim_freq.txt'
nbsimul=($(awk '{ $1}' ${dirin}/simul_data.txt | uniq))
for i in ${nbsimul[@]}
do
	echo ${i}
	prefix0=$(echo ${prefix1}${i})
	prefix=$(echo ${prefix1}${i}_1000)
	echo ${prefix}
	sbatch -W -J ${prefix} -o ${dirout}/log/${prefix}.out -e ${dirout}/log/${prefix}.err --wrap="${script}/9_run_simul_data_1000.sh ${script} ${dirin} ${dirout} ${freq} ${n_col_snpid} 1 ${n_pop} ${pop_id} ${prefix} ${prefix0}"
done
```
The steps from 2.3. to 3.2.2 where repeated for depth 10 and 100. 

## 4. Group colonies on genetic ancestries 
Using q-vectors of genetic background estimated by AM we can group colonies from our simulations into homogeneous populations (necessary to perform further analysis).
```bash
sce_list='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15' #list of scenarios kept for analysis (here we did not use all our scenarios as we made extreme scenarios for other purposes or with changing parameters)
n_col_snpid=2
type=(simul simul_data1000 simul_data)
DEPTH=(10 30 100)
for a in ${type[@]}
do
	for i in ${DEPTH[@]}
	do
		echo ${a} ${i}
		sbatch -W --mem=200G --wrap=" ${script}/10_prep_mix.sh ${dirin} ${dirout} ${script} ${a} ${i} ${sce_list} ${n_col_snpid}"
	done
done
```

## 5. Run HPM (script for HPM can be found at https://github.com/BertrandServin/beethoven)
Definition of additional parameters necessary for this part of the pipeline. 
```bash
n_col_snpid=2
type=(simul simul_data1000 simul_data)
DEPTH=(10 30 100)
```
For each homogeneous population we can now run HPM 
```bash
for a in ${type[@]}
do
	for i in ${DEPTH[@]}
	do
		n=$(ls ${dirin}/Mix_${type}/depth${DEPTH}_*.txt| wc -l)
		for j in $(seq 1 ${n})
		do
			echo ${a} ${i} ${j}
			sbatch -J hom -W --mem=200G --wrap=" python .local/bin/model_genoqueen_hom ${dirin}/Mix_${a} depth${i}_${j}.txt count${i}_${j}.txt ${n_col_snpid} 1 sim_model${i}_${j} 100000"
		done
	done
done
x=$(squeue -u seynard | grep 'hom' |wc -l)
while [ ${x} -gt 0 ]
do
	wait 10m
	x=$(squeue -u seynard | grep 'hom' |wc -l)
done
for a in ${type[@]}
do
	mv ${dirin}/Mix_${a}/sim_model* ${dirout}/Mix_${a}/sim_model* 
done
```

## 6. Post-hoc analysis
Definition of additional parameters necessary for this part of the pipeline. 
```bash
n_col_snpid=2
type=(simul simul_data1000 simul_data)
DEPTH=(10 30 100)
```
Now we estimate statistics on the output of the model HPM for each of the homogeneous population. These statistics are: genotyping error, genotyping error as a function of MAF, genotyping error as a function of genotype probability, genotyping probability as a function of MAF, three way analysis of genotyping error as function of MAF bin and genotype probability, genotype calibration and genotype probability
```bash
for a in ${type[@]}
do
	for i in ${DEPTH[@]}
	do
		n=$(ls ${dirin}/Mix_${a}/depth${i}_*.txt| wc -l)
		for j in $(seq 1 ${n})
		do
			echo ${a} ${i} ${j}
			sbatch -J sum --mem=100G --wrap="python ${script}/11_summary_simul_mix.py ${dirin} ${dirout} ${a} ${i} ${j}"
		done
	done
done
x=$(squeue -u seynard | grep 'sum' |wc -l)
while [ ${x} -gt 0 ]
do
	wait 10m
	x=$(squeue -u seynard | grep 'sum' |wc -l)
done
```
And finally we combine these summary statistics
```bash
for a in ${type[@]}
do
	for i in ${DEPTH[@]}
	do
		echo ${a} ${i}
		sbatch --mem=50G --wrap="${script}/12_combine_summary.sh ${dirin} ${dirout} ${a} ${i}"
	done
done		
```

## 7. Validation
### 7.1. Real data pool and males
For the purpose of validation of our models we had access to 34 colonies having one the one hand a pool of workers sequenced and on the other 4 males individually sequenced.\
This dataset provides a comparison for queen genotype reconstruction based on using male offspring as 'gametes' of the queen and reconstruction based on a pool of workers made using our models. 

Definition of additional parameters necessary for this part of the pipeline. 
```bash
t_recom=0.02325581 #1 recombination across 43 meiosis to estimate genetic map using Liu et al (2015)
seq_error=0.001 #sequencing error of 10^-3
prefix=mosar
freq=freq_admix${n_pop}.txt
n_col_snpid=4
```

```bash
sbatch -J prep_MOSAR -o ${dirin}/log/prep_MOSAR.out -e ${dirin}/log/prep_MOSAR.err -W --wrap="${script}/13_prep_MOSAR.sh ${script} ${dirout} ${dir_save} ${dirin} ${dir_popoolation} ${fasta} ${dirin}/vcf.vcf ${dirin}/vcf_males.vcf ${pool_size_def} ${nbjobs}"
sbatch -W --wrap="${script}/14_mean_depth.sh ${dirin}/depth_mosar.txt ${dirin}/mean_depth_mosar.txt"
```

Run sequencially AM and HPM.\
Due to the sample size we did not perform any grouping of the colonies based on genetic composition before estimating the queen genotype, this step could be added if more samples were available
```bash
sbatch -W --mem=50G -J ${prefix} -o ${dirout}/log/${prefix}.out -e ${dirout}/log/${prefix}.err --wrap="${script}/15_run_data.sh ${script} ${dirin} ${dirout} ${freq} ${n_col_snpid} ${ncpu} ${n_pop} ${pop_id} ${prefix} ${snp_list}"
```

As male offspring of the queen are equivalent to gametes (haploid), it is possible to estimate the genotype probability of the queen based on their genotypes.  
```bash
sbatch -W -J male_2_queen -o ${dirout}/log/male_2_queen.out -e ${dirout}/log/male_2_queen.err --wrap="python ${script}/0_male_proba_geno.py ${dirout}/mosar ${dirin}/geno_males_${prefix}.txt ${seq_error}"
```

Finally, we used ADMIXTURE to estimate genetic composition of male offspring and queen, reconstructed from these offspring. 
```bash
sbatch -W --mem=50G -o ${dirout}/log/admix_males_${prefix}.out -e ${dirout}/log/admix_males_${prefix}.err --wrap="${script}/16_admix_males.sh ${script} ${dirin} ${dirout} ${freq} ${n_col_snpid} ${ncpu} ${n_pop} ${pop_id} ${prefix} ${snp_list}"
```

### 7.2. Public data 
Finally, Liu et al. 2015 (10.1186/s13059-014-0566-0) made available sequence data for 46 individuals, coming from 3 colonies. With 2 colonies having 15 males and 1 queen sequenced and the last colony having 13 males and the queen sequenced.\
This dataset allows to perform resampling of male offspring (multiple times, 100 bootstrap, and different numbers, 4, 6, 8 and 10) and reconstruct the queen genotype from theirs. Therefore benchmarking the standard method for queen genotype reconstruction.


Definition of additional parameters necessary for this part of the pipeline. 
```bash
prefix=Liu
t_recom=0.02325581 #1 recombination across 43 meiosis to estimate genetic map using Liu et al (2015)
seq_error=0.001 # sequencing error of 10^-3
```

```bash
nind_test='4,6,8,10' # number of individuals sampled 
col='Colony1,Colony2,Colony3'
sbatch -W --wrap="${script}/17_make_geno_proba_Liu.sh ${script} ${dirin} ${dirout} Liu ${col} ${nind_test} ${nboot} ${seq_error} all 0"
```

