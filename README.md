# GWAS for varroa resistance in an heterogeneous honeybee population

<!-- TOC depthFrom:2 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->
- [1. Introduction](#1-Introduction)
- [2. Data preparation](#2-Data-preparation)
	- [2.1. Initialisation of parameters](#21-Initialisation-of-parameters)
	- [2.2. Paths and parameters description](#22-Paths-and-parameters-description)
	- [2.3. Prep data for reference populations](#23-Prep-data-for-reference-populations)
	- [2.4. Prep pool sequence data](#24-Prep-pool-sequence-data)
- [3. Genetic composition estimation and queen genotype reconstruction](#3Genetic-compo-geno)
	- [3.1.Estimation of genetic composition](#31-Estimation-of-genetic-composition)
	- [3.2.Group by homogeneous genetic composition](#32-Group-by-homogeneous-genetic-composition)
	- [3.3.Reconstruct queen genotype within homogeneous groups](#33-Reconstruct-queen-genotype-within-homogeneous-groups)
- [4. Phenotypes ](#4-Phenotypes)
- [5. Linkage disequilibrium](#5-LD)
- [6. Simple GWAS](#6-GWAS)
- [7. Meta GWAS](#7-MetaGWAS)
	- [7.1. Define groups](#71-Define-groups)
	- [7.2. MANTRA](#72-MANTRA)
	- [7.3. Mash](#73-mash)
- [8. Analysis on genes](#Genes)
	- [8.1. Marker vs genes correspondence](#81-Marker-vs-genes-correspondence)
	- [8.2. Statistical analysis on genes](#82-Statistical-analysis-on-genes)

<!-- /TOC -->

## 1. Introduction
This repository contains the scripts to perform the analysis described in the paper:
A large genome wide association study on a heterogeneous honeybee, Apis mellifera, population reveals the polygenic genetic architecture for resistance to the major honeybee parasite, Varroa destructor\
Sonia E Eynard, Alain Vignal, Benjamin Basso, Olivier Bouchez, Tabatha Bulach, Yves Le Conte, Benjamin Bainat, Axel Decourtye, Lucie Genestout, Matthieu Guichard, François Guillaume, Emmanuelle Labarthe, Rachid Mahla, Fanny Mondet, Markus, Neudischko, Florence Phocas, Yannick Poquet, Christina Sann, Rémi-Félix Serre, Kamila Tabet, Bertrand Servin \
bioRxiv , doi: \

The genome version used is HAV3.1 (Wallberg et al. 2019, doi: 10.1186/s12864-019-5642-0), the VCF for diversity panel from Wragg et al. 2022 (doi: 10.1111/1755-0998.13665) with 870 and 628 individuals were used.
All scripts to perform simulations and analysis are available in this GitHub repository, scripts for the statistical models to estimate genetic composition and reconstruct queen genotype are available at the GitHub repository : https://github.com/BertrandServin/beethoven.
This pipeline is available in the script RUN.sh.

## 2. Data preparation
### 2.1. Initialisation of parameters
Here is the list of the programs, with their versions, necessary for the analysis.
```bash
module load system/R-4.1.1_gcc-9.3.0
module load system/pandoc-2.1.3
module load bioinfo/bcftools-1.6
module load bioinfo/gemma-v0.97
module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.15
module load bioinfo/admixture_linux-1.3.0
module load bioinfo/plink-v1.90b5.3
module load bioinfo/shapeit.v2.904
module load bioinfo/samtools-1.8
module load bioinfo/bedtools-2.27.1
module load bioinfo/blast-2.2.26
module load system/Python-3.6.3
module load system/jdk-12.0.2
module load bioinfo/interproscan-5.44-79.0
module load system/Python-3.7.4
```

### 2.2. Paths and parameters description
Here we describe the paths and parameters that we will use for the analysis
```bash
dir="~/GWAS"
script="~/GWAS/scripts"
dir_out="~/GWAS/results"
dir_in="~/GWAS/data"
dir_save="~/save"
fasta=${dir_save}/Fasta/GCF_003254395.2_Amel_HAv3.1_genomic.fna
vcf_name="MetaGenotypesCalled870_raw_snps"
vcf_sansfiltre=${dir_in}/${vcf_name}.vcf.gz
vcf_file=${dir_in}/${vcf_name}_allfilter.vcf
pop_id="Mellifera,Caucasia,Ligustica_Carnica"
n_pop=$(echo $(IFS=","; set -f; set -- $pop_id; echo $#))
pop_id2="Mellifera,Caucasia,Ligustica_Carnica,hybrid"
pheno="Data_BS.csv"
pheno_mito='compareDepthsBeeVarroa.txt'
snp50k='HAV3_1_50000.txt'
maf_min=0.01
compo_threshold=0.8
missing_rate=0.05
ncpu=10
unif_threshold=0.99
nbjobs=10
d_min=10
d_max=50
pool_size_def=150
```

### 2.3. Prep data for reference populations
In Wragg et al. 2022 () a vcf for honeybee population representative of the genetic diversity in Europe has been developped. We use this diversity panel in this analysis to describe reference populations for the three main subspecies in Europe, and thus present in our french heterogeneous population, that are: Apis mellifera mellifera, Apis mellifera ligustica & carnica and Apis mellifera caucasia.
```bash
sbatch -W --mem=100G --job-name='prep_seqapipop' --wrap="${script}/prep_seqapipop.sh ${script} ${dir_in} ${vcf_file} ${n_pop} ${pop_id} ${unif_threshold} ${ncpu}"
```

### 2.4. Prep pool sequence data
```bash
sbatch -W --mem=200G --job-name='prep_pileup' --wrap="${script}/prep_pileup.sh ${script} ${dir_in} ${fasta} ${vcf_name} ${vcf_sansfiltre} ${vcf_file} ${snp50k} ${nbjobs} ${d_min} ${d_max} ${pool_size_def}"
```

## 3. Genetic composition estimation and queen genotype reconstruction (script for the statistical models available at https://github.com/BertrandServin/beethoven)
### 3.1. Estimation of genetic composition
```bash
Ncol=$(awk "{print NF;exit}" ${dir_in}/depth.txt)
colname=($(head -n1 ${dir_in}/depth.txt))
col_non_bs=($(echo ${colname[@]//BS*/}))
n_col_snpid=$(echo ${#col_non_bs[@]})
n_col=$((${Ncol}-${n_col_snpid}))
colname=`echo $(echo ${colname[@]}) | tr ' ' ','`
IFS=',' read -r -a array <<< "$colname"
jobnum=()
for i in $(seq $((${n_col_snpid}+1)) 1 ${Ncol})
do
	j=$((${i}-1))
	echo $j
	cname=${array[${j}]}
	job=`sbatch --mem=100G --job-name='het' --wrap="${script}/run_het.sh ${dir_in} ${cname} ${n_pop} ${snp50k} ${j}"`
	job=$(echo $job | cut -d' ' -f4 )
	jobnum+=(${job})
done
pat=$(echo ${jobnum[@]}|tr " " "|")
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
while [ ${x} -gt 1 ]
do
	sleep 10m
	x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
done
rm ${dir_in}/tmp* ${dir_in}/*_st_het.geno.bz2
mkdir -p ${dir_in}/input
mv ${dir_in}/sim_depth_count*.txt ${dir_in}/input/
bzip2 ${dir_in}/input/*
```

### 3.2. Group by homogeneous genetic composition
```bash
choice='_50k'
sbatch -W --wrap="Rscript ${script}/q_matrix.r ${dir_in} ${n_pop} ${pheno} ${pop_id} ${choice}"
sbatch -W --wrap="Rscript ${script}/group_q_matrix.r ${dir_in} . ${pop_id2} ${pop_id} ${compo_threshold} ${pheno} ${choice}"
```

### 3.3. Reconstruct queen genotype within homogeneous groups
```bash
Ncol=$(awk "{print NF;exit}" ${dir_in}/depth.txt)
colname=($(head -n1 ${dir_in}/depth.txt))
col_non_bs=($(echo ${colname[@]//BS*/}))
n_col_snpid=$(echo ${#col_non_bs[@]})
pop_idi=($(echo ${pop_id2}| tr "," "\n"))
jobnum=()
for i in ${pop_idi[@]}
do
	echo ${i}
	job=`sbatch -J hom_freq_${i} --wrap="${script}/run_hom_freq.sh ${dir_in} ${dir} ${script} ${i} ${ncpu} ${n_col_snpid} ${maf_min}"`
	job=$(echo $job | cut -d' ' -f4 )
	jobnum+=(${job})
done
pat=$(echo ${jobnum[@]}|tr " " "|")
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
while [ ${x} -gt 1 ]
do
	sleep 30m
	x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
done
awk '{print $1" "$2" "$3" "$4}' ${dir_in}/depth.txt > ${dir_in}/allele_id_final.txt
tail -n +2 ${dir_in}/allele_id_final.txt  > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/allele_id_final.txt
```

## 4. Phenotypes
```bash
sbatch -o ${dir}/log/prep_pheno.out  -e ${dir}/log/prep_pheno.err -W --wrap="Rscript ${script}/prep_pheno.r ${dir_in} ${dir_out} ${pheno} ${pop_id2} ${pheno_mito}"
```

## 5. Linkage disequilibrium
```bash
mkdir ${dir_in}/ld
sbatch -J 'list_ld' --mem=20G -W --wrap="Rscript ${script}/list_ld.r ${dir_in}/seqapipop ${dir_in}/ld ${pop_id2} ${compo_threshold}"
jobnum=()
pop_idi=($(echo ${pop_id2}| tr "," "\n"))
for i in ${pop_idi[@]}
do
	echo ${i}
	sbatch --mem=200G -W --wrap="plink --bfile ${dir_in}/seqapipop/seqapipop --keep ${dir_in}/ld/list_${i}_ld.txt --keep-allele-order --make-bed --out ${dir_in}/ld/${i}"
	chr=($(cut -f1 -d " " ${dir_in}/snp_list.txt | sort | uniq ))
	unset 'chr[${#chr[@]}-1]'
	echo ${chr[@]}
	for c in ${chr[@]}
	do
		echo $c
		job=`sbatch --mem=200G --wrap="plink --bfile ${dir_in}/ld/${i} --chr $c --make-bed --out ${dir_in}/ld/${i}_${c}; plink --bfile ${dir_in}/ld/${i}_${c} --ld-window-r2 0.1 --out ${dir_in}/ld/ld_${i}_${c} --r2 inter-chr"`
		jobnum+=(${job})
	done
done
pat=$(echo ${jobnum[@]}|tr " " "|")
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
while [ ${x} -gt 1 ]
do
	sleep 30m
	x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
done
```

## 6. Simple GWAS
```bash
list=(`echo ${pop_id2} | sed 's/,/\n/g'`)
jobnum=()
for i in ${list[@]}
do
	echo ${i}
	if [[ -f ${dir_in}/geno_hom_${i}.bgs ]]
	then
		echo 'pop to run'
		job=`sbatch -J ${i} --mem=100G -o ${dir}/log/${i}.out -e ${dir}/log/${i}.err --wrap="${script}/run_type.sh ${dir} ${dir_in} ${dir_out} ${script} ${i} ${maf_min} ${missing_rate}"`
		job=$(echo $job | cut -d' ' -f4 )
		jobnum+=(${job})
	fi
done
pat=$(echo ${jobnum[@]}|tr " " "|")
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
while [ ${x} -gt 1 ]
do
	sleep 30m
	x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
done
list=(`echo ${pop_id2} | sed 's/,/\n/g'`)
jobnum=()
for i in ${list[@]}
do
	echo ${i}
	if [[ -f ${dir_in}/geno_hom_${i}.bgs ]]
	then
		job=`sbatch -J ${i} --mem=100G -o ${dir}/log/plot${i}.out -e ${dir}/log/plot${i}.err --wrap="Rscript ${script}/check_gwas.r ${dir_out} ${i} 0.1"`
		job=$(echo $job | cut -d' ' -f4 )
		jobnum+=(${job})
	fi
done
pat=$(echo ${jobnum[@]}|tr " " "|")
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
while [ ${x} -gt 1 ]
do
	sleep 30m
	x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
done
bzip2 ${dir_out}/summary_gwas_*.txt
```

## 7. Meta GWAS
### 7.1. Define groups
```bash
pheno_id=($(ls ${dir_out}/pheno_*.txt))
pheno_list=()
for i in ${pheno_id[@]}
do
	pheno=$(echo ${i} | cut -d '/' -f9)
	pheno=$(echo ${pheno} | awk '{gsub(/pheno_/,"")}1')
	pheno=$(echo ${pheno} | awk '{gsub(/.txt/,"")}1')
	ele=($(echo ${pop_id2} |awk -F, '{for (i=1;i<=NF;i++)print $i}'))
	for j in ${ele[@]}
	do
		pheno=$(echo ${pheno} | awk '{gsub(/_'"${j}"'/,"")}1')
		pheno=$(echo ${pheno} | awk '{gsub(/'"${j}"'/,"")}1')
	done
	pheno_list[${#pheno_list[@]}]=${pheno}
done
pheno_list=($(echo "${pheno_list[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
list_pop=($(echo ${pop_id2[@]}| tr "," "\n"))
delete=Caucasia
list_pop=("${list_pop[@]/$delete}")
```

### 7.2. MANTRA
Morris. 2012
```bash
for i in ${pheno_list[@]}
do
	job=`sbatch -J 'mantra' -o ${dir}/log/mantra${i}.out -e ${dir}/log/mantra${i}.err --mem=100G --wrap="${script}/run_mantra.sh ${dir} ${script} ${dir_out} ${list_pop} ${i} egs"`
	job=$(echo $job | cut -d' ' -f4 )
	jobnum+=(${job})
done
pat=$(echo ${jobnum[@]}|tr " " "|")
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
while [ ${x} -gt 1 ]
do
	sleep 30m
	x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
done
```

### 7.3. Mash
Urbut et al. 2019
```bash
jobnum=()
for i in ${pheno_list[@]}
do
	job=`sbatch -o ${dir}/log/analysis_mash_${i}_egs1.out -e ${dir}/log/analysis_mash_${i}_egs1.err -J 'mash' --mem=200G --wrap="Rscript ${script}/run_mash.r ${dir_out} ${i} egs 0.2 ${list_pop}"`
	job=$(echo $job | cut -d' ' -f4 )
	jobnum+=(${job})
done
pat=$(echo ${jobnum[@]}|tr " " "|")
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
while [ ${x} -gt 1 ]
do
	sleep 30m
	x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
done
```

## 8. Analysis on genes
### 8.1. Marker vs genes correspondence
Variant Effect Predictor. This Ensembl's tool provides, for each marker tested it's location and distance to nearby genes, the type of variation associated (upstream, downstream, intronic) and its effect (mutation ...).
```bash
var1='/'
var2='+'
awk -v var1="${var1}" -v var2="${var2}" '{print $1" "$2" "$2" "$3""var1""$4" "var2}' ${dir_in}/allele_id_final.txt > ${dir_in}/input_vep.txt
awk '$1=="10"{$1="CM009940.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="11"{$1="CM009941.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="12"{$1="CM009942.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="13"{$1="CM009943.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="14"{$1="CM009944.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="15"{$1="CM009945.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="16"{$1="CM009946.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="1"{$1="CM009931.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="2"{$1="CM009932.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="3"{$1="CM009933.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="4"{$1="CM009934.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="5"{$1="CM009935.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="6"{$1="CM009936.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="7"{$1="CM009937.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="8"{$1="CM009938.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
awk '$1=="9"{$1="CM009939.2"}1' ${dir_in}/input_vep.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/input_vep.txt
cp ${dir_in}/input_vep.txt ${dir}/ensembl-vep/
cd ${dir}/ensembl-vep/
./vep -i input_vep.txt -o output_vep.txt --species apis_mellifera --database --genomes --force_overwrite
cd
cp ${dir}/ensembl-vep/output_vep.txt ${dir_out}/output_vep.txt
```

### 8.2. Statistical analysis on genes
mBAT-combo (GCTA suite, Li et al. 2022)
```bash
awk -F',' '{print $1" "$2" "$3" "$6}' ${dir_in}/proteins_48_403979.csv > ${dir_in}/gene_list.txt
sed -i -e 's/"linkage group LG1"/1/g' ${dir_in}/gene_list.txt
sed -i -e 's/"linkage group LG2"/2/g' ${dir_in}/gene_list.txt
sed -i -e 's/"linkage group LG3"/3/g' ${dir_in}/gene_list.txt
sed -i -e 's/"linkage group LG4"/4/g' ${dir_in}/gene_list.txt
sed -i -e 's/"linkage group LG5"/5/g' ${dir_in}/gene_list.txt
sed -i -e 's/"linkage group LG6"/6/g' ${dir_in}/gene_list.txt
sed -i -e 's/"linkage group LG7"/7/g' ${dir_in}/gene_list.txt
sed -i -e 's/"linkage group LG8"/8/g' ${dir_in}/gene_list.txt
sed -i -e 's/"linkage group LG9"/9/g' ${dir_in}/gene_list.txt
sed -i -e 's/"linkage group LG10"/10/g' ${dir_in}/gene_list.txt
sed -i -e 's/"linkage group LG11"/11/g' ${dir_in}/gene_list.txt
sed -i -e 's/"linkage group LG12"/12/g' ${dir_in}/gene_list.txt
sed -i -e 's/"linkage group LG13"/13/g' ${dir_in}/gene_list.txt
sed -i -e 's/"linkage group LG14"/14/g' ${dir_in}/gene_list.txt
sed -i -e 's/"linkage group LG15"/15/g' ${dir_in}/gene_list.txt
sed -i -e 's/"linkage group LG16"/16/g' ${dir_in}/gene_list.txt
sed -i -e 's/"Un"/-9/g' ${dir_in}/gene_list.txt
sed -i -e 's/"mitochondrion MT"/-9/g' ${dir_in}/gene_list.txt
sed -i -e 's/"//g' ${dir_in}/gene_list.txt
awk '!seen[$0]++' ${dir_in}/gene_list.txt > ${dir_in}/tmp_gene_list && mv ${dir_in}/tmp_gene_list ${dir_in}/gene_list.txt
sed -i '/^-9/d' ${dir_in}/gene_list.txt

jobnum=()
for i in ${list_pop[@]}
do
	echo ${i}
	for j in ${pheno_list[@]}
		echo ${j}
		n=$(grep 'number of analyzed individuals' ${dir_out}/gemma_${i}_${j}_freq_lmm_egs_cov.log.txt | cut -d '=' -f2)
		awk -v var="${n}" '{print $2" "$5" "$6" "$7" "$8" "$9" "$12" "var}' ${dir_out}/gemma_${i}_${j}_freq_lmm_egs_cov.assoc.txt > ${dir_out}/${j}_in_${i}.ma
		var="SNP A1 A2 freq BETA SE P N"
		sed -i "1s/.*/$var/" ${dir_out}/${j}_in_${i}.ma
		job=`sbatch --mem=8G --wrap="./gcta-1.94.1 --bfile ${i} --mBAT-combo ${dir_out}/${j}_in_${i}.ma --mBAT-gene-list ${dir_in}/gene_list.txt --out ${dir_out}/${j}_in_${i} --thread-num 1 --diff-freq 1 --mBAT-wind 5 --fastBAT-ld-cutoff 0.8 --mBAT-print-all-p"`
		jobnum+=(${job})
	done
done
pat=$(echo ${jobnum[@]}|tr " " "|")
x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
while [ ${x} -gt 1 ]
do
	sleep 30m
	x=$(squeue -u seynard | grep -Eow "$pat" |wc -l)
done
```
