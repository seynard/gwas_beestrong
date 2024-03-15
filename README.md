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
- [5. Simple GWAS](#6-GWAS)
- [6. Meta GWAS](#7-MetaGWAS)
	- [6.1. Define groups](#71-Define-groups)
	- [6.2. MANTRA](#72-MANTRA)
	- [6.3. Mash](#73-mash)
- [7. Linkage disequilibrium](#5-LD)
- [8. Analysis on genes](#Genes)
	- [8.1. Marker vs genes correspondence](#81-Marker-vs-genes-correspondence)

<!-- /TOC -->

## 1. Introduction
This repository contains the scripts to perform the analysis described in the paper: "Sequence-based genome-wide association studies reveal the polygenic architecture of <i>Varroa destructor</i> resistance in Western honey bees <i>Apis mellifera</i>"\
Sonia E Eynard, Fanny Mondet, Benjamin Basso, Olivier E Bouchez, Tabatha Bulach, Yves Le Conte, Benjamin Dainat, Axel Decourtye, Lucie Genestout, Matthieu Guichard, Francois Guillaume, Emmanuelle Labarthe, Barbara Locke, Rachid Mahla, Joachim R de Miranda, Markus Neuditschko, Florence Phocas, Yannick Poquet, Christina Sann, Remi-Felix Serre, Kamila Tabet, Alain Vignal, Bertrand Servin \
bioRxiv , doi: https://doi.org/10.1101/2024.02.16.580755\

Using pool sequence data we performed a genome wide association study for multiple traits linked to the infestation by <i>Varroa destructor</i> of the honeybee and its resistance (measured as behavioural response). The GWAS was first done for each trait, in each of the identified population. Thereafter, we performed metaGWAS, for each trait, across the populations using two different Bayesian methods: the transethnic meta-analysis method from MANTRA (Morris 2011, https://doi.org/10.1002/gepi.20630) and the statistical method mash (Urbut et al. 2019, https://doi.org/10.1038/s41588-018-0268-8). In addition, we looked into genes associated with the markers that were found significant, their function and how the variant change affected them using the Ensembl's tool Variant Effect Predictor (VEP). 

This analysis was performed using the latest genome version, AmelHAV3.1 (Wallberg et al. 2019, https://doi.org/10.1186/s12864-019-5642-0, GenBank accession number: GCA_003254395.2). Only SNPs that were previously identified with high confidence were used for the rest of the analysis. The list of these SNPs, and the associated vcfs 870 and 628 individuals, were made available by Wragg et al. (2022, https://doi.org/10.1111/1755-0998.13665). In addition, the data from Wragg et al. (2022, https://doi.org/10.1111/1755-0998.13665) were used to describe reference populations for the three predominant sub-species of Europe, <i>Apis mellifera mellifera</i>, <i>Apis mellifera ligustica & carnica</i> and <i>Apis mellifere caucasia</i>.
A step of genotype reconstruction for the honeybee queen was performed following the protocol described in Eynard et al (2022, https://doi.org/10.1111/1755-0998.13685). The statistical models to estimate genetic composition and reconstruct queen genotype are available at the GitHub repository : https://github.com/BertrandServin/beethoven.

The pipeline described below is available in the script run.sh. All scripts called by this pipeline are available in this repository.

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
dir="~/GWAS" # initial repository
script="~/GWAS/scripts" # scripts
dir_out="~/GWAS/results" # results
dir_in="~/GWAS/data" # data
dir_save="~/save" # save (containing data generally used and in storage)
fasta=${dir_save}/Fasta/GCF_003254395.2_Amel_HAv3.1_genomic.fna # fasta file for the genome assembly used, in our case AmelHAV3.1 (Wallberg et al. 2019,  https://doi.org/10.1186/s12864-019-5642-0)
vcf_name="MetaGenotypesCalled870_raw_snps" # name of the vcf for high confidence SNPs described in Wragg et al. 2022 ( https://doi.org/10.1111/1755-0998.13665)
vcf_sansfiltre=${dir_in}/${vcf_name}.vcf.gz # initial vcf, full, without filters
vcf_file=${dir_in}/${vcf_name}_allfilter.vcf # final vcf, after filters, provided by Wragg et al. 2022 ( https://doi.org/10.1111/1755-0998.13665)
pop_id="Mellifera,Caucasia,Ligustica_Carnica" # reference populations of interest for our study
n_pop=$(echo $(IFS=","; set -f; set -- $pop_id; echo $#)) # nb of reference populations of interest
ppop_id2="Mellifera,Caucasia,Ligustica_Carnica,hybrid" # populations on which the analysis is run
pheno="Data_BS.csv" # phenotype file
pheno_mito='compareDepthsBeeVarroa.txt' # sequencing summary file, containing for each pool sequencing experiment the sequencing depth for honeybee, for varroa and its ratio and standardised values
snp50k='HAV3_1_50000.txt' # list of 50K SNP selected to represent honeybee reference population and used for genetic composition estimation
maf_min=0.01 # MAF filter (>0.01)
missing_rate=0.05 # missing rate filter (<0.05)
compo_threshold=0.8 # genetic composition threshold to group pool sequencing experiments into homogeneous groups
unif_threshold=0.99 # genetic composition threshold to build reference populations
d_min=10 # minimum sequencing depth
d_max=50 # maximmum sequencing depth
pool_size_def=150 # definition of pool size
ncpu=10 # nb of cpu used
nbjobs=10 # nb of jobs run
```

### 2.3. Prep data for reference populations
As described above, we used the data available from Wragg et al. (2022, 10.1111/1755-0998.13665) to create reference populations for the major European honeybee sub-species. These reference population are used thereafter to estimate allele frequency that are usefull to estimate genetic composition for each of our pool sequencing experiments.
```bash
sbatch -W --mem=100G --job-name='prep_seqapipop' --wrap="${script}/prep_seqapipop.sh ${script} ${dir_in} ${vcf_file} ${n_pop} ${pop_id} ${unif_threshold} ${ncpu}"
```

### 2.4. Prep pool sequence data
The script prep_pileup.sh describe the protocol to go from BAM file to depth and count files in the framework of pool sequencing. Here we used tools available in popoolation2 (Kofler et al. 2011, 10.1093/bioinformatics/btr589 ), and the mpileup tool. We also applied some correction and filtering using taylor made scripts. The filnal output of this step is the acquisition of depth.txt and count.txt files providing, for each SNP (in rows) and each pool sequencing experiment (in columns) the depth of sequencing and the count of allele (reference or alternative). Using these files it is possible to estimate for each pool sequencing experiment the allele frequency observed in the pool.   
```bash
sbatch -W --mem=200G --job-name='prep_pileup' --wrap="${script}/prep_pileup.sh ${script} ${dir_in} ${fasta} ${vcf_name} ${vcf_sansfiltre} ${vcf_file} ${snp50k} ${nbjobs} ${d_min} ${d_max} ${pool_size_def}"
```

## 3. Genetic composition estimation and queen genotype reconstruction (script for the statistical models available at https://github.com/BertrandServin/beethoven)
This step of our pipeline follows the procedure described in Eynard et al. (2022,  https://doi.org/10.1111/1755-0998.13685)
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
In this study we focused on multiple phenotypes linked to the infestation by varroa and the resistance behaviours expressed by the honeybee. Four phenotypes for varroa infestation were measured: the varroa infestation in the brood, on the adult honeybee (so called phoretic varroa), using the sequencing depth of varroa mitochondrial DNA in the pool as a proxy for its infestation on adult honeybee and a complete varroa load. After statistical analysis we observed that these phenotypes are highly correlated and all represent the first axis of the principal component on phenotypes, with about 60 % of the variance explained. We therefore decided to use the first axis of the PCA as a composite phenotype representing overall varroa infestation. The two other phenotypes are linked to 'resistance' mechanisms of the honeybee towards this infestation. Colonies were scored for mite non reproduction (MNR, described in Eynard et al. (2020, https://doi.org/10.3390/insects11080492) and Mondet et al. (2020,  https://doi.org/10.3390/insects11090595)) and for recapping of infested brood cells. The phenotypes were adjusted to match a Normal distribution using Empiral Bayes method for MNR and logit transformation for the recapping.  
```bash
sbatch -o ${dir}/log/prep_pheno.out  -e ${dir}/log/prep_pheno.err -W --wrap="Rscript ${script}/prep_pheno.r ${dir_in} ${dir_out} ${pheno} ${pop_id2} ${pheno_mito}"
```

## 5. Simple GWAS
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

## 6. Meta GWAS
### 6.1. Define groups and phenotypes
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

### 6.2. MANTRA
Meta-GWAS performed across groups for each phenotypes using the MANTRA method described in Morris (2011,  https://doi.org/10.1002/gepi.20630)
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

### 6.3. Mash
Meta-GWAS performed across groups for each phenotypes using the mash method described in Urbut et al. (2019,  https://doi.org/10.1038/s41588-018-0268-8)
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

## 7. Linkage disequilibrium
Linkage disequilibrium is a crucial information when one wants to correlate significant markers for a specific trait and their underlying genes. Using genotype data from Wragg et al. (2022,  https://doi.org/10.1111/1755-0998.13665) we estimated linkage disequilibrium between SNP across the 16 honeybee chromosomes for groups of similar genetic composition than our groups for GWAS ('pure' mellifera, 'pure' ligustica & carnica and hybrides).
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
		job=`sbatch --mem=200G --wrap="plink --bfile ${dir_in}/ld/${i} --chr $c --keep-allele-order --make-bed --out ${dir_in}/ld/${i}_${c}; plink --bfile ${dir_in}/ld/${i}_${c} --keep-allele-order --ld-window-r2 0.1 --out ${dir_in}/ld/ld_${i}_${c} --r2 inter-chr"`
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

All the results of these different analysis were combined and discussed to underline the genetic mechanisms involved in varroa infestation and resistance to such infestation in a heterogeneous honeybee population. Difference between the different genetic groups were also highlighted and discuss to provide guidance for further potential development in the selection for resistant bees.
