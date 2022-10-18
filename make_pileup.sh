#! /bin/bash
###### Script to produce pileup (prerequisite to perform allele frequency estimation from poolseq with Popoolation2) ######

module load bioinfo/samtools-1.8
module load system/Java8

#dir_out='/work/project/dynagen/seynard/GWAS/BeeStrongHAV3_1'
#dir_popoolation='/usr/local/bioinfo/src/PoPoolation2/popoolation2_1201'
#fasta_seq='/genphyse/dynagen/BeeStrong/Fasta/GCF_003254395.2_Amel_HAv3.1_genomic.fna' 
#vcf=''/work/project/cytogen/Alain/seqapipopOnHAV3_1/combineGVCFs/The870vcf/MetaGenotypesCalled870_raw_snps_allfilter.vcf'
#script='/home/seynard/scriptGWAS'
dir_out=$1
dir_popoolation=$2
fasta=$3 
vcf_position=$4
script=$5
bam_file=$6

n=$7
output_pileup=${n}.pileup
sbatch -W -J pileup_${n} -o ${dir_out}/log/pileup_${n}.out -e ${dir_out}/log/pileup_${n}.err --wrap="module load bioinfo/samtools-1.8; samtools mpileup -I -l ${dir_out}/snp_pos.txt -f ${fasta} -C 50 -q 20 -Q 20 ${bam_file} -o ${dir_out}/${output_pileup}"
output_sync=${n}.sync
sbatch --mem=40G -W -J sync_${n} -o ${dir_out}/log/sync_${n}.out -e ${dir_out}/log/sync_${n}.err --wrap="java -ea -Xmx10g -jar ${dir_popoolation}/mpileup2sync.jar --fastq-type sanger --min-qual 20 --input ${dir_out}/${output_pileup} --output ${dir_out}/${output_sync}"
sed -i "s/a/A/g;s/t/T/g;s/c/C/g;s/g/G/g" ${dir_out}/${n}.sync
echo 'DONE pileup and sync'
if [ -e ${dir_out}/${n}.count ]	
then
echo 'colony ${n} already in' 
else 
module load system/Python-3.6.3
python ${script}/combine_col.py ${dir_out}/allele_id.txt ${dir_out}/${n} ${dir_out}/${n}.count
awk '$2=="'"${n}"'"{$(NF+2)=" done"} 1' ${dir_out}/BamList > ${dir_out}/tmp && mv ${dir_out}/tmp ${dir_out}/BamList
fi
rm ${dir_out}/${n}.sync
zip -j ${dir_out}/${n}.pileup.zip ${dir_out}/${n}.pileup
#tar -czvf ${dir_out}/${n}_pileup.tar.gz ${dir_out}/${n}.pileup
rm ${dir_out}/${n}.pileup


