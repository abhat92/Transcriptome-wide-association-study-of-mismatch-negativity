#$ -l h_vmem=7.5G
#$ -l tmem=7.5G
#$ -l h_rt=30:00:00
#$ -t 1-22
#$ -j y
#$ -cwd
#$ -S /bin/bash


### Script to generate frequency and info files from imputed data
pe=/SAN/neuroscience/PEIC/projects/EEG_genetics_consortium/genetic_data/Elliot/imputed_data/Elliot2.vcfs
vcftools=/share/apps/genomics/vcftools_0.1.12b/bin/vcftools
bcftools=/share/apps/genomics/bcftools-1.3.1/bin/bcftools
user=/SAN/neuroscience/PEIC/projects/EEG_genetics_consortium/genetic_data/Elliot/imputed_data/Elliot2.vcfs

## Merged pchip and pe
mkdir -p $user/info
 	$vcftools --gzvcf $pe/$SGE_TASK_ID.vcf.gz --get-INFO INFO --out $user/info/$SGE_TASK_ID
 	$vcftools --gzvcf $pe/$SGE_TASK_ID.vcf.gz --freq2 --out $user/info/$SGE_TASK_ID
 	$vcftools --gzvcf $pe/$SGE_TASK_ID.vcf.gz --hardy --out $user/info/$SGE_TASK_ID
        $bcftools query -f '%ID\n' $pe/$SGE_TASK_ID.vcf.gz -o $user/info/$SGE_TASK_ID.rsid.tmp
	echo "SNP_ID" | cat - $user/info/$SGE_TASK_ID.rsid.tmp > $user/info/$SGE_TASK_ID.rsid
	echo merging imputation info for ... $user/info/$SGE_TASK_ID.snpinfo
	echo -e "CHROM:POS:REF:ALT\tSNP_ID\tINFO\tOBS(HOM1/HET/HOM2)\tP_HWE\tFREQ_A1\tFREQ_A2" > $user/info/$SGE_TASK_ID.snpinfo
	paste <(awk 'BEGIN {OFS=":"} {print $1,$2,$3,$4}' $user/info/$SGE_TASK_ID.INFO) <(cut -f1 $user/info/$SGE_TASK_ID.rsid) <(cut -f5 $user/info/$SGE_TASK_ID.INFO) <(cut -f3,6 $user/info/$SGE_TASK_ID.hwe) <(cut -f5,6 $user/info/$SGE_TASK_ID.frq) | sed '1d' >> $user/info/$SGE_TASK_ID.snpinfo

	## Cleanup
	#rm $pe/info/*.INFO
	#rm $pe/info/*.hwe
	#rm $pe/info/*.log
	#rm $pe/info/*.frq
	#rm $pe/info/*.rsid*

# ## List all multialleic SNPs
# rm /home/molpsych/data/snp_info/all.multiallelic.txt
# for i in {1..22}
# do
#     echo finding multiallelic SNPs on chr $SGE_TASK_ID
#     awk '{if ($2 !=".") print $2}' /home/molpsych/data/snp_info/$SGE_TASK_ID.all.snpinfo | sort | uniq -c | awk '{if ($1>1) print $2}' >> /home/molpsych/data/snp_info/all.multiallelic.txt
# done


