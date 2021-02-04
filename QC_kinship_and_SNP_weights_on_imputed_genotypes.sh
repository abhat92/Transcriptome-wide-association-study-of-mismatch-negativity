#!/bin/bash
#$ -l h_vmem=15.5G
#$ -l tmem=15.5G
#$ -l h_rt=08:00:00
#$ -j y
#$ -cwd
#$ -S /bin/bash

# set path to imputed data

imputed_path=/SAN/neuroscience/PEIC/projects/EEG_genetics_consortium/genetic_data/Mei/imputed_data/imputed_plink

# set path to genotyped data

genotyped_path=/SAN/neuroscience/PEIC/projects/EEG_genetics_consortium/genetic_data/Mei

# add path to plink executable

plink_path=/share/apps/plink-1.90

#######################
## MERGE PLINK FILES ##
#######################

# merge binary plink files that contain hrc-imputed data for each chromosome

$plink_path/plink --bfile $imputed_path/1 --merge-list $imputed_path/plink_file_prefixes_by_chr --make-bed --out pe.hrc.imputed

# check if the SNPs for which more than one position has been observed are among the typed variants

# extract multi-position SNP list

grep "Warning: Multiple" pe.hrc.imputed.log | sed "s/.*'rs/rs/" | sed "s/'.//" > multi_position_snps

# sort both multi position SNP list and typed snp list

sort multi_position_snps > multi_position_snps_sorted.tmp
sort typed_snps_list > typed_snps_list_sorted.tmp

# match the multi_position_snps list with the list of typed SNPs

comm -12 multi_position_snps_sorted.tmp typed_snps_list_sorted.tmp > multi_position_typed_snps

multi_typed=`wc -l multi_position_typed_snps`

echo "$multi_typed SNPs seen to present several positions in imputed data are among the typed variants"

rm *.tmp

# first run yields an error indicating that there are 1993 SNPs with 3 or more variants, even if the "biallelic-only strict" flag was used when converting the ".vcf" files into binary plink files.
# PLINK suggests this may be due to strand inconsistencies and that a strand flip for those 1993 SNPs (-flip *.missnp) could solve the problem
# however, it is probably best to just remove these


# As the first run didn't manage to remove all the multi-allelic SNPs, we ran the removal script one more time 
# THis was on 25.09.18

for i in {1..22}
do

plink --bfile $imputed_path/$i.no.missnp --exclude pe.hrc.imputed-merge.missnp --make-bed --out $i.no.missnp.2

done

# create file with prefixes of new binary plink files for chromosomes 2 to 22

ls | grep no.missnp.2 | sed 's/missnp.2.*/missnp.2/' | uniq | sort | grep -vw "1.no.missnp.2" > plink.binaries.without.missnp.2 

# now merge the binary files from where the multi-allelic SNPs have been removed

$plink_path/plink --bfile $imputed_path/1.no.missnp.2 --merge-list $imputed_path/plink.binaries.without.missnp.2 --make-bed --out pe.hrc.imputed.2

# remove binary plink files for individual chromosomes

rm *.no.missnp*

# remove genotyped SNPs from merged imputed genotype data plink files

plink --bfile pe.hrc.imputed.2 --exclude typed_snps_list --make-bed --out no.typed.hrc.imputed

# delete binary plink files that contain all SNPs, including those that were typed

rm pe.hrc.imputed.2.bed pe.hrc.imputed.2.bim pe.hrc.imputed.2.fam

#update FAM file so it matches the Family IDs in the PEIC phenotypic table and the FAM file of the p19.clean plink fileset

# plink --bfile no.typed.hrc.imputed --update-ids hrc.impute_peic_ids --make-bed --out no.typed.hrc.imputed.updated

# remove plink files with old Family IDs (plink adds a "~" at the end)

# rm no.typed.hrc.imputed.bed no.typed.hrc.imputed.bim no.typed.hrc.imputed.fam

# merge binary plink files that contain typed genotypic data with those that contain imputed genotypic data

plink --bfile $genotyped_path/mclean_fhg19.QCed --bmerge no.typed.hrc.imputed --make-bed --out mclean.typed.and.imputed

# delete the files that contain imputed data only

rm  no.typed.hrc.imputed.updated.bed no.typed.hrc.imputed.updated.bim no.typed.hrc.imputed.updated.fam

# Generally, there is a bunch of pairs of SNPs that have the exact same position
# In most cases, this happens because some of the typed SNPs that are imputed are not assigned an rsID, but an ID based on the chromosome, the position and the reference and alternative alleles and, thus, the imputed and typed SNPs, although they are the same variant appear as two different variants with the same position. 
# The best way to go is to extract the list of IDs of the imputed versions of the SNPs ("chr:position:refA:altA")

grep "Warning: Variants" pe.typed.and.imputed.log | sed "s/.*and '//" | sed "s/' have.*//" > duplicated_imputed_snps
 
plink --bfile pe.typed.and.imputed --exclude duplicated_imputed_snps --make-bed --out mclean_fgh19.typed.and.imputed

# remove the plink files that include the duplicated imputed SNPs

rm pe.typed.and.imputed.bed pe.typed.and.imputed.bim pe.typed.and.imputed.fam

# set the name of the dataset as a variable

dataset="mclean_fgh19.typed.and.imputed"

#start with QC-in the data

#do initial filtering at both SNP and sample levels;
#for SNPs, do the filtering based on call rate, MAF, HWE departure and mendel error rate
#for samples, remove those with low call rate; the sex check and inbreeding-coeff based filtering will be applied later
#also test for batch/imputation effects on the missingness (--test-missing)

plink --bfile $dataset.2 --geno 0.05 --maf 0.01 --hwe 0.000001 --mind 0.05 --me 0.05 0.1 --test-missing mperm=1000 midp --make-bed --out $dataset.clean

# delete unQC'ed plink files

rm $dataset.2.bed $dataset.2.bim $dataset.2.fam

# identify those SNPs with significantly different missingness rates in controls and patients and remove them

awk '{if ($5 < 0.000005) { print $2 }}' pe.typed.and.imputed.clean.missing > differential.missingness.snps

diff_miss_snps=`wc -l differential.missingness.snps | awk -F " " '{print $1}'`

echo "$diff_miss_snps SNPs show significantly different call rates between cases and controls at p value < 5e-06"

# remove the SNPs that show significantly different call rates between cases and controls from the plink dataset

plink --bfile $dataset.clean --exclude differential.missingness.snps --make-bed --out $dataset.clean

# remove previous version of plink files

rm $dataset.clean*~

##infer sex from genetic data and compare to report

##A priori, we would need to infer sex from genetic data and compare it to reported sex as a sample-level QC measure
##However, this dataset was previouly QC'ed and SNPs from sex chromosomes are missing, so we cannot do the sex check in this case

##plink --bfile $dataset.clean --check-sex --out $dataset.clean

#calculate inbreeding coefficients

plink --bfile $dataset.clean --het --out $dataset.clean.inbreeding

#identify individuals with an inbreeding coefficient above 0.1

awk '{if ($6 >= 0.1 || $6 <= -0.1) { print $2 }}' $dataset.clean.inbreeding.het | sed '1d' > $dataset.inbred.samples

wc -l $dataset.inbred.samples | awk -F " " '{ print $1 }' > inbred.ind.count.tmp

inbred_count=`cat inbred.ind.count.tmp`

if [ "$inbred_count" == "0" ];then

echo "No individuals with high level of inbreeding (absolute departure from expected heterozygosity > 0.1) detected"

else

echo "$inbred_count individuals with high inbreeding coefficient identified"

#take IDs of highly inbred individuals and put them in the proper format for the egrep call

cat $dataset.inbred.samples | tr '\n' '|' | sed 's/|*$//g' > egrep_ready_codes.tmp

egrep_codes=`cat egrep_ready_codes.tmp`

#get positions for rows where the inbred id codes and the individuals ids in the .fam file match

awk -F " " '{ print $2 }' $dataset.clean.fam > all_id_list.tmp

egrep_codes=`cat egrep_ready_codes.tmp`

egrep -n $egrep_codes all_id_list.tmp  | sed 's/:.*/p/g' | tr '\n' ';' | sed 's/;*$//g' >  inbred.id.positions.tmp

subset_lines=`cat inbred.id.positions.tmp`

#use those positions to extract inbred individuals from .fam file using sed

sed -n $subset_lines $dataset.clean.fam > $dataset.inbred.samples.fam

#remove individuals with high inbreeding coefficient from the dataset

plink --bfile $dataset.clean --remove $dataset.inbred.samples.fam --make-bed --out $dataset.clean1

#delete temporary files generated so far

rm *.tmp

#delete $dataset.clean files from initial QC step

rm $dataset.clean.*

#rename $dataset.clean1 files into $dataset.clean

mv $dataset.clean1.bed $dataset.clean.bed
mv $dataset.clean1.bim $dataset.clean.bim

exit

fi

#update the Family IDs in the .fam file we have just generated using the new Family IDs Johan generated from the genealogies
#those are stored in the PEIC phenotypic table

 plink --bfile $dataset.clean --update-ids pe19_peic_ids --out $dataset.clean --make-bed

#extract the whole list of SNPs from .bim file so they can be annotated in the UCSC genome browser

awk -F "\t" '{print $2}' $dataset.clean.bim > $dataset.clean.snp.list

#for a mixed model association analysis where a kinship matrix will be included, there is no need to perform a PCA analysis because
# a) no additional adjustment is necessary to account for population stratification and b) we won't be removing ethnical outliers

#the next step is to generate a kinship matrix using LDAK to a) remove samples with very high kinship, b) use it to estimate
#trait heritabilities and c) feed it to the mixed association model

#generate thinned set of SNPs for LD

ldak --bfile $dataset.clean --thin prune --window-prune 0.2 --window-kb 1000

#generate kinship matrix

ldak --calc-kins-direct prune --bfile $dataset.clean --ignore-weights YES --power -1 --extract prune.in --kinship-raw YES

#identify duplicates and pairs of samples that show very high kinship values (> 0.95) and remove one sample for each pair

ldak --bfile $dataset.clean --filter dups --grm prune --max-rel 0.95

plink --bfile $dataset.clean --remove dups.lose --make-bed --out $dataset.QCed

#calculate SNP weights; this takes a lot of time!!

ldak --bfile $dataset.clean --cut-weights sections --keep dups.keep

ldak --bfile $dataset.clean --calc-weights-all sections --keep dups.keep

#estimate heritability of the binary case-control variable to be used in the calculation of the 
#Endophenotype Ranking Value (ERV) for each endophenotype

#call R script to extract list of case control samples that are included in the .fam file of the QC'ed dataset
#also create covariate table for case control analysis

case_con_covs="age gender center"

Rscript --vanilla extract_case_control.R $dataset $case_con_covs

#generate weighted kinship matrix for those samples

ldak --calc-kins-direct kins.case.con --bfile $dataset.clean --weights ./sections/weights.all --power -0.25 --kinship-raw YES --keep case_control_ids

#perform Haseman Elston Regression for the binary case-control trait
#first adjust kinship matrix for covariates

# ldak --adjust-grm kins.case.con2 --grm kins.case.con --covar case.con.cov.dummy

#Perform Has_Els regression

# when using the adjusted kinship matrix, the heritability estimation goes to high (0.93 +- 0.09)

ldak --he case_control_he --grm kins.case.con --covar case.con.cov.dummy --pheno case_controls --prevalence 0.05 --kinship-details NO

###############################
## MAP TO LAST snpDV RELEASE ##
###############################

# create a copy of the dataset with the rsIDs of the last release of dbSNP
# (last complete release is v.150 as of 2018-04-24)

# generate coordinate-based IDs for the dataset using the information in the BIM file
# concatenate chromosome, basepair, A1 and A2

awk -F "\t" '{print $1, ":", $4, ":", $5, ":", $6}' $dataset.QCed.bim | tr -d ' ' > $dataset.clean.posIDs

# generate file with both position-based IDs and mixed position-based and rs IDs.

paste -d "\t" $dataset.clean.bim $dataset.clean.posIDs | awk -F "\t" '{print $2, $7}' > $dataset.clean.ID.tab

# assign position-based IDs to all SNPs in the dataset feeding the $dataset.clean.ID.tab file to the --update-name flag in plink

plink --bfile $dataset.clean --update-name $dataset.clean.ID.tab 2 1 --make-bed --out $dataset.clean.posIDs

# now assign rsIDs in the last release of the dbSNP using the position-(and allele-)based IDs for matching

# set path to directory that contains the table with position-based IDs and rsIDs

id_path="/SAN/neuroscience/PEIC/ref_files/UCSC/hg19"

# the reference and alternative alleles of dbSNP will be in the same order as the A1 and A2 alleles in our dataset for some SNPs, and in the opposite order for some other SNPs; thus, we will have to do the update in two steps, each time using one of the 2 columns with position and allele-based IDs for dbSNP (chr:position:ref:alt OR chr:position:alt:ref)

# using chr:position:ref:alt for matching 

plink --memory 9500 --bfile $dataset.clean.posIDs --update-name $id_path/dbSNP.150.hg19.tab.no.dups 3 1 --make-bed --out $dataset.clean.150.1

# using chr:position:alt:ref

plink --memory 9500 --bfile $dataset.clean.150.1 --update-name $id_path/dbSNP.150.hg19.tab.no.dups 3 2 --make-bed --out $dataset.clean.150

# assign reference allele in the annotations as A2 allele in the dataset

plink --bfile $dataset.clean.150 --a2-allele $id_path/dbSNP.150.hg19.tab.no.dups 4 3 --make-bed --out $dataset.clean.a2.rs.tmp

# assign alternative allele in the annotations as A1 allele in the dataset

plink --bfile $dataset.clean.a2.rs.tmp --a1-allele $id_path/dbSNP.150.hg19.tab.no.dups 5 3 --make-bed --out $dataset.dbSNP.150.matched
