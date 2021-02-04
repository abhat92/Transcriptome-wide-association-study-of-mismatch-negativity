#!/bin/bash

#set prefix of plink-format dataset (.bed, .bim, .fam)
#if files not in working directory, provide full path

dataset_in="McLean_OmniExpress_Beadstudio.FHG19"

dataset="mclean_fhg19"

#plink files in

#prefix="/cluster/project9/PE/data/QCed_plink_files/pe19"

#start with QC-in the data

#do initial filtering at both SNP and sample levels;
#for SNPs, do the filtering based on call rate, MAF, HWE departure and mendel error rate
#for samples, remove those with low call rate; the sex check and inbreeding-coeff based filtering will be applied later

# genotypes are not sorted; sort them by just applying make-bed

plink --bfile $dataset_in --make-bed --out $dataset.sorted

# rename duplicated sample in FAM file

Rscript --vanilla rename_duplicated_sample_IDs.R

# remove old $dataset.sorted.fam file and rename the $dataset.renamed.fam file to $dataset.sorted.fam

rm $dataset.sorted.fam
mv $dataset.renamed.fam $dataset.sorted.fam

# for Mei's dataset, remove those SNPs labelled as being in chromosome 0 in the BIM file, that is, those that are unplaced.

# for now, we will keep those in chromosome 23 (X) so we can do the sex check

# don't apply a Hardy-Weinberg Equilibrium based filtering now; it can be a problem with the X chromosomes; we can do it after removing the X chromosome variants 

# there are 3 individuals (parents + daughter) that have multiple entries (18 entries in total); 
# In the "mclean_fhg19.fam" file the FIDs for all 18 entries have been renamed to "dupliFAM" and the IIDs of the 15 redudant IIDs into IID.duplicate.number; the 15 redudant entries will be removed  

# split the x-chromosome for the sexcheck afterwards

plink --bfile $dataset.sorted --geno 0.05 --maf 0.01 --mind 0.05 --me 0.05 0.1 --make-bed --out $dataset.clean.0 --chr 1-23 --remove duplicated_samples_for_removal --split-x hg19

# 29-05-2018: I have finally managed to remove the 15 redundant entries for the 3 samples from the dataset; so now I can continue with the rest of the QC.

# let's do that tomorrow

#infer sex from genetic data and compare to report

plink --bfile $dataset.clean.0 --check-sex 0.6 0.6 --out $dataset.clean

# identify those samples that show a mismatch between X-chromosome genotypes-inferred sex and reported sex and removed them

# keep those for which the reported sex was ambiguous

grep PROBLEM  mclean_fhg19.clean.sexcheck | grep -vw 0 | awk -F " " '{ print $1, $2}' > $dataset.clean.sex.mismatch.samples

plink --bfile $dataset.clean.0 --make-bed --remove $dataset.clean.sex.mismatch.samples --out $dataset.clean.1

# impute sex for all samples using genotypes (there are some samples left with ambiguous sex in the FAM file)

plink --bfile $dataset.clean.1 --impute-sex 0.6 0.6 --make-bed --out dataset.clean.1

#calculate inbreeding coefficients

plink --bfile $dataset.clean.1 --het --out $dataset.clean.1.inbreeding

#identify individuals with an inbreeding coefficient above 0.1 or below -0.1

awk '{if ($6 >= 0.1 || $6 <= -0.1) { print $0 }}' $dataset.clean.1.inbreeding.het | sed '1d' | awk -F " " '{ print $1, $2}' > $dataset.inbred.samples.0

wc -l $dataset.inbred.samples | awk -F " " '{ print $1 }' > inbred.ind.count.tmp

inbred_count=`cat inbred.ind.count.tmp`

if [ "$inbred_count" == "0" ];then

echo "No individuals with high level of inbreeding (absolute departure from expected heterozygosity > 0.1) detected"

else

echo "$inbred_count individuals with high inbreeding coefficient identified"

# remove individuals with either too-high or too-low inbreeding coefficient

# also, get rid of the X chromosome SNPs; we won't be using anymore

plink --bfile $dataset.clean.1 --chr 1-22 --remove $dataset.inbred.samples --make-bed --out $dataset.clean.2

exit

fi

# now that we have got rid of the X chromosome, we can apply the HWE threshold

plink --bfile $dataset.clean.2 --hwe 0.000001 --make-bed --out $dataset.QCed.0

# remove files generated in the intermediate steps

#update the Family IDs in the .fam file we have just generated using the new Family IDs Johan generated from the genealogies
#those are stored in the PEIC phenotypic table

# plink --bfile $dataset.QCed --update-ids pe19_peic_ids --out $dataset.QCed.1 --make-bed

#extract the whole list of SNPs from .bim file so they can be annotated in the UCSC genome browser

# awk -F "\t" '{print $2}' $dataset.clean.bim > $dataset.clean.snp.list

#for a mixed model association analysis where a kinship matrix will be included, there is no need to perform a PCA analysis because
# a) no additional adjustment is necessary to account for population stratification and b) we won't be removing ethnical outliers

#the next step is to generate a kinship matrix using LDAK to a) remove samples with very high kinship, b) use it to estimate
#trait heritabilities and c) feed it to the mixed association model

#generate thinned set of SNPs for LD

ldak --bfile $dataset.QCed --thin prune --window-prune 0.2 --window-kb 1000

#generate kinship matrix

ldak --calc-kins-direct prune --bfile $dataset.QCed --ignore-weights YES --power -1 --extract prune.in --kinship-raw YES

#identify duplicates and pairs of samples that show very high kinship values (> 0.95) and remove one sample for each pair

ldak --bfile $dataset.QCed --filter dups --grm prune --max-rel 0.95

#calculate SNP weights; this takes a lot of time!!

ldak --bfile $dataset.clean --cut-weights sections --keep dups.keep

ldak --bfile $dataset.clean --calc-weights-all sections --keep dups.keep

#estimate heritability of the binary case-control variable to be used in the calculation of the 
#Endophenotype Ranking Value (ERV) for each endophenotype

#call R script to extract list of case control samples that are included in the .fam file of the QC'ed dataset
#also create covariate table for case control analysis

case_con_covs="age gender center"

Rscript --vanilla extract_case_control.R $dataset $phenotype $case_con_covs

#generate weighted kinship matrix for those samples

ldak --calc-kins-direct kins.case.con --bfile $dataset.clean --weights ./sections/weights.all --power -0.25 --kinship-raw YES --keep case_control_ids

#perform Haseman Elston Regression for the binary case-control trait
#first adjust kinship matrix for covariates

ldak --adjust-grm kins.case.con2 --grm kins.case.con --covar case.con.cov.dummy

#The kinship matrix adjustemnt gives an error; for now, just estimate heritability with the unadjusted matrix

#Perform Has_Els regression

ldak --he case_control_he --grm kins.case.con2 --covar case.con.cov.dummy --pheno case_controls --prevalence 0.035 --kinship-details NO

