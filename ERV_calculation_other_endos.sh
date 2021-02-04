#!/bin/bash

#$ -l h_vmem=12.5G
#$ -l tmem=12.5G
#$ -l h_rt=06:00:00
#$ -t 1-11
#$ -j y
#$ -cwd
#$ -S /bin/bash


# script to estimate SNP-heritability for an endophenotype and its genetic
# correlation with schizophrenia/psychosis from several datasets and perform a
# meta-analysis

# we first need to generate the Genetic Relatedness Matrices (GRM) for each endophenotype

# set the paths and prefixes of the binary PLINK files for PEIC

peic_path="/SAN/neuroscience/PEIC/projects/EEG_genetics_consortium/genetic_data/PEIC/hard_called_genotypes"
peic_prefix="pe.dbSNP.150.matched"

# generate a vector with the names of the endophenotypes

endo_names=("digit_symbol" "digit_span_forward" "full_IQ" "BLOCK_PERCENT" "RAVLT_imm_corr" "RAVLT_del_corr" "DIGIT_PERCENT" "WBV_KL" "LVV_KL" "P300A_KL" "P300L_KL" )

# get number of cores

cores=`nproc`
cores=$(($cores-1))

# extract the name of the endophenotype for which to run the analysis

endopheno=${endo_names[$SGE_TASK_ID - 1]}

# generate a list of pruned set of SNPs, use it to create a GRM-s for each dataset 
# and run a bivariate-GREML with GCTA

# extract the firs two columns (FID and IID) of the phenotypic table for the endophenotype so they can
# be used to limit the SNP-pruning and the GRM calculation to those samples

# for i in {0..10}
# do

# endopheno=${endo_names[$i]}

awk -F "\t" '{print $1, $2}' ./ERV_calculation/endo_tables/${endopheno}_phenotype.txt > ./ERV_calculation/endo_tables/${endopheno}.samples 

# generate pruned set of SNPs

plink --bfile $peic_path/$peic_prefix --indep-pairwise 50 5 0.2 --threads $cores --keep ./ERV_calculation/endo_tables/${endopheno}.samples --out ./ERV_calculation/pruned_snps/${endopheno}  

# create GRM

gcta --bfile $peic_path/$peic_prefix --thread-num $cores --make-grm --extract ./ERV_calculation/pruned_snps/${endopheno}.prune.in --keep ./ERV_calculation/endo_tables/${endopheno}.samples --out ./ERV_calculation/grms/${endopheno}

# run bivariate-GREML to calculate the SNP-heritability of both traits (MMN Amplitude at FZ and psychosis yes/no) and the genetic correlation between them

gcta --reml-bivar --grm ./ERV_calculation/grms/${endopheno} --pheno ./ERV_calculation/endo_tables/${endopheno}_phenotype.txt --covar ./ERV_calculation/endo_tables/${endopheno}_categorical_covariates.txt --qcovar ./ERV_calculation/endo_tables/${endopheno}_quantitative_covariates.txt --out ./ERV_calculation/erv_results/${endopheno} --reml-maxit 1000 --reml-bivar-prevalence 0.003

# done
