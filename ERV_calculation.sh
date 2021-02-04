#!/bin/bash

# script to estimate SNP-heritability for an endophenotype and its genetic
# correlation with schizophrenia/psychosis from several datasets and perform a
# meta-analysis

# we first need to generate the Genetic Relatedness Matrices (GRM) for each dataset

# set the paths and prefixes of the binary PLINK files for each dataset

peic_path="/SAN/neuroscience/PEIC/projects/EEG_genetics_consortium/genetic_data/PEIC/hard_called_genotypes"
peic_prefix="pe.dbSNP.150.matched"

maryland_path="/SAN/neuroscience/PEIC/projects/EEG_genetics_consortium/genetic_data/Elliot/imputed_data/Elliot2.vcfs/imputed_plink/matched_files"
maryland_prefix="MPRC_Hong.typed.and.imputed.clean.dbSNP.150.matched"

harvard_path="/SAN/neuroscience/PEIC/projects/EEG_genetics_consortium/genetic_data/Mei/imputed_data/imputed_plink/matched_files"
harvard_prefix="mclean_fgh19.typed.and.imputed.dbSNP.150.matched"

# put those in string arrays/vectors

gene_paths=($peic_path $maryland_path $harvard_path)
gene_prfxs=($peic_prefix $maryland_prefix $harvard_prefix)

# generate another vector with the names of the dataset

ds_names=("PEIC" "Maryland" "Harvard")

# in a loop, generate a list of pruned set of SNPs, use it to create a GRM-s for each dataset 
# and run a bivariate-GREML with GCTA

for i in {0..2}
do

	# extract the path to and the prefix of the corresponding genetic dataset

	gene_path=${gene_paths[i]}
	gene_prfx=${gene_prfxs[i]}
	
	# extract the firs two columns (FID and IID) of the phenotypic table for the dataset so they can
	# be used to limit the SNP-pruning and the GRM calculation to those samples

	awk -F "\t" '{print $1, $2}' ./ERV_calculation/${ds_names[i]}_phenotypes.txt > ./ERV_calculation/${ds_names[i]}.samples 

	# generate pruned set of SNPs

	plink --bfile $gene_path/$gene_prfx --indep-pairwise 50 5 0.2 --keep ./ERV_calculation/${ds_names[i]}.samples --out ./ERV_calculation/${ds_names[i]}  

	# create GRM

	gcta --bfile $gene_path/$gene_prfx --make-grm --extract ./ERV_calculation/${ds_names[i]}.prune.in --keep ./ERV_calculation/${ds_names[i]}.samples --out ./ERV_calculation/${ds_names[i]}

	# run bivariate-GREML to calculate the SNP-heritability of both traits (MMN Amplitude at FZ and psychosis yes/no) and the genetic correlation between them

	gcta --reml-bivar --grm ./ERV_calculation/${ds_names[i]} --pheno ./ERV_calculation/${ds_names[i]}_phenotypes.txt --covar  ./ERV_calculation/${ds_names[i]}_categorial_covariates.txt --qcovar ./ERV_calculation/${ds_names[i]}_quantitative_covariates.txt --out ./ERV_calculation/${ds_names[i]} --reml-maxit 1000 --reml-bivar-prevalence 0.003

done






