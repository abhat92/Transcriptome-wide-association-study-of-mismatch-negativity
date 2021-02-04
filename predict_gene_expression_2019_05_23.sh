#$ -l h_vmem=7.5G
#$ -l tmem=7.5G
#$ -l h_rt=12:00:00
#$ -j y
#$ -cwd
#$ -S /bin/bash


### IMPUTE GENE EXPRESSION WITH PREDIXCAN ##

## PEIC ##

# dataset="PEIC"

# set path to FAM file

# fam_path="/SAN/neuroscience/PEIC/projects/EEG_genetics_consortium/genetic_data/PEIC/hard_called_genotypes/mmn_subset/pe.dbSNP.150.matched.mmn.subset.fam"

# run PrediXcan for each tissue

## Cortex
# ./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Brain_Cortex/gtex_v7_Brain_Cortex_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Brain_Cortex_Gene_Expression_Matrix

# Hypothalamus
# ./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Brain_Hypothalamus/gtex_v7_Brain_Hypothalamus_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Brain_Hypothalamus_Gene_Expression_Matrix

## Hippocampus
# ./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Brain_Hippocampus/gtex_v7_Brain_Hippocampus_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Brain_Hippocampus_Gene_Expression_Matrix

## Pituitary
# ./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Pituitary/gtex_v7_Pituitary_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Pituitary_Gene_Expression_Matrix

## Whole Blood
# ./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Whole_Blood/gtex_v7_Whole_Blood_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Whole_Blood_Gene_Expression_Matrix

# Brain ACC BA24
# ./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Brain_Anterior_cingulate_cortex_BA24/gtex_v7_Brain_Anterior_cingulate_cortex_BA24_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Brain_Anterior_cingulate_cortex_BA24_Gene_Expression_Matrix

## Amygdala
# ./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Brain_Amygdala/gtex_v7_Brain_Amygdala_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Brain_Amygdala_Gene_Expression_Matrix

## Frontal Cortex
# ./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Brain_Frontal_Cortex_BA9/gtex_v7_Brain_Frontal_Cortex_BA9_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Brain_Frontal_Cortex_BA9_Gene_Expression_Matrix

## DLPFC
# ./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights ./common_mind/CommonMindDB.v2/DLPFC_newMetax.db --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/DLPFC_Gene_Expression_Matrix


## MEI ##

dataset="Mei"

# set path to FAM file

fam_path="/SAN/neuroscience/PEIC/projects/EEG_genetics_consortium/genetic_data/Mei/imputed_data/imputed_plink/matched_files/mclean_fgh19.typed.and.imputed.dbSNP.150.matched.fam"

# run PrediXcan for each tissue

## Cortex
./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Brain_Cortex/gtex_v7_Brain_Cortex_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Brain_Cortex_Gene_Expression_Matrix

# Hypothalamus
./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Brain_Hypothalamus/gtex_v7_Brain_Hypothalamus_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Brain_Hypothalamus_Gene_Expression_Matrix

## Hippocampus
./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Brain_Hippocampus/gtex_v7_Brain_Hippocampus_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Brain_Hippocampus_Gene_Expression_Matrix

## Pituitary
./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Pituitary/gtex_v7_Pituitary_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Pituitary_Gene_Expression_Matrix

## Whole Blood
./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Whole_Blood/gtex_v7_Whole_Blood_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Whole_Blood_Gene_Expression_Matrix

# Brain ACC BA24
./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Brain_Anterior_cingulate_cortex_BA24/gtex_v7_Brain_Anterior_cingulate_cortex_BA24_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Brain_Anterior_cingulate_cortex_BA24_Gene_Expression_Matrix

## Amygdala
./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Brain_Amygdala/gtex_v7_Brain_Amygdala_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Brain_Amygdala_Gene_Expression_Matrix

## Frontal Cortex
./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Brain_Frontal_Cortex_BA9/gtex_v7_Brain_Frontal_Cortex_BA9_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Brain_Frontal_Cortex_BA9_Gene_Expression_Matrix

## DLPFC
./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights ./common_mind/CommonMindDB.v2/DLPFC_newMetax.db --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/DLPFC_Gene_Expression_Matrix

## ELLIOT ##

# dataset="Elliot"

# set path to FAM file

# fam_path="/SAN/neuroscience/PEIC/projects/EEG_genetics_consortium/genetic_data/Elliot/imputed_data/Elliot2.vcfs/imputed_plink/matched_files/MPRC_Hong.typed.and.imputed.clean.dbSNP.150.matched.fam"

# run PrediXcan for each tissue

## Cortex
# ./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Brain_Cortex/gtex_v7_Brain_Cortex_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Brain_Cortex_Gene_Expression_Matrix

# Hypothalamus
# ./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Brain_Hypothalamus/gtex_v7_Brain_Hypothalamus_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Brain_Hypothalamus_Gene_Expression_Matrix

## Hippocampus
# ./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Brain_Hippocampus/gtex_v7_Brain_Hippocampus_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Brain_Hippocampus_Gene_Expression_Matrix

## Pituitary
# ./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Pituitary/gtex_v7_Pituitary_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Pituitary_Gene_Expression_Matrix

## Whole Blood
# ./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Whole_Blood/gtex_v7_Whole_Blood_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Whole_Blood_Gene_Expression_Matrix

# Brain ACC BA24
# ./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Brain_Anterior_cingulate_cortex_BA24/gtex_v7_Brain_Anterior_cingulate_cortex_BA24_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Brain_Anterior_cingulate_cortex_BA24_Gene_Expression_Matrix

## Amygdala
# ./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Brain_Amygdala/gtex_v7_Brain_Amygdala_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Brain_Amygdala_Gene_Expression_Matrix

## Frontal Cortex
# ./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights /SAN/neuroscience/PEIC/projects/TWAS/PredictDB/Brain_Frontal_Cortex_BA9/gtex_v7_Brain_Frontal_Cortex_BA9_imputed_europeans_tw_0.5_signif.db  --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/Brain_Frontal_Cortex_BA9_Gene_Expression_Matrix

## DLPFC
# ./PrediXcan.py --predict --dosages /SAN/neuroscience/PEIC/projects/TWAS/$dataset/dosage_mats_2 --dosages_prefix chrom --samples $fam_path --weights ./common_mind/CommonMindDB.v2/DLPFC_newMetax.db --output_prefix /SAN/neuroscience/PEIC/projects/TWAS/$dataset/prediction_results_2/DLPFC_Gene_Expression_Matrix

