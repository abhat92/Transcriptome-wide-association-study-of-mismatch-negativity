#!/bin/bash

## PEIC ##

dataset="PEIC"

# Amygdala
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Brain_Amygdala_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Brain_Amygdala

# ACC
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Brain_Anterior_cingulate_cortex_BA24_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Brain_Anterior_cingulate_cortex_BA24

# Cortex
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Brain_Cortex_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Brain_Cortex

# Frontal Cortex
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Brain_Frontal_Cortex_BA9_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Brain_Frontal_Cortex

# Hippocampus
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Brain_Hippocampus_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Brain_Hippocampus

# Hypothalamus
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Brain_Hypothalamus_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Brain_Hypothalamus

# Pituitary
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Pituitary_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Pituitary

# Whole Blood
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Whole_Blood_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Whole_Blood

#DLPFC
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/DLPFC_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_DLPFC

## MEI ##

dataset="Mei"

# Amygdala
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Brain_Amygdala_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Brain_Amygdala

# ACC
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Brain_Anterior_cingulate_cortex_BA24_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Brain_Anterior_cingulate_cortex_BA24

# Cortex
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Brain_Cortex_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Brain_Cortex

# Frontal Cortex
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Brain_Frontal_Cortex_BA9_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Brain_Frontal_Cortex

# Hippocampus
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Brain_Hippocampus_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Brain_Hippocampus

# Hypothalamus
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Brain_Hypothalamus_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Brain_Hypothalamus

# Pituitary
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Pituitary_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Pituitary

# Whole Blood
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Whole_Blood_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Whole_Blood

#DLPFC
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/DLPFC_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_DLPFC

## ELLIOT ##

dataset="Elliot"

# Amygdala
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Brain_Amygdala_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Brain_Amygdala

# ACC
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Brain_Anterior_cingulate_cortex_BA24_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Brain_Anterior_cingulate_cortex_BA24

# Cortex
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Brain_Cortex_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Brain_Cortex

# Frontal Cortex
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Brain_Frontal_Cortex_BA9_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Brain_Frontal_Cortex

# Hippocampus
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Brain_Hippocampus_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Brain_Hippocampus

# Hypothalamus
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Brain_Hypothalamus_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Brain_Hypothalamus

# Pituitary
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Pituitary_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Pituitary

# Whole Blood
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/Whole_Blood_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_Whole_Blood

#DLPFC
./PrediXcan.py --assoc --pred_exp ./$dataset/prediction_results_2/DLPFC_Gene_Expression_Matrix_predicted_expression.txt --pheno ./$dataset/fzamp_phenos_for_predixcan.txt --linear --output_prefix ./$dataset/association_results_2/FZAmplitude_DLPFC

