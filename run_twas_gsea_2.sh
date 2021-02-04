#$ -l h_vmem=12.5G
#$ -l tmem=12.5G
#$ -l h_rt=02:00:00
#$ -t 460
#$ -j y
#$ -cwd
#$ -S /bin/bash

# set necessary variables

dataset="meta_2"

# for this, we need to run my local, conda-based installation of R; set the path

r_path="/cluster/project9/PE/user/aritz/software/anaconda2/lib/R/bin"

# in a loop, perform enrichment analysis for each tissue

# run 10 versions combining:
# 	- either TWAS t-statistic or TWAS log10p_value as outcome
#       - either both Ollie's and the CNS-GO genesets, ollie's gene-sets only, the GO-CNS genesets only,
#         the actin-related ones or the immune-related one
# don't run self-contained analyses

# first, generate 10 folders in which to store the 10 types of enrichment analysis

mkdir -p ./pathway_enrichment/$dataset/ollie.cns.tstats.results
mkdir -p ./pathway_enrichment/$dataset/ollie.cns.logpval.results
mkdir -p ./pathway_enrichment/$dataset/ollie.tstats.results
mkdir -p ./pathway_enrichment/$dataset/ollie.logpval.results
mkdir -p ./pathway_enrichment/$dataset/actin.tstats.results
mkdir -p ./pathway_enrichment/$dataset/actin.logpval.results
mkdir -p ./pathway_enrichment/$dataset/cns.tstats.results
mkdir -p ./pathway_enrichment/$dataset/cns.logpval.results
mkdir -p ./pathway_enrichment/$dataset/immune.tstats.results
mkdir -p ./pathway_enrichment/$dataset/immune.logpval.results

# select the tissue for which to perform the enrichment analysis

tissue="Brain_Frontal_Cortex"

# run TWAS-GSEA script

# -log10Pvalue as outcome / Ollie's pathways + CNS pathways from GO combined

$r_path/Rscript TWAS-GSEA.V1.3.R \
        --twas_results ./pathway_enrichment/$dataset/twas_results/$tissue.logpval.corrected.results_table.txt \
        --pos ./pathway_enrichment/$dataset/position_files/$tissue.position_file.txt \
        --gmt_file ./pathway_enrichment/ollie_and_cns_go_gene_sets.gmt \
        --self_contained F \
        --competitive T \
        --expression_ref ./pathway_enrichment/$dataset/expression_mats/$tissue.expression_mat.txt \
        --probit_P_as_Z F \
        --output ./pathway_enrichment/$dataset/ollie.cns.logpval.results/$tissue.ollie.cns.logpval.gsea.corrected.results

# t-statistic as outcome / Ollie's pathways + CNS pathways from GO combined

$r_path/Rscript TWAS-GSEA.V1.3.R \
        --twas_results ./pathway_enrichment/$dataset/twas_results/$tissue.tstats.corrected.results_table.txt \
        --pos ./pathway_enrichment/$dataset/position_files/$tissue.position_file.txt \
        --gmt_file ./pathway_enrichment/ollie_and_cns_go_gene_sets.gmt \
        --self_contained F \
        --competitive T \
        --expression_ref ./pathway_enrichment/$dataset/expression_mats/$tissue.expression_mat.txt \
        --probit_P_as_Z F \
        --output ./pathway_enrichment/$dataset/ollie.cns.tstats.results/$tissue.ollie.cns.tstats.gsea.corrected.results

# -log10Pvalue as outcome / Ollie's pathways only

$r_path/Rscript TWAS-GSEA.V1.3.R \
        --twas_results ./pathway_enrichment/$dataset/twas_results/$tissue.logpval.corrected.results_table.txt \
        --pos ./pathway_enrichment/$dataset/position_files/$tissue.position_file.txt \
        --gmt_file ./pathway_enrichment/Pocklington2015_134sets_LoFi.gmt \
        --self_contained F \
        --competitive T \
        --expression_ref ./pathway_enrichment/$dataset/expression_mats/$tissue.expression_mat.txt \
        --probit_P_as_Z F \
        --output ./pathway_enrichment/$dataset/ollie.logpval.results/$tissue.ollie.logpval.gsea.corrected.results

# t-statistic as outcome / Ollie's pathways only

$r_path/Rscript TWAS-GSEA.V1.3.R \
        --twas_results ./pathway_enrichment/$dataset/twas_results/$tissue.tstats.corrected.results_table.txt \
        --pos ./pathway_enrichment/$dataset/position_files/$tissue.position_file.txt \
        --gmt_file ./pathway_enrichment/Pocklington2015_134sets_LoFi.gmt \
        --self_contained F \
        --competitive T \
        --expression_ref ./pathway_enrichment/$dataset/expression_mats/$tissue.expression_mat.txt \
        --probit_P_as_Z F \
        --output ./pathway_enrichment/$dataset/ollie.tstats.results/$tissue.ollie.tstats.gsea.corrected.results

# -log10Pvalue as outcome / Actin pathways only

$r_path/Rscript TWAS-GSEA.V1.3.R \
        --twas_results ./pathway_enrichment/$dataset/twas_results/$tissue.logpval.corrected.results_table.txt \
        --pos ./pathway_enrichment/$dataset/position_files/$tissue.position_file.txt \
        --gmt_file ./pathway_enrichment/actin_gene_sets.gmt \
        --self_contained F \
        --competitive T \
        --expression_ref ./pathway_enrichment/$dataset/expression_mats/$tissue.expression_mat.txt \
        --probit_P_as_Z F \
        --output ./pathway_enrichment/$dataset/actin.logpval.results/$tissue.actin.logpval.gsea.corrected.results

# t-statistic as outcome / Ollie's pathways only

$r_path/Rscript TWAS-GSEA.V1.3.R \
        --twas_results ./pathway_enrichment/$dataset/twas_results/$tissue.tstats.corrected.results_table.txt \
        --pos ./pathway_enrichment/$dataset/position_files/$tissue.position_file.txt \
        --gmt_file ./pathway_enrichment/actin_gene_sets.gmt \
        --self_contained F \
        --competitive T \
        --expression_ref ./pathway_enrichment/$dataset/expression_mats/$tissue.expression_mat.txt \
        --probit_P_as_Z F \
        --output ./pathway_enrichment/$dataset/actin.tstats.results/$tissue.actin.tstats.gsea.corrected.results

# -log10Pvalue as outcome / CNS pathways only

$r_path/Rscript TWAS-GSEA.V1.3.R \
        --twas_results ./pathway_enrichment/$dataset/twas_results/$tissue.logpval.corrected.results_table.txt \
        --pos ./pathway_enrichment/$dataset/position_files/$tissue.position_file.txt \
        --gmt_file ./pathway_enrichment/cns_gene_sets.gmt \
        --self_contained F \
        --competitive T \
        --expression_ref ./pathway_enrichment/$dataset/expression_mats/$tissue.expression_mat.txt \
        --probit_P_as_Z F \
        --output ./pathway_enrichment/$dataset/cns.logpval.results/$tissue.cns.logpval.gsea.corrected.results

# t-statistic as outcome / CNS pathways only

$r_path/Rscript TWAS-GSEA.V1.3.R \
        --twas_results ./pathway_enrichment/$dataset/twas_results/$tissue.tstats.corrected.results_table.txt \
        --pos ./pathway_enrichment/$dataset/position_files/$tissue.position_file.txt \
        --gmt_file ./pathway_enrichment/cns_gene_sets.gmt \
        --self_contained F \
        --competitive T \
        --expression_ref ./pathway_enrichment/$dataset/expression_mats/$tissue.expression_mat.txt \
        --probit_P_as_Z F \
        --output ./pathway_enrichment/$dataset/cns.tstats.results/$tissue.cns.tstats.gsea.corrected.results

# -log10Pvalue as outcome / IMMUNE pathways only

$r_path/Rscript TWAS-GSEA.V1.3.R \
        --twas_results ./pathway_enrichment/$dataset/twas_results/$tissue.logpval.corrected.results_table.txt \
        --pos ./pathway_enrichment/$dataset/position_files/$tissue.position_file.txt \
        --gmt_file ./pathway_enrichment/immune_gene_sets.gmt \
        --self_contained F \
        --competitive T \
        --expression_ref ./pathway_enrichment/$dataset/expression_mats/$tissue.expression_mat.txt \
        --probit_P_as_Z F \
        --output ./pathway_enrichment/$dataset/immune.logpval.results/$tissue.immune.logpval.gsea.corrected.results

# t-statistic as outcome / IMMUNE pathways only

$r_path/Rscript TWAS-GSEA.V1.3.R \
        --twas_results ./pathway_enrichment/$dataset/twas_results/$tissue.tstats.corrected.results_table.txt \
        --pos ./pathway_enrichment/$dataset/position_files/$tissue.position_file.txt \
        --gmt_file ./pathway_enrichment/immune_gene_sets.gmt \
        --self_contained F \
        --competitive T \
        --expression_ref ./pathway_enrichment/$dataset/expression_mats/$tissue.expression_mat.txt \
        --probit_P_as_Z F \
        --output ./pathway_enrichment/$dataset/immune.tstats.results/$tissue.immune.tstats.gsea.corrected.results

