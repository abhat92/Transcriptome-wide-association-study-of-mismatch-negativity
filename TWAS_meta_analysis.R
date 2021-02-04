# META-ANALYSIS OF TWAS RESULTS

# clear up working environment

rm(list = ls())

# this script needs to be run with 3.3.2 version of the root R installed in the cluster

# load necessary libraries

library(bacon)
library(data.table)

# load .RData files that contain the association results for each dataset,
# remove the unnecessary files and leave only the list with the results of the bacon command
# , the name of the tests (phenotype-tissue) and the number of effective predictors
# and rename them with the name of the dataset

datasets <- c("PEIC", "Mei", "Elliot")

for(j in 1:length(datasets)){
  
# load .RData file

load(paste("./", datasets[j], "/association_results_2/explore_twas_association_results.RData", sep = ""))

# remove all objects except "corrected_assoc", "assoc_names" and "eff_pred_counts"
  
# objects to keep
  
keep_1 <- c("assoc_list", "corrected_assoc", "assoc_names", "datasets", "j")

keep_2 <- vector("list", length = length(datasets))

for(k in 1:length(keep_2)){
  
  keep_2[[k]] <- ls()[grep(datasets[k], ls())]
  
}

keep_2 <- unlist(keep_2)

rm(list = setdiff(ls(), c("sample_sizes", keep_1, keep_2)))

# rename remaining objects by adding the name of the dataset

assign(paste(datasets[j], "assoc_list", sep = "_"), assoc_list)
assign(paste(datasets[j], "corrected_assoc", sep = "_"), corrected_assoc)
assign(paste(datasets[j], "assoc_names", sep = "_"), assoc_names)
# assign(paste(datasets[j], "eff_pred_counts", sep = "_"), eff_pred_counts)


# remove the objects with the generic name

rm(list = c("assoc_list", "corrected_assoc", "assoc_names"))

}

# do the meta-analysis using 2 methods:

# 2 ways of performing the meta-analysis with the "meta" functions:
#   - Feed non-corrected p-values and set "corrected = TRUE", in "meta"
#   - Feed corrected p-values and set "corrected = FALSE" in "meta"

# Try the 4 ways (2 thresholding methods * 2 "meta" modes) and compare the results

# for each tissue and each dataset, generate a "bonferroni-corrected" p-value variable by multiplying the
# p-values by the number of effective gene expression traits

# do it on both the uncorrected and corrected association results

# first create lists of lists for uncorrected association results, corrected association results
# names of association tests and number of effective predictors

assoc_ls_list <- vector("list", length = length(datasets))
correc_ls_list <- vector("list", length = length(datasets))
assoc_name_list <- vector("list", length = length(datasets))
# eff_pred_list <- vector("list", length = length(datasets))

for(i in 1:length(datasets)){
  
  assoc_ls_list[[i]] <- get(paste(datasets[i], "assoc_list", sep = "_"))
  correc_ls_list[[i]] <- get(paste(datasets[i], "corrected_assoc", sep = "_"))
  assoc_name_list[[i]] <- get(paste(datasets[i], "assoc_names", sep = "_"))
#   eff_pred_list[[i]] <- get(paste(datasets[i], "eff_pred_counts", sep = "_"))
  
}

# perform inverse-variance weighted fixed effects meta-analysis using the "meta" function in bacon
# do it on both uncorrected and corrected data

# get number of associations to do meta-analysis for (the number is the same across all datasets
# so getting it from one of them is ok)

assoc_count <- length(assoc_name_list[[1]])

merged_uncorrect <- vector("list", length = length(datasets))
merged_correct <- vector("list", length = length(datasets))

for(i in 1:assoc_count){
  
  # merge association tables from all datasets for the corresponding test (phenotype-tissue)
  # both with uncorrected and corrected results
  
  uncorrect_seed <- assoc_ls_list[[1]][[i]]
  correct_seed <- correc_ls_list[[1]][[i]]
  colnames(uncorrect_seed)[-1] <- paste(datasets[1], colnames(uncorrect_seed)[2:ncol(uncorrect_seed)], sep = "_")
  colnames(correct_seed)[-1] <- paste(datasets[1], colnames(correct_seed)[2:ncol(correct_seed)], sep = "_")
  
  for(j in 2:length(datasets)){
    
    # extract the association tables for the corresponding dataset and test
    
    uncorrect_tab <- assoc_ls_list[[j]][[i]]
    correct_tab <- correc_ls_list[[j]][[i]]
    colnames(uncorrect_tab)[-1] <- paste(datasets[j], colnames(uncorrect_tab)[2:ncol(uncorrect_tab)], sep = "_")
    colnames(correct_tab)[-1] <- paste(datasets[j], colnames(correct_tab)[2:ncol(correct_tab)], sep = "_")
    
    # merge association tables with the seeds
    
    uncorrect_seed <- merge(uncorrect_seed, uncorrect_tab, by = "gene", all = TRUE)
    correct_seed <- merge(correct_seed, correct_tab, by = "gene", all = TRUE)
    
    }
 
 # remove entries with NA-s in 2 of the datasets (cannot perform a meta-analysis on them!)
 
 uncor_nas <- apply(uncorrect_seed, 1, function(x) length(x[is.na(x)]))
 cor_nas <- apply(correct_seed, 1, function(x) length(x[is.na(x)]))
 
 uncorrect_no_nas <- uncorrect_seed[uncor_nas == 0,]
 correct_no_nas <- correct_seed[cor_nas == 0,]
 
 # store merged association tables in corresponding slot of the purposedly created empty lists
  
  merged_uncorrect[[i]] <- uncorrect_no_nas
  merged_correct[[i]] <- correct_no_nas
  
}

# now use the merged tables to generate bacon objects and perform a meta-analysis

# generate empty lists where to store the results for multiple tissues

uncor_meta_res <- vector("list", length = assoc_count)
cor_meta_res <- vector("list", length = assoc_count)

for(i in 1:assoc_count){
  
  # extract tstats, effect sizes and standard errors
  
  # uncorrected results
  
  uncor_tstat_mat <- as.matrix(merged_uncorrect[[i]][,grep("_t", colnames(merged_uncorrect[[i]]))])
  uncor_es_mat <- as.matrix(merged_uncorrect[[i]][,grep("_beta", colnames(merged_uncorrect[[i]]))])
  uncor_se_mat <- as.matrix(merged_uncorrect[[i]][,grep("_se.beta", colnames(merged_uncorrect[[i]]))])
  
  # corrected results
  
  cor_tstat_mat <- as.matrix(merged_correct[[i]][,grep("_t", colnames(merged_correct[[i]]))])
  cor_es_mat <- as.matrix(merged_correct[[i]][,grep("_beta", colnames(merged_correct[[i]]))])
  cor_se_mat <- as.matrix(merged_correct[[i]][,grep("_se.beta", colnames(merged_correct[[i]]))])
  
  # generate bacon objects
  
  uncor_bacon <- bacon(teststatistics = uncor_tstat_mat, effectsizes = uncor_es_mat, standarderrors = uncor_se_mat)
  cor_bacon <- bacon(teststatistics = cor_tstat_mat, effectsizes = cor_es_mat, standarderrors = cor_se_mat)
  
  # again, it looks like bacon doesn't properly store the effect-sizes and standard errors, so we will
  # have to add them manually
  
  uncor_bacon@effectsizes <- uncor_es_mat
  uncor_bacon@standarderrors <- uncor_se_mat
  
  cor_bacon@effectsizes <- cor_es_mat
  cor_bacon@standarderrors <- cor_se_mat
  
  # now perform meta-analysis using the "meta" function
  # set "corrected=TRUE" for uncorrected data and "corrected=FALSE" for corrected data
  
  uncor_meta <- meta(uncor_bacon, corrected = FALSE)
  cor_meta <- meta(cor_bacon, corrected = FALSE)
  
  # build meta_analysis result data frame from bacon objects
  
  # meta analysis results corrected at meta-analysis stage
  uncor_meta_df <- data.frame(gene = merged_uncorrect[[i]]$gene, beta = es(uncor_meta), tstats = tstat(uncor_meta), p_value = pval(uncor_meta), se.beta = se(uncor_meta))
  uncor_meta_df <- uncor_meta_df[,c(1, grep("meta", colnames(uncor_meta_df)))]
  colnames(uncor_meta_df) <- gsub("[.]meta", "", colnames(uncor_meta_df))
  
  # meta-analysis results from previously corrected association results
  
  cor_meta_df <- data.frame(gene = merged_correct[[i]]$gene, beta = es(cor_meta), tstats = tstat(cor_meta), p_value = pval(cor_meta), se.beta = se(cor_meta))
  cor_meta_df <- cor_meta_df[,c(1, grep("meta", colnames(cor_meta_df)))]
  colnames(cor_meta_df) <- gsub("[.]meta", "", colnames(cor_meta_df))
  
  # compare uncorrected vs corrected results-based meta-analysis results for beta, tstats, p-value and 
  # standard errors
  
  # open .pdf device 
  
  pdf(paste("./meta_analysis_results_2/", assoc_name_list[[1]][[i]], "_uncor_vs_cor_meta_results_comparison.pdf", sep = ""), width = 10, height = 10)
  
  # set graphic parameters to 2 rows and 2 columns
  
  par(mfrow = c(2,2))
  
  # set paramaters to be plotted
  
  plot_pars <- colnames(uncor_meta_df)[-1]
  
  # generated scatterplots for both types of meta-analysis results for each parameter in a loop
  
  for(j in 1:length(plot_pars)){
    
    plot(uncor_meta_df[,plot_pars[j]], cor_meta_df[,plot_pars[j]], xlab = paste(plot_pars[j], "(meta-analysis on uncorrected values)"), ylab = paste(plot_pars[j], "(meta-analysis on corrected values)"))
    abline(a = 0, b = 1, col = "red")
  }
  
  # close pdf device
  
  dev.off()
  
  # add tissue name
  
  uncor_meta_df$tissue <- gsub("FZAmplitude_", "", assoc_name_list[[1]][i])
  cor_meta_df$tissue <- gsub("FZAmplitude_", "", assoc_name_list[[1]][i])
  
  # calculate FDR
  
  uncor_meta_df$tissue_FDR <- p.adjust(p = uncor_meta_df$p_value, method = "fdr")
  cor_meta_df$tissue_FDR <- p.adjust(p = cor_meta_df$p_value, method = "fdr")
  
  # store the results of the meta-analyses in the corresponding lists
  
  uncor_meta_res[[i]] <- uncor_meta_df
  cor_meta_res[[i]] <- cor_meta_df
}

# collapse the results for each method across tissues into an only data.frame

all_uncor_meta <- as.data.frame(do.call("rbind", uncor_meta_res))
all_cor_meta <- as.data.frame(do.call("rbind", cor_meta_res))

# generate a table that contains the meta-analysis results for all genes across all tissues

# read in gene annotations

gene_annot <- as.data.frame(fread("formatted.gencode.v19.gene.annotations.txt", header = T, sep = "\t"))

# generate new gene-id by remove the number added to the Ensembl IDs

all_uncor_meta$gene_2 <- gsub("[.].*", "", all_uncor_meta$gene)
all_cor_meta$gene_2 <- gsub("[.].*", "", all_cor_meta$gene)

# do the same on the gene annotations

gene_annot$gene_id_2 <- gsub("[.].*", "", gene_annot$gene_id)

# merge gene annotations results table using Ensembl Gene ID as key

all_uncor_annot_res <- merge(gene_annot[,c("gene_id_2", "chr", "start", "end", "gene_name", "gene_type", "gene_status")], all_uncor_meta, by.x = "gene_id_2", by.y = "gene_2", all.y = TRUE)
all_cor_annot_res <- merge(gene_annot[,c("gene_id_2", "chr", "start", "end", "gene_name", "gene_type", "gene_status")], all_cor_meta, by.x = "gene_id_2", by.y = "gene_2", all.y = TRUE)

# re-order table by p-value

all_uncor_annot_res <- all_uncor_annot_res[order(all_uncor_annot_res$p_value),]
all_cor_annot_res <- all_cor_annot_res[order(all_cor_annot_res$p_value),]

# add phenotype column

all_uncor_annot_res$MMN_phenotype <- "FZ_amplitude"
all_cor_annot_res$MMN_phenotype <- "FZ_amplitude"

# re-organize columns

all_uncor_annot_res <- all_uncor_annot_res[,c("MMN_phenotype", "tissue", "gene_name", "beta", "tstats", "p_value", "se.beta", "gene","gene_id_2", "chr", "start", "end", "gene_type", "gene_status", "tissue_FDR")]
colnames(all_uncor_annot_res)[colnames(all_uncor_annot_res) == "gene"] <- "gene_id"

all_cor_annot_res <- all_cor_annot_res[,c("MMN_phenotype", "tissue", "gene_name", "beta", "tstats", "p_value", "se.beta", "gene","gene_id_2", "chr", "start", "end", "gene_type", "gene_status", "tissue_FDR")]
colnames(all_cor_annot_res)[colnames(all_cor_annot_res) == "gene"] <- "gene_id"

###############
## FDR BASED ##
###############

# calculate overall FDRs

all_uncor_annot_res$overall_FDR <- p.adjust(all_uncor_annot_res$p_value, method = "fdr")
all_cor_annot_res$overall_FDR <- p.adjust(all_cor_annot_res$p_value, method = "fdr")

# write table out

write.csv(all_uncor_annot_res, "./meta_analysis_results_2/uncorrected_meta_analysis_results_all_tissues_with_fdr.csv", row.names = F, quote = F)
write.csv(all_cor_annot_res, "./meta_analysis_results_2/corrected_meta_analysis_results_all_tissues_with_fdr.csv", row.names = F, quote = F)

# subset using a lenient multiple-testing correction

uncor_fdr_sig_res <- all_uncor_annot_res[all_uncor_annot_res$overall_FDR < 0.15, ]
cor_fdr_sig_res <- all_cor_annot_res[all_cor_annot_res$overall_FDR < 0.15, ]

# write table out

write.csv(uncor_fdr_sig_res, "./meta_analysis_results_2/uncorrected_meta_analysis_significant_results_at_fdr_0.15.csv", row.names = F, quote = F)
write.csv(cor_fdr_sig_res, "./meta_analysis_results_2/corrected_meta_analysis_significant_results_at_fdr_0.15.csv", row.names = F, quote = F)

# select the genes that have an in-tissue FDR < 0.05 (only the uncorrected results)

uncor_fdr_sig_2 <- all_uncor_annot_res[all_uncor_annot_res$tissue_FDR < 0.05,]

# extract the list of tissues that have, at least, 1 significant gene

sig_tissues <- unique(uncor_fdr_sig_2$tissue)

# get all the combinations of numbers from 1 to the number of "significant tissues"

tissue_n <- length(sig_tissues)

combs <- unlist(lapply(c(1:tissue_n), function(x) as.list(as.data.frame(combn(c(1:tissue_n), m = as.numeric(x))))), recursive = FALSE)

# for all combinations, extract the corresponding tissue(s)'s results, calculate the overall FDR and
# get the number of significant hits at FDR < 0.05

# create empty vector in which to store number of significant hits obtained with each combination

hit_vec <- vector("numeric", length = length(combs))

for(i in 1:length(combs)){
  
  # extract tissue(s) to be used
  
  tissue_sub <- sig_tissues[combs[[i]]]
  
  # subset the uncorrected TWAS meta-analysis results for those tissues
  
  twas_res_subset <- all_uncor_annot_res[all_uncor_annot_res$tissue %in% tissue_sub, ]
  
  # re-calculate overall FDR
  
  twas_res_subset$overall_FDR <- p.adjust(p = twas_res_subset$p_value, method = "fdr")
  
  # extract significant hits
  
  sig_twas_res <- twas_res_subset[twas_res_subset$overall_FDR < 0.05, ]
  
  # store number of significant hits in vector
  
  hit_vec[i] <- nrow(sig_twas_res)
  
  # print out number of singificant hits and the tissues involved
  
  print(paste(nrow(sig_twas_res), "signiticant hits with", paste(tissue_sub, collapse = " and ")))
  
}

# get the combination of tissues that provides the highest number of hits

top_tissues <- sig_tissues[combs[[which(hit_vec == max(hit_vec))]]]

# subset uncorrected twas results table for those tissues

twas_res <- all_uncor_annot_res[all_uncor_annot_res$tissue %in% top_tissues,]

# re-calculate overall FDR

twas_res$overall_FDR <- p.adjust(p = twas_res$p_value, method = "fdr")

# write table out

write.csv(twas_res, "./meta_analysis_results_2/twas_results.csv", row.names = F, quote = F)

# select significant hits

sig_twas_res <- twas_res[twas_res$overall_FDR < 0.05,]

# write table out

write.csv(twas_res, "./meta_analysis_results_2/significant_twas_results.csv", row.names = F, quote = F)

# write out the results for each tissue separately, too

for(i in 1:length(top_tissues)) {
  
  # subset TWAS results table for corresponding tissue
  
  tissue_res <- twas_res[twas_res$tissue == top_tissues[i],]
  
  # write table out
  
  write.csv(tissue_res, paste("./meta_analysis_results_2/", top_tissues[i], "_twas_results.csv", sep = ""), row.names = F, quote = F)
  
}

# save object in working environment into .RData file

save.image("./meta_analysis_results_2/TWAS_meta_analysis.RData")





