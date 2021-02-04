# visualization of TWAS results

# clear up working environment

rm(list = ls())

# load necessary libraries

library(ggplot2)
library(qusage)
library(org.Hs.eg.db)
library(annotate)
library(NMF)
library(tidyr)
library(RColorBrewer)
library(tidyr)
library(reshape2)

#####################################
## GENE-ASSOCIATION MANHATTAN PLOT ##
#####################################

# read in PrediXcan meta-analysis results

meta_res <- read.csv("./meta_analysis_results_2/twas_results.csv")

# replace "Brain Frontal Cortex" by "Frontal Brain Cortex"

meta_res$tissue <- as.character(meta_res$tissue)
meta_res$tissue[meta_res$tissue == "Brain_Frontal_Cortex"] <- "Frontal_Brain_Cortex"

# get list of tissues for which we have PrediXcan results

tissues <- as.character(unique(meta_res$tissue))

# get a data frame with columns: "GENE", "CHR", "START", "COEF", "FDR", "TISSUE"

man_df <- meta_res[,c("gene_name", "chr", "start", "beta", "overall_FDR", "tissue")]
colnames(man_df) <- c("GENE", "CHR", "START", "COEF", "FDR", "TISSUE")

# remove chr-s from chromosome values

man_df$CHR <- as.numeric(gsub("chr", "", man_df$CHR))

# generate logFDR value

man_df$logFDR <- sign(man_df$COEF)*-log(man_df$FDR, 10)

# generate chromosome-wise coloring pallete

mypalette <- c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B") # chr color palette

# set threshold for suggestive evidence of association

sugg <- 0.1

# import function gene.manhattan which is just the function gg.manhattan to generate ggplot-based Manhattan plots
# but tweaked a little bit to work with gene-level results

# the function was taken from the GitHub account of "pcgoddard"

source("gene_manhattan_function.R")

# in a loop, subset the visualization data.frame for each tissue, get the corresponding highlighting genes
# and generate the manhattan plot

for(i in 1:length(tissues)){
 
  # open .pdf file
  
  pdf(paste("./results_visualization/", tissues[i], "_predixcan_manhattan_plot.pdf", sep = ""), width = 11)
  
  # subset visualization data for the tissue
    
  tissue_man <- man_df[man_df$TISSUE == tissues[i], ]
  
  # get gene-names to highlit
  
  hlight_genes <- as.character(tissue_man[tissue_man$FDR < 0.1, "GENE"])
  
  # generate Manhattan plot
  
  gene.manhattan(df = tissue_man, threshold = 0.05, sugg = sugg, hlight = hlight_genes, col = mypalette, ylims = c(-max(man_df$logFDR)-0.1, max(man_df$logFDR)+0.1), title = gsub("_", " ", tissues[i]))

  # close connection to PDF
  
  dev.off()
}

#################################
## GENE-SET ENRICHMENT BARPLOT ##
#################################

# read in table of significant gene-set enrichment results (for the same subset, of course!!)

gs_enrich <- read.csv("./pathway_enrichment/meta_2/ollie.logpval.gsea.results.csv")

# prepare data.frame for ggplot visualization

gs_df <- gs_enrich[,c("GeneSet", "overall_FDR")]
colnames(gs_df) <- c("GeneSet", "FDR")

# remove the prefixes from the enriched pathways

gs_df$GeneSet <- gsub(".*GO_", "", gs_df$GeneSet)
gs_df$GeneSet <- gsub("OLLIE_", "", gs_df$GeneSet)

# order the dataset from smallest to largest FDR

gs_df <- gs_df[order(gs_df$FDR),]

# remove underscores from GeneSet names and put them in lower case

gs_df$GeneSet <- gsub("_", " ", gs_df$GeneSet)
gs_df$GeneSet <- tolower(gs_df$GeneSet)

# create horizontal bar-plot

pdf("./results_visualization/gene_set_enrichment_results_bar_plot_2.pdf", height = 4, width = 10)

g <- ggplot(data = gs_df, aes(x = GeneSet, y = -log(FDR, 10)))
g <- g + geom_bar(stat = "identity", aes(fill = -log(FDR, 10)), color = "black")
g <- g + geom_label(aes(x = GeneSet, y = 0.025, label = GeneSet), fill = "white", color = "black", hjust = 0, vjust = 0.5, size = 6)
g <- g + geom_hline(yintercept = -log(0.05, 10), color = "red")
g <- g + geom_hline(yintercept = -log(0.1, 10), color = "red", linetype = "dashed")
g <- g + scale_x_discrete(limits=rev(as.character(gs_df$GeneSet)))
g <- g + xlab("")
g <- g + ylab("-log10(FDR)")
g <- g + labs(fill = "-log10(FDR)")
g <- g + scale_fill_gradient(low = "lightskyblue", high = "navy")
g <- g + theme_bw()
g <- g + theme(legend.position = "none", axis.title.x=element_text(size = 16, face = "bold"), axis.text.x=element_text(size = 14), axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(), legend.title = element_text(face = "bold", size = 10))
g <- g + coord_flip()

print(g)

dev.off()

##########################################
## TISSUE-SPECIFIC ASSOCIATIONS HEATMAP ##
##########################################

# do it separately for each gene-set and also include for the top associated genes

# first remove the genes for which we don't have Gene Symbol

meta_res <- meta_res[!is.na(meta_res$gene_name),]

# set number of top genes to visualize

top_g <- 10

# create list in which to store the data.frames containing the tissue-specific TWAS results
# for each gene-set

gs_assoc_list <- vector("list", length = nrow(gs_enrich)+1)

# the first one will contain the results for the top associated genes

# collapse table by selecting the minimum p-value for each gene

gene_pvals <- aggregate(p_value ~ gene_name, data = meta_res, FUN = "min")

# re-order table of genes with their smalles p-values by increasing p-value

gene_pvals <- gene_pvals[order(gene_pvals$p_value),]

# extract list of genes with the smalles p-values

top_genes <- as.character(gene_pvals$gene_name[1:10])

# extract statistics of top genes

top_genes_df <- meta_res[meta_res$gene_name %in% top_genes,c("tissue", "gene_name", "beta", "tstats", "p_value")]

# add column indicating what gene-set it refers to

top_genes_df$GeneSet <- "TOP TWAS GENES"

# add that table to the first position of the list and do the same for the enriched gene-sets

gs_assoc_list[[1]] <- top_genes_df

# read in the GMT files that contain the pathways

# Ollie's gene-sets

ollie_gmt <- read.gmt("./pathway_enrichment/Pocklington2015_134sets_LoFi.gmt")

# for each significantly enriched pathway, extract list of Entrez IDs, find corresponding Gene Symbols
# and intersect it with the list of the Gene Symbols in the multixcan meta-analysis results table,
# where the tissue-specific results are also included


for(i in 1:nrow(gs_enrich)){
  
  # get GeneSet name
  
  gs_name <- as.character(gs_enrich[i, "GeneSet"])
  
  # extract list of Entrez IDs
  
  entrez_vec <- ollie_gmt[gs_name][[1]]
  
  # print out number of entrez ids identified
  
  print(paste(length(entrez_vec), "Entrez IDs identified for", gs_name))
  
  # get list of Gene Symbols for those entrez id-s
  
  gene_names <- getSYMBOL(entrez_vec, data='org.Hs.eg')
  
  # print out number of Gene Symbols obtained
  
  print(paste(length(gene_names), "Gene Symbols for", gs_name))
  
  # subset table that contains gene names and their lowest p-values across tissues
  
  gs_res <- gene_pvals[gene_pvals$gene_name %in% gene_names, ]
  
  # re-order that data.frame by increasing p-value 
  
  gs_res <- gs_res[order(gs_res$p_value),]
  
  # get list of top genes from the gene-set
  
  gs_top_genes <- as.character(gs_res$gene_name[1:min(top_g, nrow(gs_res))])
  
  # extract tissue-specific association statistics for those genes
  
  gs_top_res <- meta_res[meta_res$gene_name %in% gs_top_genes, c("tissue", "gene_name", "beta", "tstats", "p_value")]
  
  # re-format name of GeneSet
  
  gs_name_1 <- gsub(".*GO_", "", gsub("OLLIE_", "", gs_name))
  gs_name_1 <- gsub("_" , " ", gs_name_1)
  gs_name_1 <- toupper(gs_name_1)
  
  # add column with GeneSet name
  
  gs_top_res$GeneSet <- gs_name_1
  
  # put table in corresponding slot of list
  
  gs_assoc_list[[i+1]] <- gs_top_res 
  
}

# collapse the gene-set specific table into an only data frame for visualization

gs_assoc_df <- as.data.frame(do.call("rbind", gs_assoc_list))

# rename columns

colnames(gs_assoc_df)[colnames(gs_assoc_df) == "gene_name"] <- "GeneSymbol"
colnames(gs_assoc_df)[colnames(gs_assoc_df) == "tissue"] <- "Tissue"

# change "Brain Frontal Cortex" to "Frontal Brain Cortex"

gs_assoc_df$Tissue <- as.character(gs_assoc_df$Tissue)
gs_assoc_df$Tissue[gs_assoc_df$Tissue == "Brain_Frontal_Cortex"] <- "Frontal_Brain_Cortex"

# remove underscores from tissue names

gs_assoc_df$Tissue <- gsub("_", " ", gs_assoc_df$Tissue)


# now, using ggplot, plot tissue vs GeneSymbol, filling the tiles by color mapping the tstatistic

pdf("./results_visualization/tissue_specific_twas_results_for_top_genes_heatmap.pdf", width = 16, height = 12)

g <- ggplot(data = gs_assoc_df, aes(x = Tissue, y = GeneSymbol))
g <- g + geom_tile(aes(fill = tstats))
g <- g + scale_fill_gradient2(low = "blue", mid = "lightyellow", high = "red", midpoint = 0)
g <- g + theme_bw()
g <- g + facet_wrap(~GeneSet, scales = "free", nrow = 3, ncol = 3)
print(g)

dev.off()

# try aheatmap

gs_assoc_df[is.na(gs_assoc_df)] <- 0

# re-order table

# add a column that assigns an order priority to the GeneSets
# the top TWAS genes have priority and, then, the enriched Geneset by TWAS-GSEA FDR

gs_priorities <- data.frame(GeneSet = unique(gs_assoc_df$GeneSet), priority = c(1:length(unique(gs_assoc_df$GeneSet))))

gs_assoc_df <- merge(gs_assoc_df, gs_priorities, by = "GeneSet")

# get the most distinctive colors

n <- 9
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector_2 <- sample(col_vector, n)

# sort the table by GeneSet priority and the signed_mat_t

gs_assoc_df <- gs_assoc_df[order(gs_assoc_df$priority),]

# put table in wide format

gs_assoc_wide <- dcast(GeneSymbol + GeneSet ~ Tissue, data = gs_assoc_df, value.var = "beta")

# fill up the NA-s with 0-s

gs_assoc_wide[is.na(gs_assoc_wide)] <- 0

# generate heatmap with aheatmap

pdf("./results_visualization/tissue_specific_twas_results_for_top_genes_heatmap_2.pdf", width = 12, height = 10, onefile = FALSE)

aheatmap(x = as.matrix(gs_assoc_wide[,c(3,4)]), border_color = "black", Colv = NA, cexRow = 5, annColors = list(GeneSet = col_vector_2), annLegend = TRUE, width = 12, height = 12, annRow = gs_assoc_wide[,"GeneSet"], labRow = as.character(gs_assoc_wide$GeneSymbol), labCol = gsub("_", " ", colnames(gs_assoc_wide[,c(3,4)])))

dev.off()

# try plotting each GeneSet separately

genesets <- unique(gs_assoc_df$GeneSet)

# open pdf device

pdf("./results_visualization/tissue_specific_twas_results_for_top_genes_heatmap_3.pdf", width = 22, height = 16, onefile = FALSE)

# set graphic paramaters to generate a 3x3 plot figure

par(mfrow = c(3,3))

# open for loop

for(i in 1:length(genesets)){
  
  # extract data for the corresponding geneset
  
  gs_data <- gs_assoc_wide[gs_assoc_wide$GeneSet == genesets[i],]
  
  # plot if out in a heatmap
  
  # remove x-axis labels if the plot doesn't go in the last row
  
  # if(i < 7){
  
  # aheatmap(x = as.matrix(gs_data[,c(3:5)]), fontsize  = 10, main = genesets[i], border_color = "grey75", Rowv = TRUE,  Colv = NA, labRow = as.character(gs_data$GeneSymbol), labCol = c("", "", ""))

  # } else if (i > 6) {
    
  aheatmap(x = as.matrix(gs_data[,c(3,4)]), fontsize = 10, main = genesets[i], border_color = "grey75", Rowv = TRUE,  Colv = NA, labRow = as.character(gs_data$GeneSymbol), labCol = gsub("_", " ", colnames(gs_data[,c(3,4)])))
    
  # }    
}

# close pdf device

dev.off()

## Using the coefficients generated by MultiXcan

# read in the table with the gene coefficients generate by Multixcan for the corresponding
# tissue combinations for all datasets

datasets <- c("PEIC", "Mei", "Elliot")

# in a loop, read in the coefficients files and store the results in a list

coeff_list <- vector("list", length = length(datasets))

for(i in 1:length(datasets)){
  
  # read in coefficient table
  
  coef_tab <- read.table(paste("./", datasets[i], "/multixcan_results_2/multixcan_coeffs.88", sep = ""), header = T, sep = "\t")
  
  # create an across-dataset key combining tissue (called variable) and gene
  
  coef_tab$gene_tissue <- paste(coef_tab$gene, coef_tab$variable, sep = "-")
  
  # add dataset name to the "param" column
  
  colnames(coef_tab)[grep("param", colnames(coef_tab))] <- paste(datasets[i], "coeff", sep = "_")
  
  # add table to the list
  
  coeff_list[[i]] <- coef_tab
    
}

# merge the tables of coefficients for each gene and tissue across datasets

# extract first as seed

seed_tab <- coeff_list[[1]]

# we only need the coefficient and the key column

seed_tab <- seed_tab[,c(1,4)]

# add the data for the other two datasets by merging iteratively

for(i in 2:length(datasets)){
  
  # extract table and select necessary columns
  
  coef_tab <- coeff_list[[i]]
  coef_tab <- coef_tab[,c(1,4)]
  
  # merge with seed
  
  merged_tab <- merge(seed_tab, coef_tab, by = "gene_tissue", all = TRUE)
  seed_tab <- merged_tab
}

# generate independent tissue and gene variables from the key column and remove it

seed_tab$gene <- gsub("-.*", "", seed_tab$gene_tissue)
seed_tab$tissue <- gsub(".*-", "", seed_tab$gene_tissue)
seed_tab <- seed_tab[,-grep("gene_tissue", colnames(seed_tab))]

# re-order the columns

coeff_df <- seed_tab[,c("gene", "tissue", colnames(seed_tab)[grep("coeff", colnames(seed_tab))])]

# for those gene-tissue combinations with coefficients from more than 1 dataset, combine them by
# calculating the average weighted by the square-root of the sample size

# get sample sizes

sample_n <- c(254, 71, 403)

# for each gene/tissue calculate the combined "coefficient"

# create empty vector in which to store the results

meta_coeffs <- vector("numeric", length = length(sample_n))

# open for loop

for(i in 1:nrow(coeff_df)){
  
  # extract vector of coefficients and remove NA-s
  
  coeff_vec <- coeff_df[i, grep("coeff", colnames(coeff_df))]
  coeff_vec <- coeff_vec[!is.na(coeff_vec)]
  
  # if coefficients from more than 1 dataset, combine them in a weighted average
  # otherwise, just add the only coefficient as the "combined coefficient"
  
  if(length(coeff_vec) == 1){
    
    meta_coeffs[i] <- coeff_vec
    
  } else if (length(coeff_vec) > 1){
    
    # calculate square-root-of-sample-size-weighted average
    
    wgt_avg <- sum(coeff_vec*sqrt(sample_n))/sum(sqrt(sample_n))
    
    # assign weighted average to meta-coeffs vector
    
    meta_coeffs[i] <- wgt_avg
  }
  
  # print out iteration
  
  print(paste("Coefficients for", i, "of", nrow(coeff_df), "gene-tissue pairs combined"))
}

# add that meta-coeff vector to the coeff_df table

coeff_df$Meta_coeff <- meta_coeffs

# select only the gene, the tissue and the correspoding aggregated coefficient

coeff_df_2 <- coeff_df[,c("gene","tissue", "Meta_coeff")]

# reshape the coeff_df_2 data.frame into a wide format so the coefficients for each tissue are in separate
# columns

coeffs_wide <- spread(coeff_df_2, tissue, Meta_coeff)
colnames(coeffs_wide)[-1] <- paste(colnames(coeffs_wide)[-1], "_beta", sep = "")

# add gene symbols to this by merging with the appropriate columns of the meta-analysis results table

# get columns that contain the PrediXcan-assoc results: we want to remove them

pred_cols <- unlist(lapply(c("_beta", "_se.beta", "_tstats"), function(x) grep(x, colnames(meta_res))))

# remove them

meta_res_2 <- meta_res[,-pred_cols]

# merge this meta-analysis results table with the meta-coeffs table in wide format

meta_res_coeffs <- merge(meta_res_2, coeffs_wide, by = "gene", all.x = TRUE)

# write table out in .csv format

write.csv(meta_res_coeffs, "./multixcan_meta_analysis_2/multixcan_meta_analysis_results_subset_88.csv", row.names = F, quote = F)

# add GeneSet membership and geneset priority data for plotting

gs_assoc_2 <- merge(gs_assoc_df[,c("GeneSet","GeneSymbol", "priority")], meta_res_coeffs[,c("gene", "gene_name", "multixcan_pval", colnames(meta_res_coeffs)[grep("_beta", colnames(meta_res_coeffs))])], by.x = "GeneSymbol", by.y = "gene_name")

# now generate signed -log10Pvalues for each gene and tissue

sign_logPvals <- -log(gs_assoc_2$multixcan_pval,10)*sign(gs_assoc_2[,grep("_beta", colnames(gs_assoc_2))])

# give appropriate column names

colnames(sign_logPvals) <- gsub("_beta", "_logPval", colnames(sign_logPvals))

# add the signed -log10Pvalues to the main plotting data.frame

gs_assoc_3 <- as.data.frame(cbind(gs_assoc_2, sign_logPvals))

# replace the NA-s with 0-s

gs_assoc_3[is.na(gs_assoc_3)] <- 0

# # get the maximum absolute logPvalue and its sign
# 
# max_logpval <- apply(gs_assoc_3[,grep("_logPval", colnames(gs_assoc_3))], 1, function(x) x[max(abs(x))])
# sign_max_logpval <- sign(max_logpval)

# remove the "_logPval" addition from the tissue colnames

colnames(gs_assoc_3) <- gsub("_logPval", "", colnames(gs_assoc_3))

# # add the maximum logpval and its sign to the data.frame
# 
# gs_assoc_2$max_beta <- max_beta
# gs_assoc_2$sign_max_beta <- sign_max_beta

# plot a heatmap for the top genes of each GeneSet (including the TOP TWAS genes) separately

genesets <- unique(gs_assoc_df$GeneSet)

## logPVALUES

# calculate signed log-pvalue

gs_assoc_df$sign_logPvals <- sign(gs_assoc_df$beta)*-log(gs_assoc_df$p_value, 10)

# put the data.frame in wide format for heatmap plotting

gs_assoc_wide_2 <- dcast(GeneSymbol + GeneSet + priority ~ Tissue, data = gs_assoc_df, value.var = "sign_logPvals")

# replace NA-s by 0-s

gs_assoc_wide_2[is.na(gs_assoc_wide_2)] <- 0

# generate breaks to map colors to numbers in the heatmap

# get maximum absolute logPvalue

max_abs_logpval <- max(abs(gs_assoc_df$sign_logPvals[!is.na(gs_assoc_df$sign_logPvals)]))

# create a symmetric numeric vector centered on 0 that goes from minus to plus largest absolute logPvalue

col_breaks <- seq(from = -max_abs_logpval, to = max_abs_logpval, by = 2*max_abs_logpval/100)

# open pdf device

pdf("./results_visualization/tissue_specific_twas_results_for_top_genes_heatmap_logpvals.pdf", width = 15, height = 12, onefile = FALSE)

# set graphic paramaters to generate a 3x3 plot figure

par(mfrow = c(2,2))

# open for loop

for(i in 1:length(genesets)){
  
  # extract data for the corresponding geneset
  
  gs_data <- gs_assoc_wide_2[gs_assoc_wide_2$priority == i,]
  
  # extract name of Geneset
  
  geneset <- unique(gs_data$GeneSet)
  
  # generate breaks
  
  aheatmap(x = as.matrix(gs_data[,c(4,5)]), breaks = col_breaks, scale = "none", fontsize = 10, main = geneset, border_color = "white", Rowv = TRUE,  Colv = NA, labRow = as.character(gs_data$GeneSymbol), labCol = gsub("_", " ", colnames(gs_data[,c(4,5)])))
  
  # }    
}

# close pdf device

dev.off()

## TSTATS

# put the data.frame in wide format for heatmap plotting

gs_assoc_wide_3 <- dcast(GeneSymbol + GeneSet + priority ~ Tissue, data = gs_assoc_df, value.var = "tstats")

# replace NA-s by 0-s

gs_assoc_wide_3[is.na(gs_assoc_wide_3)] <- 0

# generate breaks to map colors to numbers in the heatmap

# get maximum absolute logPvalue

max_abs_tstat <- max(abs(gs_assoc_df$tstats[!is.na(gs_assoc_df$tstats)]))

# create a symmetric numeric vector centered on 0 that goes from minus to plus largest absolute logPvalue

col_breaks <- seq(from = -max_abs_tstat, to = max_abs_tstat, by = 2*max_abs_tstat/100)

# turn the "Frontal_Brain_Cortex" column name into "Frontal Cortex"

colnames(gs_assoc_wide_3)[colnames(gs_assoc_wide_3) == "Frontal Brain Cortex"] <- "Frontal Cortex"


# open pdf device

# pdf("./results_visualization/tissue_specific_twas_results_for_top_genes_heatmap_tstats.pdf", width = 15, height = 12, onefile = FALSE)

# set graphic paramaters to generate a 3x3 plot figure

# par(mfrow = c(2,2))

# open for loop

for(i in 1:length(genesets)){
  
  # open pdf file
  
  pdf(paste("./results_visualization/", paste(genesets[i], collapse = "_"), "2_heatmap.pdf", sep = ""), onefile = FALSE, width = 4, height = 5)
  
  # extract data for the corresponding geneset
  
  gs_data <- gs_assoc_wide_3[gs_assoc_wide_3$priority == i,]
  
  # extract name of Geneset
  
  geneset <- unique(gs_data$GeneSet)
  
  # generate breaks
  
  aheatmap(x = as.matrix(gs_data[,c(4,5)]), legend = FALSE, breaks = col_breaks, scale = "none", fontsize = 12, main = geneset, border_color = "white", Rowv = TRUE,  Colv = NA, labRow = as.character(gs_data$GeneSymbol), labCol = gsub("_", " ", colnames(gs_data[,c(4,5)])))
  
  # close pdf device
  
  dev.off()
  
  # }    
}

# close pdf device

# dev.off()

##################################################
## ACROSS TISSUES/DATASETS RESULTS CORRELATIONS ##
##################################################

# put the table that contains the coefficients for each gene in each tissue and dataset in wide format

# define function to be able to use "spread" by feeding a vector of column names as the "value" argument

myspread <- function(df, key, value) {
  # quote key
  keyq <- rlang::enquo(key)
  # break value vector into quotes
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% gather(variable, value, !!!s) %>%
    unite(temp, !!keyq, variable) %>%
    spread(temp, value)
}

# applied the newly defined "myspread" function 

coeffs_wide_2 <- myspread(coeff_df, tissue, c("PEIC_coeff", "Mei_coeff", "Elliot_coeff"))

# generate correlation matrix for the 9 tissue-dataset-specific Multixcan-generate coefficient vectors
# the coefficients vary in a very wide range (specially for those from DLPFC) so we should better use
# "Spearman's Rho" to avoid outliers distorting the results too much.

# get the numbers of the columns that contain that data

coef_cols <- setdiff(grep("_coeff", colnames(coeffs_wide_2)), grep("Meta", colnames(coeffs_wide_2)))

cormat <- cor(x = as.matrix(coeffs_wide_2[coef_cols]), method = "spearman", use = "pairwise.complete.obs")

# change the colnames and rownames so they are suitable for plotting

colnames(cormat) <- gsub("_coeff", "", colnames(cormat))
colnames(cormat) <- gsub("BA9_", "", colnames(cormat))
colnames(cormat) <- gsub("Mei", "Harvard", colnames(cormat))
colnames(cormat) <- gsub("Elliot", "Maryland", colnames(cormat))

rownames(cormat) <- colnames(cormat)

# generate annotations vectors for dataset and tissue

tissue_vec <- gsub("Cortex_.*", "Cortex", colnames(cormat))
tissue_vec <- gsub("DLPFC_.*", "DLPFC", tissue_vec)

dataset_vec <- gsub(".*_", "", colnames(cormat))

# apply last touch on the rownames and colnames of the correlation matrix

colnames(cormat) <- gsub("Cortex_", "Cortex-", colnames(cormat))
colnames(cormat) <- gsub("DLPFC_", "DLPFC-", colnames(cormat))
colnames(cormat) <- gsub("_", " ", colnames(cormat))
rownames(cormat) <- colnames(cormat)


pdf("./results_visualization/multixcan_gene_coefficients_correlations_heatmap.pdf", onefile = FALSE, width = 12)

aheatmap(cormat, main = "TWAS result correlations: coefficients", cexCol = 0, annCol = list(tissue = tissue_vec, dataset = dataset_vec))

dev.off()

# save objects in working directory into .RData file

save.image("./results_visualization/visualize_TWAS_results.RData")
