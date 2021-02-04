# script to prepare the PrediXcan outputs for the TWAS-GSEA script prepared by Oliver Pain (available from
# his GitHub page)

# clean working environment

rm(list = ls())

# load necessary libraries

library(data.table)
library(qusage)
library(GSEABase)
library(WGCNA)

#####################################
## PREDIXCAN META-ANALYSIS RESULTS ##
#####################################

# read in the csv file that contains the PrediXcan meta-analsysis results

# set the path to the folder that contains that file

meta_path <- "/SAN/neuroscience/PEIC/projects/TWAS/meta_analysis_results_2" 

meta_res <- read.csv(paste(meta_path, "/corrected_meta_analysis_results_all_tissues_with_fdr.csv", sep = ""), header = T)

# extract the list of tissues we have TWAS results for (well, the tissues for which we are going
# to present the TWAS results)

tissues <- unique(meta_res$tissue)

# define names of datasets included

datasets <- c("PEIC", "Mei", "Elliot")

# create folders in which to store TWAS-GSEA-ready files from the meta-analsys results

# dataset directory

dir.create("./pathway_enrichment/meta_2", showWarnings = FALSE)

# subdirectories to store formatted twas result files, position files and
# formatted imputed expression matrices

dir.create("./pathway_enrichment/meta_2/twas_results", showWarnings = FALSE)
dir.create("./pathway_enrichment/meta_2/position_files", showWarnings = FALSE)
dir.create("./pathway_enrichment/meta_2/expression_mats", showWarnings = FALSE)

# now, in a loop, prepare the files necessary for the TWAS-GSEA for each tissue

for(i in 1:length(tissues)){
  
  # extract the results for the corresponding tissue
  
  tissue_res <- meta_res[meta_res$tissue == tissues[i],]

  # identify duplicated gene names

  gene_freqs <- table(tissue_res$gene_name)

  dupli_genes <- names(gene_freqs)[gene_freqs > 1]

  # if duplicated genes identified, only keep the entry with the largest absolute t-score
  
  if(length(dupli_genes) > 0){
    
    # for each duplicated gene, identify the entry with the largest absolute t-score and remove the rest
    
    # generate empty list in which to store rows to be removed
    
    rem_rows <- vector("list", length = length(dupli_genes))
    
    for(l in 1:length(dupli_genes)){
      
      print(paste(length(dupli_genes), "duplicated genes; selecting one entry for each"))
      
      # get the rows that correspond to result for that gene
      
      dupli_rows <- which(!is.na(match(tissue_res$gene_name, dupli_genes[l])))
      
      # of those, identify the one that contains the smalles p-value
      
      smallest_p <- intersect(dupli_rows, which(tissue_res$p_value == min(tissue_res$p_value[dupli_rows])))
      
      # get the rows to remove and store them in list
      
      rem_rows[[l]] <- setdiff(dupli_rows, smallest_p)
      
    }
    
    # remove row of duplicated genes with the non-maximum absolute t-value
    
    tissue_res_2 <- tissue_res[-unlist(rem_rows),]
    
  } else {
    
    print(paste(length(dupli_genes), "duplicated genes; no need to remove repeated entries"))
    
    
    tissue_res_2 <- tissue_res
  }


  # add a "FILE" column indicating the name of the file the weights were taken from
  
  # give a unique name in each row as, otherwise, the TWAS_GSEA.R script will identify them
  # as duplicated features and remove them
  
  # set full path to the folder that contains the several files
  
  predictDB_path <- "/SAN/neuroscience/PEIC/projects/TWAS/PredictDB//"
  
  tissue_res_2$FILE <- paste(predictDB_path, make.names(tissue_res_2$gene_name), "_expression_weights", sep = "")
  
  # create 2 versions of the TWAS results file: one using the t-statistic as the TWAS.Z variable
  # and the other using the 
  
  tissue_res_2$logPval <- -log(tissue_res_2$p_value, 10)
  
  # remove the "chr" from the chromosome values
  
  tissue_res_2$chr <- gsub("chr", "", tissue_res_2$chr)
  
  # order table by chromosome and by start position
  
  tissue_res_2 <- tissue_res_2[order(as.numeric(tissue_res_2$chr), as.numeric(tissue_res_2$start)),]
    
  # extract the columns for the reformatted TWAS results in the proper order and re-name them
  
  tissue_res_3 <- tissue_res_2[,c("FILE","gene_name", "chr", "start", "end", "tstats", "p_value")]
  colnames(tissue_res_3) <- c("FILE", "ID", "CHR", "P0", "P1", "TWAS.Z", "TWAS.P")
  
  tissue_res_4 <- tissue_res_2[,c("FILE","gene_name", "chr", "start", "end", "logPval", "p_value")]
  colnames(tissue_res_4) <- c("FILE", "ID", "CHR", "P0", "P1", "TWAS.Z", "TWAS.P")
  
  # write tables out
  
  write.table(tissue_res_3, paste("./pathway_enrichment/meta_2/twas_results/", tissues[i], ".tstats.corrected.results_table.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)
  write.table(tissue_res_4, paste("./pathway_enrichment/meta_2/twas_results/", tissues[i], ".logpval.corrected.results_table.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)
  
  # generate position file, too
  
  pos_file <- tissue_res_3[,c("FILE", "ID", "CHR", "P0", "P1"),]
  
  # remove the path up to "PredictDB//" from the "FILE" column and rename it to "WGT"
  
  pos_file$FILE <- gsub(".*[/]",  "" ,pos_file$FILE)
  colnames(pos_file)[1] <- "WGT"
  
  # write table out
  
  write.table(pos_file, paste("./pathway_enrichment/meta_2/position_files/", tissues[i], ".position_file.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)

  # in a loop, merge the expression matrices from all datasets

  # read in expression matrix for the first dataset to use it as seed
  
  # first get path to corresponding folder
  
  expr_path <- paste("./", datasets[1], "/prediction_results_2/", sep = "")
  
  # get names of files that contain expression matrices
  
  expr_fnames <- list.files(expr_path)[grep("expression.txt", list.files(expr_path))]
  
  # grep tissue name on that list of file names to extract that corresponds to the tissue
  
  tissue_fname <- expr_fnames[grep(tissues[i], expr_fnames)]
  
  # read in expression matrix
  
  expr_seed <- read.table(paste(expr_path, tissue_fname, sep = ""), header = T, sep ="\t")
  
  # read in corresponding phenotypes table
  
  phenos <- read.table(paste("./", datasets[1], "/fzamp_phenos_for_predixcan.txt", sep = ""), header = T, sep = "\t")
  
  # select only the samples for which we have phenotypic data
  
  expr_seed <- expr_seed[expr_seed$IID %in% phenos$IID, ]
  
  # remove the transcript id number from the gene-ids (just because they are missing in the DLPFC ID-s
  # and we wouldn't be able to merge that table with the rest)
  
  colnames(expr_seed) <- gsub("[.].*", "", colnames(expr_seed))
  
  # add IID as rownames
  
  rownames(expr_seed) <- expr_seed$IID
  
  # transpose and return it into a data.frame format
  
  expr_seed_t <- as.data.frame(t(expr_seed))
  
  
  # now, in another loop, merge the expression matrices from the rest of the datasets and for the
  # corresponding tissue to the seed matrix
  
  for(j in 2:length(datasets)){
    
    # first get path to corresponding folder
    
    expr_path <- paste("./", datasets[j], "/prediction_results_2/", sep = "")
    
    # get names of files that contain expression matrices
    
    expr_fnames <- list.files(expr_path)[grep("expression.txt", list.files(expr_path))]
    
    # grep tissue name on that list of file names to extract that corresponds to the tissue
    
    tissue_fname <- expr_fnames[grep(tissues[i], expr_fnames)]
    
    # read in expression matrix
    
    expr_mat <- read.table(paste(expr_path, tissue_fname, sep = ""), header = T, sep ="\t")
    
    # read in corresponding phenotypes table
    
    phenos <- read.table(paste("./", datasets[j], "/fzamp_phenos_for_predixcan.txt", sep = ""), header = T, sep = "\t")
    
    # subset expression matrix so it only contains IID-s for which we have phenotypic data
    
    expr_mat <- expr_mat[expr_mat$IID %in% phenos$IID, ]
    
    # remove the transcript id number from the gene-ids (just because they are missing in the DLPFC ID-s
    # and we wouldn't be able to merge that table with the rest)
    
    colnames(expr_mat) <- gsub("[.].*", "", colnames(expr_mat))
    
    # add IID as rownames
    
    rownames(expr_mat) <- expr_mat$IID
    
    # transpose and return it into a data.frame format
    
    expr_mat_t <- as.data.frame(t(expr_mat))
    
    # merge both expression matrices
    
    merged_expr <- merge(expr_seed_t, expr_mat_t, by = 0, all = TRUE)
    
    # add column "Row.names" as actual rownames, remove it and rename the data.frame into expr_seed_t
    
    rownames(merged_expr) <- merged_expr$Row.names
    merged_expr <- merged_expr[,-1]
    expr_seed_t <- merged_expr
    
  }

  # we need the gene names rather than the gene-ids to match this matrix to the multixcan results (although
  # the latter are much better identifiers than the former)
  
  # merge the necessary annotations in the meta-analysis results with the expression matrix
  
  annots_expr <- merge(tissue_res[,c("gene_id_2", "gene_name")], expr_seed_t, by.x = 1, by.y = 0, all.y = TRUE)

  # add the "_expression_weights" string to the gene_names
  
  annots_expr$gene_name <- as.character(annots_expr$gene_name)
  annots_expr$gene_name[!is.na(annots_expr$gene_name)] <- paste(annots_expr$gene_name[!is.na(annots_expr$gene_name)], "_expression_weights", sep = "")
  
  # apply make.names on the new gene-names
  
  annots_expr$gene_name <- make.names(annots_expr$gene_name)
  
  # add "FID" and "IID" in the corresponding rows of the "gene_name" column
  
  annots_expr$gene_name[match(c("FID", "IID"), annots_expr$gene_id_2)] <- c("FID", "IID")

  # select the rows for the genes for which we have PrediXcan meta-analysis data
  
  # get IDs of genes for which we have data
  
  if(length(dupli_genes) > 0){
  
  twas_ids <- tissue_res$gene_id_2[-unlist(rem_rows)]
  
  } else {
    
    twas_ids <- tissue_res$gene_id_2
  }

  # get row positions for those genes
  
  row_pos <- which(!is.na(match(annots_expr$gene_id_2, twas_ids)))
  
  # and row positions for the FID and IID
  
  id_pos <- which(!is.na(match(annots_expr$gene_name, c("FID", "IID"))))
  
  annots_expr <- annots_expr[c(id_pos, row_pos), ]
  
  # add "gene_name" column as rownames and deleted both the "gene_id_2" and the "gene_name" columns
  
  rownames(annots_expr) <- annots_expr$gene_name
  annots_expr <- annots_expr[,-grep("gene", colnames(annots_expr))]

  # transpose and re-format into a data frame
  
  annots_expr_t <- as.data.frame(t(annots_expr))
  
  # make genes be in the same order as in the meta-analysis results / position file
  
  annots_expr_t <- annots_expr_t[,c(1,2, match(pos_file$WGT, colnames(annots_expr_t)))]

  # now compare gene order with that in the position file as an additional check that there is no
  # problem with gene order
  
  print(paste("Are the genes in the same order in expression matrix and position file?", identical(colnames(annots_expr_t)[-c(1,2)], pos_file$WGT)))
  
  # also make sure genes are in the same order in the meta-analysis results table and in the position file
  
  print(paste("Are the genes in the same order in the meta-analysis results table and the position file?", identical(tissue_res_3$ID, pos_file$ID)))

  # now, we can write the expression matrix out
  
  write.table(annots_expr_t, paste("./pathway_enrichment/meta_2/expression_mats/", tissues[i], ".expression_mat.txt", sep = ""), col.names = T, row.names = F, sep = " ", quote = F)

}

# print notice indicating that the preparation of files for Ollie's TWAS-GSEA pipeline has finished

print("Preparation of files for Ollie's TWAS-GSEA script done")

# save environment objects into an RData file

# save.image(paste("./pathway_enrichment/multixcan_meta_2/prepare_files_for_twas_gsea_2_subset_", subset, ".RData", sep = ""))

