# combine GSEA results performed on meta-analysis results of gene-level analysis

# clean working environment

rm(list = ls())

# define vector of gsea parameters

gsea_pars <- "ollie.logpval"

# get the number of tests performed in each modality

gsea_tests <- 210

# generate empty list in which to store the GSEA results tables for each approach

gsea_ls <- vector("list", length = length(gsea_pars))

# and an empty vector in which to store the number of hits obtained with each parameter combination

hit_vec <- vector("numeric", length = length(gsea_pars))

# in a loop, read the results from the competitive test for each GSEA parameter combination, combine those and
# calculate overall, across-tissue FDR using the total number of tests

for(i in 1:length(gsea_pars)){
  
  # get list of filenames that contain results
  
  # path to folder
  
  gsea_path <- paste("./pathway_enrichment/meta_2/", gsea_pars[i], ".results", sep = "")
  
  gsea_fnames_1 <- list.files(gsea_path)[grep("competitive.txt", list.files(gsea_path))]
  
  gsea_fnames <- gsea_fnames_1
  
  if(length(gsea_fnames) > 0){
  
  # create empty list in which to store each results table
  
  gsea_list <- vector("list", length = length(gsea_fnames))
  
  for(j in 1:length(gsea_fnames)){
  
   # read results table
    
    gsea_res <- read.table(paste(gsea_path, "/", gsea_fnames[j], sep = ""), header = T, sep = " ")
    
    # extract tissue name from file name
    
    tissue <- gsub("[.].*", "", gsea_fnames[j])
    
    # add tissue-name to table
    
    gsea_res$tissue <- tissue
    
    # store table in corresponding slot of list
    
    gsea_list[[j]] <- gsea_res
  }
  # read in all files and "rbind" them
  
  gsea_tab <- as.data.frame(do.call("rbind", gsea_list))
  
  
  # calculate overall FDR using the corresponding number of standard linear regression tests across tissues
  
  gsea_tab$overall_FDR <- p.adjust(gsea_tab$P, method = "fdr", n = gsea_tests[i])

  # put table in corresponding slot of list
  
  gsea_ls[[i]] <- gsea_tab
  
  # get number of significant hits at FDR < 10%
  
  sig_hits <- nrow(gsea_tab[gsea_tab$overall_FDR < 0.1,])
  
  # store that number in the corresponding slot in vector
  
  hit_vec[i] <- sig_hits
  
  }
  else {
    
    print("No gene-set made it into the competitive test")
  }
}

# get the position in the gsea-parameter list of the parameter combination that yields the maximum number of
# hits

top_pars <- which(hit_vec == max(hit_vec))

# extract corresponding table of gene-set enrichment results

top_gsea <- gsea_ls[[top_pars]]

# subset it so only the significant hits remain

sig_res <- top_gsea[top_gsea$overall_FDR < 0.1, ]

# write table out

write.csv(sig_res, "./pathway_enrichment/meta_2/ollie.logpval.gsea.results.csv", )
