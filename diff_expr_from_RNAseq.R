# simple script to perform differential expression analysis from RNA-seq gene-count data

# clean working environment

rm(list = ls())

# load necessary libraries

library(car)
library(limma)
library(fastDummies)
library(NMF)
library(ggplot2)
library(ggpubr)

## READ IN AND PREPARE DATA ##

# set path to folder where we have the data

data_path <- "./BrainSpan/rnaseq_gene_expression"

# read in RPKM matrix

rpkms <- read.csv(paste(data_path, "/expression_matrix.csv", sep = ""), header = F, row.names = 1)

# read in gene-annotations

gannots <- read.csv(paste(data_path, "/rows_metadata.csv", sep = ""), header = T, row.names = 1)

# read in phenotypes

phenotypes <- read.csv(paste(data_path, "/columns_metadata.csv", sep = ""), header = T, row.names = 1)

# add Ensembl-IDs as rownames of the gene_count matrix

rownames(rpkms) <- gannots$ensembl_gene_id

# create unique ID for each sample combining the donor-name and the brain structure acronym

phenotypes$donor_structure <- paste(phenotypes$donor_name, phenotypes$structure_acronym)

# add that as column names to gene-count matrix

colnames(rpkms) <- phenotypes$donor_structure

# do the same with gene-annotations and phenotypes

rownames(gannots) <- gannots$ensembl_gene_id
rownames(phenotypes) <- phenotypes$donor_structure

#check if genes and samples are in the same order

print(paste("Are the samples in the same order?", identical(colnames(rpkms), rownames(phenotypes)), sep = " "))
print(paste("Are the genes in the same order?", identical(rownames(rpkms), rownames(gannots)), sep = " "))

# generate plot with density distribution of RPKM values for 60 randomly selected genes

gene_sample <- sample(x = rownames(rpkms), size = 60)

# raw RPKM distributions

pdf("./BrainSpan/RPKM_distributions.pdf", width = 12, height = 8)

# set graphic parameters to fit several plots in each page

par(mfrow=c(2,3))

# in a loop, plot density distributions for each gene

for(i in 1:length(gene_sample)){

plot(density(as.numeric(as.matrix(rpkms[gene_sample[i],]))), xlab = "RPKM", main = gene_sample[i])

# close loop

}

# close connection to pdf file

dev.off()

# log-RPKM distributions

pdf("./BrainSpan/log_RPKM_distributions.pdf", width = 12, height = 8)

# set graphic parameters to fit several plots in each page

par(mfrow=c(2,3))

# in a loop, plot density distributions for each gene

for(i in 1:length(gene_sample)){
  
  plot(density(log(as.numeric(as.matrix(rpkms[gene_sample[i],]))+2^-18)), xlab = "log(RPKM + 2^-18)", main = gene_sample[i])
  
  # close loop
  
}

# close connection to pdf file

dev.off()

# also generate overall distributions

# raw RPKM

pdf("./BrainSpan/overall_RPKM_distribution.pdf")

plot(density(as.numeric(as.matrix(rpkms))), xlab = "RPKM")

dev.off()

# log-RPKM

pdf("./BrainSpan/overall_log_RPKM_distribution.pdf")

plot(density(log(as.numeric(as.matrix(rpkms)))+2^-18), xlab = "log(RPKM + 2^-18)")
abline(v = -4, col = "red")
abline(v = -4.5, col = "red")
abline(v = -5, col = "red")

dev.off()

# generate log2-RPKM matrix

log_rpkm <- log(rpkms, 2)

# 2^-5 looks like a good threshold to filter lowly expressed genes out.

log_thr <- -5

# remove those genes with a expression below the threshold in 90% of more of the samples

frac_thr <- 0.9

# calculate lowly expressed sample fraction per gene

low_counts <- apply(log_rpkm, 1, function(x) length(x[x < log_thr]))
low_fracs <- low_counts/ncol(log_rpkm)

# remove lowly expressed genes

filt_expr <- rpkms[low_fracs < frac_thr,]

# from the "age" variable extract the number and the unit and add them as variables in
# the phenotypic table

no_unit <- do.call("rbind", strsplit(as.character(phenotypes$age), split = " "))
colnames(no_unit) <- c("age_no", "age_unit")
phenotypes <- as.data.frame(cbind(phenotypes, no_unit))
phenotypes$age_no <- as.numeric(as.character(phenotypes$age_no))

## I have manually checked the ages and corresponding genders and I have manually
## defined 9 groups of 4 - 6 individuals.

## create that variable (age_group) and add it to the phenotypic table

# group definitions

group_defs <- list(c(8, "pcw", 12,"pcw"), c(13, "pcw", 16, "pcw"), c(17, "pcw", 24, "pcw"), c(25, "pcw", 37, "pcw"), c(4, "mos", 1, "yrs"), c(2, "yrs", 4, "yrs"), c(8, "yrs", 15, "yrs"), c(18, "yrs", 23, "yrs"), c(30, "yrs", 40, "yrs"))
group_names <- c("8 - 12 pcw", "13 - 16 pcw", "17 - 24 pcw", "25 - 37 pcw", "4 mos - 1 yr", "2 - 4 yrs", "8 - 15 yrs", "18 - 23 yrs", "30 - 40 yrs")

# add empty variable to phenotypic table

phenotypes$age_group <- NA

# in a loop, assign each group name to the samples that meet the corresponding age criteria

for(i in 1:length(group_defs)){
  
  # extract age-group criteria
  
  group_def <- group_defs[[i]]
  
  # get positions of sample whose age is equal or greater of the lower age limit
  
  low_age_rows <- intersect(which(phenotypes$age_no >= as.numeric(group_def[1])), which(phenotypes$age_unit == group_def[2]))
  
  # do the same with the upper age limit
  
  high_age_rows <- intersect(which(phenotypes$age_no <= as.numeric(group_def[3])), which(phenotypes$age_unit == group_def[4]))
  
  # get rows that fullfill both conditions
  # for all ranges except for the "4 mos - 1 yr" one, that's the intersection
  
  if(group_def[2] == group_def[4]){
  
  age_range_rows <- intersect(low_age_rows, high_age_rows)
  
  } else {
    
    age_range_rows <- union(low_age_rows, high_age_rows)
  }
  # in the positions give by the age_range_rows vector, add the corresponding age group name to the empty age_group variable
  
  phenotypes$age_group[age_range_rows] <- group_names[i]
  
}
# 
# define reference groups in the variables we are including in the regressions
# use the groups with the largest sample-sizes as reference

# turn the age-group variable into a factor

phenotypes$age_group <- as.factor(phenotypes$age_group)

# convert the age-group variable into a set of dummy variables and add it to the table

phenotypes <- dummy_cols(phenotypes, select_columns = "age_group")

# read in manually-generated table of brain strucuture names and the corresponding broader regions

brain_regions <- read.csv(paste(data_path, "/brain_structure_groupings.csv", sep = ""))

# add brain regions to the phenotypes table

phenotypes <- merge(phenotypes, brain_regions, by = "structure_name", all.x = TRUE)

# indentify and remove the brain structures with too small a sample size

small_strs <- names(table(phenotypes$structure_acronym))[table(phenotypes$structure_acronym) < 2]

phenotypes <- phenotypes[!phenotypes$structure_acronym %in% small_strs,]
filt_expr <- filt_expr[,phenotypes$donor_structure]

# make sure samples are still in the same order

print(paste("Are the samples in the same order in phenotypes and expression tables?", identical(colnames(filt_expr), phenotypes$donor_structure)))

# remove empty levels from brain structures variables

phenotypes$structure_acronym <- droplevels(phenotypes$structure_acronym)

# identify the group(s) with largest sample size for the covariates that will go into the model

gender_ref <- names(table(phenotypes$gender))[table(phenotypes$gender) == max(table(phenotypes$gender))]
str_ref <-  names(table(phenotypes$structure_acronym))[table(phenotypes$structure_acronym) == max(table(phenotypes$structure_acronym))]

# set those groups as reference (in case there are more than 1 group selected in a varible
# just choose the first in the vector)

phenotypes$gender <- relevel(phenotypes$gender, ref = gender_ref[1])
phenotypes$structure_acronym <- relevel(phenotypes$structure_acronym, ref = str_ref[1])

# build model

reg_formula <- "~ age_group + gender + structure_acronym"
reg_model <- model.matrix(as.formula(reg_formula), data = phenotypes)

#check potential colinearity problems in the model

diag_phenos <- phenotypes
diag_phenos$gene_expr <- as.numeric(filt_expr[1,])

diag_lm <- lm(as.formula(paste("gene_expr", reg_formula)), data = diag_phenos)

diag_vif <- vif(diag_lm)

print(paste("The DF-adjusted GVIF value for the variable of interest is", round(diag_vif[1,3], digits = 2), sep = " "))

if(diag_vif[1,3] > 2){
  print("Variance of coefficients inflated beyond acceptable levels")
} else if(diag_vif[1,3] < 2){
  print("No coefficient variance inflation problems with the model")
}

# calculate the correlation between samples from the same patient by using it as blocking
# variable

# run duplicateCorrelation

id_cors <- duplicateCorrelation(object = filt_expr, design = reg_model, block = phenotypes$donor_name)

# extract consensus intra-sample correlation value

cons_cor <- id_cors$consensus.correlation

# in a loop, run linear regressions for each group for both brain cortex and frontal brain cortex

# make age-group names R compatible and get list of unique age-groups

colnames(phenotypes) <- make.names(colnames(phenotypes))
phenotypes$age_group <- 

age_groups <- as.character(unique(phenotypes$age_group))

# open loop

for(i in 1:length(age_groups){
  
  # extract the name of the age group for which we will generate signatures
  
  age_group <- age_groups[i]

  # build model
  
  sig_formula <- paste("~ ", age_group, " + gender + structure_acronym")
  sig_model <- model.matrix(as.formula(sig_formula), data = phenotypes)
  
  
  # 
})

####################################
## IDENTIFY BRAIN-REGION CLUSTERS ##
####################################

# fit a model with age, gender and structure-acronym and remove the effect of the first two

# identify the group(s) with largest sample size for the variables that will go into the model

age_ref <- names(table(phenotypes$age))[table(phenotypes$age) == max(table(phenotypes$age))]
gender_ref <- names(table(phenotypes$gender))[table(phenotypes$gender) == max(table(phenotypes$gender))]
str_ref <-  names(table(phenotypes$structure_acronym))[table(phenotypes$structure_acronym) == max(table(phenotypes$structure_acronym))]

# set those groups as reference (in case there are more than 1 group selected in a varible
# just choose the first in the vector)

phenotypes$age <- relevel(phenotypes$age, ref = age_ref[1])
phenotypes$gender <- relevel(phenotypes$gender, ref = gender_ref[1])
phenotypes$str_ref <- relevel(phenotypes$structure_acronym, ref = str_ref[1])

# build model

clust_formula <- "~ age + gender + structure_acronym"
clust_model <- model.matrix(as.formula(clust_formula), data = phenotypes)

#check potential colinearity problems in the model

diag_phenos <- phenotypes
diag_phenos$gene_expr <- as.numeric(filt_expr[1,])

diag_lm <- lm(as.formula(paste("gene_expr", clust_formula)), data = diag_phenos)

diag_vif <- vif(diag_lm)

print(paste("The DF-adjusted GVIF value for the variable of interest is", round(diag_vif[1,3], digits = 2), sep = " "))

if(diag_vif[1,3] > 2){
  print("Variance of coefficients inflated beyond acceptable levels")
} else if(diag_vif[1,3] < 2){
  print("No coefficient variance inflation problems with the model")
}


# calculate the correlation between samples from the same patient by using it as blocking
# variable

# run duplicateCorrelation

id_cors <- duplicateCorrelation(object = filt_expr, design = clust_model, block = phenotypes$donor_name)

# extract consensus intra-sample correlation value

cons_cor <- id_cors$consensus.correlation

# run linear regressions on all genes using lmFit; add the donor-name as blocking variable
# and the consensus intra-donor correlation

gene_fit <- lmFit(object = filt_expr, design = clust_model, block = phenotypes$donor_name, correlation = cons_cor)

# extract coefficients from the regression fit

gene_coefs <- gene_fit$coef

# get columns of coefficients we want to remove

rem_eff_cols <- unlist(lapply(c("age", "gender"), function(x) grep(x, colnames(gene_coefs))))

# calculate aggregate effect per each gene and sample

aggr_eff <- gene_coefs[,rem_eff_cols]%*%t(clust_model[,rem_eff_cols])

# remove the aggregated effect of age and gender from the raw log-rpkms

corr_expr <- filt_expr - aggr_eff

# calculate correlation matrix between samples

cormat <- cor(corr_expr, method = "pearson")

# plot correlation matrix as a heatmap using "aheatmap"

pdf("./BrainSpan/age_and_gender_corrected_rnaseq_cormatrix.pdf", onefile = FALSE)

aheatmap(cormat, annCol = phenotypes$structure_acronym)

dev.off()

# impossible to see anything in that heatmap

# try with a PCA

pca_res <- prcomp(t(corr_expr), scale = TRUE)

# extract sample coordinates and add the to the phenotypes

samp_coords <- pca_res$x
phenos_coords <- merge(phenotypes, samp_coords, by.x = "donor_structure", by.y = 0)
phenos_coords <- as.data.frame(phenos_coords)

# calculate fraction of variance explained by each Principal Component

# extract standard deviations, turn them into variances and divide them by the total variants

pca_sdevs <- pca_res$sdev 
pca_vars <- pca_sdevs^2
pca_var_fracs <- pca_vars / sum(pca_vars)

# plot all combinations of the first 4 dimensions

# create ggplot objects for each plot and then arrange them in an only page

# combinations of PC-s

pc_combs <- combn(x = c(1:4), m = 2)

# create empty list in which to store ggplot objects

g_list <- vector("list", length = ncol(pc_combs))

# open for loop

for(i in 1:ncol(pc_combs)){
  
  # get PC names
  
  pc1 <- paste("PC", pc_combs[1,i], sep = "")
  pc2 <- paste("PC", pc_combs[2,i], sep = "")
  
  # create vector of columns we need to plot and subset the data
  
  plot_cols <- c("structure_acronym", pc1, pc2)
  plot_phenos <- phenos_coords[,plot_cols]
  
  # rename them
  
  colnames(plot_phenos) <- c("structure_acronym", "PC1", "PC2")
  
  # generate ggplot object
  
  g <- ggplot(data = plot_phenos, aes(x = PC1, y = PC2))
  g <- g + geom_point(aes(color = structure_acronym, shape = structure_acronym))
  g <- g + stat_ellipse(aes(fill = structure_acronym), alpha = 0.25, geom = "polygon")  
  g <- g + xlab(paste(pc1, " (", round(pca_var_fracs[pc_combs[1,i]]*100, digits = 2), "%)", sep = ""))
  g <- g + ylab(paste(pc2, " (", round(pca_var_fracs[pc_combs[2,i]]*100, digits = 2), "%)", sep = ""))
  g <- g + theme_bw()
  

  # store ggplot object in the list
  
  g_list[[i]] <- g
  
}

# open pdf file

pdf("./BrainSpan/age_and_gender_corrected_rnaseq_pca_plots.pdf", height = 12, width = 22)

# arrange all ggplot-objects using ggarrange

ggarrange(plotlist = g_list, ncol = 3, nrow = 2)

# close connection to pdf file

dev.off()

# according to the PCA-s, the only relatively sensible grouping of the brain structures would be
# "cerebellar cortex" vs the rest

# but, of course, that is not useful at all; so we will do it manually following anatomical
# location

# read in manually-generated table of brain strucuture names and the corresponding broader regions

brain_regions <- read.csv(paste(data_path, "/brain_structure_groupings.csv", sep = ""))

# add brain regions to the phenotypes table

phenotypes <- merge(phenotypes, brain_regions, by = "structure_name", all.x = TRUE)

#############################
## DIFFERENTIAL EXPRESSION ##
#############################

# generate and add two new variables to the phenotypic table:
# Brain cortex vs. the rest
# Frontal cortex vs. the rest

# get the rows that contain that for each

cortex_rows <- setdiff(grep("cortex", phenotypes$structure_name), grep("cerebellar", phenotypes$structure_name))

frontal_terms <- c("dorsolateral", "orbital", "primary motor ", "ventrolateral")
frontal_rows <- unlist(lapply(frontal_terms, function(x) grep(x, phenotypes$structure_name)))

# generate empty variables and fill them with "yes/no" values using the row-numbers in the vectors above

phenotypes$brain_cortex <- NA
phenotypes$brain_cortex[cortex_rows] <- "yes"
phenotypes$brain_cortex[-cortex_rows] <- "no"

phenotypes$frontal_cortex <- NA
phenotypes$frontal_cortex[frontal_rows] <- "yes"
phenotypes$frontal_cortex[-frontal_rows] <- "no"

# turn these variables into factors

phenotypes$brain_cortex <- as.factor(phenotypes$brain_cortex)
phenotypes$frontal_cortex <- as.factor(phenotypes$frontal_cortex)

# set reference group for each variable that goes into the models

# for brain and frontal cortex variables: the "no" group

phenotypes$brain_cortex <- relevel(phenotypes$brain_cortex, ref = "no")
phenotypes$frontal_cortex <- relevel(phenotypes$frontal_cortex, ref = "no")

## BRAIN CORTEX SIGNATURES ##

# subset the dataset for brain cortex

cortex_phenos <- phenotypes[phenotypes$brain_cortex == "yes",]
# 
# remove the "embryonic" regions as they are confounded with the youngest age group

# cortex_phenos <- cortex_phenos[cortex_phenos$brain_region != "embryonic",]
cortex_expr <- filt_expr[,colnames(filt_expr) %in% cortex_phenos$donor_structure]

# check if samples are in the same order

print(paste("Are samples in the same order in expression and phenotypic tables for Brain Cortex?", identical(colnames(cortex_expr), cortex_phenos$donor_structure)))

# remove the levels that are empty from the age group and brain region variables

cortex_phenos$age_group <- droplevels(cortex_phenos$age_group)
cortex_phenos$brain_region <- droplevels(cortex_phenos$brain_region)

# convert the age-group variable into a set of dummy variables and add it to the table

cortex_phenos <- dummy_cols(cortex_phenos, select_columns = "age_group")

# for gender and brain region,
# identify the group(s) with largest sample size for the variables that will go into the model

gender_ref <- names(table(cortex_phenos$gender))[table(cortex_phenos$gender) == max(table(cortex_phenos$gender))]
region_ref <-  names(table(cortex_phenos$brain_region))[table(cortex_phenos$brain_region) == max(table(cortex_phenos$brain_region))]
  
# set those groups as reference (in case there are more than 1 group selected in a varible
# just choose the first in the vector)

cortex_phenos$gender <- relevel(cortex_phenos$gender, ref = gender_ref[1])
cortex_phenos$brain_region <- relevel(cortex_phenos$brain_region, ref = region_ref[1])

# make sure column names are R-compatible

colnames(cortex_phenos) <- make.names(colnames(cortex_phenos))

# build regression models

# create formula

age_terms <- colnames(cortex_phenos)[grep("age_group_", colnames(cortex_phenos))]

reg_formula <- paste("~", paste(age_terms, collapse = " + "), "+ gender")

cortex_model <- model.matrix(as.formula(reg_formula), data = cortex_phenos)

#check potential colinearity problems in the model

cortex_diag <- cortex_phenos
cortex_diag$gene_expr <- as.numeric(cortex_expr[1,])

cortex_diag_lm <- lm(as.formula(paste("gene_expr", reg_formula)), data = cortex_diag)

cortex_vif <- vif(cortex_diag_lm)

print(paste("The DF-adjusted GVIF value for the variable of interest is", round(cortex_vif[1,3], digits = 2), sep = " "))

if(cortex_vif[1,3] > 2){
  print("Variance of coefficients inflated beyond acceptable levels")
} else if(cortex_vif[1,3] < 2){
  print("No coefficient variance inflation problems with the model")
}

# calculate the correlation between samples from the same patient by using it as blocking
# variable

# run duplicateCorrelation

cortex_cors <- duplicateCorrelation(object = cortex_expr, design = cortex_model, block = cortex_phenos$donor_name)

# extract consensus intra-sample correlation value

cons_cortex_cor <- cortex_cors$consensus.correlation

# run linear regressions on all genes using lmFit; add the donor-name as blocking variable
# and the consensus intra-donor correlation

gene_fit <- lmFit(object = cortex_expr, design = cortex_model, block = cortex_phenos$donor_name, correlation = cons_cortex_cor)

# save working environment into .RData file

save.image("./BrainSpan/diff_expr_from_RNAseq.RData")
