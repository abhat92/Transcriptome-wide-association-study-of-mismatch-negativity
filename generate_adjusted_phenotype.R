# simple script to generated adjusted phenotypes by fitting them to a linear mixed model and extracting the 
# effects of undesired covariates

# load necessary libraries

library(data.table)
library(lme4qtl)
library(lme4)

# set variables

# select phenotype to be adjusted

pheno <- "MMNCZArea"

# path where genetic data (kinship matrix and fam file) can be found

geno_path <- "/cluster/project9/PE/user/aritz/re_QCed_PEIC_genotypes/typed_genotypes_only"

# read in PEIC dataset

peic <- read.table("PEIC_dataset_2018_07_16.txt", header = T, sep = "\t")

# read in kinship matrix

kin_mat <- as.data.frame(fread(paste(geno_path, "/prune.grm.raw", sep = ""), header = F, sep = " "))

# read in fam file

fam <- read.table(paste(geno_path, "/pe19.clean.fam", sep = ""), header = F, sep = " ")

# add numbers as column names in kinship matrix (rownames are automatically assigned numbers as names)

colnames(kin_mat) <- c(1:ncol(kin_mat))

# add additional id column to fam file with numeric id-s

fam$id_2 <- c(1:nrow(fam))

# add numeric ID (id_2) to phenotype table by merging it with fam file (original ID and numeric ID only)

peic_2 <- merge(peic, fam[,c("V2", "id_2")], by.x = "id", by.y = "V2", all.x = TRUE)

# add variable that codes gender with "male" and "female" values

cat_gender <- vector("character", length = nrow(peic_2))

cat_gender[peic_2$gender == "1"] <- "male"
cat_gender[peic_2$gender == "2"] <- "female"

peic_2$cat_gender <- as.factor(cat_gender)

# define the fixed effect covariates

covars <- c("group", "age", "cat_gender", "MMN_lab")

# subset phenotypic table so only the samples that have non-NA values for phenotype, the covariates and 
# the numeric ID are included

# count number of NA-s per sample

na_count <- apply(peic_2[,c("id_2", pheno, covars)], 1, function(x) length(x[is.na(x)]))

# subset table so it contains only the samples with an NA count of 0

fit_phenos <- peic_2[na_count == 0,]

# removed unused levels from factors of interest

fit_phenos$group <- droplevels(fit_phenos$group)
fit_phenos$cat_gender <- droplevels(fit_phenos$cat_gender)
fit_phenos$MMN_lab <- droplevels(fit_phenos$MMN_lab)

# subset the kinship matrix so it contains only the samples that will go into the regression

fit_kinmat <- kin_mat[which(!is.na(match(rownames(kin_mat), fit_phenos$id_2))), which(!is.na(match(colnames(kin_mat), fit_phenos$id_2)))]

# reorder both the kinship matrix and the phenotypic table (the subsetted ones) according to the numeric ID-s

fit_phenos <- fit_phenos[order(fit_phenos$id_2),]

fit_kinmat <- fit_kinmat[order(as.numeric(rownames(fit_kinmat))), order(as.numeric(colnames(fit_kinmat)))]

# check if the row names and column names of the subsetted kinship matrix fit with the numeric ID-s of the
# subsetted phenotypic table

print(paste("Are the rownames in the kinship matrix and the numeric id-s in the phenotypic table identical?", identical(rownames(fit_kinmat), as.character(fit_phenos$id_2))))
print(paste("Are the colnames in the kinship matrix and the numeric id-s in the phenotypic table identical?", identical(colnames(fit_kinmat), as.character(fit_phenos$id_2))))

# make sure that gender is coded as 

# fit linear mixed model using reference coding

# create string to be fed as formula

ref_formula <- paste(pheno, paste(paste(covars, collapse = " + "), "(1|id_2)", sep = " + "), sep = " ~ ")

# perform mixed model regression feeding the formula, the phenotype table and the kinship matrix

lmm_fit_ref <- relmatLmer(formula = as.formula(ref_formula), data = fit_phenos, relmat = list(id_2 = as.matrix(fit_kinmat)))

# extract fixed effects

fixed_effects <- fixef(lmm_fit_ref)

# create a model matrix with dummy variables so we can multiply it by the fixed effects

# generate formula

lm_formula <- paste(pheno, paste(covars, collapse = " + "), sep = " ~ ")

# use it to generate model matrix

mod_mat <- model.matrix(as.formula(lm_formula), data = fit_phenos)

# get vector of accumulated fixed effects per sample

fix_ef_vec <- fixed_effects%*%t(mod_mat)

# remove accumulated fixed effects from raw outcome

adjust_pheno <- fit_phenos[,pheno] - fix_ef_vec

# give it a try at fitting a simple linear model without mixed effects

lm_fit <- lm(formula = as.formula(lm_formula), data = fit_phenos)

# extract residuals

lm_res <- residuals(lm_fit)

# plot lmm-based adjusted phenotype vs lm residuals and the raw phenotype

pdf(paste(pheno, "raw_vs_adjusted_phenos.pdf", sep = "_"))
plot(x = lm_res, y = adjust_pheno, main = "linear vs. mixed model adjusted phenotype", xlab = paste("adjusted", pheno, "(linear model)"), ylab = paste("adjusted", pheno, "(mixed model)"))
plot(x = fit_phenos[,pheno], y = lm_res, main = "raw vs linear model adjusted phenotype", ylab = paste("adjusted", pheno, "(linear model)"), xlab = paste("raw", pheno))
plot(x = fit_phenos[,pheno], y = adjust_pheno, main = "raw vs mixed model adjusted phenotype",  xlab = paste("raw", pheno), ylab = paste("adjusted", pheno, "(mixed model)"))
dev.off()

# the names of the entries in the linear-model-adjusted phenotype corresponde to the rownames in fit_pheno
# so we can just add the adjusted variable to the fit phenos table

fit_phenos$adjusted_pheno <- lm_res

# now merge the fam table (only the FID and the IID) and the fit_phenos table (only the id and the adjusted_pheno)

adj_pheno_tab <- merge(x = fam[,c("V1", "V2")], y = fit_phenos[,c("id", "adjusted_pheno")], by.x = 2, by.y = 1)

# re-order columns

adj_pheno_tab <- adj_pheno_tab[,c("V1", "V2", "adjusted_pheno")]

# PrediXcan expects the phenotype file to have headers, so add proper ones

colnames(adj_pheno_tab) <- c("FID", "IID", paste("adjusted", pheno, sep = "_"))

# write table out

write.table(adj_pheno_tab, paste("adjusted", pheno, "table.txt", sep = "_"), col.names = T, row.names = F, sep = "\t", quote = F)

# save objects in environment into .RData file

save.image(paste("adjusted", pheno, "generation.RData", sep = "_"))
