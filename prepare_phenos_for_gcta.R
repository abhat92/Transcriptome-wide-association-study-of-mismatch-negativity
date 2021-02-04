# simple script to prepare the phenotypic tables for the bivariate-GREML analysis with GCTA

# path to the folder that contains the merged phenotypes for all datasets

pheno_path <- "/SAN/neuroscience/PEIC/projects/EEG_genetics_consortium/phenotypic_data/merged_phenotypes"

# read merged phenotypic table

pheno_tab <- read.table(paste(pheno_path, "/MMN_Merged_2018_12_06.txt", sep = ""), header = T, sep = "\t")

# we are not going to use the relatives for this analysis, so get rid of them

# pheno_tab <- pheno_tab[pheno_tab$group != "relative", ]

# recode group variable into "0" and "1"

pheno_tab$group_num <- ifelse(pheno_tab$group == "control", 0, 1)

# in a loop, subset the phenotypes for each dataset and merge them with the corresponding FAM file to
# add the Family-IDs (FIDs)

# in different stages of the analysis, we have used different names for the datasets

datasets_1 <- c("PEIC", "Mei", "Elliot")

datasets_2 <- c("PEIC", "Harvard", "Maryland")

# open for loop

for(i in 1:length(datasets_1)){

	# subset phenotypes table

	ds_phenos <- pheno_tab[pheno_tab$Sample == datasets_2[i],]

	# read in phenotypic table prepared for PrediXcan as it contains FIDs

	fid_phenos <- read.table(paste("./", datasets_1[i], "/fzamp_phenos_for_predixcan.txt", sep = ""), header = T, sep = "\t")

	# add the FIDs to the rest of the phenotypes
	# remember that, for the Harvard dataset from Mei, the ids of the main table are the FID-s of the FAM file, not the IIDs	

	if(datasets_1[i] != "Mei") {

	all_phenos <- merge(ds_phenos, fid_phenos[,c("IID", "FID")], by.x = "id", by.y = "IID")

	# re-organize columns and re-name them

        all_phenos <- all_phenos[,c("FID", "id", "MMNFZAmplitude", "group_num", "age", "gender", "MMN_lab")]
        colnames(all_phenos)[colnames(all_phenos) == "id"] <- "IID"

	} else {

	all_phenos <- merge(ds_phenos, fid_phenos[,c("FID", "IID")], by.x = "id", by.y = "FID")
	
	# re-organize columns and re-name them

        all_phenos <- all_phenos[,c("id", "IID", "MMNFZAmplitude", "group_num", "age", "gender", "MMN_lab")]
        colnames(all_phenos)[colnames(all_phenos) == "id"] <- "FID"

	}
	
	# also remember that all Elliot's samples from Maryland were scanned in the same machine and, thus,
	# the "MMN_lab" variable shouldn't be included in the table of categorial covariates

	# get number of MMN-labs present
	# if only one, just remove that variable
	# otherwise, run a regression on the phenotype with age, gender and MMN_lab, remove the effect of the latter
	# from the phenotype and use the adjusted phenotype for the analysis

	# remove missing levels from MMN_lab

	all_phenos$MMN_lab <- droplevels(all_phenos$MMN_lab)
	
	# get unique value count

	mmn_lab_count <- length(unique(all_phenos$MMN_lab))

	if(mmn_lab_count > 1){

	adj_reg <- lm(formula = MMNFZAmplitude ~ age + gender + group_num + MMN_lab, data = all_phenos)
	
	# extract the coefficients estimated in the regression

	reg_coefs <- coef(adj_reg)

	# create a model matrix with formula above

	adj_modmat <- model.matrix(MMNFZAmplitude ~ age + gender + group_num + MMN_lab, data = all_phenos)
	
	# calculate aggregate effect of MMN-labs and extract it from the phenotype

	aggr_eff <- reg_coefs[grep("MMN_lab", names(reg_coefs))]%*%t(adj_modmat[,grep("MMN_lab", colnames(adj_modmat))]) 
	adj_pheno <- all_phenos$MMNFZAmplitude - aggr_eff
	
	# add the adjusted phenotype to the table and set it as the phenotype to be used in the analysis

	all_phenos$MMNFZAmplitude_adj <- adj_pheno[1,]

	just_phenos <- all_phenos[,c("FID", "IID", "MMNFZAmplitude_adj", "group_num")]
	
	} else {
        
        just_phenos <- all_phenos[,c("FID", "IID", "MMNFZAmplitude", "group_num")]
	}
	

	# now, create three separate tables:
	# - The phenotypes of interest
	# - The categorical covariates
	# - The quantitative covarites

#	just_phenos <- all_phenos[,c("FID", "IID", "MMNFZAmplitude", "group_num")]
	cat_covars <- all_phenos[,c("FID", "IID", "gender")]
	quant_covars <- all_phenos[,c("FID", "IID", "age")]

	# write tables out

	write.table(just_phenos, paste("./ERV_calculation/", datasets_2[i], "_phenotypes.txt", sep = ""), col.names = F, row.names = F, sep = "\t", quote = F)
	write.table(cat_covars, paste("./ERV_calculation/", datasets_2[i], "_categorial_covariates.txt", sep = ""), col.names = F, row.names = F, sep = "\t", quote = F)
	write.table(quant_covars, paste("./ERV_calculation/", datasets_2[i], "_quantitative_covariates.txt", sep = ""), col.names = F, row.names = F, sep = "\t", quote = F)

}
	




