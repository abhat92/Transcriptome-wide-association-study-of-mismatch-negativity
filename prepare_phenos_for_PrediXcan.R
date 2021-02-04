# set necessary variables

# path to phenotypes

pheno_path  <- "/SAN/neuroscience/PEIC/projects/EEG_genetics_consortium/phenotypic_data/merged_phenotypes"

# dataset

dataset="Mei"

# read in table that contains phenotypes for all datasets

phenos <- read.table(paste(pheno_path, "MMN_Merged_2018_12_06.txt", sep = "/"), header = T, sep = "\t")

# read in FAM file and generate table with FID, IID and the adjusted MMN Phenotype

# depending on the dataset the FAM file is in different location

if(dataset == "Elliot"){
  
  # read in FAM file
  
  fam <- read.table("/SAN/neuroscience/PEIC/projects/EEG_genetics_consortium/genetic_data/Elliot/imputed_data/Elliot2.vcfs/imputed_plink/matched_files/MPRC_Hong.typed.and.imputed.clean.dbSNP.150.matched.fam", header = F, sep = " ")
  
  # merge FID and IID with phenotypic table
  
  fam_phenos <- merge(fam[,c(1,2)], phenos[c("id", "FZAmp_adjusted")], by.x = 2, by.y = 1)
  colnames(fam_phenos) <- c("FID", "IID", "FZAmp_adjusted")
  
  # write table out
  
  write.table(fam_phenos, paste("./", dataset, "/fzamp_phenos_for_predixcan.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)
  
} else if (dataset == "Mei") {
  
  # read in FAM file
  
  fam <- read.table("/SAN/neuroscience/PEIC/projects/EEG_genetics_consortium/genetic_data/Mei/imputed_data/imputed_plink/matched_files/mclean_fgh19.typed.and.imputed.dbSNP.150.matched.fam", header = F, sep = " ")
  
  # merge FID and IID with phenotypic table
  
  fam_phenos <- merge(fam[,c(1,2)], phenos[c("id", "FZAmp_adjusted")], by.x = 1, by.y = 1)
  colnames(fam_phenos) <- c("FID", "IID", "FZAmp_adjusted")
  
  # write table out
  
  write.table(fam_phenos, paste("./", dataset, "/fzamp_phenos_for_predixcan.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)
  
  
} else if (dataset == "PEIC"){
  
  # read in FAM file
  
  fam <- read.table("/cluster/project9/PE/user/aritz/re_QCed_PEIC_genotypes/typed_and_imputed_genotypes/pe.dbSNP.150.matched.fam", header = F, sep = " ")
  
  # merge FID and IID with phenotypic table
  
  fam_phenos <- merge(fam[,c(1,2)], phenos[c("id", "FZAmp_adjusted")], by.x = 2, by.y = 1)
  fam_phenos <- fam_phenos[,c(2,1,3)]
  colnames(fam_phenos) <- c("FID", "IID", "FZAmp_adjusted")
  
  # write table out
  
  write.table(fam_phenos, paste("./", dataset, "/fzamp_phenos_for_predixcan.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)
  
}
