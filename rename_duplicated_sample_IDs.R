# read in FAM file

fam <- read.table("mclean_fhg19.sorted.fam", header = F, sep = " ")

# get number of entries per Individual-ID

iid_freqs <- table(fam$V2)

# get IIDs with more than 1 entry

multi_iids <- names(iid_freqs[iid_freqs > 1])

# turn fam object into a matrix

fam_mat <- as.matrix(fam)

# in a loop, rename the FIDs and IIDs of the rows that correspond to each duplicated IID

for (i in 1:length(multi_iids)) {
  
  # get rows in FAM file
  
  fam_rows <- which(fam_mat[,2] == multi_iids[i])

  # update FID

  fam_mat[fam_rows, 1] <- rep("dupliFAM", times = length(fam_rows))

  # update IID 
  
  fam_mat[fam_rows, 2] <- c(fam_mat[fam_rows[1], 2], paste(fam_mat[fam_rows[2:length(fam_rows)], 2], seq(from = 1, to = (length(fam_rows)-1), by = 1), sep = "_"))

}

# check if the renaming worked well

print(fam_mat[unlist(lapply(multi_iids, function(x) grep(x, fam_mat[,2]))),])

# write.out renamed FAM file

write.table(fam_mat, "mclean_fhg19.renamed.fam", col.names = F, row.names = F, sep = " ", quote = F)

# select duplicated entries to be removed from the dataset and output the table

del_samps <- fam_mat[intersect(grep("dupliFAM", fam_mat[,1]), grep("_", fam_mat[,2])), c(1,2)]

write.table(del_samps, "duplicated_samples_for_removal", col.names = F, row.names = F, sep = " ", quote = F)