if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install("dada2")
BiocManager::install("interp") 
BiocManager::install("Biostrings")
install.packages("writexl", repos = "http://cran.us.r-project.org")
install.packages("data.table", repos = "http://cran.us.r-project.org")

#load packages 
library(BiocManager)
library(dada2)
library(Biostrings)
library(readr)
packageVersion("readr")
library("writexl")
library(data.table)


print("libraries loaded")

# Read the CSV file
seqtab_nochim_imported <- fread("/home/estagaman/PacBio_05_29_24/seqtab_nochim.csv", header = TRUE)

# Set the first column as row names
rownames(seqtab_nochim_imported) <- seqtab_nochim_imported$V1

# Remove the first column (which are now row names)
seqtab_nochim <- seqtab_nochim_imported[, -1, with = FALSE]

rownames(seqtab_nochim) = seqtab_nochim_imported$V1

print("seqtab import done")

seqtab_matrix = as.matrix(seqtab_nochim)

rownames(seqtab_matrix) = seqtab_nochim_imported$V1

print("converted to matrix")

taxa <- assignTaxonomy(seqtab_matrix, "/home/estagaman/silva_nr99_v138.1_wSpecies_train_set.fa", multithread=TRUE)

print("taxonomy assigned. saving. . . ")

write.csv(taxa, file = "/home/estagaman/PacBio_05_29_24/taxa.csv")

save.image("tax_assign_workspace.RData")
