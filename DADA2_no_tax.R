#basic package installation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install("dada2")
BiocManager::install("interp") 
BiocManager::install("Biostrings")
install.packages("writexl", repos = "http://cran.us.r-project.org")
install.packages("readr", repos = "http://cran.us.r-project.org")
#load packages 
library(BiocManager)
library(dada2)
library(Biostrings)
library(readr)
packageVersion("readr")
library("writexl")

#load fastq files 
path <- "PacBio_05_29_24/demux_fastq" #CHANGE TO LOCATION OF CCS FILES 
list.files(path)

#reading in names of fastq files 
fnFs <- sort(list.files(path, pattern=".fastq.gz", full.names = TRUE)) #no reverse reads

sample.names <-basename(fnFs)

#REMOVE PRIMERS using CUTADAPT-------------------------------------------------------------------------------------

#cutadapt <- path.expand("/usr/bin/cutadapt")
#system2(cutadapt, args = "--version") 

cut_dir <- file.path("/home/estagaman/PacBio_05_29_24/demux_fastq", "cutadapt")
#if (!dir.exists(cut_dir)) dir.create(cut_dir)

cut <- file.path(cut_dir, basename(fnFs))

names(cut) <- sample.names

# It's good practice to keep some log files so let's create some
# file names that we can use for those 
#cut_logs <- path.expand(file.path(cut_dir, paste0(sample.names, ".log")))

F27 <- "TTTCTGTTGGTGCTGATATTGCAGRGTTYGATYMTGGCTCAG"
R1492 <- "ACTTGCCTGTCGCTCTATCTTCRGYTACCTTGTTACGACTT"
rc <- dada2:::rc

#cutadapt_args <- c("-g", F27, "-a",dada2:::rc(R1492),
                   #"--revcomp",
                   #"-n", 2)

# Loop over the list of files, running cutadapt on each file.  If you don't have a vector of sample names or 
# don't want to keep the log files you can set stdout = "" to output to the console or stdout = NULL to discard
#for (i in seq_along(fnFs)) {
  #system2(cutadapt, 
          #args = c(cutadapt_args,
                   #"-o", cut[i], 
                   #fnFs[i]),
          #stdout = cut_logs[i])  
#}

# quick check that we got something
head(list.files(cut_dir))
#---------------------------------------------------------------

#inspect length distribution
#lens.fn <- lapply(cut, function(fn) nchar(getSequences(fn))) 
#lens <- do.call(c, lens.fn)
#hist(lens, breaks =100, xlim=c(1000, 2000))   #should all be slightly less than 1500

#inspect read quality 
#plotQualityProfile(cut[1:2]) #reads visualization

#filtering 
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs)<-sample.names

#filtFs
packageVersion("")

out<- filterAndTrim(cut, filtFs, minQ = 3, minLen = 100, maxLen = 1600,
                    maxN=0, maxEE=2, rm.phix=FALSE, 
                    compress=TRUE, multithread = TRUE)  
head(out)

#estimating error rates 
errF<- learnErrors(filtFs, errorEstimationFunction = PacBioErrfun, BAND_SIZE=32, multithread = TRUE) #increased band size because of high indel rate
plotErrors(errF, nominalQ = TRUE)

#sample interference 
dd1 <- dada(filtFs, err=errF, BAND_SIZE = 32, multithread = TRUE) #adjusted band size due to higher rate of indels

#dd1[[1]] #tells number of sequence variants found

#saveRDS(dd1, file.path(path.rds, "dada2.rds"))

#read tracking
track <- cbind(reads = out[,1], filtered=out[,2], denoised=sapply(dd1, function(x) sum(x$denoised)))
trackdf = as.data.frame(track)
file_path_track = "/home/estagaman/PacBio_05_29_24/track.csv"
write.csv(trackdf, file = file_path_track, col.names = TRUE)

#construct ASV table
seqtab<- makeSequenceTable(dd1)
#seqtab
#View(seqtab)
#dim(seqtab)
table(nchar(getSequences(seqtab))) #table to inspect distribution of sequence lengths

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method= "consensus", minFoldParentOverAbundance=3.5, multithread=TRUE, verbose=TRUE)
#dim(seqtab.nochim)
chim_removed = sum(seqtab.nochim/sum(seqtab))
write(chim_removed, file = "/home/estagaman/PacBio_05_29_24/chim_removal_percent.txt")
#View(seqtab.nochim)

#taxonomy assignment with AssignTaxonomy
class(seqtab.nochim)

write.csv(seqtab.nochim, file = "/home/estagaman/PacBio_05_29_24/seqtab_nochim.csv", row.names = TRUE, col.names= TRUE, na = "NA")

seqtab_nochim_imported <- read.csv("/home/estagaman/PacBio_05_29_24/seqtab_nochim.csv", header = TRUE, row.names=1)

seqtab_matrix = as.matrix(seqtab_nochim_imported)

taxa <- assignTaxonomy(seqtab_matrix, "/home/estagaman/silva_nr99_v138.1_wSpecies_train_set.fa", multithread=TRUE)

write.csv(taxa, file = "/home/estagaman/PacBio_05_29_24/taxa.csv")

#------------------------------------------------
#remove rownames for easier inspection
# taxa.print <- taxa # Removing sequence rownames for display only
# rownames(taxa.print) <- NULL
# head(taxa.print)

#convert taxa to matrix 
# taxonomydf <- as.data.frame(taxa) #converts taxa to data frame if useful 
# View(taxonomydf)
# 
# taxa_matrix = as.matrix(taxonomydf)

#------------------------------------------------------------------------------------------------
#CREATE PHYLOSEQ OBJECT
# library(phyloseq)
# 
# ps_species <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
#                        tax_table(taxa_matrix))
# 
# # Creating a DNAStringSet object from taxa names and assigning names to it
# dna <- Biostrings::DNAStringSet(taxa_names(ps_species))
# names(dna) <- taxa_names(ps_species)
# 
# # Merging the DNAStringSet with the phyloseq object
# ps_species <- merge_phyloseq(ps_species, dna)
# 
# # Assigning new names to the taxa in the phyloseq object
# taxa_names(ps_species) <- paste0("ASV", seq(ntaxa(ps_species)))
# 
# OTU1 = as(otu_table(ps_species), "matrix") #change ps to coll_taxa if tax_glom has been used 
# TAXA1 = as(tax_table(ps_species), "matrix")
# OTUdf = as.data.frame(OTU1)
# TAXAdf = as.data.frame(TAXA1)
# #META1 = as(sample_data(coll_taxa), "matrix")
# #METAdf = as.data.frame(META1)
# 
# #transpose OTU table for microbiomeanalyst
# #View(OTUdf)
# OTU_trans <- t(OTUdf)
# 
# #file paths for microbiomeanalyst files 
# file_path_OTU <- "/Users/elise/Downloads/PacBio_1/1985_fastq/test_set/OTU.csv"  #change path 
# #file_path_META <- "/Users/elise/Downloads/DADA2_PacBio_tutorial/META.ma.csv"
# file_path_TAXA <- "/Users/elise/Downloads/PacBio_1/1985_fastq/test_set/TAXA.csv"   #change path 
# 
# # Use write.table to save the data table as a tab-delimited text file
# write.table(OTU_trans, file = file_path_OTU, sep = "\t", row.names = TRUE, na = "NA")
# #write.csv(METAdf, file = file_path_META, row.names = TRUE)
# write.csv(TAXAdf, file = file_path_TAXA, row.names = TRUE, na = "NA")
