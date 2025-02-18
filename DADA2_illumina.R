if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")
BiocManager::install("interp") 

library(dada2)
remove.packages("dada2")
library(dada2); packageVersion("dada2")

getwd()
system("mkdir dada2tutorial")
setwd("/Users/elise/dada2tutorial/")
system("curl -O https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip")
system("unzip miseqsopdata.zip")

path <- "/Users/elise/Downloads/RacialDataHealthyControls/FASTQ" #containing the fastq files after unzipping.
list.files(path)

#reading in names of fastq files 
fnFs <- sort(list.files(path, pattern="_R1_001.trimmed.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.trimmed.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

#inspect read quality 
plotQualityProfile(fnFs[1:2]) #forward reads visualization
plotQualityProfile(fnRs[1:2]) #reverse reads 

#filtering 
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs)<-sample.names
names(filtRs)<-sample.names

filtFs

out<- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                    maxN=0, maxEE=c(2,2), truncQ = 2, rm.phix=TRUE, 
                    compress=TRUE, multithread=TRUE)
head(out)

#estimating error rates 
errF<- learnErrors(filtFs, multithread = TRUE)
errR<- learnErrors(filtRs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)

#sample interference 
dadaFs <- dada(filtFs, err=errF, multithread = TRUE)
dadaRs <- dada(filtRs, err=errR, multithread = TRUE)

dadaFs[[1]] #tells number of sequence variants found 

#MERGING: forward and reverse 
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 3, verbose=TRUE) 
head(mergers[[1]])

#construct ASV table 
seqtab<- makeSequenceTable(mergers)
seqtab
View(seqtab)
dim(seqtab)

table(nchar(getSequences(seqtab))) #table to inspect distribution of sequence lengths

#remove chimeras 
seqtab.nochim <- removeBimeraDenovo(seqtab, method= "consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim/sum(seqtab))
View(seqtab.nochim)

#tracking reads through entire pipeline 
#look at individual steps for drops and overall drop in reads
getN<- function(x) sum(getUniques(x))
track<- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track)<- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
track


#installing DECIPHER
BiocManager::install("DECIPHER") 
library(DECIPHER); packageVersion("DECIPHER")

dna <- DNAStringSet(getSequences(seqtab.nochim)) #makes string of ASVs 
load("/Users/elise/Downloads/SILVA_SSU_r138_2019.RData")
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) #took kind of a long time
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

taxa<-taxid

taxonomydf<- as.data.frame(taxa)
View(taxonomydf)




# LOAD PHYLOSEQ 
BiocManager::install("phyloseq") 
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw())





#IMPORT METADATA TABLE "samdf"

library(readxl)
samdf <- read_excel("Downloads/RacialDataHealthyControls/RacialDataNewHealthyControls.xlsx")
samdf<- as.data.frame(samdf)

samples.out <- rownames(seqtab.nochim)
rownames(samdf) <- samples.out
samdf<-samdf[, -1]
View(samdf)
class(samdf)


#CONSTRUCTING PHYLOSEQ OBJECT: otu table and taxonomy table 

# Creating phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

# Creating a DNAStringSet object from taxa names and assigning names to it
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
# Merging the DNAStringSet with the phyloseq object
ps <- merge_phyloseq(ps, dna)

# Assigning new names to the taxa in the phyloseq object
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

ps







#ALPHA DIVERSITY AND BETA DIVERSITY ANALYSIS 

# Creating an alpha diversity plot
plot_richness(ps, x="Age", measures=c("Shannon", "Simpson"), color="Race")

# Performing Bray Curtis beta diversity and NMDS ordination
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="PCoA", distance="bray")

# Creating a plot of the NMDS ordination with color based on the "When" variable
plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray PCoA")

# Creating a bar plot of the top 20 taxa based on their sums in the phyloseq object
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="family") + facet_wrap(~When, scales="free_x")






#ALPHA diversity statistical analysis 
library("vegan")
data_shannon <- diversity(otu_table(ps), index = "shannon") 
data_shannon



#----------FILTERING out non-bacteria and mycoplasma 

#remove non-bacteria
ps<-subset_taxa(ps, tax_table(ps)[, "domain"] == "Bacteria")
View(tax_table(ps))

#remove mycoplasma
ps<-subset_taxa(ps, tax_table(ps)[, "genus"] != "Mycoplasma")



#------COLLAPSING TAXONOMY according to level desired

coll_taxa<- tax_glom(ps, taxrank = rank_names(ps)[6], NArm=FALSE) #NArm=TRUE if i want to remove NAs for that level


#------MAKE DATA FRAMES of otu table, taxonomy table, and metadata table
#check OTU and tax tables
coll_taxa
View(otu_table(coll_taxa))
View(tax_table(coll_taxa))
colnames(otu_table(coll_taxa))==rownames(tax_table(coll_taxa))

#making OTU table, taxa, and metadata into data frames 
OTU1 = as(otu_table(coll_taxa), "matrix") #change ps to coll_taxa if tax_glom has been used 
TAXA1 = as(tax_table(coll_taxa), "matrix")
OTUdf = as.data.frame(OTU1)
TAXAdf = as.data.frame(TAXA1)
META1 = as(sample_data(coll_taxa), "matrix")
METAdf = as.data.frame(META1)

#------MICROBIOMEANALYST FORMATTING: microbiomeanalyst is a website that makes cool graphs of alpha and beta diversity for you

#transpose OTU table for microbiomeanalyst
View(OTUdf)
OTU_trans <- t(OTUdf)
View(OTU_trans)

#file paths for microbiomeanalyst files 
file_path_OTU <- "/Users/elise/Downloads/RacialDataHealthyControls/OTU.ma.txt" #change to your path
file_path_META <- "/Users/elise/Downloads/RacialDataHealthyControls/META.ma.csv" #change to your path
file_path_TAXA <- "/Users/elise/Downloads/RacialDataHealthyControls/TAXA.ma.csv" #change to your path 


# Use write.table to save the data table as a tab-delimited text file
write.table(OTU_trans, file = file_path_OTU, sep = "\t", row.names = TRUE, na = "NA")
write.csv(METAdf, file = file_path_META, row.names = TRUE)
write.csv(TAXAdf, file = file_path_TAXA, row.names = TRUE, na = "NA")




#--------SETUP FOR MAASLIN

#DECIDE WHAT TAX LEVELS TO INCLUDE: domain, phylum, class, order, family, and genus separated by "_" 
TAXAdf$taxa_all<-paste(TAXAdf$domain, TAXAdf$phylum, TAXAdf$class, TAXAdf$order, TAXAdf$family, TAXAdf$genus, sep = "_") #if I want full taxonomy listed
TAXAdf$taxa_fg<-paste(TAXAdf$family, TAXAdf$genus, sep = "_") #if I just want to name by family and genus 
View(TAXAdf)


#checking if column names in abundances data match rownames for taxa (sequences in same order)
colnames(OTUdf)==rownames(TAXAdf)

#IF ALL TRUE --> YES, GO AHEAD 

#assign taxonomy of each ASV to column names for OTU table
for (n in 1:ncol(OTUdf)) {
  if (TAXAdf$taxa_fg[n] != "NA_NA") {
    colnames(OTUdf)[n]<- TAXAdf$taxa_fg[n]
  } else {
    colnames(OTUdf)[n]<- paste("unclassified", "_", n)
  }
}

colnames(OTUdf)<-TAXAdf$taxa_fg #choose which column from taxa depending on which taxa i want to include
View(OTUdf)

#------for filtering out unclassified taxa

#INSTALL and LOAD TIDYVERSE/DPLYR
install.packages("tidyverse")
library("tidyverse")

#makes 1 dataframe without unclassified, other dataframe of only unclassified 
OTUdf_NAonly <- OTUdf %>% dplyr::select(starts_with("unclassified",ignore.case = TRUE))
OTUdf_NArm <- OTUdf %>% select(-contains("unclassified"))
View(OTUdf_NAonly)
View(OTUdf)
#compute sum unclassified for each sample, add row to OTU_NArm
NAsums<-c()

for (x in 1:nrow(OTUdf_NAonly)) {
  sum<- sum(OTUdf_NAonly[x, ])
  print(paste(x, ":", sum))
  NAsums[x]<- sum
  }

OTUdf_NAonly$total<-NAsums
OTUdf_NArm$other<-OTUdf_NAonly$total

#---------------------

#SAVING FILES at the end for MaAslin

# Load the "writexl" package/library
library("writexl")

#save rownames
#if feature table and metadata have different orders, do this twice 
sample_fasta_names<-c(rownames(OTUdf))
write(sample_fasta_names, "Downloads/RacialDataHealthyControls/DADA2_samplenames")

# Write abundances/taxonomy table to an Excel file located at "/Downloads/DADA2_features.xlsx"
write_xlsx(OTUdf, "Downloads/RacialDataHealthyControls/DADA2_features.xlsx") #set to keep classified separate

# Write metadata to an Excel file located at "/Downloads/DADA2_metadata.xlsx"
write_xlsx(samdf, "Downloads/RacialDataHealthyControls/DADA2_metadata.xlsx")
