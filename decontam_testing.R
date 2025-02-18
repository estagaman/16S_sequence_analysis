#making data tables for decontam 

#import ASV table and taxa table if necessary 
#seqtab_matrix and taxa_matrix 

#import metadata table 
#PacBio_metadata_all 

library(readxl)
PacBio_metadata_all <- read_excel("PacBio_05_29_24/PacBio_metadata_all.xlsx", 
                                  sheet = "PS metadata - all ")
View(PacBio_metadata_all)

project_fastq = PacBio_metadata_all["fastq"]

project_fastq = unlist(c(project_fastq))
class(project_fastq)

#remove samples and taxa from other projects --> project_ps  
project_ps <- prune_samples(project_fastq, ps_species)
project_ps <- prune_taxa(taxa_sums(project_ps) > 0, project_ps)

#remove non-bacteria 
project_ps <-subset_taxa(project_ps, tax_table(project_ps)[, "Kingdom"] == "Bacteria")

#check things out 
sample_names(project_ps)

#need two separate phyloseq objects, one for each batch 
batch1 <- PacBio_metadata_all[PacBio_metadata_all$batch == 1, "fastq"]
batch2 <- PacBio_metadata_all[PacBio_metadata_all$batch == 2, "fastq"]
allfastq <- PacBio_metadata_all["fastq"]

batch1 = unlist(c(batch1))
batch2 = unlist(c(batch2))
allfastq = unlist(c(allfastq))

batch1_ps <- prune_samples(batch1, project_ps) #batch1_ps 
batch2_ps <- prune_samples(batch2, project_ps) #batch2_ps

#make metadata tables for each 
batch1_md <- subset(PacBio_metadata_all, batch == 1)
batch2_md <- subset(PacBio_metadata_all, batch == 2)

#batch1_md <- as.matrix(batch1_md)
#batch2_md <- as.matrix(batch2_md)
#PacBio_metadata_all <- as.matrix(PacBio_metadata_all)

#format metadata so sample name is the row name 
row.names(batch1_md) = batch1
row.names(batch2_md) = batch2
row.names(PacBio_metadata_all) = allfastq

#batch1_md <- sample_data(batch1_md)
#batch2_md <- sample_data(batch2_md)

#add the metadata tables to the phyloseq objects
#sample_data(batch1_ps) <- batch1_md
#sample_data(batch2_ps) <- batch2_md
otu = as(otu_table(project_ps), "matrix")
tax = as(tax_table(project_ps), "matrix")

otu1 = as(otu_table(batch1_ps), "matrix")
tax1 = as(tax_table(batch1_ps), "matrix")

otu2 = as(otu_table(batch2_ps), "matrix")
tax2 = as(tax_table(batch2_ps), "matrix")

# Create an index vector based on the desired order
desired_order = unlist(c(rownames(otu_table(project_ps))))
index <- match(rownames(PacBio_metadata_all), desired_order)

desired_order_1 = unlist(c(rownames(otu_table(batch1_ps))))
index_1 <- match(rownames(batch1_md), desired_order_1)

desired_order_2 = unlist(c(rownames(otu_table(batch2_ps))))
index_2 <- match(rownames(batch2_md), desired_order_2)

# Reorder the rows of the matrix based on the index vector
all_rows = rownames(PacBio_metadata_all)[order(index)]
PacBio_metadata_all <- PacBio_metadata_all[order(index), ]
PacBio_metadata_all = as.data.frame(PacBio_metadata_all)
rownames(PacBio_metadata_all) <- all_rows

batch1_rows = rownames(batch1_md)[order(index_1)]
batch1_md <- batch1_md[order(index_1), ]
batch1_md = as.data.frame(batch1_md)
rownames(batch1_md) <- batch1_rows

batch2_rows = rownames(batch2_md)[order(index_2)]
batch2_md <- batch2_md[order(index_2), ]
batch2_md = as.data.frame(batch2_md)
rownames(batch2_md) <- batch2_rows

#batch1_md <- as.data.frame(batch1_md)
#batch2_md <- as.data.frame(batch2_md)

#check for matching rownames
#for (name in rownames(batch1_md)){
 # here = name %in% rownames(otu_table(batch1_ps))
  #if (!here){
    #print()
 # }
#}

class(otu)
class(PacBio_metadata_all)
class(batch2_md)

#make phyloseq object including the metadata 
batch1_ps <- phyloseq(otu_table(otu1, taxa_are_rows=FALSE), 
                                   tax_table(tax1),
                                    sample_data(batch1_md))
                                   
batch2_ps <- phyloseq(otu_table(otu2, taxa_are_rows=FALSE), 
                      tax_table(tax2),
                      sample_data(batch2_md))

project_ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
                      tax_table(tax),
                      sample_data(PacBio_metadata_all))

#-------------------------
#now we can run decontam 
BiocManager::install("decontam")
library("decontam")

install.packages("ggplot2")
library("ggplot2")

sample_data(batch2_ps)$is.neg <- sample_data(batch2_ps)$Sample_Control == "Control"
#contamdf.prev <- isContaminant(batch2_ps, method="prevalence", neg="is.neg")
#table(contamdf.prev$contaminant)

#more strict threshold 
contamdf2.prev <- isContaminant(batch2_ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf2.prev$contaminant)

head(which(contamdf2.prev$contaminant))


#presence/absence in negative controls and true samples 
ps.pa <- transform_sample_counts(batch2_ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_Control == "Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf2.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#-----------------------making new, decontaminated phyloseqs for each batch 

head(contamdf.prev)

length(rownames(contam2))

contam <- subset(contamdf.prev, contaminant == TRUE)
contam1 <- subset(contamdf1.prev, contaminant == TRUE)
contam2 <- subset(contamdf2.prev, contaminant == TRUE)

noncontam <- subset(contamdf.prev, contaminant == FALSE)
noncontam1 <- subset(contamdf1.prev, contaminant == FALSE)
noncontam2 <- subset(contamdf2.prev, contaminant == FALSE)

#pull out just ASV names from contam and noncontam 
contaminant_asvs = unlist(c(rownames(contam)))
contaminant_asvs_1 = unlist(c(rownames(contam1)))
contaminant_asvs_2 = unlist(c(rownames(contam2)))

# Prune the original phyloseq object to exclude contaminants
decontam_ps <- prune_taxa(non_contaminant_asvs, project_ps)
decontam1_ps <- prune_taxa(non_contaminant_asvs_1, batch1_ps)
decontam2_ps <- prune_taxa(non_contaminant_asvs_2, batch2_ps)

#ps of just the contaminants for each batch 
contam_only_ps <- prune_taxa(contaminant_asvs, project_ps)
contam_only_1_ps <- prune_taxa(contaminant_asvs_1, batch1_ps)
contam_only_2_ps <- prune_taxa(contaminant_asvs_2, batch2_ps)


#---------making taxa tables with prevalence included
taxa = as(tax_table(contam_only_ps), "matrix")
taxa = cbind(taxa,contam$p.prev)
colnames(taxa)[8] <- "prevalence"

taxa1 = as(tax_table(contam_only_1_ps), "matrix")
taxa1 = cbind(taxa1,contam1$p.prev)
colnames(taxa1)[8] <- "prevalence"

taxa2 = as(tax_table(contam_only_2_ps), "matrix")
taxa2 = cbind(taxa2,contam2$p.prev)
colnames(taxa2)[8] <- "prevalence"

taxa = as.data.frame(taxa)
taxa1 = as.data.frame(taxa1)
taxa2 = as.data.frame(taxa2)

write.csv(taxa, file = "/home/estagaman/PacBio_05_29_24/contaminants/taxa_both", row.names = TRUE, na = "NA")
write.csv(taxa1, file = "/home/estagaman/PacBio_05_29_24/contaminants/taxa_1", row.names = TRUE, na = "NA")
write.csv(taxa2, file = "/home/estagaman/PacBio_05_29_24/contaminants/taxa_2", row.names = TRUE, na = "NA")
