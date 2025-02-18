#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 17:10:48 2024

@author: elise
"""
import pandas as pd


otu = pd.read_csv("/Users/elise/Downloads/PacBio_05_29_24/OTU_ps.csv", sep ='\t')
taxa = pd.read_csv("/Users/elise/Downloads/PacBio_05_29_24/TAXA_ps.csv", sep = ",")

otu.head()
taxa.head()

#to make indexes into ASV# instead of sequence
indexes = []
for i in range(1, len(taxa)+ 1): 
    name = "ASV" + str(i)
    indexes.append(name) 

taxa.index = indexes 
otu.index = indexes

#reformat OTU table 
ASV_list = otu.index

otu.set_index("#NAME")

read_counts = pd.DataFrame(otu.sum(axis='columns'))

read_counts = read_counts.set_index(ASV_list)

read_counts.head()

read_counts.columns = ["count"]

read_counts["Genus"] = taxa["Genus"]
read_counts["Species"] = taxa["Species"]
read_counts.head()

#-------------------------------------------------------------------
#remove non-bacteria 

taxa.head()

read_counts["Kingdom"] = taxa["Kingdom"]
read_counts["Phylum"] = taxa["Phylum"]

bact_only = read_counts.loc[read_counts["Kingdom"]=="Bacteria"]
phylum_adj = bact_only.dropna(subset=["Phylum"])


read_counts = phylum_adj


#------------------------------------------------------------------


#calculate total number of reads and ASVs 

total_reads = read_counts["count"].sum()
print(total_reads, "total reads")

total_ASVs = len(read_counts["Genus"])
print(total_ASVs, "total ASVs")


#percent of reads and ASVs classified to species level 
noNAs = read_counts.dropna(subset=["Species"])
spec_reads = noNAs["count"].sum()
spec_ASVs = len(noNAs["Genus"])

#for staph specifically
staph = read_counts.loc[read_counts["Genus"]=="Staphylococcus"]
staph_reads = staph["count"].sum()
staph_ASVs = len(staph["count"])

staph_noNAs = staph.dropna(subset=["Species"])
staph_spec_ASVs = len(staph_noNAs["Genus"])
staph_spec_reads = staph_noNAs["count"].sum()


print(spec_reads, "out of", total_reads, "reads classified at species level")
print(round(spec_reads/total_reads * 100, 2), "%")

print(spec_ASVs, "out of", total_ASVs, "ASVs classified at species level")
print(round(spec_ASVs/total_ASVs * 100, 2), "%")

print(staph_spec_reads, "out of", staph_reads, "staph reads classified at species level")
print(round(staph_spec_reads/staph_reads * 100, 2), "%")

print(staph_spec_ASVs, "out of", staph_ASVs, "staph ASVs classified at species level")
print(round(staph_spec_ASVs/staph_ASVs * 100, 2), "%")
