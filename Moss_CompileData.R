################################################################################
################################################################################
################################################################################
################################################################################
### Peatland mosses: read taxa (ASV), add metadata and more
### Authors: Rachel Korn
### korn@cumulonimbus.at University of Fribourg 2021
################################################################################


library("phyloseq")


rm(list = ls())


################################################################################
### Arguments

source("Moss_readMetadataAndASV.R")


################################################################################
### Read pitcher taxa

asvDir <- "../smp_results"


### Prokaryotes
prok.l <- readTaxa(primer = "16S")
prok.l


list2env(prok.l, globalenv())

prok.a <- moss.a
rm(moss, moss.a, moss.s, moss.r, prok.l)


### Eukaryotes
euk.l <- readTaxa(primer = "18S")
euk.l

list2env(euk.l, globalenv())


euk.a <- moss.a
rm(moss, moss.a, moss.s, moss.r, euk.l)


## Save
saveRDS(prok.a, "rds/Moss_prok.a.RDS")
saveRDS(euk.a, "rds/Moss_euk.a.RDS")


### Merge both
## Renumber taxa
taxa_names(prok.a) <- paste0("ASV", seq(ntaxa(prok.a)))
taxa_names(prok.a)


taxa_names(euk.a) <- paste0("ASV", seq(ntaxa(prok.a) + 1,
                                       ntaxa(prok.a) + ntaxa(euk.a)))
taxa_names(euk.a)


moss <- merge_phyloseq(prok.a, euk.a)
tax_table(moss)
moss


### Save
saveRDS(moss, "rds/Moss_moss.RDS")


################################################################################
################################################################################
################################################################################
################################################################################
