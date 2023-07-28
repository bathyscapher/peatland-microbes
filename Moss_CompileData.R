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
# rm(moss, moss.a, moss.s, moss.r, prok.l)


### Eukaryotes
euk.l <- readTaxa(primer = "18S")
euk.l

list2env(euk.l, globalenv())


euk.a <- moss.a
rm(moss.a, euk.l)


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
moss


### Save
saveRDS(moss, "rds/Moss_moss.RDS")


### Export taxa table
taxa <- as.data.frame(tax_table(moss)@.Data)
otus <- as.data.frame(t(otu_table(moss)@.Data))
taxa_otus <- merge(taxa,
                   otus,
                   by = 0)
names(taxa_otus)[1] <- "ASV"
taxa_otus[, -c(1:7)] <- round(taxa_otus[, -c(1:7)],
                              digits = 7)
taxa_otus <- taxa_otus[do.call(order, taxa_otus[, c(2:7)]), ]



# write.table(moss_taxa_otus, "Korn_et_al_Taxa-table.csv",row.names = FALSE)



################################################################################
################################################################################
################################################################################
################################################################################
