################################################################################
### Peatland microbes network analysis
### Code to perform Mantel and partial Mantel tests between co-occurrence
### network and EV-environmental variables (light, temp, hum, pH).
### Analysis performed for all species and separately for Eu- and Prokaryotes
### Sub-Authors: Louis-Félix Bersier and Sarah M Gray
################################################################################


# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("phyloseq")

library("reshape2")
library("igraph")
library("Hmisc")
library("ecodist")


rm(list = ls())
gc()


################################################################################
### 1) From Rachel's code: matrix of otus taxonomy and "abundance" in the
#      40 sampling sites and pruning otus with less than 8 occurrences

moss <- readRDS("rds/Moss_moss.RDS")
moss

## Extract ASVs, taxonomy and merge
otus <- data.frame(t(otu_table(moss)))
tax <- data.frame(tax_table(moss))

otus.tax <- merge(tax, otus, by = 0)
names(otus.tax)[1] <- "ASV"
head(otus.tax)


## Fill unannotated genera with higher-level taxa
otus.tax$Genus <- ifelse(is.na(otus.tax$Genus), otus.tax$Family, otus.tax$Genus)
otus.tax$Genus <- ifelse(is.na(otus.tax$Genus), otus.tax$Order, otus.tax$Genus)
otus.tax$Genus <- ifelse(is.na(otus.tax$Genus), otus.tax$Class, otus.tax$Genus)
otus.tax$Genus <- ifelse(is.na(otus.tax$Genus), otus.tax$Phylum, otus.tax$Genus)


## Abbreviate Incertae_Sedis to IS and merge it with higher taxonomy
otus.tax$Genus[otus.tax$Genus == "Incertae_Sedis"] <- "IS"
otus.tax[4:6][otus.tax[4:6] == "Incertae_Sedis"] <- NA

otus.tax$Genus[!is.na(otus.tax$Family) & otus.tax$Genus == "IS"] <- paste(
  otus.tax$Family[!is.na(otus.tax$Family) & otus.tax$Genus == "IS"],
  otus.tax$Genus[!is.na(otus.tax$Family) & otus.tax$Genus == "IS"],
  sep = ":"
)
otus.tax$Genus[!is.na(otus.tax$Order) & otus.tax$Genus == "IS"] <- paste(
  otus.tax$Order[!is.na(otus.tax$Order) & otus.tax$Genus == "IS"],
  otus.tax$Genus[!is.na(otus.tax$Order) & otus.tax$Genus == "IS"],
  sep = ":"
)
otus.tax$Genus[!is.na(otus.tax$Class) & otus.tax$Genus == "IS"] <- paste(
  otus.tax$Class[!is.na(otus.tax$Class) & otus.tax$Genus == "IS"],
  otus.tax$Genus[!is.na(otus.tax$Class) & otus.tax$Genus == "IS"],
  sep = ":"
)
otus.tax$Genus[!is.na(otus.tax$Phylum) & otus.tax$Genus == "IS"] <- paste(
  otus.tax$Phylum[!is.na(otus.tax$Phylum) & otus.tax$Genus == "IS"],
  otus.tax$Genus[!is.na(otus.tax$Phylum) & otus.tax$Genus == "IS"],
  sep = ":"
)

grep("IS", otus.tax$Genus, value = TRUE)


otus.tax$Genus <- make.unique(otus.tax$Genus)
rownames(otus.tax) <- otus.tax$ASV


## Subset to taxa occuring in at least 8 samples (= at least 20 %)
otus.pa <- otus
otus.pa[otus.pa > 0] <- 1
otus.pa <- otus.pa[rowSums(otus.pa) > 7, ]
summary(rowSums(otus.pa))

otus.tax <- otus.tax[otus.tax$ASV %in% rownames(otus.pa), ]
summary(rowSums(otus.pa))


# create separate dataframe for bacteria and eukaryota
otus.tax.euka <- otus.tax[otus.tax$Domain=="Eukaryota",]
otus.tax.prok <- otus.tax[otus.tax$Domain=="Bacteria",]


## Reformat
otus.tax <- otus.tax[, -c(1:7)]
otus.tax <- data.frame(otus.tax)
str(otus.tax)

otus.tax.euka <- otus.tax.euka[, -c(1:7)]
otus.tax.euka <- data.frame(otus.tax.euka)
str(otus.tax.euka)

otus.tax.prok <- otus.tax.prok[, -c(1:7)]
otus.tax.prok <- data.frame(otus.tax.prok)
str(otus.tax.prok)


################################################################################
### 2) For each species, compute its "optimal" environmental condition
#      (based on a weighted average of the values of the 4 EV environmental variables
#      on the abundance)
#   Rem: the 4 EV variables correspond to the ones of the SEM analysis (Habitat condition)

# Load EV value:
EV.sites <- read.csv("manteltest/EVMeansLFB.csv", header=T, row.names = 1)
str(EV.sites)

pairs(t(EV.sites[1:5,]))
# --> we see that Soil_reaction is highly correlated to Nitrogen

# Compute "optimum" and store in otus.tax

# add empty variables for optima in otus.tax
otus.tax$opt.light <- 0
otus.tax$opt.temp <- 0
otus.tax$opt.hum <- 0
otus.tax$opt.pH <- 0

otus.tax.prok$opt.light <- 0
otus.tax.prok$opt.temp <- 0
otus.tax.prok$opt.hum <- 0
otus.tax.prok$opt.pH <- 0

otus.tax.euka$opt.light <- 0
otus.tax.euka$opt.temp <- 0
otus.tax.euka$opt.hum <- 0
otus.tax.euka$opt.pH <- 0



# Compute optima:
# all data
for(i in 1:4) {
  for(j in 1:dim(otus.tax)[1]) {
    otus.tax[j,40+i] <- crossprod(log(as.numeric(otus.tax[j,1:40])+1),
                                  as.numeric(EV.sites[i,])) / sum(log(otus.tax[j,1:40]+1))
  }
}

str(otus.tax)


# prokaryota
for(i in 1:4) {
  for(j in 1:dim(otus.tax.prok)[1]) {
    otus.tax.prok[j,40+i] <- crossprod(log(as.numeric(otus.tax.prok[j,1:40])+1),
                                       as.numeric(EV.sites[i,])) / sum(log(otus.tax.prok[j,1:40]+1))
  }
}

str(otus.tax.prok)


# eukaryota
for(i in 1:4) {
  for(j in 1:dim(otus.tax.euka)[1]) {
    otus.tax.euka[j,40+i] <- crossprod(log(as.numeric(otus.tax.euka[j,1:40])+1),
                                       as.numeric(EV.sites[i,])) / sum(log(otus.tax.euka[j,1:40]+1))
  }
}

str(otus.tax.euka)


################################################################################
### 3) Correlation matrix with Spearman's rank correlation (as in Rachel's
#      code, but without pruning of correlations with r>0.7 and p<0.005). All
#      data are kept.
#      Compute distance matrices for each EV environmental variable.

# Matrix of correlation (1 if both taxa occur in same sites and have the same
# ordering in abundances in these sites; note that r can be negative)

# All taxa
otus.dist <- rcorr(t(otus.tax[,1:40]), type = "spearman")$r
otus.dist
str(otus.dist)
# transform the correlation matrix in a distance matrix (dist class)
otus.dist <- as.dist(otus.dist)
str(otus.dist)

# Prokaryota
otus.dist.prok <- rcorr(t(otus.tax.prok[,1:40]), type = "spearman")$r
otus.dist.prok
str(otus.dist.prok)
# transform the correlation matrix in a distance matrix (dist class)
otus.dist.prok <- as.dist(otus.dist.prok)
str(otus.dist.prok)

# Eukarya
otus.dist.euka <- rcorr(t(otus.tax.euka[,1:40]), type = "spearman")$r
otus.dist.euka
str(otus.dist.euka)
# transform the correlation matrix in a distance matrix (dist class)
otus.dist.euka <- as.dist(otus.dist.euka)
str(otus.dist.euka)


# Matrices of Euclidean distances for the EV values (Nitrogen excluded)
# all taxa
light.dist <- dist(otus.tax[,41])
temp.dist <- dist(otus.tax[,42])
hum.dist <- dist(otus.tax[,43])
pH.dist <- dist(otus.tax[,44])

# prokaryota
light.dist.prok <- dist(otus.tax.prok[,41])
temp.dist.prok <- dist(otus.tax.prok[,42])
hum.dist.prok <- dist(otus.tax.prok[,43])
pH.dist.prok <- dist(otus.tax.prok[,44])

# eukaryota
light.dist.euka <- dist(otus.tax.euka[,41])
temp.dist.euka <- dist(otus.tax.euka[,42])
hum.dist.euka <- dist(otus.tax.euka[,43])
pH.dist.euka <- dist(otus.tax.euka[,44])


################################################################################
### 4) Simple Mantel tests

# set the number of permutations to estimate the p-value of the Mantel tests (10'000 in the manuscript)
np <- 1000

### 4.1) Simple Mantel tests - All taxa
m.light <- mantel(otus.dist~light.dist,mrank=T,nperm=np)
m.temp <- mantel(otus.dist~temp.dist,mrank=T,nperm=np)
m.hum <- mantel(otus.dist~hum.dist,mrank=T,nperm=np)
m.pH <- mantel(otus.dist~pH.dist,mrank=T,nperm=np)

### 4.2) Simple Mantel tests - Prokaryota
m.light.prok <- mantel(otus.dist.prok~light.dist.prok,mrank=T,nperm=np)
m.temp.prok <- mantel(otus.dist.prok~temp.dist.prok,mrank=T,nperm=np)
m.hum.prok <- mantel(otus.dist.prok~hum.dist.prok,mrank=T,nperm=np)
m.pH.prok <- mantel(otus.dist.prok~pH.dist.prok,mrank=T,nperm=np)

### 4.3) Simple Mantel tests - Eukaryota
m.light.euka <- mantel(otus.dist.euka~light.dist.euka,mrank=T,nperm=np)
m.temp.euka <- mantel(otus.dist.euka~temp.dist.euka,mrank=T,nperm=np)
m.hum.euka <- mantel(otus.dist.euka~hum.dist.euka,mrank=T,nperm=np)
m.pH.euka <- mantel(otus.dist.euka~pH.dist.euka,mrank=T,nperm=np)


################################################################################
### 5) Partial Mantel tests

### 5.1) Partial Mantel tests - All taxa
mp.light <- mantel(otus.dist~light.dist+temp.dist+hum.dist+pH.dist,mrank=T,
                   nperm=np)
mp.temp <- mantel(otus.dist~temp.dist+light.dist+hum.dist+pH.dist,mrank=T,
                  nperm=np)
mp.hum <- mantel(otus.dist~hum.dist+light.dist+temp.dist+pH.dist,mrank=T,
                 nperm=np)
mp.pH <- mantel(otus.dist~pH.dist+light.dist+temp.dist+hum.dist,mrank=T,
                nperm=np)


### 5.2) Partial Mantel tests - Prokaryota
mp.light.prok <- mantel(otus.dist.prok~light.dist.prok+temp.dist.prok+hum.dist.prok+pH.dist.prok,
                        mrank=T,nperm=np)
mp.temp.prok <- mantel(otus.dist.prok~temp.dist.prok+light.dist.prok+hum.dist.prok+pH.dist.prok,
                       mrank=T,nperm=np)
mp.hum.prok <- mantel(otus.dist.prok~hum.dist.prok+light.dist.prok+temp.dist.prok+pH.dist.prok,
                      mrank=T,nperm=np)
mp.pH.prok <- mantel(otus.dist.prok~pH.dist.prok+light.dist.prok+temp.dist.prok+hum.dist.prok,
                     mrank=T,nperm=np)


### 5.3) Partial Mantel tests - Eukaryota
mp.light.euka <- mantel(otus.dist.euka~light.dist.euka+temp.dist.euka+hum.dist.euka+pH.dist.euka,
                        mrank=T,nperm=np)
mp.temp.euka <- mantel(otus.dist.euka~temp.dist.euka+light.dist.euka+hum.dist.euka+pH.dist.euka,
                       mrank=T,nperm=np)
mp.hum.euka <- mantel(otus.dist.euka~hum.dist.euka+light.dist.euka+temp.dist.euka+pH.dist.euka,
                      mrank=T,nperm=np)
mp.pH.euka <- mantel(otus.dist.euka~pH.dist.euka+light.dist.euka+temp.dist.euka+hum.dist.euka,
                     mrank=T,nperm=np)


### Results:
m.light
mp.light
m.temp
mp.temp
m.hum
mp.hum
m.pH
mp.pH

m.light.prok
mp.light.prok
m.temp.prok
mp.temp.prok
m.hum.prok
mp.hum.prok
m.pH.prok
mp.pH.prok

m.light.euka
mp.light.euka
m.temp.euka
mp.temp.euka
m.hum.euka
mp.hum.euka
m.pH.euka
mp.pH.euka


################################################################################
### 6) Are co-occurring taxa similar in moss taxa preference ? (see SEM result)

otus.tax$opt.mossNMDS <- 0
otus.tax$opt.Altitude <- 0
otus.tax$opt.CanopyCover <- 0
# Compute optima:

for(j in 1:dim(otus.tax)[1]) {
  otus.tax[j,45] <- crossprod(log(as.numeric(otus.tax[j,1:40])+1),
                              as.numeric(EV.sites[8,])) / sum(log(otus.tax[j,1:40]+1))
}
for(j in 1:dim(otus.tax)[1]) {
  otus.tax[j,46 ] <- crossprod(log(as.numeric(otus.tax[j,1:40])+1),
                               as.numeric(EV.sites[7,])) / sum(log(otus.tax[j,1:40]+1))
}
for(j in 1:dim(otus.tax)[1]) {
  otus.tax[j,47] <- crossprod(log(as.numeric(otus.tax[j,1:40])+1),
                              as.numeric(EV.sites[6,])) / sum(log(otus.tax[j,1:40]+1))
}


# Distance matrices
moss.dist <- dist(otus.tax[,45])
alt.dist <- dist(otus.tax[,46])
ccov.dist <- dist(otus.tax[,47])

m.moss <- mantel(otus.dist~moss.dist,
                 mrank=T,nperm=np)
pm.moss1 <- mantel(otus.dist~moss.dist+alt.dist+ccov.dist,
                   mrank=T,nperm=np)

m.moss
pm.moss1
