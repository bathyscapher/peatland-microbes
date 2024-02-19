################################################################################
### Peatland microbes network analysis
### Authors: Rachel Korn
### korn@cumulonimbus.at University of Fribourg 2021 - 2024
################################################################################

library("phyloseq")
library("reshape2")
library("igraph")
library("Hmisc")
library("indicspecies")


rm(list = ls())
gc()

set.seed(7113)


################################################################################
### ASV data ###################################################################
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


## Abbreviate Incertae_Sedis to IS and merge it with its higher taxonomy
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


## Make a copy for indicator species analysis
otus.tax.is <- otus.tax


## Subset to taxa occuring in at least 8 samples (= at least 20 %)
otus.pa <- otus
otus.pa[otus.pa > 0] <- 1
otus.pa <- otus.pa[rowSums(otus.pa) > 7, ]
summary(rowSums(otus.pa))

otus.tax.is <- otus.tax[otus.tax.is$ASV %in% rownames(otus.pa), ]
summary(rowSums(otus.pa))


## Reformat
otus.tax.is <- otus.tax.is[, -c(1:7)]
otus.tax.is <- data.frame(t(otus.tax.is))


## Metadata
mossMeta <- read.table("csv/Mosses_Metadata_unscaled.csv", sep = "\t",
                       header = TRUE)
summary(rownames(otus.tax.is) == mossMeta$FullID)


### Indicator species for soil reaction
range(mossMeta$Soil.reaction)
mossMeta$Soil.reactionC <- mossMeta$Soil.reaction
mossMeta$Soil.reactionC[mossMeta$Soil.reactionC < 2] <- "Strong acidic"
mossMeta$Soil.reactionC[mossMeta$Soil.reactionC < 4] <- "Acidic"
mossMeta$Soil.reactionC[mossMeta$Soil.reactionC < 5] <- "Moderate acidic"
mossMeta$Soil.reactionC
table(mossMeta$Soil.reactionC)


indval <- multipatt(otus.tax.is, mossMeta$Soil.reactionC,
                    control = how(nperm = 999), duleg = TRUE)
summary(indval)


sr <- data.frame(indval$str)


### Export indicator species table
sr.otus <- sr
sr.otus$ASV <- rownames(sr.otus)
sr.otus <- merge(sr.otus, otus.tax[, c(1:7)], by = "ASV")
sr.otus[, 2:4] <- format(round(sr.otus[, 2:4], 2), nsmall = 3)


# write.table(sr.otus, "csv/Mosses_IndicatorSpecies.csv", sep = "\t",
#             row.names = FALSE)


## Arrange for network analysis
sr[sr < 0.7] <- NA
sr$Acidic[!is.na(sr$Acidic)] <- "Acidic"
sr$Moderate.acidic[!is.na(sr$Moderate.acidic)] <- "Moderate acidic"
sr$Strong.acidic[!is.na(sr$Strong.acidic)] <- "Strong acidic"


sr$SoilReaction <- apply(sr, 1, function(x) x[!is.na(x)][1])
sr <- sr["SoilReaction"]
sr$ASV <- rownames(sr)


sr <- sr[complete.cases(sr), ]
table(sr)
table(sr$SoilReaction)


################################################################################
# Co-occurrence network with rank correlation ##################################
## Presence/absence and taxa that occur at least in 8 (so that they can all
## occur within a site) samples and then, subset


## Total abundance
otus.tax$Abundance <- rowSums(otus.tax[, 8:47])


## Correlation analysis based on Spearman's coefficient
otus.dist <- rcorr(t(otus), type = "spearman")
otus.cor <- otus.dist$r
otus.cor.p <- otus.dist$P


## Multiple testing correction using Benjamini-Hochberg standard FDR correction
otus.cor.p <- p.adjust(otus.cor.p, method = "BH")


## Positive and netagive cooccurence at given coefficient and p-value cutoff
cor.cutoff <- 0.7
p.cutoff <- 0.005


otus.cor[which(otus.cor >= (- cor.cutoff) & otus.cor <= cor.cutoff)] <- 0
otus.cor[which(otus.cor.p > p.cutoff)] <- 0
diag(otus.cor) <- 0


## Delete rows and columns with sum = 0
otus.cor <- otus.cor[which(rowSums(otus.cor) != 0), ]
otus.cor <- otus.cor[, which(colSums(otus.cor) != 0)]
dim(otus.cor)


### Create a graph
g3 <- graph.adjacency(otus.cor,
                      weight = TRUE, # weight = correlation
                      mode = "undirected")
g3 <- simplify(g3) # remove duplicate and loop edges


### Add taxonomy
## Get taxa in network and merge with taxonomy
target <- data.frame(ASV = V(g3)$name,
                     to_sort = seq(1:length(V(g3)$name)))


otu.target <- merge(target, otus.tax, by = "ASV")
otu.target <- otu.target[order(otu.target$to_sort), ]


table(sr$SoilReaction)
intersect(otu.target$ASV, sr$ASV)
length(intersect(otu.target$ASV, sr$ASV))


## Add indicator species
otu.target <- merge(otu.target, sr, by = "ASV", all.x = TRUE)
# table(otu.target$SoilReaction)
otu.target <- otu.target[order(otu.target$to_sort), ]


E(g3)
E(g3)$weight
V(g3)
V(g3)$name


### Add edge attributes
E(g3)$scale <- ifelse(E(g3)$weight < 0,
                      E(g3)$weight * -1,
                      E(g3)$weight)
E(g3)$sign <- ifelse(E(g3)$weight < 0,
                     "Negative",
                     "Positive")


### Add vertex attributes
g3 <- set_vertex_attr(g3, "degree", value = degree(g3))

g3 <- set_vertex_attr(g3, "ASV", value = otu.target$ASV)
V(g3)$ASV

g3 <- set_vertex_attr(g3, "Domain", value = otu.target$Domain)
V(g3)$Domain

g3 <- set_vertex_attr(g3, "Phylum", value = otu.target$Phylum)
V(g3)$Phylum

g3 <- set_vertex_attr(g3, "Genus", value = otu.target$Genus)
V(g3)$Genus

g3 <- set_vertex_attr(g3, "Abundance", value = sqrt(otu.target$Abundance))
V(g3)$Abundance


### Indicator species
V(g3)$SoilReaction <- otu.target$SoilReaction


## Export graph
write_graph(g3,
            paste("graphml/Co-occurrence", cor.cutoff, p.cutoff, ".graphml",
                  sep = "_"),
            format = "graphml")


## Network statistics
gsize(g3) # = ecount(g3)
vcount(g3)
vertex_connectivity(g3)

table(E(g3)$sign)
table(V(g3)$Domain)
table(V(g3)$SoilReaction)

table(degree(g3))
g3.degree <- as.data.frame(table(degree(g3)))
summary(degree(g3))
barplot(sort(degree_distribution(g3)))


## Subnetworks
components(g3)
