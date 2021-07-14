################################################################################
################################################################################
################################################################################
################################################################################
### SMP Read taxa (OTU or ASV), add metadata and much more
### Authors: Rachel Korn
### korn@cumulonimbus.at University of Fribourg 2019/2020
################################################################################


library("phyloseq")
library("reshape2")
library("igraph")
# library("vegan")
library("Hmisc")
library("indicspecies")


setwd("~/Seafile/_FNS_2016-NicheDynamCoexist/SarraceniaMicrobiomeProject/Mosses")


rm(list = ls())


set.seed(7113)


################################################################################
### Pitcher data
# prok.a <- readRDS("SMP_prok.a.RDS")
# euk.a <- readRDS("SMP_euk.a.RDS")
moss <- readRDS("SMP_moss40.RDS")


################################################################################
otus <- data.frame(t(otu_table(moss)))


## Extract taxonomy and merge
tax <- data.frame(tax_table(moss))
otus.tax <- merge(tax, otus, by = 0)
names(otus.tax)[1] <- "ASV"


## Fill unannotated genera with higher-level taxa
otus.tax$Genus <- ifelse(is.na(otus.tax$Genus), otus.tax$Family, otus.tax$Genus)
otus.tax$Genus <- ifelse(is.na(otus.tax$Genus), otus.tax$Order, otus.tax$Genus)
otus.tax$Genus <- ifelse(is.na(otus.tax$Genus), otus.tax$Class, otus.tax$Genus)
otus.tax$Genus <- ifelse(is.na(otus.tax$Genus), otus.tax$Phylum, otus.tax$Genus)


## Abbreviate Incertae_Sedis to IS and merge it with higher taxonomy
otus.tax$Genus[otus.tax$Genus == "Incertae_Sedis"] <- "IS"

subset(otus.tax[1:7], subset = (Phylum == "Incertae_Sedis" |
                                  Order == "Incertae_Sedis" |
                                  Class == "Incertae_Sedis" |
                                  Family == "Incertae_Sedis"))
otus.tax[1:7][otus.tax[1:7] == "Incertae_Sedis"] <- NA

otus.tax$Genus <- ifelse(!is.na(otus.tax$Family) & otus.tax$Genus == "IS",
                         paste(otus.tax$Family, otus.tax$Genus, sep = ":"),
                         otus.tax$Genus)
otus.tax$Genus <- ifelse(!is.na(otus.tax$Order) & otus.tax$Genus == "IS",
                         paste(otus.tax$Order, otus.tax$Genus, sep = ":"),
                         otus.tax$Genus)
otus.tax$Genus <- ifelse(!is.na(otus.tax$Class) & otus.tax$Genus == "IS",
                         paste(otus.tax$Class, otus.tax$Genus, sep = ":"),
                         otus.tax$Genus)
otus.tax$Genus <- ifelse(!is.na(otus.tax$Phylum) & otus.tax$Genus == "IS",
                         paste(otus.tax$Phylum, otus.tax$Genus, sep = ":"),
                         otus.tax$Genus)
grep("IS", otus.tax$Genus, value = TRUE)


otus.tax$Genus <- make.unique(otus.tax$Genus)
rownames(otus.tax) <- otus.tax$Genus


## Make a copy for indicator species analysis
otus.tax.is <- otus.tax


## Only taxa that occur at least in 8 samples (= at least 20 %) and then, subset
otus.pa <- otus
otus.pa[otus.pa > 0] <- 1
otus.pa <- otus.pa[rowSums(otus.pa) > 7, ]
summary(rowSums(otus.pa))
subset.inc <- rownames(otus.pa)


otus.tax.is <- otus.tax[otus.tax.is$ASV %in% subset.inc, ]


## Reformat
otus.tax.is <- otus.tax.is[, -c(1:7)]
otus.tax.is <- data.frame(t(otus.tax.is))


## Metadata
mossMeta <- read.table("Mosses_Metadata_unscaled.csv", sep = "\t",
                       header = TRUE)
rownames(otus.tax.is) == mossMeta$FullID


### Indicator species
## For sites
indval <- multipatt(otus.tax.is, mossMeta$Site,
                    control = how(nperm = 999), duleg = TRUE)
summary(indval)


## For altitude
range(mossMeta$Altitude)
mossMeta$AltitudeC <- mossMeta$Altitude
mossMeta$AltitudeC[mossMeta$AltitudeC < 800] <- "Lowland"
mossMeta$AltitudeC[mossMeta$AltitudeC < 1100] <- "Montane"
mossMeta$AltitudeC[mossMeta$AltitudeC < 1500] <- "Subalpine"
mossMeta$AltitudeC
table(mossMeta$AltitudeC)


indval <- multipatt(otus.tax.is, mossMeta$AltitudeC,
                    control = how(nperm = 999), duleg = TRUE)
summary(indval)


alt <- data.frame(indval$str)
alt[alt < 0.7] <- NA
alt$Lowland[!is.na(alt$Lowland)] <- "Lowland"
alt$Montane[!is.na(alt$Montane)] <- "Montane"
alt$Subalpine[!is.na(alt$Subalpine)] <- "Subalpine"


alt$Altitude <- apply(alt, 1, function(x) x[!is.na(x)][1])
alt <- alt[-c(1:3)]
alt$Genus <- rownames(alt)

alt[complete.cases(alt),]
table(alt)


## For soil reaction
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
sr[sr < 0.7] <- NA
sr$Acidic[!is.na(sr$Acidic)] <- "Acidic"
sr$Moderate.acidic[!is.na(sr$Moderate.acidic)] <- "Moderate acidic"
sr$Strong.acidic[!is.na(sr$Strong.acidic)] <- "Strong acidic"


sr$SoilReaction <- apply(sr, 1, function(x) x[!is.na(x)][1])
sr <- sr[-c(1:3)]
sr$Genus <- rownames(sr)

sr <- sr[complete.cases(sr), ]
table(sr)
table(sr$SoilReaction)

## For moss composition: quite an arbitrary partition
# range(mossMeta$MossHabitat)
# mossMeta$MossHabitatC <- mossMeta$MossHabitat
# mossMeta$MossHabitatC[mossMeta$MossHabitatC < -0.033865] <- "Raised bog"
# mossMeta$MossHabitatC[mossMeta$MossHabitatC < 1.356875] <- "Transient bog"
# mossMeta$MossHabitatC[mossMeta$MossHabitatC < 2.747615] <- "Fen"
# mossMeta$MossHabitatC
#
#
# indval <- multipatt(otus.tax.is, mossMeta$MossHabitatC,
#                     control = how(nperm = 999), duleg = TRUE)
# summary(indval)


## For canopy cover
range(mossMeta$CanopyCover)
mossMeta$CanopyCoverC <- mossMeta$CanopyCover
mossMeta$CanopyCoverC[mossMeta$CanopyCoverC < 40] <- "Low cover"
mossMeta$CanopyCoverC[mossMeta$CanopyCoverC < 70] <- "Intermediate cover"
mossMeta$CanopyCoverC[mossMeta$CanopyCoverC < 90] <- "High cover"
mossMeta$CanopyCoverC


indval <- multipatt(otus.tax.is, mossMeta$CanopyCoverC,
                    control = how(nperm = 999), duleg = TRUE)
summary(indval)


cc <- data.frame(indval$str)
cc[cc < 0.7] <- NA
cc$High.cover[!is.na(cc$High.cover)] <- "High Cover"
cc$Intermediate.cover[!is.na(cc$Intermediate.cover)] <- "Intermediate cover"
cc$Low.cover[!is.na(cc$Low.cover)] <- "Low cover"


cc$CanopyCover <- apply(cc, 1, function(x) x[!is.na(x)][1])
cc <- cc[-c(1:3)]
cc$Genus <- rownames(cc)


## Merge them
indic <- cbind(sr, cc, alt)
indic <- indic[, -c(2, 4)]


################################################################################
### Co-occurrence network with rank correlation
## Presence/absence and taxa that occur at least in 8 (so that they can all
## occurr within a site) samples and then, subset
otus <- data.frame(t(otu_table(moss)))

otus.pa <- data.frame(t(otu_table(moss)))
otus.pa[otus.pa > 0] <- 1
otus.pa <- otus.pa[rowSums(otus.pa) > 7, ]
summary(rowSums(otus.pa))
subset.inc <- rownames(otus.pa)

otus <- otus[rownames(otus) %in% subset.inc, ]


## Total abundance
otus.tax$Abundance <- rowSums(otus.tax[, 8:47])


## Correlation analysis based on spearman's co-efficient
otus.dist <- rcorr(t(otus), type = "spearman")
otus.cor <- otus.dist$r
otus.cor.p <- otus.dist$P


## Multiple testing correction using Benjamini-Hochberg standard FDR correction
otus.cor.p <- p.adjust(otus.cor.p, method = "BH")


## Positive and netagive cooccurence at given coefficient and p-value cutoff
cor.cutoff <- 0.4
p.cutoff <- 0.001


otus.cor[which(otus.cor >= (- cor.cutoff) & otus.cor <= cor.cutoff)] <- 0
otus.cor[which(otus.cor.p > p.cutoff)] <- 0


## Delete rows and columns with sum = 0
otus.cor <- otus.cor[which(rowSums(otus.cor) != 1), ]
otus.cor <- otus.cor[, which(colSums(otus.cor) != 0)]


### Create a graph
g3 <- graph.adjacency(otus.cor, weight = TRUE, mode = "undirected")
g3 <- simplify(g3) # remove duplicate and loop edges


### Add taxonomy
## Get taxa in network and merge with taxonomy
target <- data.frame(ASV = V(g3)$name,
                     to_sort = seq(1:length(V(g3)$name)))


otu.target <- merge(target, otus.tax, by = "ASV")
otu.target <- otu.target[order(otu.target$to_sort), ]


## Add indicator species
otu.target <- merge(otu.target, sr, by = "Genus", all.x = TRUE)
otu.target <- otu.target[order(otu.target$to_sort), ]


##
# otu.target$DomainInt <- as.numeric(factor(otu.target$Domain))
# otu.target[, c(3, 50)]


E(g3)
V(g3)
V(g3)$name
V(g3)$degree <- degree(g3) # connectivity


### Add edge attributes
## Pairwise correlations as weight
otus.cor[lower.tri(otus.cor, diag = TRUE)] <- NA
otus.cor.m <- melt(otus.cor, value.name = "rho")
otus.cor.m$rho[otus.cor.m$rho == 0] <- NA
otus.cor.m <- otus.cor.m[complete.cases(otus.cor.m), ]


E(g3)$width <- otus.cor.m$rho


# E(g3)$color <- ifelse(E(g3)$width > 0, '#009E73', '#D55E00')
E(g3)$scale <- ifelse(E(g3)$width < 0, E(g3)$width * -1, E(g3)$width)
E(g3)$sign <- ifelse(E(g3)$width < 0, "Negative", "Positive")


### Add vertex attributes
g3 <- set_vertex_attr(g3, "ASV", value = otu.target$ASV)
V(g3)$ASV

g3 <- set_vertex_attr(g3, "Genus", value = otu.target$Genus)
V(g3)$Genus
V(g3)$name <- V(g3)$Genus

g3 <- set_vertex_attr(g3, "Domain", value = otu.target$Domain)
V(g3)$Domain


g3 <- set_vertex_attr(g3, "Abundance", value = sqrt(otu.target$Abundance))
V(g3)$Abundance

V(g3)$size <- V(g3)$Abundance


### Louvain community detection
## Set negative edges to zero
E(g3)$width <- ifelse(E(g3)$width < 0, 0, E(g3)$width)
E(g3)$width

g3.louv <- cluster_louvain(g3, weights = E(g3)$width)
length(g3.louv)
sizes(g3.louv)
table(sizes(g3.louv))
membership(g3.louv)
modularity(g3.louv)


V(g3)$louvain.member <- g3.louv$membership


### Indicator species
# V(g3)$CanopyCover <- otu.target$CanopyCover
V(g3)$SoilReaction <- otu.target$SoilReaction
# V(g3)$Altitude <- otu.target$Altitude



## Export graph
# write_graph(g3, "~/Documents/Seafile/_FNS_2016-NicheDynamCoexist/SarraceniaMicrobiomeProject/Mosses/Cytoscape/Moss_Co-occurrence_8.graphml",
            # format = "graphml")
write_graph(g3, "~/Documents/Seafile/_FNS_2016-NicheDynamCoexist/SarraceniaMicrobiomeProject/Mosses/Cytoscape/Moss_Co-occurrence_8_pub.graphml",
            format = "graphml")


## Network statistics
vcount(g3)
ecount(g3)


table(V(g3)$Domain)
table(degree(g3))
mean(degree(g3))
summary(degree(g3))


## Negative and positive edges
summary(otus.cor.m$rho > 0)


## Subnetworks
components(g3)


################################################################################
################################################################################
################################################################################
################################################################################
