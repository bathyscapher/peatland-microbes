################################################################################
################################################################################
################################################################################
################################################################################
### SMP Moss Microbiome
### Author: korn@cumulonimbus.at University of Fribourg 2021
################################################################################


library("phyloseq")
library("ggplot2")
theme_set(theme_bw(base_size = 20) +
            theme(rect = element_rect(fill = "transparent")))
library("gridExtra")
library("vegan")
library("gplots")
library("reshape2")


rm(list = ls())


################################################################################
### Moss data
moss.pa <- readRDS("rds/Moss_prok.a.RDS")
moss.ea <- readRDS("rds/Moss_euk.a.RDS")
moss <- readRDS("rds/Moss_smp.RDS")



sample_data(moss.pa)$Site <- ordered(sample_data(moss.pa)$Site,
                                     levels = c("LT", "LM", "LV", "LE", "CB"))
sample_data(moss.ea)$Site <- ordered(sample_data(moss.ea)$Site,
                                     levels = c("LT", "LM", "LV", "LE", "CB"))


################################################################################
### Convert ps to df
tax.m <- as.data.frame(moss@tax_table@.Data)
otu.m <- as.data.frame(t(otu_table(moss)))

tax.otu.m <- merge(tax.m, otu.m, by = 0, all = TRUE)
rownames(tax.otu.m) <- tax.otu.m$Row.names
names(tax.otu.m)[1] <- "Taxa"


rm(tax.m, otu.m)


################################################################################
### Plot taxa
tax.otu.m[, -c(1:7)] <- log1p(tax.otu.m[, -c(1:7)])


tax.otu.m.m <- reshape2::melt(tax.otu.m, id.vars = c("Taxa", "Domain",
                                                     "Phylum", "Class",
                                                     "Order", "Family",
                                                     "Genus"),
                              variable.name = "FullID", value.name = "Count")
tax.otu.m.m$Site <- substr(tax.otu.m.m$FullID, 1, 2)
tax.otu.m.m$Site <- ordered(tax.otu.m.m$Site,
                            levels = c("LT", "LM", "LV", "LE", "CB"))
tax.otu.m.m$Domain[tax.otu.m.m$Domain == "Archaea"] <- "A."


## Reorder taxonomically
tax.otu.m.m[, c(2:7)] <- lapply(tax.otu.m.m[, c(2:7)], as.factor)


tax <- c("Euryarchaeota", "Nanoarchaeota", ## Archaea
         ## Proteobacteria
         "Bdellovibrionota", "Campilobacterota", "Desulfobacterota",
         "Proteobacteria",
         "Myxococcota", ## Deltaproteobacteria
         "Planctomycetota", "Verrucomicrobiota", ## PVC group
         "Gemmatimonadota", ## FCB group
         ## Singles
         "Acidobacteriota", "Bacteroidota", "Chloroflexi", "Dependentiae",
         ## Terrabacteria
         "Actinobacteriota", "Armatimonadota", "Cyanobacteria", "Firmicutes",
         ## ?
         "WPS-2", "FCPU426", "Patescibacteria", "MBNT15",
         "Elusimicrobiota", "Spirochaetota", "NB1-j", "Abditibacteriota",
         "RCP2-54", "Latescibacterota", "Sumerlaeota", "Fibrobacterota", "WS4",
         "Deinococcota", "SAR324_clade(Marine_group_B)",
         ### Eukaryotes ###
         ## Amoeba
         "Amoebozoa_ph", "Dictyostelia", "Gracilipodida", "Schizoplasmodiida",
         "Protosteliida",
         ## SAR:Heterokonta/Stramenopiles
         "Bicosoecida", "Ochrophyta_ph", "Diatomea", "Hyphochytriomycetes",
         "MAST-3", "MAST-12", "Peronosporomycetes",
         "Cercozoa", ## SAR: Rhizaria
         ## Harosa:Alveolata
         "Ciliophora", "Dinoflagellata", "Protalveolata",
         "Euglenozoa", "Heterolobosea", "Parabasalia", ## Excavata
         ## Algae
         "Chlorophyta_ph", "Cryptophyceae_ph", "Phragmoplastophyta",
         "Klebsormidiophyceae",
         "Incertae_Sedis",
         "Nucleariidae_and_Fonticula_group", ## Holomycota
         ## Fungi
         "Ascomycota", "Basidiomycota", "Blastocladiomycota", "Chytridiomycota",
         "Cryptomycota", "LKM15", "Zoopagomycota",
         "Holozoa_ph", ## Choanozoa
         "Aphelidea", ## Sister group of fungi
         ## Metazoa
         "Annelida", "Platyhelminthes", "Nematozoa", "Rotifera", "Tardigrada",
         "Gastrotricha", "Arthropoda")


## Pronounce small values even more
tax.otu.m.m$Count[tax.otu.m.m$Count < 0.01 & tax.otu.m.m$Count > 0] <- 0.01
summary(tax.otu.m.m$Count)


## Plot
ggplot(tax.otu.m.m, aes(x = factor(Phylum, level = rev(tax)), y = Count,
                        fill = Site)) +
  geom_bar(stat = "identity") +
  facet_grid(Domain ~ Site, scales = "free", space = "free_y") +
  ylab(expression(paste(ln[e], "(Abundance + 1) [%]"))) +
  xlab("") +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1))
# ggsave("Moss_Phyla.pdf", width = 8.27, height = 11.69)


### Inspect taxa and frequencies
df <- tax.otu.m[tax.otu.m$Phylum == "Ascomycota", ]
summary(colSums(df[, -c(1:7)]) == 0)

unique(df$Genus)


################################################################################
### Alpha-diversity
## Prokaryotes
otu.p <- data.frame(
  otu_table(transform_sample_counts(
    moss.pa, function(abund) {1 * (abund > 0)})
  ))
otu.p$alpha <- rowSums(otu.p)
otu.p <- otu.p[c("alpha")]
otu.p$FullID <- rownames(otu.p)
otu.p$Domain <- "Prokaryotes"


## Eukaryotes
otu.e <- data.frame(otu_table(
  transform_sample_counts(moss.ea,
                          function(abund) {1 * (abund > 0)})
))
otu.e$alpha <- rowSums(otu.e)
otu.e <- otu.e[c("alpha")]
otu.e$FullID <- rownames(otu.e)
otu.e$Domain <- "Eukaryotes"


## Merge
otus <- rbind(otu.p, otu.e)
rm(otu.p, otu.e)


## Merge with metadata
meta.m <- data.frame(sample_data(moss.pa))


otus <- merge(otus, meta.m, by = "FullID", all = TRUE)
# otus$FullID == rownames(meta.m)
otus$Domain <- ordered(otus$Domain, levels = c("Prokaryotes", "Eukaryotes"))


## Alpha-diversity
ggplot(otus, aes(x = Site, y = alpha)) +
  geom_boxplot(color = "black", size = 0.2, outlier.shape = NA) +
  geom_jitter(aes(fill = Site), height = 0, width = 0.3, size = 2, shape = 21) +
  stat_summary(aes(group = Site), fun = mean,
               colour = "black", geom = "point", fill = "gray",
               shape = 8, size = 3.5, show.legend = FALSE) +
  facet_grid(Domain ~ Site, scales = "free") +
  xlab("") +
  ylab(expression(alpha-diversity)) +
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_blank())
# ggsave("Moss_Alpha.pdf", width = 11.69, height = 5)


## Mean alpha-diversity by site and domain
aggregate(alpha ~ Site + Domain, data = otus, mean)


## Hypothesis test for difference of mean alpha between fen and raised bog sites
### Test for normality
# pdf("Moss_QQNormAlpha.pdf", width = 11.69, height = 6)
par(mfrow = c(1, 2))
qqnorm(otus$alpha[otus$Domain == "Prokaryotes"], main = "Prokaryotes")
qqline(otus$alpha[otus$Domain == "Prokaryotes"], lty = 2)
qqnorm(otus$alpha[otus$Domain == "Eukaryotes"], main = "Eukaryotes", ylab = "")
qqline(otus$alpha[otus$Domain == "Eukaryotes"], lty = 2)
# dev.off()


## One-way ANOVA followed by Tukey's HSD
### Prokaryotes
res.aov <- aov(alpha ~ Site, data = otus[otus$Domain == "Prokaryotes", ])
summary(res.aov)

TukeyHSD(res.aov)


### Eukaryotes
res.aov <- aov(alpha ~ Site, data = otus[otus$Domain == "Eukaryotes", ])
summary(res.aov)

TukeyHSD(res.aov)


################################################################################
### Ordination
## Prokaryotes
set.seed(128252)

moss.pa.log <- transform_sample_counts(moss.pa,
                                       function(otu) {log1p(otu)})
moss.nmds <- ordinate(moss.pa.log,
                      method = "NMDS", distance = "bray", k = 2,
                      trymax = 50)
moss.nmds


nmds.prok <- plot_ordination(moss.pa.log, moss.nmds, shape = "Sector",
                color = "Site", title = NULL, label = "FullID") +
  stat_ellipse(aes(group = Site), type = "t", linetype = 2, size = 0.2) +
  geom_point(size = 3) +
  coord_fixed(ratio = 1) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.box = "vertical", legend.margin = margin()) +
  xlim(-1, 0.9) +
  ylim(-1, 0.5) +
  xlab("nMDS1") +
  ylab("nMDS2")


## Eukaryotes
set.seed(128252)

moss.ea.log <- transform_sample_counts(moss.ea, function(otu) {log1p(otu)})
moss.nmds <- ordinate(moss.ea.log, method = "NMDS", distance = "bray", k = 2,
                      trymax = 50)
moss.nmds


nmds.euk <- plot_ordination(moss.ea.log, moss.nmds, shape = "Sector",
                             color = "Site", title = NULL, label = "FullID") +
  stat_ellipse(aes(group = Site), type = "t", linetype = 2, size = 0.2) +
  geom_point(size = 3) +
  coord_fixed(ratio = 1) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.box = "vertical", legend.margin = margin()) +
  xlim(-1, 0.9) +
  ylim(-1, 0.5) +
  xlab("nMDS1") +
  ylab("")


## Arrange plots side by side and export
nmds.both <- arrangeGrob(nmds.prok, nmds.euk, nrow = 1)


 # ggsave("Moss_nMDS.pdf", nmds.both, width = 11.69, height = 8.27)


################################################################################
################################################################################
### Prepare metadata for SEM
## Read metadata
mossMeta <- read.table("csv/MossesMetadata.csv", sep = "\t", header = TRUE)


mossMeta$FullID <- gsub("_", "", mossMeta$FullID)
mossMeta$Sector <- substr(mossMeta$FullID, 3, 3)


## Ellenberg values
ev <- read.table("csv/EVMeans.csv", sep = "\t", header = TRUE)


## Merge
mossMeta <- merge(mossMeta, ev, by = "Site")


## Canopy cover
cover <- read.table("csv/CanopyCover.csv", header = TRUE, sep = "\t",
                    fill = FALSE, check.names = TRUE)
leaves <- read.table("csv/LeafSamples.csv",
                     header = TRUE, sep = ",", fill = FALSE)


cover <- merge(leaves, cover[, c(1:2, 8)],
               by.x = c("Site", 'PlantOld'), by.y = c("Site", 'PlantOld'))
cover$Sector <- substr(cover$FullID, 4, 4)


## Calculate mean canopy cover sector-wise
cover.mean <- aggregate(cover$CanopyCover, list(cover$Site, cover$Sector), mean,
                        na.rm = TRUE)
names(cover.mean) <- c("Site", "Sector", "CanopyCover")


mossMeta <- merge(mossMeta, cover.mean, by = c("Site", "Sector"))


## Dry weight of sampled mosses as proxy of habitat size and sampling intensity
weight <- read.table("csv/MossSamples.csv", header = TRUE, sep = ",",
                     fill = FALSE, check.names = TRUE)
weight$FullID <- gsub("_", "", weight$FullID)


mossMeta <- merge(mossMeta, weight[, c(1, 5)], by = "FullID")


################################################################################
### Reduce coordinates to 1D with MDS
rownames(mossMeta) <- mossMeta$FullID


xy.mds <- cmdscale(dist(mossMeta[, c(14:15)]), eig = TRUE, k = 1)
xy.mds$GOF # should be > 0.8
distance <- as.data.frame(xy.mds$points)
colnames(distance) <- "Distance"
distance$FullID <- row.names(distance)


xy.mds.df <- as.data.frame(xy.mds$points)
colnames(xy.mds.df) <- "Points"
xy.mds.df$Plant <- rownames(xy.mds.df)


ggplot(xy.mds.df) +
  geom_text(aes(x = 0, y = Points, label = Plant), na.rm = TRUE) +
  scale_x_continuous(breaks = seq(-1, 1, 1)) +
  xlab("") +
  ylab("MDS 1")
# ggsave("Moss_Distance_1D.pdf", width = 3, height = 6)


mossMeta <- merge(mossMeta, distance, by = "FullID")


### Moss as habitat
det <- read.table("csv/MossComposition.csv", sep = "\t", header = TRUE,
                  row.names = 1)

det[det == "x"] <- 1
det[det == ""] <- 0
det[] <- lapply(det, as.integer)
det <- log1p(t(det))

set.seed(128252)

det.nmds <- metaMDS(det, distance = "bray", k = 1, trymax = 75,
                    autotransform = FALSE)


det.nmds
stressplot(det.nmds)


# pdf("Moss_DetnMDS.pdf", height = 8.27, width = 8.27)
margin(5, 5, 10, 5)
plot(det.nmds)
# dev.off()


det.nmds.sc <- data.frame(scores(det.nmds))
names(det.nmds.sc) <- "MossComposition"
det.nmds.sc$Site <- substr(rownames(det.nmds.sc), 1, 2)
det.nmds.sc$Sector <- substr(rownames(det.nmds.sc), 3, 3)


mossMeta <- merge(mossMeta, det.nmds.sc, by = c("Site", "Sector"))


################################################################################
### Microbiomes
## Prokaryotes
moss.prok.mb <- as.data.frame(otu_table(moss.pa))
moss.prok.mb <- log1p(moss.prok.mb)

set.seed(128252)


moss.prok.nmds <- vegan:::metaMDS(moss.prok.mb, distance = "bray",
                                  trymax = 50, k = 2, autotransform = FALSE)
moss.prok.nmds

stressplot(moss.prok.nmds)


moss.prok.nmds.sc <- as.data.frame(scores(moss.prok.nmds)[, 1])
names(moss.prok.nmds.sc) <- "mb.moss.prok"
moss.prok.nmds.sc$FullID <- rownames(moss.prok.nmds.sc)


mossMeta <- merge(mossMeta, moss.prok.nmds.sc, by = "FullID", all = TRUE)


## Eukaryotes
moss.euk.mb <- as.data.frame(otu_table(moss.ea))
moss.euk.mb <- log1p(moss.euk.mb)

set.seed(128252)

moss.euk.nmds <- vegan:::metaMDS(moss.euk.mb, distance = "bray",
                                 trymax = 50, k = 2, autotransform = FALSE)
moss.euk.nmds

stressplot(moss.euk.nmds)


moss.euk.nmds.sc <- as.data.frame(scores(moss.euk.nmds)[, 1])
names(moss.euk.nmds.sc) <- "mb.moss.euk"
moss.euk.nmds.sc$FullID <- rownames(moss.euk.nmds.sc)


mossMeta <- merge(mossMeta, moss.euk.nmds.sc, by = "FullID", all = TRUE)


## Both
moss.mb <- as.data.frame(otu_table(moss))
moss.mb <- log1p(moss.mb)

set.seed(128252)

moss.nmds <- vegan:::metaMDS(moss.mb, distance = "bray",
                             trymax = 50, k = 2, autotransform = FALSE)
moss.nmds

stressplot(moss.nmds)


moss.nmds.sc <- as.data.frame(scores(moss.nmds)[, 1])
names(moss.nmds.sc) <- "mb.moss"
moss.nmds.sc$FullID <- rownames(moss.nmds.sc)


mossMeta <- merge(mossMeta, moss.nmds.sc, by = "FullID", all = TRUE)


### Export data
# write.table(mossMeta, "csv/Mosses_Metadata_unscaled.csv", sep = "\t",
#             row.names = FALSE)


################################################################################
### Scale numeric variables that will be used in the SEM
str(mossMeta[, -c(1:9, 11:15)])

par(mfrow = c(2, 1))
boxplot(mossMeta[, -c(1:9, 11:15)])
abline(h = 0, lty = 2, col = "gray")


mossMeta[, -c(1:9, 11:15)] <- apply(mossMeta[, -c(1:9, 11:15)], 2, scale)


boxplot(mossMeta[, -c(1:9, 11:15)])
abline(h = 0, lty = 2, col = "gray")
dev.off()


### Export data
# write.table(mossMeta, "csv/Mosses_Metadata_scaled.csv", sep = "\t",
#             row.names = FALSE)


################################################################################
### Correlation and covariance
corPalette <- colorRampPalette(c("#D55E00", "white",  "#56B4E9"))(n = 19)


## Sort covariance data.frame
covs <- data.frame(cov(mossMeta[, c(10, 16:17, 19:21, 23:26)],
                       use = "complete"))
covs <- covs[c("Altitude", "Light", "Temperature", "Humidity", "Soil.reaction",
               "Nitrogen", "CanopyCover", "Moss_g", "Distance",
               "MossComposition")]
covs <- covs[colnames(covs), ]


pdf("Moss_Metadata_Covariance.pdf", height = 8.27, width = 8.27)
heatmap.2(as.matrix(covs), dendrogram = "none", key = TRUE, key.title = "",
          cellnote = round(covs, digits = 2), notecol = "black", notecex = 0.6,
          trace = "none", distfun = function(x) {x}, cexCol = 1.5, cexRow = 1.5,
          symm = TRUE, margins = c(16, 16), col = corPalette,
          tracecol = "black", Rowv = FALSE, Colv = FALSE)
dev.off()


################################################################################
theme_set(theme_bw(base_size = 9) +
            theme(rect = element_rect(fill = "transparent")))


## Alpha-diversity vs environmental variables
mossMeta <- read.table("csv/Mosses_Metadata_unscaled.csv", sep = "\t",
                       header = TRUE)
mossMeta$Site <- ordered(mossMeta$Site,
                         levels = c("LT", "LM", "LV", "LE", "CB"))

# names(otus)
# names(mossMeta)
mossMeta <- merge(mossMeta, otus[, c(1:3)])
names(mossMeta)

mossMeta.s <- mossMeta[, c(2, 10, 16:17, 19:21, 23:31)]
names(mossMeta.s)


mossMeta.s.m <- melt(mossMeta.s, id.vars = c("Domain", "alpha", "Site"),
                     variable.name = "Variable", value.name = "Value")
names(mossMeta.s.m)


ggplot(mossMeta.s.m, aes(x = Value, y = alpha, color = Site)) +
  geom_smooth(aes(group = 1), color = "gray",
              method = "lm", formula = y ~ x, se = FALSE) +
  geom_point(pch = 21) +
  facet_grid(Domain ~ Variable, scales = "free", space = "free_y") +
  ylab(expression(alpha-diversity)) +
  xlab("m | - | - | - | - | - | % | g | m | - | - | - | -") +
  scale_color_viridis_d() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave("Moss_AlphaVsEnvironment.pdf", height = 6, width = 11.69)


################################################################################
################################################################################
################################################################################
################################################################################
