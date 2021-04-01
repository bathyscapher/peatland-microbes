################################################################################
################################################################################
################################################################################
################################################################################
### Mosses: Match Ellenberg's indicator values with a species list
### Authors: Rachel Korn
### korn@cumulonimbus.at University of Fribourg 2020
################################################################################


library("plyr")
library("reshape2")
library("ggplot2")
theme_set(theme_bw(base_size = 20) +
            theme(rect = element_rect(fill = "transparent")))


rm(list = ls())


setwd("~/EVA/EllenbergsValuesAutomaton/")


################################################################################
### Read Ellenberg's values
ev <- lapply(list.files(pattern = "Ellenberg.*.csv"), read.delim, sep = "\t")


## Convert list of data.frames to data.frame
ev <- rbind.fill(ev)


## Resort columns
ev <- ev[c("Name", "L", "T", "K", "F", "R", "N", "S", "LF", "LF_B", "SUB")]


## Reformat columns
ev[, 2:11] <- as.data.frame(apply(ev[, 2:11], 2,
                                  function(x) gsub("\\?|x", NA, x)))
ev[, 2:11] <- as.data.frame(apply(ev[, 2:11], 2,
                                  function(x) gsub("-|B|b|~|=|\\(|)", "", x)))
ev[, c(2:8)] <- as.data.frame(apply(ev[, c(2:8)], 2, as.integer))


################################################################################
### Read vegetation list
setwd("~/peatland-microbes/")


vascu <- read.table("csv/VascularPlants.csv", sep = "\t", header = TRUE,
                    quote = "")

## Delete certain entries (trees, non-bog species)
vascu <- vascu[!grepl("Lycopus europaeus L.|Angelica sylvestris L.|Pinus sylvestris L.|Epilobium hirsutum L.|Anthoxanthum odoratum L.|Betula pubescens Ehrh.|Festuca rubra aggr.|Menyanthes trifoliata L.|Melampyrum pratense L.|Picea abies \\(L.\\) H. Karst.|Pinus mugo subsp. uncinata \\(DC.\\) Domin|Filipendula ulmaria \\(L.\\) Maxim.|Calluna vulgaris \\(L.\\) Hull",
                    vascu$Taxon), ]


vascu <- vascu[!(vascu$Taxon == "Potentilla erecta (L.) Raeusch." &
                   vascu$Site == "LM" |
                   vascu$Taxon == "Phragmites australis (Cav.) Steud." &
                   vascu$Site == "LM" |
                   vascu$Taxon == "Silene flos-cuculi (L.) Clairv." &
                   vascu$Site == "LM" |
                   vascu$Taxon == "Scirpus sylvaticus L." & vascu$Site == "LM"),
               ]


bryo <- read.table("csv/Bryophyta.csv", sep = "\t", header = TRUE, quote = "")


bryo2 <- read.table("csv/Bryophyta_2018.csv",
                   sep = "\t", header = TRUE, quote = "")


plants <- rbind(vascu, bryo, bryo2)
# sort(names(bryo)) == sort(names(vascu))
rm(vascu, bryo, bryo2)


## Remove duplicates on site-level
plants <- plants[!duplicated(plants[ , c("Taxon", "Site")]), ]


### Format species names
## Replace accented letters
plants$SpeciesShort <- chartr("áéíóúňäöü", "aeiounaou", plants$Taxon)


## Delete abbreviated author names
plants$SpeciesShort <- gsub("[A-Z]{1,2}[a-z]*\\.", "", plants$SpeciesShort,
                            perl = TRUE)


## Delete f.
plants$SpeciesShort <- gsub("f\\.", "", plants$SpeciesShort,
                            perl = TRUE)


## Delete author names
plants$SpeciesShort <- gsub(" [A-Z]{1,2}[a-z]*", "",
                            plants$SpeciesShort, perl = TRUE)


## Delete ampersands
plants$SpeciesShort <- gsub("\\&", "", plants$SpeciesShort, perl = TRUE)


## Delete author names in brackets
plants$SpeciesShort <- gsub("\\([A-Z]{1,2}[a-z]* *\\)", "",
                            plants$SpeciesShort, perl = TRUE)


## Delete brackets and trailing white space
plants$SpeciesShort <- gsub(" *$", "", gsub("\\(\\)", "", plants$SpeciesShort))


## Reduce multiple spaces
plants$SpeciesShort <- gsub(" {2,}", " ", plants$SpeciesShort)


# levels(as.factor(plants$SpeciesShort))


## Get column number + 1
cols <- dim(plants)[2] + 1


## Add empty columns for Ellenberg's values to species list
plants$CheckPoint <- NA
plants[names(ev[2:11])] <- NA


################################################################################
### Fuzzy match by species. Set cost to conserve names (no insertions and
### substitutions, deletions allowed)
gottem <- lapply(plants$SpeciesShort, agrep, x = ev$Name, value = FALSE,
                 max.distance = c(all = 1,
                                  deletions = 2,
                                  insertions = 2,
                                  substitutions = 0))


################################################################################
### Fill Ellenberg's values into plants
EVA <- function(matches, index){
  if(length(matches) == 1)
    {plants[index, cols:(cols + 10)] <- ev[matches, ]}
  else
  {plants[index, cols:(cols + 10)] <- c(paste(length(matches), "match(es)"),
                                          rep(NA, 10))}
  }


## Update target data.frame
plants[, cols:(cols + 10)] <- do.call(rbind,
                                      Map(EVA, gottem, seq_along(gottem)))


## See mismatches
mismatches <- plants[grep("match", plants$CheckPoint), ]


levels(as.factor(mismatches$SpeciesShort))


## Revise taxonomic changes = 0 matches
# Carex viridula = Carex flava agg. ???
# Sarracenia purpurea = n.e.
# Sphagnum inundatum = ???
# Sphagnum recurvum aggr. = Sphagnum fallax (Klinggr.) Klinggr. (mucronatum, recurvum) ???


plants[plants$SpeciesShort == "Drosera ×obovata", cols:(cols + 10)] <-
  ev[ev$Name == "Drosera x obovata", ]


plants[plants$SpeciesShort == "Cephalozia connivens subsp. connivens",
       cols:(cols + 10)] <-
  ev[ev$Name == "Cephalozia connivens (Dicks.) Lindb.", ]

plants[plants$SpeciesShort == "Odontoschisma fluitans", cols:(cols + 10)] <-
  ev[ev$Name == "Cladopodiella fluitans (Nees) Buch", ]

plants[plants$SpeciesShort == "Pinus mugo subsp. uncinata", cols:(cols + 10)] <-
  ev[ev$Name == "Pinus mugo rotundata", ]

plants[plants$SpeciesShort == "Polytrichum formosum aggr.", cols:(cols + 10)] <-
  ev[ev$Name == "Polytrichum formosum Hedw.", ]

plants[plants$SpeciesShort == "Pseudoscleropodium purum", cols:(cols + 10)] <-
  ev[ev$Name == "Scleropodium purum (Hedw.) Limpr.", ]

plants[plants$SpeciesShort == "Ranunculus acris subsp. friesianus",
       cols:(cols + 10)] <-
  ev[ev$Name == "Ranunculus acris agg.", ]

plants[plants$SpeciesShort == "Silene flos-cuculi", cols:(cols + 10)] <-
  ev[ev$Name == "Lychnis flos-cuculi", ]

plants[plants$SpeciesShort == "Sphagnum magellanicum aggr.",
       cols:(cols + 10)] <-
  ev[ev$Name == "Sphagnum magellanicum Brid. (medium)", ]

plants[plants$SpeciesShort == "Sphagnum majus subsp. majus",
       cols:(cols + 10)] <-
  ev[ev$Name == "Sphagnum majus (Russ.) C. Jens. (dusenii)", ]

plants[plants$SpeciesShort == "Sphagnum palustre subsp. palustre",
       cols:(cols + 10)] <-
  ev[ev$Name == "Sphagnum palustre L. (cymbifolium)", ]

plants[plants$SpeciesShort == "Sphagnum subnitens subsp. subnitens",
       cols:(cols + 10)] <-
  ev[ev$Name == "Sphagnum subnitens Russ. & Warnst. ex Warnst. (plumosum)", ]

plants[plants$SpeciesShort == "Straminergon stramineum", cols:(cols + 10)] <-
  ev[ev$Name == "Calliergon stramineum (Brid.) Kindb.", ]

plants[plants$SpeciesShort == "Trichophorum cespitosum subsp. cespitosum",
       cols:(cols + 10)] <-
  ev[ev$Name == "Trichophorum cespitosum cespitosum (austr.)", ]

plants[plants$SpeciesShort == "Ulota crispa s.l.", cols:(cols + 10)] <-
  ev[ev$Name == "Ulota crispa (Hedw.) Brid. var. crispa", ]

plants[plants$SpeciesShort == "Vaccinium microcarpum", cols:(cols + 10)] <-
  ev[ev$Name == "Vaccinium oxycoccos microcarpum", ]


## Multiple matches
plants[plants$SpeciesShort == "Atrichum undulatum", cols:(cols + 10)] <-
  ev[ev$Name == "Atrichum undulatum (Hedw.) P. Beauv. var. undulatum", ]

plants[plants$SpeciesShort == "Aulacomnium palustre", cols:(cols + 10)] <-
  ev[ev$Name == "Aulacomnium palustre (Hedw.) Schwaegr. var. palustre", ]

plants[plants$SpeciesShort == "Betula pubescens", cols:(cols + 10)] <-
  ev[ev$Name == "Betula pubescens", ]

plants[plants$SpeciesShort == "Carex canescens", cols:(cols + 10)] <-
  ev[ev$Name == "Carex canescens canescens", ]

plants[plants$SpeciesShort == "Cephalozia bicuspidata", cols:(cols + 10)] <-
  ev[ev$Name == "Cephalozia bicuspidata (L.) Dum. var. bicuspidata", ]

plants[plants$SpeciesShort == "Dactylorhiza maculata", cols:(cols + 10)] <-
  ev[ev$Name == "Dactylorhiza maculata", ]

plants[plants$SpeciesShort == "Dactylorhiza majalis", cols:(cols + 10)] <-
  ev[ev$Name == "Dactylorhiza majalis (O. latifolia)", ]

plants[plants$SpeciesShort == "Galium palustre", cols:(cols + 10)] <-
  ev[ev$Name == "Galium palustre palustre", ]

plants[plants$SpeciesShort == "Melampyrum pratense", cols:(cols + 10)] <-
  ev[ev$Name == "Melampyrum pratense ssp. pratense", ]

plants[plants$SpeciesShort == "Molinia caerulea", cols:(cols + 10)] <-
  ev[ev$Name == "Molinia caerulea caerulea", ]

plants[plants$SpeciesShort == "Potentilla erecta", cols:(cols + 10)] <-
  ev[ev$Name == "Potentilla erecta (Tormentilla erecta)", ]

plants[plants$SpeciesShort == "Ranunculus flammula", cols:(cols + 10)] <-
  ev[ev$Name == "Ranunculus flammula flammula", ]

plants[plants$SpeciesShort == "Salix repens", cols:(cols + 10)] <-
  ev[ev$Name == "Salix repens", ]

plants[plants$SpeciesShort == "Sphagnum capillifolium", cols:(cols + 10)] <-
  ev[ev$Name == "Sphagnum capillifolium (L.) Hedw. var. capillifolium (nemoreum)", ]

plants[plants$SpeciesShort == "Trichophorum cespitosum", cols:(cols + 10)] <-
  ev[ev$Name == "Trichophorum cespitosum cespitosum (austr.)", ]

plants[plants$SpeciesShort == "Vaccinium oxycoccos", cols:(cols + 10)] <-
  ev[ev$Name == "Vaccinium oxycoccos oxycoccos (Oxyc. palustr.)", ]


## See remaining mismatches
mismatches <- plants[grep("match", plants$CheckPoint), ]


levels(as.factor(mismatches$SpeciesShort))


## Summarize by site
plants[, c(12:18)] <- as.data.frame(apply(plants[, c(12:18)], 2, as.integer))


colnames(plants)[12:21] <- c("Light", "Temperature", "Continentality",
                             "Humidity", "Soil reaction", "Nitrogen", "Salt",
                             "Lifeform", "Leaves", "Substrate")


plants.m <- melt(plants[, c(8, 12:18)], variable.name = "EV",
                 value.name = "Value")

plants.m$EV <- as.factor(plants.m$EV)



### Plot
plants.m$Site <- ordered(plants.m$Site,
                         levels = c("CB", "LE", "LV", "LT", "LM"))


## Exclude Continentality and Salt
plants.m <- plants.m[plants.m$EV != "Continentality" & plants.m$EV != "Salt", ]


ggplot(plants.m, aes(x = Site, y = Value)) +
  geom_boxplot(size = 0.2, outlier.shape = NA) +
  geom_jitter(aes(color = Site), na.rm = TRUE, pch = 21,
              height = 0, width = 0.3, size = 2) +
  facet_grid( ~ EV, scales = "free") +
  theme(legend.position = "top", legend.direction = "horizontal",
        axis.text.x = element_blank()) +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  xlab("") +
  ylab("Ellenberg's indicator value")
# ggsave("Moss_EV.pdf", width = 11.69, height = 5)


### Summarize EV by site
ev.mean <- aggregate(plants[, 12:18], list(plants$Site), mean, na.rm = TRUE)
names(ev.mean)[1] <- "Site"


write.table(ev.mean, "csv/Moss_EVMeans.csv", sep = "\t", row.names = FALSE)


################################################################################
################################################################################
################################################################################
################################################################################

