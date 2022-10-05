################################################################################
################################################################################
################################################################################
################################################################################
### Mosses SEM
### Author: korn@cumulonimbus.at University of Fribourg 2021
################################################################################


library("lavaan")
library("semPlot")
# library("blavaan")
library("gplots")
library("GGally")


# sessionInfo()


rm(list = ls())


setwd("~/repos/peatland-microbes/csv/")


mossMeta <- read.table("Mosses_Metadata_scaled.csv", sep = "\t", header = TRUE)


dim(mossMeta)
colnames(mossMeta)
# mMeta.df$FullID


################################################################################
### Testing 123
check1.vars <- c("Site", "Altitude", "Light", "Temperature", "Humidity",
                 "Soil.reaction", "Nitrogen", "CanopyCover", "Moss_g",
                 "Distance", "MossComposition")
newdata2 <- mossMeta[check1.vars]
print(round(cor(newdata2[, -1], use = "complete"), digits = 1))


newdata2$Site <- ordered(newdata2$Site,
                         levels = c("CB", "LE", "LV", "LT", "LM"))


# ggpairs(newdata2, aes(colour = Site, alpha = 0.4),
#         upper = list(continuous = wrap("cor", size = 3)),
#         lower = list(continuous = wrap("points", alpha = 0.3, size = 2),
#                      combo = wrap("dot", alpha = 0.4, size = 0.2))) +
#   theme_bw(base_size = 10) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave("Moss_MetadataPairs.pdf", width = 11.69, height = 11.69)


### Correlation and covariance
my_palette <- colorRampPalette(c("red", "white",  "blue"))(n = 19)


## Sort (out) correlation data.frame
# cors <- data.frame(cor(mossMeta[, -c(1:9, 11:15, 18, 22)], use = "complete"))
# cors <- cors[c("Altitude", "Distance", "CanopyCover", "Light", "Moss_g",
#                "MossHabitat", "Temperature", "Humidity", "Soil.reaction",
#                "Nitrogen", "mb.moss.prok", "mb.moss.euk", "mb.moss")]
# cors <- cors[colnames(cors), ]
#
#
# pdf("Moss_Metadata_Correlation.pdf", height = 8.27, width = 8.27)
# heatmap.2(as.matrix(cors), dendrogram = "none", key = TRUE, key.title = "",
#           cellnote = round(cors, digits = 2), notecol = "black", notecex = 0.7,
#           trace = "none", distfun = function(x) {x}, cexCol = 1.4, cexRow = 1.4,
#           symm = TRUE, margins = c(16, 16), col = my_palette,
#           tracecol = "black", Rowv = FALSE, Colv = FALSE)
# dev.off()


################################################################################
### SEM
## EV composite
ev <-
'mb.moss.prok ~ Light + Temperature + Humidity + Soil.reaction
mb.moss.euk ~ Light + Temperature + Humidity + Soil.reaction
mb.moss ~ Light + Temperature + Humidity + Soil.reaction'


fit.ev <- sem(ev, data = mossMeta)
summary(fit.ev, rsq = TRUE)


## Extract coefficients
Light.p <- coef(fit.ev)[[1]]
Temperature.p <- coef(fit.ev)[[2]]
Humidity.p <- coef(fit.ev)[[3]]
Soil.reaction.p <- coef(fit.ev)[[4]]


coef(fit.ev)[1:4]; Light.p; Temperature.p; Humidity.p; Soil.reaction.p


Light.e <- coef(fit.ev)[[5]]
Temperature.e <- coef(fit.ev)[[6]]
Humidity.e <- coef(fit.ev)[[7]]
Soil.reaction.e <- coef(fit.ev)[[8]]


coef(fit.ev)[5:8]; Light.e; Temperature.e; Humidity.e; Soil.reaction.e


Light.m <- coef(fit.ev)[[9]]
Temperature.m <- coef(fit.ev)[[10]]
Humidity.m <- coef(fit.ev)[[11]]
Soil.reaction.m <- coef(fit.ev)[[12]]


coef(fit.ev)[9:12]; Light.m; Temperature.m; Humidity.m; Soil.reaction.m


## Compute composite "ev"
mossMeta$ev.p <- Light.p * mossMeta$Light +
  Temperature.p * mossMeta$Temperature +
  Humidity.p * mossMeta$Humidity +
  Soil.reaction.p * mossMeta$Soil.reaction

mossMeta$ev.e <- Light.e * mossMeta$Light +
  Temperature.e * mossMeta$Temperature +
  Humidity.e  * mossMeta$Humidity +
  Soil.reaction.e * mossMeta$Soil.reaction

mossMeta$ev.m <- Light.m * mossMeta$Light +
  Temperature.m * mossMeta$Temperature +
  Humidity.m * mossMeta$Humidity +
  Soil.reaction.m * mossMeta$Soil.reaction


### Run SEM with the manually computed composite
ev.comp <-
'mb.moss.prok ~ ev.p
mb.moss.euk ~ ev.e
mb.moss ~ ev.m'

fit.ev.comp <- sem(ev.comp, data = mossMeta)


## Control if the standard coefficient and z-value (standardised coefficient)
## are similar
summary(fit.ev, rsq = TRUE, standardized = TRUE)
summary(fit.ev.comp, rsq = TRUE, standardized = TRUE)


rm(fit.ev.comp, ev, ev.comp, Light.p, Temperature.p, Humidity.p,
   Soil.reaction.p, Light.e, Temperature.e, Humidity.e, Soil.reaction.e,
   Light.m, Temperature.m, Humidity.m, Soil.reaction.m, moss,
   newdata2, check1.vars)


################################################################################
### SEM
set.seed(12614)


moss <-
'MossComposition ~ Altitude + CanopyCover + Distance
# ev.p ~ Altitude + CanopyCover + Distance
# ev.e ~ Altitude + CanopyCover + Distance
ev.m ~ Altitude + CanopyCover + Distance

mb.moss ~ Distance + Altitude + CanopyCover + ev.m + MossComposition + Moss_g
# mb.moss.prok ~ Distance + Altitude + CanopyCover + ev.p + MossComposition + Moss_g + mb.moss.euk
# mb.moss.euk  ~ Distance + Altitude + CanopyCover + ev.e + MossComposition + Moss_g + mb.moss.prok

ev.m ~~ MossComposition

# ev.p ~~ MossComposition
# ev.e ~~ MossComposition
# ev.e ~~ ev.p
# mb.moss.prok ~~ mb.moss.euk
'


fit.moss <- sem(moss, data = mossMeta, estimator = "MLM")
summary(fit.moss, rsq = TRUE, fit.measures = TRUE)


modificationindices(fit.moss, minimum.value = 3)


## Update model
fit.moss.up <- update(fit.moss, add = "ev.p ~~ ev.e")
summary(fit.moss.up, fit.measures = TRUE, rsq = TRUE)

modificationindices(fit.moss, minimum.value = 3, op = "~")


fit.moss.b <- bsem(moss, data = mossMeta, n.chains = 4,
                   bcontrol = list(cores = 6),
                   control = list(max_treedepth = 10))
summary(fit.moss.b, rsq = TRUE, fit.measures = TRUE)


# Value                                       -156.593       0.473 # wo mh
# Value                                       -163.824       0.413 # w mh
# Value                                       -188.741       0.363 # w mh + exp
# Value                                       -163.831       0.412


### Plot
pdf(file = "Moss_SEM.pdf", height = 8.27, width = 8.27)
semPaths(fit.moss.b, what = "est", whatLabels = "est", residuals = FALSE,
         intercepts = FALSE, sizeMan = 5, sizeMan2 = 3, edge.label.cex = 0.5,
         fade = FALSE, layout = "circle2", style = "mx", nCharNodes = 0,
         posCol = "#009e73ff", negCol = "#d55e00ff", edge.label.color = "black",
         layoutSplit = TRUE, curve = 2, curvature = 2, fixedStyle = 1,
         exoCov = FALSE, rotation = 1)
dev.off()


pdf(file = "Moss_SEM_composites.pdf", height = 8.27, width = 8.27)
semPaths(fit.ev, what = "est", whatLabels = "est", residuals = FALSE,
         intercepts = FALSE, sizeMan = 5, sizeMan2 = 3, edge.label.cex = 0.5,
         fade = FALSE, layout = "tree", style = "mx", nCharNodes = 0,
         posCol = "#009e73ff", negCol = "#d55e00ff", edge.label.color = "black",
         layoutSplit = TRUE, curve = 2, curvature = 2, fixedStyle = 1,
         exoCov = FALSE, rotation = 4)
dev.off()



################################################################################
### Simpler approach
moss.v2 <-
'
MossHabitat ~ Altitude + CanopyCover + Distance
# ev.p ~ Altitude + CanopyCover + Distance
# ev.e ~ Altitude + CanopyCover + Distance


mb.moss.prok ~ Moss_g + Distance + Altitude + CanopyCover + MossHabitat #+ ev.p
mb.moss.euk ~ Moss_g + Distance + Altitude + CanopyCover + MossHabitat #+ ev.e


# ev.p ~~ MossHabitat
# ev.e ~~ MossHabitat
# ev.p ~~ ev.e

mb.moss.prok ~~ mb.moss.euk'


fit.moss.v2 <- sem(moss.v2, data = mossMeta)
summary(fit.moss.v2, rsq = TRUE, fit.measures = TRUE)


modificationindices(fit.moss.v2, minimum.value = 3)


## Update model
# fit.moss.up <- update(fit.moss, add = "ev.p ~ Distance
#                       ev.e ~ Distance")
# summary(fit.moss.up, fit.measures = TRUE, rsq = TRUE)
#
# modificationindices(fit.moss, minimum.value = 3, op = "~")



fit.moss.v2.b <- bsem(moss.v2, data = mossMeta, n.chains = 4,
                   bcontrol = list(cores = 6),
                   control = list(max_treedepth = 10, adapt_delta = 0.9))
summary(fit.moss.v2.b, rsq = TRUE, fit.measures = TRUE)



### Plot
pdf(file = "Moss_SEM_v2_tst.pdf", height = 8.27, width = 8.27)
semPaths(fit.moss.v2.b, what = "est", whatLabels = "est", residuals = FALSE,
         intercepts = FALSE, sizeMan = 5, sizeMan2 = 3, edge.label.cex = 0.5,
         fade = FALSE, layout = "tree", style = "mx", nCharNodes = 0,
         posCol = "#009e73ff", negCol = "#d55e00ff", edge.label.color = "black",
         layoutSplit = TRUE, curve = 2, curvature = 2, fixedStyle = 1,
         exoCov = FALSE, rotation = 1)
dev.off()



################################################################################
################################################################################
################################################################################
################################################################################
