################################################################################
################################################################################
################################################################################
################################################################################
### Peatland microbes: read ASV
### Authors: Rachel Korn
### korn@cumulonimbus.at University of Fribourg 2021
################################################################################
### Returned items are (processed in the following order):
## * moss = full dataset
## * moss.g = taxa merged on genus level
## * moss.s = pruned from noise taxa
## * moss.a = abundance filtered


# rm(list = ls())
#
#
# asvDir <- "/home/r/Seafile/SMP_results"
# primer <- "18S"
# setwd(asvDir)


################################################################################
readTaxa <- function (primer = c("16S", "18S"),
merge.mosses = c(TRUE, FALSE)) {


  ##############################################################################
  ### Read taxa (ASV)
  if (primer == "16S")
    {
    seqtab.nochim <- readRDS(paste(asvDir,
                                   "SMP_16S_reseq/seqtab.nochim_SMP16S.rds",
                                   sep = "/"))

    taxa <- readRDS(paste(asvDir, "SMP_16S_reseq/taxa_SMP16S.rds", sep = "/"))
    }


  if (primer == "18S")
    {
    seqtab.nochim <- readRDS(paste(asvDir,
                                   "SMP_18S_reseq/seqtab.nochim_SMP18S.rds",
                                   sep = "/"))
    taxa <- readRDS(paste(asvDir, "SMP_18S_reseq/taxa_SMP18S.rds", sep = "/"))
    }


  moss <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
                   tax_table(taxa))


  dna <- Biostrings::DNAStringSet(taxa_names(moss))
  names(dna) <- taxa_names(moss)
  moss <- merge_phyloseq(moss, dna)
  taxa_names(moss) <- paste0("ASV", seq(ntaxa(moss)))


  ### Rename taxonomic ranks
  colnames(tax_table(moss)) <- c("Domain", "Phylum", "Class", "Order",
                                 "Family", "Genus")


  ##############################################################################
  ### Remove non-moss samples
  smpMeta <- read.table("~/repos/peatland-microbes/csv/smpMeta.csv",
                         sep = "\t", header = TRUE, row.names = 1)
  smpMeta$FullID <- rownames(smpMeta)


  smpMeta <- sample_data(smpMeta)
  moss <- merge_phyloseq(moss, smpMeta)
  # head(sample_data(moss))


  moss <- subset_samples(moss, Succession == "Moss")


  ### Merge resequenced samples
  sample_data(moss)$FullID <- gsub("r", "", sample_data(moss)$FullID)
  moss <- merge_samples(moss, "FullID")
  sample_data(moss) <- smpMeta


  ### Merge resequenced moss samples (samples *M and *Q) and restore metadata
  sample_data(moss)$FullID <- gsub("Q", "M", sample_data(moss)$FullID)
  moss <- merge_samples(moss, "FullID")


  ############################################################################
  ### Aggregate taxa by genus
  moss.g  <- tax_glom(moss, taxrank = rank_names(moss)[6],
                      NArm = FALSE, bad_empty = c("", " ", "\t"))


  ##############################################################################
  ### Read metadata
  mossMeta <- read.table("~/repos/peatland-microbes/csv/MossesMetadata.csv",
                         sep = "\t", header = TRUE)
  mossMeta$FullID <- gsub("_", "", mossMeta$FullID)
  rownames(mossMeta) <- mossMeta$FullID


  ## Merge ps object with metadata
  mossMeta <- sample_data(mossMeta)
  moss.g <- merge_phyloseq(moss.g, mossMeta)
  head(sample_data(moss.g))
  # sample_names(moss.g) == mossMeta$FullID


  ##############################################################################
  ### Remove spurious taxa
  if(primer == "16S")
    {
    moss.s <- subset_taxa(moss.g, !(Domain %in% c("unknown", "Eukaryota", NA) |
                                      Phylum %in% c("Eukaryota_unclassified",
                                                    NA) |
                                      Order %in% c("Chloroplast") |
                                      Family %in% c("Mitochondria")))
    }


  if(primer == "18S")
    {
    moss.s <- subset_taxa(moss.g, !(Domain %in% c("Bacteria", "unknown") |
                                      Phylum %in% c("Eukaryota_unclassified",
                                                    "Myxogastria",
                                                    "Apicomplexa",
                                                    "Neocallimastigomycota",
                                                    "Mollusca", "Vertebrata",
                                                    "Microsporidia",
                                                    "Mucoromycota",
                                                    "Archaeorhizomycetes", NA) |
                                      Class %in% c("Insecta", "Ellipura",
                                                   "Embryophyta", "Arachnida",
                                                   "Heterophyidae",
                                                   "Ichthyophonae",
                                                   "Arthropoda_unclassified",
                                                   "unclassified_Hexapoda",
                                                   "Ascomycota_unclassified",
                                                   "Agaricomycetes",
                                                   "Basidiomycota_unclassified",
                                                   "Exobasidiomycetes",
                                                   "Microbotryomycetes",
                                                   "Dothideomycetes",
                                                   "Pucciniomycetes",
                                                   "Spiculogloeomycetes",
                                                   "Ustilaginomycetes",
                                                   "Taphrinomycetes",
                                                   "Lecanoromycetes",
                                                   "Basidiomycota_unclassified",
                                                   "Sordariomycetes",
                                                   "Pezizomycetes",
                                                   "Entorrhizomycetes",
                                                   "Agaricostilbomycetes",
                                                   "Orbiliomycetes",
                                                   "Eurotiomycetes",
                                                   "Leotiomycetes",
                                                   "Entomophthoromycetes",
                                                   "Dacrymycetes") |
                                      Order %in% c("Filobasidiales",
                                                   "Spizellomycetales",
                                                   "Rhytismatales",
                                                   "Dothideomycetes",
                                                   "Mytilinidiales",
                                                   "Naohideales",
                                                   "Pleosporales",
                                                   "Kickxellales",
                                                   "Archaeorhizomycetes",
                                                   "Atractiellomycetes",
                                                   "Trichosporonales") |
                                      Family %in% c("Pleosporaceae",
                                                    "Mrakiaceae",
                                                    "Teratosphaeriaceae",
                                                    "Leotiaceae",
                                                    "Pleosporales_fa",
                                                    "Didymellaceae",
                                                    "Phaeosphaeriaceae",
                                                    "Rhynchogastremataceae",
                                                    "Phaeotremellaceae") |
                                      Genus %in% c("Lecophagus",
                                                   "Genolevuria",
                                                   "Trimorphomyces") |
                                      Phylum == "Arthropoda" & Class %in% NA |
                                      Phylum == "Ascomycota" & Class %in% NA |
                                      Phylum == "Basidiomycota" &
                                      Class %in% NA))
    }


  ## Remove empty taxa
  any(taxa_sums(moss.s) == 0)
  sum(taxa_sums(moss.s) == 0)
  moss.s <- prune_taxa(taxa_sums(moss.s) > 0, moss.s)


  ## Relative abundance
  moss.r <- transform_sample_counts(moss.s, function(otu) {otu / sum(otu)})


  ## Abundance filtering
  moss.a <- filter_taxa(moss.r, function(x) {mean(x) > 0.000001}, TRUE)


  ## Remove empty taxa
  moss.a <- prune_taxa(taxa_sums(moss.a) > 0, moss.a)


  ##############################################################################
  ### Return variables
  moss.l <- list(mossMeta = mossMeta, moss = moss,
                 moss.s = moss.s,  moss.r = moss.r,
                 moss.a = moss.a)

  return(moss.l)
  }


################################################################################
################################################################################
################################################################################
################################################################################
