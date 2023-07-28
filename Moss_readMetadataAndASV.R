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


################################################################################
readTaxa <- function (primer = c("16S", "18S"),
                      merge.mosses = c(TRUE, FALSE)) {


  ##############################################################################
  ### Read taxa (ASV)
  if (primer == "16S") {
    seqtab.nochim <- readRDS(paste(asvDir,
                                   "SMP_16S_reseq/seqtab.nochim_SMP16S.rds",
                                   sep = "/"))

    taxa <- readRDS(paste(asvDir, "SMP_16S_reseq/taxa_SMP16S.rds", sep = "/"))
  }


  if (primer == "18S") {
    seqtab.nochim <- readRDS(paste(asvDir,
                                   "SMP_18S_reseq/seqtab.nochim_SMP18S.rds",
                                   sep = "/"))
    taxa <- readRDS(paste(asvDir, "SMP_18S_reseq/taxa_SMP18S.rds", sep = "/"))
  }


  moss <- phyloseq(otu_table(seqtab.nochim,
                             taxa_are_rows = FALSE),
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
  smpMeta <- read.table("csv/smpMeta.csv",
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


  ## Remove empty taxa
  any(taxa_sums(moss) == 0)
  sum(taxa_sums(moss) == 0)
  moss <- prune_taxa(taxa_sums(moss) > 0, moss)


  ##############################################################################
  ### Remove unidentified taxa
  if(primer == "16S") {
    moss.id <- subset_taxa(moss,
                           !(Domain %in% c("unknown", "Eukaryota", NA) |
                               Phylum %in% c("Eukaryota_unclassified",
                                             NA)))
  }

  if(primer == "18S") {
    moss.id <- subset_taxa(moss, !(Domain %in% c("Bacteria", "unknown") |
                                     Phylum %in% c("Eukaryota_unclassified",
                                                   NA)))
  }

  ## Remove empty taxa (there should be none)
  # sum(taxa_sums(moss.s) == 0)
  moss.id <- prune_taxa(taxa_sums(moss.id) > 0, moss.id)

  message(paste("Note:",
                sum(taxa_sums(moss.id) == 0),
                "empty taxa were removed after removing unidentified taxa."))


  message(paste(round(100 - 100 / ntaxa(moss) * ntaxa(moss.id),
                      digits = 2),
                "% unidentified",
                ifelse(primer == "16S",
                       "prokaryote",
                       "eukaryote"),
                "taxa were removed from the samples."))


  ### Remove spurious taxa
  if(primer == "16S") {
    moss.s <- subset_taxa(moss.id,
                          !(Order %in% c("Chloroplast") |
                              Family %in% c("Mitochondria")))
  }

  if(primer == "18S") {
    moss.s <- subset_taxa(moss.id,
                          !(Phylum %in% c("Myxogastria",
                                          "Apicomplexa",
                                          "Neocallimastigomycota",
                                          "Mollusca", "Vertebrata",
                                          "Microsporidia",
                                          "Mucoromycota",
                                          "Archaeorhizomycetes") |
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
                              Phylum == "Basidiomycota"))
  }



  ## Remove empty taxa (there should be none)
  if (sum(taxa_sums(moss.s) == 0)) {
    message(paste("Warning:",
                  sum(taxa_sums(moss.s) == 0), "empty taxa were removed."))
    moss.s <- prune_taxa(taxa_sums(moss.s) > 0, moss.s)
  }



  message(paste(round(100 - 100 / ntaxa(moss.id) * ntaxa(moss.s),
                      digits = 2),
                "% unidentified",
                ifelse(primer == "16S",
                       "prokaryote",
                       "eukaryote"),
                "taxa were removed from the samples."))


  ##############################################################################
  ### Aggregate taxa by genus
  moss.g  <- tax_glom(moss.s,
                      taxrank = rank_names(moss.s)[6],
                      NArm = FALSE, bad_empty = c("", " ", "\t"))

  message(paste(ntaxa(moss.g),
                "of",
                ntaxa(moss.s),
                "ASVs (",
                round(100 / ntaxa(moss.s) * ntaxa(moss.g),
                      digits = 2),
                "%) remained after aggregating",
                ifelse(primer == "16S",
                       "prokaryote",
                       "eukaryote"),
                "ASVs on genus level"))


  ##############################################################################
  ### Read metadata
  mossMeta <- read.table("csv/MossesMetadata.csv",
                         sep = "\t", header = TRUE)
  mossMeta$FullID <- gsub("_", "", mossMeta$FullID)
  rownames(mossMeta) <- mossMeta$FullID


  ## Merge ps object with metadata
  mossMeta <- sample_data(mossMeta)
  moss.g <- merge_phyloseq(moss.g, mossMeta)

  # sample_names(moss.g) == mossMeta$FullID

  sample_data(moss.g)$FullID <- rownames(sample_data(moss.g))
  sample_data(moss.g)$Sector <- substr(sample_data(moss.g)$FullID, 3, 3)
  sample_data(moss.g)$Succession <- NULL
  head(sample_data(moss.g))


  ## Relative abundance
  moss.r <- transform_sample_counts(moss.g, function(otu) {otu / sum(otu)})


  ## Abundance filtering
  moss.a <- filter_taxa(moss.r, function(x) {mean(x) > 0.000001}, TRUE)


  ## Remove empty taxa
  moss.a <- prune_taxa(taxa_sums(moss.a) > 0, moss.a)

  ## Remove empty taxa (there should be none)
  if (sum(taxa_sums(moss.a) == 0)) {
    message(paste("Warning:",
                  sum(taxa_sums(moss.a) == 0), "empty taxa were removed."))
    moss.a <- prune_taxa(taxa_sums(moss.a) > 0, moss.a)
  }


  ##############################################################################
  ### Return variables
  moss.l <- list(mossMeta = mossMeta,
                 # moss = moss,
                 # moss.s = moss.s,
                 # moss.r = moss.r,
                 moss.a = moss.a)


  return(moss.l)
}


################################################################################
################################################################################
################################################################################
################################################################################
