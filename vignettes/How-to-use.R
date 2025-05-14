## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(Variant)
library(GenomicRanges)


## -----------------------------------------------------------------------------
gr <- GRanges("chr1", IRanges(100, 100), strand = "+")
snv <- SNV(
  VariantID = "SNV001",
  GenomicPos = gr,
  ReferenceAllele = "A",
  MutatedAllele = "G",
  transition = TRUE,
  Exon = "1",
  GeneSymbol = "TP53",
  GeneType = "protein_coding",
  MutationConsequence = "missense_variant",
  ClinicalRel = "Pathogenic",
  EnsemblID = "ENST00000269305",
  ProteinVariation = "p.Pro34Arg"
)

## -----------------------------------------------------------------------------
VariantID(snv)
GeneSymbol(snv)
MutationConsequence(snv)

