% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Constructors.R
\name{DNV}
\alias{DNV}
\title{Constructor for class DNV}
\usage{
DNV(
  VariantID,
  GenomicPos,
  ReferenceAllele,
  MutatedAllele,
  validated = FALSE,
  Exon = NA_character_,
  GeneSymbol = NA_character_,
  GeneType = NA_character_,
  MutationConsequence = NA_character_,
  ClinicalRel = NA_character_,
  EnsemblID = NA_character_,
  ProteinVariation = NA_character_
)
}
\arguments{
\item{VariantID}{Character. Unique ID of the variant.}

\item{GenomicPos}{\code{GRanges} object for the position of the variant.}

\item{ReferenceAllele}{Character. The allele in the reference genome.}

\item{MutatedAllele}{Character. The alternative or mutated allele.}

\item{validated}{Logical. TRUE if the de novo is validated.}

\item{Exon}{Character. The exon where the variant is located.}

\item{GeneSymbol}{Character. Gene name (optional).}

\item{GeneType}{Character. Type of gene e.g. "protein_coding" (optional).}

\item{MutationConsequence}{Character. Effect in the phenotype (optional).}

\item{ClinicalRel}{Character. Clinical annotation (optional).}

\item{EnsemblID}{Character. Transcript from Ensembl (optional).}

\item{ProteinVariation}{Character. HGVS protein change (optional).}
}
\value{
class \code{DNV} object
}
\description{
Constructor to create an object of class DNV
}
\examples{
gr <- GenomicRanges::GRanges("chr13", IRanges::IRanges(32914438, 32914438),
 strand = "+")
dnv <- DNV(
  VariantID = "DNV1",
  GenomicPos = gr,
  ReferenceAllele = "G",
  MutatedAllele = "T",
  validated = FALSE
)
}
