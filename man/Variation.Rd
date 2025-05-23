% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Methods.R
\name{Variation}
\alias{Variation}
\alias{Variation,Variant-method}
\title{Get Reference allele and mutated}
\usage{
Variation(object)

\S4method{Variation}{Variant}(object)
}
\arguments{
\item{object}{An object of class Variant}
}
\value{
A character string
}
\description{
Returns the reference allele and mutated
}
\examples{
gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
 strand = "+")
snv <- SNV(
  VariantID = "SNV001",
  GenomicPos = gr,
  ReferenceAllele = "A",
  MutatedAllele = "G",
  transition = TRUE
)
Variation(snv)
}
