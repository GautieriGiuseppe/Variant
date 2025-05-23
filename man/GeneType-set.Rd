% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Methods.R
\name{GeneType<-}
\alias{GeneType<-}
\alias{GeneType<-,Variant,character-method}
\title{Set Gene type}
\usage{
GeneType(object) <- value

\S4method{GeneType}{Variant,character}(object) <- value
}
\arguments{
\item{object}{An object of class Variant}

\item{value}{gene type to set}
}
\value{
The modifed object
}
\description{
Set the gene type of the variant
}
\examples{
gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
 strand = "+")
snv <- SNV(
  VariantID = "old_id",
  GenomicPos = gr,
  ReferenceAllele = "A",
  MutatedAllele = "G",
  GeneType = "type",
  transition = TRUE
)
GeneType(snv) <- "new_type"
GeneType(snv)
}
