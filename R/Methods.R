# Get description of Variant
#' @title describe
#' @description Returns a description of the variant
#' @param object An object of class Variant
#' @return A character string
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100), strand = "+")
#' snv <- SNV(
#'   VariantID = "SNV001",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   transition = TRUE
#' )
#' describe(snv)
#'
#' dnv <- DNV(
#'   VariantID = "DNV001",
#'   GenomicPos = gr,
#'   ReferenceAllele = "C",
#'   MutatedAllele = "T",
#'   validated = TRUE
#' )
#' describe(dnv)
#' @rdname describe
#' @export
setGeneric("describe", function(object) standardGeneric("describe"))
#' @rdname describe
#' @export
setMethod("describe", "Variant", function(object){
  paste("Generic Variant at", object@GenomicPos)
})
#' @rdname describe
#' @export
setMethod("describe", "SNV", function(object){
  genomepos <- object@GenomicPos
  chr <- as.character(seqnames(genomepos))
  pos <- start(genomepos)
  paste("Single Nucleotide Variant at", chr, pos,
        if (object@transition) "transition" else "transversion")
})
#' @rdname describe
#' @export
setMethod("describe", "DNV", function(object){
  genomepos <- object@GenomicPos
  chr <- as.character(seqnames(genomepos))
  pos <- start(genomepos)
  paste("De Novo Variant at", chr, pos,
        if (object@validated) "validated" else "non validated")
})
#' @rdname describe
#' @export
setMethod("describe", "ONV", function(object){
  genomepos <- object@GenomicPos
  chr <- as.character(seqnames(genomepos))
  pos <- start(genomepos)
  paste("Other Nucleotide Variant at", chr, pos, "type:", object@complexity
  )})
#' @rdname describe
#' @export
setMethod("describe", "Insertion", function(object){
  genomepos <- object@GenomicPos
  chr <- as.character(seqnames(genomepos))
  pos <- start(genomepos)
  paste("Insertion Variant at", chr, pos, "length:", object@lengthInserted
  )})
#' @rdname describe
#' @export
setMethod("describe", "Deletion", function(object){
  genomepos <- object@GenomicPos
  chr <- as.character(seqnames(genomepos))
  pos <- start(genomepos)
  paste("Deletion Variant at", chr, pos, "length:", object@lengthDeleted
  )})
#' Get Variant ID
#' @title Get Variant ID
#' @description Returns the unique ID of the variant
#' @param object An object of class Variant
#' @return A character string
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100), strand = "+")
#' snv <- SNV(
#'   VariantID = "SNV001",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   transition = TRUE
#' )
#' VariantID(snv)
#' @rdname VariantID
#' @export
setGeneric("VariantID", function(object) standardGeneric("VariantID"))
#' @rdname VariantID
#' @export
setMethod("VariantID", "Variant", function(object){object@VariantID})

# Get Ensembl ID
#' @title Get Ensembl ID
#' @description Returns the ensembl ID of the variant
#' @param object An object of class Variant
#' @return A character string
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100), strand = "+")
#' snv <- SNV(
#'   VariantID = "SNV001",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   transition = TRUE
#' )
#' EnsemblID(snv)
#' @rdname EnsemblID
#' @export
setGeneric("EnsemblID", function(object) standardGeneric("EnsemblID"))
#' @rdname EnsemblID
#' @export
setMethod("EnsemblID", "Variant", function(object){object@EnsemblID})

# Get the Genomic Position of the Variant
#' @title Get Position of variant
#' @description Returns the position of the variant
#' @param object An object of class GRanges
#' @return A character string
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100), strand = "+")
#' snv <- SNV(
#'   VariantID = "SNV001",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   transition = TRUE
#' )
#' Position(snv)
#' @rdname GenomicPos
#' @export
setGeneric("Position", function(object) standardGeneric("Position"))
#' @rdname GenomicPos
#' @export
setMethod("Position", "Variant", function(object){object@GenomicPos})

# Get the reference allele and the mutated one.
#' @title Get Reference allele and mutated
#' @description Returns the reference allele and mutated
#' @param object An object of class Variant
#' @return A character string
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
#'  strand = "+")
#' snv <- SNV(
#'   VariantID = "SNV001",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   transition = TRUE
#' )
#' Variation(snv)
#' @rdname Variation
#' @export
setGeneric("Variation", function(object) standardGeneric("Variation"))
#' @rdname Variation
#' @export
setMethod("Variation", "Variant", function(object){
  paste("The reference allele:", object@ReferenceAllele,
        " | The variant allele  :", object@MutatedAllele)
})

# Get the Mutated allele
#' @title Get Mutated Allele
#' @description Returns the mutated allele the variant
#' @param object An object of class Variant
#' @return A character string
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
#'  strand = "+")
#' snv <- SNV(
#'   VariantID = "SNV001",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   transition = TRUE
#' )
#' MutatedAllele(snv)
#' @rdname MutatedAllele
#' @export
setGeneric("MutatedAllele", function(object) standardGeneric("MutatedAllele"))
#' @rdname MutatedAllele
#' @export
setMethod("MutatedAllele", "Variant", function(object) {
  object@MutatedAllele
})
# Get the Reference allele
#' @title Get Reference Allele
#' @description Returns the reference allele the variant
#' @param object An object of class Variant
#' @return A character string
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
#'  strand = "+")
#' snv <- SNV(
#'   VariantID = "SNV001",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   transition = TRUE
#' )
#' ReferenceAllele(snv)
#' @rdname ReferenceAllele
#' @export
setGeneric("ReferenceAllele",
           function(object) standardGeneric("ReferenceAllele"))
#' @rdname ReferenceAllele
#' @export
setMethod("ReferenceAllele", "Variant", function(object) {
  object@ReferenceAllele
})
# Get the gene symbol
#' @title Get Gene Symbol
#' @description Returns the gene symbol of the variant
#' @param object An object of class Variant
#' @return A character string
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
#'  strand = "+")
#' snv <- SNV(
#'   VariantID = "SNV001",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   transition = TRUE
#' )
#' GeneSymbol(snv)
#' @rdname GeneSymbol
#' @export
setGeneric("GeneSymbol", function(object) standardGeneric("GeneSymbol"))
#' @rdname GeneSymbol
#' @export
setMethod("GeneSymbol", "Variant", function(object){object@GeneSymbol})

# Get the exon in which is located the Variant
#' @title Get Exon
#' @description Returns the exon of the variant
#' @param object An object of class Variant
#' @return A character string
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
#'  strand = "+")
#' snv <- SNV(
#'   VariantID = "SNV001",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   transition = TRUE
#' )
#' Exon(snv)
#' @rdname Exon
#' @export
setGeneric("Exon", function(object) standardGeneric("Exon"))
#' @rdname Exon
#' @export
setMethod("Exon", "Variant", function(object){object@Exon})

# Get the gene type
#' @title Get Gene Type
#' @description Returns the gene type of the variant
#' @param object An object of class Variant
#' @return A character string
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
#'  strand = "+")
#' snv <- SNV(
#'   VariantID = "SNV001",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   transition = TRUE
#' )
#' GeneType(snv)
#' @rdname GeneType
#' @export
setGeneric("GeneType", function(object) standardGeneric("GeneType"))
#' @rdname GeneType
#' @export
setMethod("GeneType", "Variant", function(object){object@GeneType})

# Get the consequence of mutation
#' @title Get Consequence
#' @description Returns the consequence of the variant
#' @param object An object of class Variant
#' @return A character string
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
#'  strand = "+")
#' snv <- SNV(
#'   VariantID = "SNV001",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   transition = TRUE
#' )
#' MutationConsequence(snv)
#' @rdname MutationConsequence
#' @export
setGeneric("MutationConsequence", 
           function(object) standardGeneric("MutationConsequence"))
#' @rdname MutationConsequence
#' @export
setMethod("MutationConsequence", "Variant",
          function(object){object@MutationConsequence})

# Get the clinical relevance
#' @title Get Clinical Relevance
#' @description Returns the clinical relevance of the variant
#' @param object An object of class Variant
#' @return A character string
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
#'  strand = "+")
#' snv <- SNV(
#'   VariantID = "SNV001",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   transition = TRUE
#' )
#' Relevance(snv)
#' @rdname Relevance
#' @export
setGeneric("Relevance", function(object) standardGeneric("Relevance"))
#' @rdname Relevance
#' @export
setMethod("Relevance", "Variant", function(object){object@ClinicalRel})

# Get the protein variation
#' @title Get Protein Variation
#' @description Returns the protein variation of the variant
#' @param object An object of class Variant
#' @return A character string
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
#'  strand = "+")
#' snv <- SNV(
#'   VariantID = "SNV001",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   transition = TRUE
#' )
#' ProteinVariation(snv)
#' @rdname ProteinVariation
#' @export
setGeneric("ProteinVariation",
           function(object) standardGeneric("ProteinVariation"))
#' @rdname ProteinVariation
#' @export
setMethod("ProteinVariation", "Variant",
          function(object){object@ProteinVariation})

# setters

# Set the Variant ID of a Variant
#' @title Set Variant ID
#' @description Set the unique ID of the variant
#' @param object An object of class Variant
#' @param value The ID to insert
#' @return The modified object
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
#'  strand = "+")
#' snv <- SNV(
#'   VariantID = "old_id",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   transition = TRUE
#' )
#' VariantID(snv) <- "new_id"
#' VariantID(snv)
#' @rdname VariantID-set
#' @export
setGeneric("VariantID<-",
           function(object, value) standardGeneric("VariantID<-"))
#' @rdname VariantID-set
#' @export
setReplaceMethod("VariantID",
                 signature(object = "Variant", value = "character"),
                 function(object, value){
                   object@VariantID <- as.character(value)
                   validObject(object)
                   return(object)
                 })

# Set the exon of a Variant
#' @title Set Exon
#' @description Set the exon of the variant
#' @param object An object of class Variant
#' @param value The exon to set
#' @return The modifed object
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
#'  strand = "+")
#' snv <- SNV(
#'   VariantID = "variantid",
#'   Exon = "old_exon",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   transition = TRUE
#' )
#' Exon(snv) <- "new_exon"
#' Exon(snv)
#' @rdname Exon-set
#' @export
setGeneric("Exon<-", function(object, value) standardGeneric("Exon<-"))
#' @rdname Exon-set
#' @export
setReplaceMethod("Exon",
                 signature(object = "Variant", value = "character"),
                 function(object, value){
                   object@Exon <- as.character(value)
                   validObject(object)
                   return(object)
                 })

# Set the reference allele of a Variant
#' @title Set reference allele
#' @description Set the reference allele of the variant
#' @param object An object of class Variant
#' @param value Reference allele
#' @return The modifed object
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
#'  strand = "+")
#' snv <- SNV(
#'   VariantID = "variantid",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   transition = TRUE
#' )
#' ReferenceAllele(snv) <- "new_reference"
#' ReferenceAllele(snv)
#' @rdname ReferenceAllele-set
#' @export
setGeneric("ReferenceAllele<-",
           function(object, value) standardGeneric("ReferenceAllele<-"))
#' @rdname ReferenceAllele-set
#' @export
setReplaceMethod("ReferenceAllele",
                 signature(object = "Variant", value = "character"),
                 function(object, value){
                   object@ReferenceAllele <- as.character(value)
                   validObject(object)
                   return(object)
                 })

# Set the mutated allele of a Variant
#' @title Set mutated allele
#' @description Set the mutated allele of the variant
#' @param object An object of class Variant
#' @param value The mutated allele
#' @return The modified object
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
#'  strand = "+")
#' snv <- SNV(
#'   VariantID = "variantid",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   transition = TRUE
#' )
#' MutatedAllele(snv) <- "new_mutated"
#' MutatedAllele(snv)
#' @rdname MutatedAllele-set
#' @export
setGeneric("MutatedAllele<-",
           function(object, value) standardGeneric("MutatedAllele<-"))
#' @rdname MutatedAllele-set
#' @export
setReplaceMethod("MutatedAllele",
                 signature(object = "Variant", value = "character"),
                 function(object, value){
                   object@MutatedAllele <- as.character(value)
                   validObject(object)
                   return(object)
                 })

# Set the gene symbol of a Variant
#' @title Set Gene symbol
#' @description Set the gene symbol of the variant
#' @param object An object of class Variant
#' @param value Gene symbol to set
#' @return The modifed object
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
#'  strand = "+")
#' snv <- SNV(
#'   VariantID = "old_id",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   GeneSymbol = "symbol",
#'   transition = TRUE
#' )
#' GeneSymbol(snv) <- "new_symbol"
#' GeneSymbol(snv)
#' @rdname GeneSymbol-set
#' @export
setGeneric("GeneSymbol<-",
           function(object, value) standardGeneric("GeneSymbol<-"))
#' @rdname GeneSymbol-set
#' @export
setReplaceMethod("GeneSymbol",
                 signature(object = "Variant", value = "character"),
                 function(object, value){
                   object@GeneSymbol <- as.character(value)
                   validObject(object)
                   return(object)
                 })

# Set the gene type of a Variant
#' @title Set Gene type
#' @description Set the gene type of the variant
#' @param object An object of class Variant
#' @param value gene type to set
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
#'  strand = "+")
#' snv <- SNV(
#'   VariantID = "old_id",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   GeneType = "type",
#'   transition = TRUE
#' )
#' GeneType(snv) <- "new_type"
#' GeneType(snv)
#' @return The modifed object
#' @rdname GeneType-set
#' @export
setGeneric("GeneType<-", function(object, value) standardGeneric("GeneType<-"))
#' @rdname GeneType-set
#' @export
setReplaceMethod("GeneType",
                 signature(object = "Variant", value = "character"),
                 function(object, value){
                   object@GeneType <- as.character(value)
                   validObject(object)
                   return(object)
                 })

# Set the mutation consequence of a Variant
#' @title Set mutation consequence
#' @description Set the mutation consequence of the variant
#' @param object An object of class Variant
#' @param value The mutation consequence
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
#'  strand = "+")
#' snv <- SNV(
#'   VariantID = "old_id",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   transition = TRUE
#' )
#' MutationConsequence(snv) <- "new_consequence"
#' MutationConsequence(snv)
#' @return The modified object
#' @rdname MutationConsequence-set
#' @export
setGeneric("MutationConsequence<-",
           function(object, value) standardGeneric("MutationConsequence<-"))
#' @rdname MutationConsequence-set
#' @export
setReplaceMethod("MutationConsequence",
                 signature(object = "Variant", value = "character"),
                 function(object, value){
                   object@MutationConsequence <- as.character(value)
                   validObject(object)
                   return(object)
                 })

# Set clinical relevance of a Variant
#' @title Set clinical relevance
#' @description Set the clinical relevance of the variant
#' @param object An object of class Variant
#' @param value The clinical relevance
#' @return The modified object
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
#'  strand = "+")
#' snv <- SNV(
#'   VariantID = "old_id",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   transition = TRUE
#' )
#' Relevance(snv) <- "new_rel"
#' Relevance(snv)
#' @rdname ClinicalRel-set
#' @export
setGeneric("Relevance<-",
           function(object, value) standardGeneric("Relevance<-"))
#' @rdname ClinicalRel-set
#' @export
setReplaceMethod("Relevance",
                 signature(object = "Variant", value = "character"),
                 function(object, value){
                   object@ClinicalRel <- as.character(value)
                   validObject(object)
                   return(object)
                 })

# Set the Ensembl ID of a Variant
#' @title Set Ensembl ID
#' @description Set the Ensembl ID of the variant
#' @param object An object of class Variant
#' @param value The Ensembl ID
#' @return The modified object
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
#'  strand = "+")
#' snv <- SNV(
#'   VariantID = "old_id",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   transition = TRUE
#' )
#' EnsemblID(snv) <- "new_id"
#' EnsemblID(snv)
#' @rdname EnsemblID-set
#' @export
setGeneric("EnsemblID<-",
           function(object, value) standardGeneric("EnsemblID<-"))
#' @rdname EnsemblID-set
#' @export
setReplaceMethod("EnsemblID",
                 signature(object = "Variant", value = "character"),
                 function(object, value){
                   object@EnsemblID <- as.character(value)
                   validObject(object)
                   return(object)
                 })

#' Set the protein variation of a Variant
#' @title Set protein variation
#' @description Set the protein variation of the variant
#' @param object An object of class Variant
#' @param value The protein variation
#' @return The modified object
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100), strand = "+")
#' snv <- SNV(
#'   VariantID = "old_id",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   transition = TRUE
#' )
#' ProteinVariation(snv) <- "new_variation"
#' ProteinVariation(snv)
#' @rdname ProteinVariation-set
#' @export
setGeneric("ProteinVariation<-", function(object, value) standardGeneric("ProteinVariation<-"))

#' @rdname ProteinVariation-set
#' @export
setReplaceMethod("ProteinVariation",
                 signature(object = "Variant", value = "character"),
                 function(object, value){
                   object@ProteinVariation <- as.character(value)
                   validObject(object)
                   return(object)
                 })
