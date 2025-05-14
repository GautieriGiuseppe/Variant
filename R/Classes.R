# Scientific Programming Project
# R package for variant objects representing different mutation types 

# Set of S4 classes for different mutation types

# S4 class representing a generic variant
#' @slot VariantID Character. ID of variant.
#' @slot GenomicPos GRanges object. Genomic position. 
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom GenomicRanges GRanges seqnames start strand
#' @importFrom methods new validObject
#' @importFrom IRanges IRanges
#' @slot Exon Character. Exons involved.
#' @slot ReferenceAllele Character. Allele Reference.
#' @slot MutatedAllele Character. Mutated allele.
#' @slot GeneSymbol Character. Gene Symbol.
#' @slot GeneType Character. The type of gene.
#' @slot MutationConsequence Character. The biological consequence of variant.
#' @slot ClinicalRel Character. The clinical relevance.
#' @slot EnsemblID Character. The Ensembl ID.
#' @slot ProteinVariation Character. The variation of protein.
#' @title Class Variant: Generic Variant
#' @description S4 class for a generic variant
#' @name Variant-class
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", 
#' IRanges::IRanges(100, 100), strand = "+")
#' v <- Variant(
#'   VariantID = "v1",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "T"
#' )
#' VariantID(v)
#' @export
setClass(
    "Variant",
    slots = list(
        VariantID = "character",
        GenomicPos = "GRanges",
        Exon = "character",
        ReferenceAllele = "character",
        MutatedAllele = "character",
        GeneSymbol = "character",
        GeneType = "character",
        MutationConsequence = "character",
        ClinicalRel = "character",
        EnsemblID = "character",
        ProteinVariation = "character"
      )
    )

setValidity("Variant", function(object){
    position <- object@GenomicPos
    if (length(position) != 1) return("must be a single GRanges object")
    if (start(position) <= 0) return("Start must be positive integer")
    if (!as.character(seqnames(position)) 
        %in% paste0("chr", c(1:22, "X", "Y")))
      return("Chromosome must be chr1-22, chrX or chrY")
    if (!as.character(strand(position)) %in% c("+", "-"))
      return("Strand must be '+' or '-' .")
    TRUE
})

# Class for Single Nucleotide Variant
#' @title Class SNV: Single Nucleotide Variant
#' @description S4 class for a single nucleotide variant
#' @name SNV-class
#' @slot transition Logical. Substitution is a transition.
#' @examples
#' gr <- GenomicRanges::GRanges("chr13", IRanges::IRanges(32914438, 32914438),
#'  strand = "+")
#' snv <- SNV(
#'   VariantID = "SNV1",
#'   GenomicPos = gr,
#'   ReferenceAllele = "G",
#'   MutatedAllele = "A",
#' )
#' @export
setClass(
    "SNV",
    contains = "Variant",
    slots = list(
        transition = "logical" # TRUE when A->G or C->T
    )
)

#' @title Is Transition
#' @description Returns TRUE if transition
#' @param object An object of class SNV 
#' @return A logical value
#' @rdname isTransition
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
#'  strand = "+")
#' snv <- SNV(VariantID = "id", GenomicPos = gr, ReferenceAllele = "A",
#'  MutatedAllele = "G", transition = TRUE)
#' isTransition(snv)
#' @export
setGeneric("isTransition", function(object) standardGeneric("isTransition"))
#' @rdname isTransition
#' @export
setMethod("isTransition", "SNV", function(object) object@transition)
#' @title Class DNV: De Novo Variant
#' @description S4 class for a de novo variant
#' @name DNV-class
#' @slot validated Logical. DNV has been confirmed
#' @examples
#' gr <- GenomicRanges::GRanges("chr13", IRanges::IRanges(32914438, 32914438),
#'  strand = "+")
#' dnv <- DNV(
#'   VariantID = "DNV1",
#'   GenomicPos = gr,
#'   ReferenceAllele = "G",
#'   MutatedAllele = "T",
#'   validated = FALSE
#' )
#' @export
setClass(
    "DNV",
    contains = "Variant",
    slots = list(
        validated = "logical"
    )
)
#' @title Class ONV: Other Nucleotide Variant
#' @description S4 class for an other nucleotide variant
#' @name ONV-class
#' @slot complexity Character. Type of other variant
#' @examples
#' gr <- GenomicRanges::GRanges("chr22", IRanges::IRanges(29091840, 29091840),
#'  strand = "+")
#' onv <- ONV(
#'   VariantID = "ONV001",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "G",
#'   complexity = "translocation"
#' )
#' describe(onv)
#' VariantID(onv)
#' @export
setClass(
    "ONV",
    contains = "Variant",
    slots = list(
        complexity = "character"
    )
)
#' @title Class Insertion: Insertion Variant
#' @description S4 class for an insertion variant
#' @name Insertion-class
#' @slot lengthInserted Integer. The length of inserted sequence.
#' @examples
#' gr <- GenomicRanges::GRanges("chr13", IRanges::IRanges(32914438, 32914438),
#'  strand = "+")
#' ins <- Insertion(
#'   VariantID = "INS1",
#'   GenomicPos = gr,
#'   ReferenceAllele = "-",
#'   MutatedAllele = "G",
#'   lengthInserted = as.integer(1)
#' )
#' isFrameshift(ins)
#' @export
setClass(
    "Insertion",
    contains = "Variant",
    slots = list(
        lengthInserted = "integer"
    )
)
#' @title Is Frameshift
#' @description Returns TRUE if is a frameshift mutation
#' @param object An object of class Insertion
#' @return A logical value
#' @rdname isFrameshift
#' @examples
#' gr <- GenomicRanges::GRanges("chr13", IRanges::IRanges(32914438, 32914438),
#'  strand = "+")
#' ins <- Insertion(
#'   VariantID = "INS1",
#'   GenomicPos = gr,
#'   ReferenceAllele = "-",
#'   MutatedAllele = "T",
#'   lengthInserted = as.integer(1)
#' )
#' isFrameshift(ins)
#' @export
setGeneric("isFrameshift", function(object) standardGeneric("isFrameshift"))
#' @rdname isFrameshift
#' @export
setMethod("isFrameshift", "Insertion", function(object){
  object@lengthInserted %% 3 != 0
})
#' @title Class Deletion: Deletion Variant
#' @description S4 class for a deletion variant
#' @name Deletion-class
#' @slot lengthDeleted Integer. The length of deleted sequence.
#' @examples
#' gr <- GenomicRanges::GRanges("chr13", IRanges::IRanges(32914438, 32914438),
#'  strand = "+")
#' del <- Deletion(
#'   VariantID = "DEL1",
#'   GenomicPos = gr,
#'   ReferenceAllele = "G",
#'   MutatedAllele = "-",
#'   lengthDeleted = as.integer(1)
#' )
#' isFrameshift(del)
#' @export
setClass(
    "Deletion",
    contains = "Variant",
    slots = list(
        lengthDeleted = "integer"
      )
  )
#' @rdname isFrameshift
#' @export
#' @examples
#' gr <- GenomicRanges::GRanges("chr13", IRanges::IRanges(32914438, 32914438),
#'  strand = "+")
#' del <- Deletion(
#'   VariantID = "DEL1",
#'   GenomicPos = gr,
#'   ReferenceAllele = "G",
#'   MutatedAllele = "-",
#'   lengthDeleted = as.integer(1)
#' )
#' isFrameshift(del)
setMethod("isFrameshift", "Deletion", function(object){
    object@lengthDeleted %% 3 != 0
})
