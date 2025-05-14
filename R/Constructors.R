# Constructors 

# Constructor for Variant class
#' @title Constructor for Class Variant
#' @description Function to create object of class Variant.
#' @param VariantID Character. Unique ID of the variant.
#' @param GenomicPos \code{GRanges} object for the position of the variant.
#' @param ReferenceAllele Character. The allele in the reference genome.
#' @param MutatedAllele Character. The alternative or mutated allele.
#' @param Exon Character. The exon where the variant is located.
#' @param GeneSymbol Character. Gene name (optional).
#' @param GeneType Character. Type of gene e.g. "protein_coding" (optional).
#' @param MutationConsequence Character. Effect in the phenotype (optional).
#' @param ClinicalRel Character. Clinical annotation (optional).
#' @param EnsemblID Character. Transcript from Ensembl (optional).
#' @param ProteinVariation Character. HGVS protein annotation (optional).
#' @return class \code{Variant} object
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100),
#'  strand = "+")
#' v <- Variant(
#'   VariantID = "v1",
#'   GenomicPos = gr,
#'   ReferenceAllele = "A",
#'   MutatedAllele = "T"
#' )
#' VariantID(v)
#' @export
Variant <- function(VariantID,
                    GenomicPos,
                    ReferenceAllele,
                    MutatedAllele,
                    Exon = NA_character_,
                    GeneSymbol = NA_character_,
                    GeneType = NA_character_,
                    MutationConsequence = NA_character_,
                    ClinicalRel = NA_character_,
                    EnsemblID = NA_character_,
                    ProteinVariation = NA_character_) {
  new("Variant",
      VariantID = VariantID,
      GenomicPos = GenomicPos,
      ReferenceAllele = ReferenceAllele,
      MutatedAllele = MutatedAllele,
      Exon = Exon,
      GeneSymbol = GeneSymbol,
      GeneType = GeneType,
      MutationConsequence = MutationConsequence,
      ClinicalRel = ClinicalRel,
      EnsemblID = EnsemblID,
      ProteinVariation = ProteinVariation)
}

# Constructor for SNV class
#' @title Constructor for class SNV
#' @description Constructor to create object of class SNV
#' @param VariantID Character. Unique ID of the variant.
#' @param GenomicPos \code{GRanges} object for the position of the variant.
#' @param ReferenceAllele Character. The allele in the reference genome.
#' @param MutatedAllele Character. The alternative or mutated allele.
#' @param transition Logical. TRUE when transition (A<->G, C<->T).
#' @param Exon Character. The exon where the variant is located.
#' @param GeneSymbol Character. Gene name (optional).
#' @param GeneType Character. Type of gene e.g. "protein_coding" (optional).
#' @param MutationConsequence Character. Effect in the phenotype (optional).
#' @param ClinicalRel Character. Clinical annotation (optional).
#' @param EnsemblID Character. Transcript from Ensembl (optional).
#' @param ProteinVariation Character. HGVS protein annotation (optional).
#' @return class \code{SNV} object
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
SNV <- function(VariantID,
                GenomicPos,
                ReferenceAllele,
                MutatedAllele,
                transition = NA,
                Exon = NA_character_,
                GeneSymbol = NA_character_,
                GeneType = NA_character_,
                MutationConsequence = NA_character_,
                ClinicalRel = NA_character_,
                EnsemblID = NA_character_,
                ProteinVariation = NA_character_) {
  
  new("SNV",
      VariantID = VariantID,
      GenomicPos = GenomicPos,
      ReferenceAllele = ReferenceAllele,
      MutatedAllele = MutatedAllele,
      Exon = Exon,
      GeneSymbol = GeneSymbol,
      GeneType = GeneType,
      MutationConsequence = MutationConsequence,
      ClinicalRel = ClinicalRel,
      EnsemblID = EnsemblID,
      ProteinVariation = ProteinVariation,
      transition = transition)
}

# Constructor for DNV class
#' @title Constructor for class DNV
#' @description Constructor to create an object of class DNV
#' @param VariantID Character. Unique ID of the variant.
#' @param GenomicPos \code{GRanges} object for the position of the variant.
#' @param ReferenceAllele Character. The allele in the reference genome.
#' @param MutatedAllele Character. The alternative or mutated allele.
#' @param validated Logical. TRUE if the de novo is validated.
#' @param Exon Character. The exon where the variant is located.
#' @param GeneSymbol Character. Gene name (optional).
#' @param GeneType Character. Type of gene e.g. "protein_coding" (optional).
#' @param MutationConsequence Character. Effect in the phenotype (optional).
#' @param ClinicalRel Character. Clinical annotation (optional).
#' @param EnsemblID Character. Transcript from Ensembl (optional).
#' @param ProteinVariation Character. HGVS protein change (optional).
#' @return class \code{DNV} object
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
DNV <- function(VariantID,
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
                ProteinVariation = NA_character_) {
  
  new("DNV",
      VariantID = VariantID,
      GenomicPos = GenomicPos,
      ReferenceAllele = ReferenceAllele,
      MutatedAllele = MutatedAllele,
      validated = validated,
      Exon = Exon,
      GeneSymbol = GeneSymbol,
      GeneType = GeneType,
      MutationConsequence = MutationConsequence,
      ClinicalRel = ClinicalRel,
      EnsemblID = EnsemblID,
      ProteinVariation = ProteinVariation)
}

# Constructor for ONV class
#' @title Constructor for class ONV
#' @description Constructor to create an object of class ONV
#' @param VariantID Character. Unique ID of the variant.
#' @param GenomicPos \code{GRanges} object for the position of the variant.
#' @param ReferenceAllele Character. The allele in the reference genome.
#' @param MutatedAllele Character. The alternative or mutated allele.
#' @param complexity Character. The type of complex mutation.
#' @param Exon Character. The exon where the variant is located.
#' @param GeneSymbol Character. Gene name (optional).
#' @param GeneType Character. Type of gene e.g. "protein_coding" (optional).
#' @param MutationConsequence Character. Effect in the phenotype (optional).
#' @param ClinicalRel Character. Clinical annotation (optional).
#' @param EnsemblID Character. Transcript from Ensembl (optional).
#' @param ProteinVariation Character. HGVS protein annotation (optional).
#' @return class \code{ONV} object
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
ONV <- function(VariantID,
                GenomicPos,
                ReferenceAllele,
                MutatedAllele,
                complexity = NA_character_,
                Exon = NA_character_,
                GeneSymbol = NA_character_,
                GeneType = NA_character_,
                MutationConsequence = NA_character_,
                ClinicalRel = NA_character_,
                EnsemblID = NA_character_,
                ProteinVariation = NA_character_) {
  
  new("ONV",
      VariantID = VariantID,
      GenomicPos = GenomicPos,
      ReferenceAllele = ReferenceAllele,
      MutatedAllele = MutatedAllele,
      complexity = complexity,
      Exon = Exon,
      GeneSymbol = GeneSymbol,
      GeneType = GeneType,
      MutationConsequence = MutationConsequence,
      ClinicalRel = ClinicalRel,
      EnsemblID = EnsemblID,
      ProteinVariation = ProteinVariation)
}

# Constructor for Insertion class
#' @title Constructor for class Insertion
#' @description Constructor to create an object of class Insertion
#' @param VariantID Character. Unique ID of the variant.
#' @param GenomicPos \code{GRanges} object for the position of the variant.
#' @param ReferenceAllele Character. The allele in the reference genome.
#' @param MutatedAllele Character. The alternative or mutated allele.
#' @param lengthInserted Integer. The length of inserted sequence.
#' @param Exon Character. The exon where the variant is located.
#' @param GeneSymbol Character. Gene name (optional).
#' @param GeneType Character. Type of gene e.g. "protein_coding" (optional).
#' @param MutationConsequence Character. Effect in the phenotype (optional).
#' @param ClinicalRel Character. Clinical annotation (optional).
#' @param EnsemblID Character. Transcript from Ensembl (optional).
#' @param ProteinVariation Character. HGVS protein annotation (optional).
#' @return class \code{Insertion} object
#' @examples
#' gr <- GenomicRanges::GRanges("chr13", IRanges::IRanges(32914438, 32914438),
#'  strand = "+")
#' ins <- Insertion(
#'   VariantID = "INS1",
#'   GenomicPos = gr,
#'   MutatedAllele = "G",
#'   lengthInserted = as.integer(1)
#' )
#' isFrameshift(ins)
#' @export
Insertion <- function(VariantID,
                      GenomicPos,
                      ReferenceAllele = "-",
                      MutatedAllele,
                      lengthInserted = NA_integer_,
                      Exon = NA_character_,
                      GeneSymbol = NA_character_,
                      GeneType = NA_character_,
                      MutationConsequence = NA_character_,
                      ClinicalRel = NA_character_,
                      EnsemblID = NA_character_,
                      ProteinVariation = NA_character_) {
  if (is.null(lengthInserted)){
      lengthInserted <- nchar(MutatedAllele)
      lengthInserted <- as.integer(lengthInserted)
  }
  
  new("Insertion",
      VariantID = VariantID,
      GenomicPos = GenomicPos,
      ReferenceAllele = ReferenceAllele,
      MutatedAllele = MutatedAllele,
      lengthInserted = lengthInserted,
      Exon = Exon,
      GeneSymbol = GeneSymbol,
      GeneType = GeneType,
      MutationConsequence = MutationConsequence,
      ClinicalRel = ClinicalRel,
      EnsemblID = EnsemblID,
      ProteinVariation = ProteinVariation)
}

# Constructor for Deletion class
#' @title Constructor for class Deletion
#' @description Constructor to create an object of class Deletion
#' @param VariantID Character. Unique ID of the variant.
#' @param GenomicPos \code{GRanges} object for the position of the variant.
#' @param ReferenceAllele Character. The allele in the reference genome.
#' @param MutatedAllele Character. The alternative or mutated allele.
#' @param lengthDeleted Integer. The length of deleted sequence.
#' @param Exon Character. The exon where the variant is located.
#' @param GeneSymbol Character. Gene name (optional).
#' @param GeneType Character. Type of gene e.g. "protein_coding" (optional).
#' @param MutationConsequence Character. Effect in the phenotype (optional).
#' @param ClinicalRel Character. Clinical annotation (optional).
#' @param EnsemblID Character. Transcript from Ensembl (optional).
#' @param ProteinVariation Character. HGVS protein annotation (optional).
#' @return class \code{Deletion} object
#' gr <- GenomicRanges::GRanges("chr13", IRanges::IRanges(32914438, 32914438),
#'  strand = "+")
#' del <- Deletion(
#'   VariantID = "DEL1",
#'   GenomicPos = gr,
#'   ReferenceAllele = "G",
#'   lengthDeleted = as.integer(1)
#' )
#' isFrameshift(del)
#' @export
Deletion <- function(VariantID,
                     GenomicPos,
                     ReferenceAllele,
                     MutatedAllele = "-",
                     lengthDeleted = NA_integer_,
                     Exon = NA_character_,
                     GeneSymbol = NA_character_,
                     GeneType = NA_character_,
                     MutationConsequence = NA_character_,
                     ClinicalRel = NA_character_,
                     EnsemblID = NA_character_,
                     ProteinVariation = NA_character_) {
  if (is.null(lengthDeleted)){
    lengthDeleted <- nchar(ReferenceAllele)
    lengthDeleted <- as.integer(lengthDeleted)
  }
  
  new("Deletion",
      VariantID = VariantID,
      GenomicPos = GenomicPos,
      ReferenceAllele = ReferenceAllele,
      MutatedAllele = MutatedAllele,
      lengthDeleted = lengthDeleted,
      Exon = Exon,
      GeneSymbol = GeneSymbol,
      GeneType = GeneType,
      MutationConsequence = MutationConsequence,
      ClinicalRel = ClinicalRel,
      EnsemblID = EnsemblID,
      ProteinVariation = ProteinVariation)
}