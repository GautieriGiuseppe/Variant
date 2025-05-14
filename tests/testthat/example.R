# example data
library(GenomicRanges)
snv1 <- SNV(
  VariantID = "SNV001",
  GenomicPos = GRanges("chr7", IRanges(117199644, 117199644), strand = "+"),
  ReferenceAllele = "G",
  MutatedAllele = "A",
  transition = TRUE,
  Exon = "15",
  GeneSymbol = "CFTR",
  GeneType = "protein_coding",
  MutationConsequence = "missense_variant",
  ClinicalRel = "Pathogenic",
  EnsemblID = "ENST00000003084",
  ProteinVariation = "p.Gly551Asp"
)
snv2 <- SNV(
  VariantID = "SNV002",
  GenomicPos = GRanges("chr17", IRanges(41276045, 41276046), strand = "+"),
  ReferenceAllele = "AG",
  MutatedAllele = "G",
  transition = FALSE,
  Exon = "2",
  GeneSymbol = "BRCA1",
  GeneType = "protein_coding",
  MutationConsequence = "frameshift_variant",
  ClinicalRel = "Pathogenic",
  EnsemblID = "ENST00000357654",
  ProteinVariation = "p.Glu23Valfs"
)
snv3 <- SNV(
  VariantID = "SNV003",
  GenomicPos = GRanges("chr17", IRanges(7676159, 7676159), strand = "+"),
  ReferenceAllele = "C",
  MutatedAllele = "G",
  transition = FALSE,
  Exon = "4",
  GeneSymbol = "TP53",
  GeneType = "protein_coding",
  MutationConsequence = "missense_variant",
  ClinicalRel = "Uncertain significance",
  EnsemblID = "ENST00000269305",
  ProteinVariation = "p.Pro72Arg"
)
dnv1 <- DNV(
  VariantID = "DNV001",
  GenomicPos = GRanges("chr17", IRanges(41276045, 41276045), strand = "+"),
  ReferenceAllele = "C",
  MutatedAllele = "T",
  validated = TRUE,
  Exon = "2",
  GeneSymbol = "TP53",
  GeneType = "protein_coding",
  MutationConsequence = "nonsense_variant",
  ClinicalRel = "Likely pathogenic",
  EnsemblID = "ENST00000269305",
  ProteinVariation = "p.Arg175Ter"
)
dnv2 <- DNV(
  VariantID = "DNV002", 
  GenomicPos = GRanges("chr2", IRanges(166233456, 166233456), strand = "+"),
  ReferenceAllele = "C",
  MutatedAllele = "T",
  validated = TRUE,
  Exon = "26",
  GeneSymbol = "SCN2A",
  GeneType = "protein_coding",
  MutationConsequence = "nonsense_variant",
  ClinicalRel = "Pathogenic",
  EnsemblID = "ENST00000409713",
  ProteinVariation = "p.Arg1587*"
)
onv1 <- ONV(
  VariantID = "ONV001",
  GenomicPos = GRanges("chr22", IRanges(29091840, 29091840), strand = "+"),
  ReferenceAllele = "A",
  MutatedAllele = "G",
  complexity = "translocation",
  Exon = NA_character_,
  GeneSymbol = "BCR",
  GeneType = "protein_coding",
  MutationConsequence = "fusion_gene",
  ClinicalRel = "Pathogenic",
  EnsemblID = "ENST00000305877",
  ProteinVariation = NA_character_
)
onv2 <- ONV(
  VariantID = "ONV002",
  GenomicPos = GRanges("chr10", IRanges(12345678, 12378901), strand = "+"),
  ReferenceAllele = "-",
  MutatedAllele = "-",
  complexity = "inversion",
  Exon = "3",
  GeneSymbol = "FGFR2",
  GeneType = "protein_coding",
  MutationConsequence = "inversion",
  ClinicalRel = "Pathogenic",
  EnsemblID = "ENST00000358487",
  ProteinVariation = NA_character_
)
insertion1 <- Insertion(
  VariantID = "INS001",
  GenomicPos = GRanges("chr19", IRanges(44908822, 44908822), strand = "+"),
  ReferenceAllele = "-",
  MutatedAllele = "T",
  lengthInserted = as.integer(1),
  Exon = "3",
  GeneSymbol = "APOE",
  GeneType = "protein_coding",
  MutationConsequence = "frameshift_variant",
  ClinicalRel = "Uncertain significance",
  EnsemblID = "ENST00000252486",
  ProteinVariation = "p.Leu28fs"
)
insertion2 <- Insertion(
  VariantID = "INS002",
  GenomicPos = GRanges("chr4", IRanges(3076600, 3076600), strand = "+"),
  MutatedAllele = "CAG",
  lengthInserted = as.integer(3),
  Exon = "1",
  GeneSymbol = "HTT",
  GeneType = "protein_coding",
  MutationConsequence = "frameshift_variant",
  ClinicalRel = "Pathogenic",
  EnsemblID = "ENST00000355072",
  ProteinVariation = "p.Gln18fs"
)
deletion1 <- Deletion(
  VariantID = "DEL001",
  GenomicPos = GRanges("chr13", IRanges(32914438, 32914438), strand = "+"),
  ReferenceAllele = "G",
  MutatedAllele = "-",
  lengthDeleted = as.integer(1),
  Exon = "7",
  GeneSymbol = "BRCA2",
  GeneType = "protein_coding",
  MutationConsequence = "frameshift_variant",
  ClinicalRel = "Pathogenic",
  EnsemblID = "ENST00000380152",
  ProteinVariation = "p.Glu23fs"
)
deletion2 <- Deletion(
  VariantID = "DEL002",
  GenomicPos = GRanges("chr13", IRanges(32914438, 32914438), strand = "+"),
  ReferenceAllele = "G",
  lengthDeleted = as.integer(1),
  Exon = "2",
  GeneSymbol = "BRCA2",
  GeneType = "protein_coding",
  MutationConsequence = "frameshift_variant",
  ClinicalRel = "Pathogenic",
  EnsemblID = "ENST00000544455",
  ProteinVariation = "p.Glu23fs"
)
example <- list(snv1 = snv1, snv2 = snv2, snv3 = snv3, dnv1 = dnv1,
                         dnv2 = dnv2, onv1 = onv1, onv2 = onv2, 
                         insertion1 = insertion1, insertion2 = insertion2,
                         deletion1 = deletion1, deletion2 = deletion2)
