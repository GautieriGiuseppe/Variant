library(testthat)
library(Variant)

# Load the example data
data("example")

test_that("example variants are valid SNV/DNV/ONV/Insertion/Deletion objects", {
  data("example", package = "Variant")
  list2env(example, envir = environment())
  
  expect_s4_class(snv1, "SNV")
  expect_s4_class(dnv1, "DNV")
  expect_s4_class(onv1, "ONV")
  expect_s4_class(insertion1, "Insertion")
  expect_s4_class(deletion1, "Deletion")
  
  expect_true(isTransition(snv1))
  expect_false(isTransition(snv2))
  expect_equal(VariantID(dnv2), "DNV002")
  expect_match(Variation(snv3), "reference allele")
  expect_true(isFrameshift(insertion1))
  expect_false(isFrameshift(insertion2))
  expect_equal(GeneSymbol(onv1), "BCR")
  expect_equal(GeneType(deletion1), "protein_coding")
})
