## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(minSNPs)
library(Biostrings) # optional
library(BiocParallel) # optional, but needed for parallel processing

## -----------------------------------------------------------------------------
isolates_from_biostring <- readBStringSet(
  system.file("extdata", "Chlamydia_mapped.fasta", package = "minSNPs"))
isolates_from_default <- read_fasta(
  system.file("extdata", "Chlamydia_mapped.fasta", package = "minSNPs"))
processed_from_biostring <- process_allele(isolates_from_biostring)
processed_from_default <- process_allele(isolates_from_default)

## -----------------------------------------------------------------------------
high_d_snps <- find_optimised_snps(seqc = processed_from_biostring,
  metric = "simpson", number_of_result = 1, max_depth = 1,
  included_positions = c(), excluded_positions = c())

## -----------------------------------------------------------------------------
discriminating_snps <- find_optimised_snps(seqc = processed_from_biostring,
  metric = "percent", number_of_result = 1, max_depth = 1,
  included_positions = c(), excluded_positions = c(),
  goi = c("A_D213", "H_S1432"))

## -----------------------------------------------------------------------------
cat("High D SNPs\n")
output_result(high_d_snps)
cat("SNPs discriminating against A_D213, H_S1432\n")
output_result(discriminating_snps)

## ---- eval=FALSE--------------------------------------------------------------
#  output_result(high_d_snps, view = "csv",
#    file_name = "high_d_snps.csv")
#  output_result(discriminating_snps, view = "csv",
#    file_name = "discriminating_snps.csv")

