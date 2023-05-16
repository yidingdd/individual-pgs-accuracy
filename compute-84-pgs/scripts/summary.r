suppressPackageStartupMessages({
  library(readr)
  library(glue)
  library(tidyverse)
  library(optparse)
  library(matrixStats)
})

# Create an option parser
parser <- OptionParser()
parser <- add_option(parser, "--pgs", dest = "pgs", help = "Path to the raw pgs file output by dapgen")
parser <- add_option(parser, "--param", dest = "param", help = "Path to heritability parameter files")
parser <- add_option(parser, "--output", dest = "output", help = "Path to output file")

# Parse the command-line arguments
parser <- parse_args(parser)

# Read in the parameters
param <- read_tsv(parser$param, col_types = 'dd')

# Read the posterior sampling of PGS
pgs_df <- read_tsv(parser$pgs, col_types = cols())
pgs_mat <- pgs_df %>% select(starts_with('SAMPLE_')) %>% as.matrix()

# Compute pgs estimate, pgs variance, accuracy and 90% credible interval
smry_df <- data.frame(
  FID_IID = pgs_df$indiv,
  pgs = rowMeans(pgs_mat),
  post_var = rowVars(pgs_mat),
  CI_90_lower = apply(pgs_mat, 1, quantile, probs = 0.05),
  CI_90_upper = apply(pgs_mat, 1, quantile, probs = 0.95)
) %>% mutate(accuracy = 1 - post_var / (param$h2 * param$pheno_var))

# Output
write_tsv(smry_df, parser$output)