suppressPackageStartupMessages({
  library(readr)
  library(tidyverse)
  library(optparse)
  library(matrixStats)
})

# Create an option parser
parser <- OptionParser()
parser <- add_option(parser, "--pca", dest = "pca", help = "Path to the PCA projection")
parser <- add_option(parser, "--output", dest = "output", help = "Path to output distance file")

# Parse the command-line arguments
parser <- parse_args(parser)

# Read in the pcs
pca_df <- read_tsv(parser$pca, col_types = cols())

pca_mat <- pca_df %>%
  select(starts_with('PC')) %>%
  as.matrix()

# Compute genetic distance 
dist_df <- data.frame(
  FID_IID = paste0(pca_df$FID, '_', pca_df$IID),
  dist = sqrt(rowSums(pca_mat ** 2))
)

# Output
write_tsv(dist_df, parser$output)