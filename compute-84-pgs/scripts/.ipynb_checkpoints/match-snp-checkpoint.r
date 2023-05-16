suppressPackageStartupMessages({
  library(readr)
  library(glue)
  library(tidyverse)
  library(optparse)
})

# Create an option parser
parser <- OptionParser()
parser <- add_option(parser, "--snp_info", dest = "snp_info", help = "Path to snp-info.txt")
parser <- add_option(parser, "--target_genome_build", dest = "target_genome_build", help = "Target genome build, 37 or 38")
parser <- add_option(parser, "--input_bed_temp", dest = "input_bed_temp", help = "Path to your target per chrom bed file, use {CHROM} to replace actual chromosome number and please don't add suffix")
parser <- add_option(parser, "--output_dir", dest = "output_dir", help = "Output directory for the prepared genotype")

# Parse the command-line arguments
parser <- parse_args(parser)

# Get target genome
target_genome_build <- paste0('POS_GRCh', parser$target_genome_build)

# Read SNP info file
snp_info <- read_tsv(parser$snp_info, col_types = 'ccddcc')
head(snp_info)

# Read in your bed file and intersect the SNPs from two files
mapped_snps_list <- list()

for (CHROM in 1:22) {
  mapped_snps_list[[CHROM]] <- read_tsv(glue(parser$input_bed_temp, ".bim"),
                                        col_types = 'ccddcc',
                                        col_names = c('CHR', 'SNP', 'CM', 'BP', 'OTHER_ALLELE', 'EFFECT_ALLELE')) %>%
    rename(!!target_genome_build := BP) %>%
    select(CHR, SNP, !!target_genome_build) %>%
    left_join(snp_info,
              by = c('CHR', target_genome_build),
              suffix = c('_target', '_weight')) %>%
    drop_na() %>%
    select('SNP_target', 'SNP_weight')
}

mapped_snps <- bind_rows(mapped_snps_list)

# Output mapped SNPs
write_tsv(mapped_snps, file.path(parser$output_dir, 'matched_snps.tsv'))