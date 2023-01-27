library(optparse)
library(bigsnpr)
library(tibble)
library(readr)
library(dplyr)
source('helper.r')

parser <- OptionParser()
parser <- add_option(parser,
    "--ld_dir",
    action = "store", type = "character"
)
parser <- parse_args(parser)
print(parser)

file.remove(file.path(parser$ld_dir, paste0("chr", 'all', ".corr.sbk")))

for (chrom in seq(22)) {
    print(chrom)
    df_snp0 <- read_tsv(
        file.path(parser$ld_dir, paste0("chr", chrom, ".snp.tsv")),
        col_types = "icdiccd"
    )
  corr0 <- readRDS(file.path(parser$ld_dir, paste0("chr", chrom, ".corr.rds")))
  ldscore0 <- readRDS(file.path(parser$ld_dir, paste0("chr", chrom, ".ldsc.rds")))
  if (chrom == 1) {
      print(chrom)
    ld_snp <- df_snp0
    ldscore <- ldscore0
    corr <- as_SFBM(corr0, compact = TRUE, backingfile = file.path(parser$ld_dir, paste0("chr", 'all', ".corr")))
  } else {
      print(chrom)
    ld_snp <- rbind(ld_snp, df_snp0)
    ldscore <- c(ldscore, ldscore0)
    corr$add_columns(corr0, nrow(corr))
  }
}

saveRDS(
    corr,
    file.path(parser$ld_dir, paste0("chr", 'all', ".corr.rds"))
)

saveRDS(
    ldscore,
    file.path(parser$ld_dir, paste0("chr", 'all', ".ldsc.rds"))
)

write_tsv(
    ld_snp,
    file.path(parser$ld_dir, paste0("chr", 'all', ".snp.tsv"))
)



