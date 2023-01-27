library(optparse)
library(bigsnpr)
library(tibble)
library(readr)
library(dplyr)
source('helper.r')

parser <- OptionParser()
parser <- add_option(parser,
    "--train_bfile",
    action = "store", type = "character"
)
parser <- add_option(parser,
    "--train_list",
    action = "store", type = "character", default = NA
)
parser <- add_option(parser,
    "--ld_dir",
    action = "store", type = "character"
)
parser <- add_option(parser,
    "--chrom",
    action = "store", type = "integer"
)
parser <- add_option(parser,
    "--n_core",
    action = "store", type = "integer"
)


parser <- parse_args(parser)
print(parser)


# select individuals
train_fam <- read_fam(parser$train_bfile)
if(is.na(parser$train_list)){
    train_ind <- NULL
}else{
    train_ind <- as.vector(read_tsv(parser$train_list, col_names = c('IID', 'FID')) %>% pull(IID))
    train_ind <- match(train_ind, train_fam$IID)
    train_ind <- train_ind[!is.na(train_ind)]
    
}


# select SNPs
train_snps <- read_bim(parser$train_bfile) 
snp_ind <- which(train_snps$CHR == parser$chrom)

# prepare bigSNP object
train_bigSNP <- prepare_bigSNP(parser$train_bfile, inds = train_ind, snps = snp_ind, ncores = parser$n_core)

POS2 <- snp_asGeneticPos(
    train_bigSNP$map$chromosome,
    train_bigSNP$map$physical.pos
)

corr0 <- snp_cor(
    train_bigSNP$genotypes,
    infos.pos = POS2,
    size = 3 / 1000,
    ncores = parser$n_core
)

ldsc <- Matrix::colSums(corr0^2)

saveRDS(
    corr0,
    file.path(parser$ld_dir, paste0("chr", parser$chrom, ".corr.rds"))
)

saveRDS(
    ldsc,
    file.path(parser$ld_dir, paste0("chr", parser$chrom, ".ldsc.rds"))
)

write_tsv(
    bind_cols(train_snps %>% slice(snp_ind), ldsc = ldsc),
    file.path(parser$ld_dir, paste0("chr", parser$chrom, ".snp.tsv"))
)
