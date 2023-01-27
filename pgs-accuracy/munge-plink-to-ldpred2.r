library(optparse)
library(readr)
library(dplyr)
library(tidyverse)

parser<- OptionParser()
parser <- add_option(parser,
    '--sumstats',
     action = 'store', type = 'character'
)
parser <- add_option(parser,
    '--type',
     action = 'store', type = 'character', default = 'continuous'
)
parser <- parse_args(parser)
print(parser)


if (parser$type == 'continuous'){
    plink_sumstats <- read_tsv(parser$sumstats,
               col_names = c("#CHROM", "POS", "ID", "REF", "ALT", "A1", "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE"),
               col_types = 'iiccccciddddc',
               skip = 1
              )

    ldpred2_sumstats <- plink_sumstats %>% select(
                                            chr = `#CHROM`,
                                            rsid = ID,
                                            pos = POS,
                                            a0 = ALT,
                                            a1 = REF,
                                            beta = BETA,
                                            beta_se = SE,
                                            n_eff = OBS_CT,
                                            p = P
    )

    ldpred2_sumstats %>% write_tsv(paste0(parser$sumstats, '.ldpred2'))
    
}

if (parser$type == 'binary'){
    plink_sumstats <- read_tsv(parser$sumstats,
                               col_names = c("#CHROM", "POS", "ID", "REF", "ALT", "A1", "FIRTH?", "TEST", "OBS_CT", "OR", "LOG(OR)_SE", "Z_STAT", "P", "ERRCODE"),
                               col_types = 'iicccccciddddc',
                               skip = 1
              )  %>% mutate(BETA = log(OR))

    ldpred2_sumstats <- plink_sumstats %>% select(
                                            chr = `#CHROM`,
                                            rsid = ID,
                                            pos = POS,
                                            a0 = ALT,
                                            a1 = REF,
                                            beta = BETA,
                                            beta_se = `LOG(OR)_SE`,
                                            n_eff = OBS_CT,
                                            p = P
    )

    ldpred2_sumstats %>% write_tsv(paste0(parser$sumstats, '.ldpred2'))
}



