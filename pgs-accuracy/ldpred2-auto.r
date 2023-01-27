library(optparse)
library(bigsnpr)
library(tibble)
library(readr)
library(dplyr)
source('helper.r')

parser <- OptionParser()

# parameters shared by auto and grid
parser <- add_option(parser,
  "--model",
  action = "store", type = "character"
)
parser <- add_option(parser,
  "--ld_dir",
  action = "store", type = "character"
)
parser <- add_option(parser,
  "--train_bfile",
  action = "store", type = "character"
)
parser <- add_option(parser,
  "--train_sumstats",
  action = "store", type = "character"
)
parser <- add_option(parser, "--n_core",
  action = "store", type = "integer"
)

parser <- add_option(parser, "--n_burn_in",
  action = "store", type = "integer", default = 100,
  help = "iterations to be discarded in the beginning of MCMC chains"
)
parser <- add_option(parser, "--n_iter",
  action = "store", type = "integer", default = 500,
  help = "iterations to be kept in the remaining of MCMC chains"
)
parser <- add_option(parser, "--out_prefix",
  action = "store", type = "character"
)

# parameters specific to auto
parser <- add_option(parser, "--report_step",
  action = "store", type = "integer", default = 5,
  help = "thinning invertals of MCMC chains in snp_ldpred2_auto"
)
parser <- add_option(parser, "--n_chain",
  action = "store", type = "integer", default = 5,
  help = "number of mini chains in snp_ldpred2_auto"
)


parser <- parse_args(parser)
print(parser)

# ------------------------------------------------------------
# train the model
# ------------------------------------------------------------

ld_snp <- read_tsv(
    file.path(parser$ld_dir, "chrall.snp.tsv"),
    col_types = "icdiccd"
  )

ldscore <- readRDS(file.path(parser$ld_dir, "chrall.ldsc.rds"))

corr <- readRDS(file.path(parser$ld_dir, "chrall.corr.rds"))

# load train sumstats
train_sumstats <- read_tsv(parser$train_sumstats, col_types = "iciccdddd") 


# match LD and train_sumstats
ld_map <- ld_snp %>% select(chr = CHR, pos = BP, a0 = A1, a1 = A2)
df_beta <- snp_match(train_sumstats, ld_map, strand_flip = F)

# run ldsc to get h2 estimation and parameters grid for tuning
ldsc <- with(df_beta, snp_ldsc(ldscore, length(ldscore),
  chi2 = (beta / beta_se)^2,
  sample_size = n_eff, blocks = NULL
))

h2_est <- ldsc[["h2"]]


# ---------------------------------------------------------------------
# LDpred2 auto
# ---------------------------------------------------------------------

# train model with ldpred2 auto
multi_auto <- snp_ldpred2_auto(
  corr, 
  df_beta, 
  h2_init = h2_est,
  burn_in = parser$n_burn_in, 
  num_iter = parser$n_iter,
  report_step = parser$report_step,
  vec_p_init = seq_log(1e-4, 0.9, parser$n_chain),
  sparse = TRUE, 
  ncores = parser$n_core,
)

# rescale weights to allelic scale, the default scale of snp_ldpred2_auto output is:
# if sumstats is on std scale, then beta is std scale, beta_sample is std scale 
# if sumstats is on allelic scale, then beta is allelic scale, beta_sample is std scale
effect_scales <- with(df_beta, sqrt(n_eff * beta_se^2 + beta^2))  

for (chain_i in 1:length(multi_auto)){
    multi_auto[[chain_i]]$sample_beta <- sweep(multi_auto[[chain_i]]$sample_beta, 1, effect_scales, '*')
}

# mini chain quality control
all_h2 <- sapply(multi_auto, function(auto) auto$h2_est)
h2 <- median(all_h2)
keep <- between(all_h2, 0.7 * h2, 1.4 * h2)
all_p <- sapply(multi_auto, function(auto) auto$p_est)
p <- median(all_p[keep])
keep <- keep & between(all_p, 0.5 * p, 2 * p)
beta_sample <- do.call( 'cbind', lapply(multi_auto[keep], function(auto) auto$sample_beta))
beta <- rowMeans(beta_sample)


# model info
model <- list(
  h2_keep = all_h2[keep],
  p_keep = all_p[keep],
  h2_est = h2_est,
  df_beta = df_beta,
  keep = keep,
  beta_sample = beta_sample,
  scale = effect_scales,
  multi_auto = multi_auto
)

                                        
# save model                                     
saveRDS(model, paste0(parser$out_prefix, ".rds"))
rm(model)                                      
                                        

# write weights
weights <- beta_sample
colnames(weights) <- paste0("SAMPLE_", seq(1, ncol(weights)))
                                        
liftOver <- runonce::download_file("http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver", "tmp-data")
bigsnpr:::make_executable(liftOver)
df_beta$pos_hg38 <- snp_modifyBuild(df_beta, liftOver, from = "hg19", to = "hg38")$pos
weights <- bind_cols(
  df_beta %>% select(
    CHR = chr, POS = pos, A1 = a0, A2 = a1, SNP = rsid, pos_hg38
  ),
  as_tibble(weights)
)
write_tsv(weights, paste0(parser$out_prefix, ".post.weight.tsv.gz"))