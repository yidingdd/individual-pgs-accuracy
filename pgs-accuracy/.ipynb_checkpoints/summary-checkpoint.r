suppressPackageStartupMessages(
    {
        library(matrixStats)
        library(dplyr)
        library(readr)
        library(bigreadr)
        library(optparse)
        library(tidyverse)
        source('/u/project/pasaniuc/yiding/projects/individual-pgs-accuracy/scripts/pipeline/helper.r')
    }
)


parser <- OptionParser()
parser <- add_option(parser,
  "--pheno_name",
  action = "store", type = "character"
)
parser <- add_option(parser,
  "--pheno_file",
  action = "store", type = "character"
)

parser <- add_option(parser,
  "--test_prs",
  action = "store", type = "character"
)


parser <- add_option(parser,
    '--train_list',
    action = 'store', type = 'character', default = NULL
)
parser <- add_option(parser,
    '--test_list',
    action = 'store', type = 'character', default = NULL
)
parser <- add_option(parser,
    '--model',
    action = 'store', type = 'character', default = NULL
)
parser <- add_option(parser,
    '--covariates',
    action = 'store', type = 'character', default = NULL
)

parser <- add_option(parser,
    '--ancestry',
    action = 'store', type = 'character', default = NULL
)
parser <- add_option(parser,
  "--out_prefix",
  action = "store", type = "character"
)
parser <- parse_args(parser)
print(parser)


# get model paramters 
model <- readRDS(parser$model)
h2_seq <- sapply(model$multi_auto, function(chain) chain$h2_est)
p_seq <- sapply(model$multi_auto, function(chain) chain$p_est)
keep <- model$keep
n_sample <- ncol(model$beta_sample)
h2_est <- mean(h2_seq[keep])
p_est <- mean(p_seq[keep])
rm(model)  
                
# get residual phenotype 
train_idx <- read_tsv(parser$train_list, col_names = c('FID', 'IID'), col_types = cols()) %>% pull(IID)
test_idx <- read_tsv(parser$test_list, col_names = c('FID', 'IID'), col_types = cols()) %>% pull(IID)

pheno <- read_tsv(parser$pheno_file, col_types = cols()) 
covaraites <- read_tsv(parser$covariates, col_types = cols()) 
reg_df <- left_join(pheno, covaraites, by = c('IID' = 'IID', 'FID' = 'FID')) %>% drop_na()
reg_df_train <- reg_df %>% slice(match(train_idx, IID)) 
reg_df_test <- reg_df %>% slice(match(test_idx, IID))
                
fit <- lm(PHENO ~ .,  reg_df_train %>% select(-FID, -IID))
pheno_res_train <- residuals(fit)
scale_pheno_res <- var(pheno_res_train)
pheno_res_test <- data.frame(
    FID_IID = paste0(reg_df_test$FID, '_', reg_df_test$IID),
    pheno = reg_df_test$PHENO,
    pheno_res = reg_df_test$PHENO - predict(fit, newdata = reg_df_test)
    
)
if ('PHENO_G' %in% colnames(reg_df_test)){
    pheno_res_test$pheno_g <- reg_df_test$PHENO_G
}
        
# load ancestry
ancestry <- read_tsv(parser$ancestry, col_types = cols()) %>% 
                mutate(FID_IID = paste0(FID, '_', IID)) %>% select(-FID, -IID)

                


# summarize prs raw 
prs_df <- read_tsv(parser$test_prs, col_types = cols())
prs_mat <- prs_df %>% select(starts_with('SAMPLE_')) %>% as.matrix()


smry_df <- data.frame(
    FID_IID = prs_df$indiv,
    mean = rowMeans(prs_mat), 
    post_var = rowVars(prs_mat),
    accuracy = 
)%>%mutate(accuracy = 1 - post_var/scale_pheno_res/h2_est)
    left_join(ancestry, by = c('FID_IID' = 'FID_IID')) %>% 
    left_join(pheno_res_test, by = c('FID_IID' = 'FID_IID')) %>%
    relocate(ancestry, .after = FID_IID) 
saveRDS(smry_df, file = paste0(parser$out_prefix, '.smry.rds'))
                                
param_list <- list(
    h2_seq = h2_seq,
    p_seq = p_seq,
    keep = keep,
    n_sample = n_sample,
    h2_est = h2_est,
    p_est  = p_est,
    scale_pheno_res = scale_pheno_res,
    scale_prior_var = scale_prior_var,
    sample_size = nrow(reg_df_train)
)
saveRDS(param_list, file = paste0(parser$out_prefix, '.param.rds'))