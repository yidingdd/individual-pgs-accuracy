#--------------------------#
# read data 
#--------------------------#
library(bigsnpr)
library(readr)
library(dplyr)

read_bim <-function(bfile){
    bim <- tryCatch(
        {
            read_delim(
                paste0(bfile, '.bim'),
                col_names = c("CHR", "SNP", "CM", "BP", "A1", "A2"),
                col_types = "icdicc",
                delim = '\t'
            )
        }, error=function(cond) {
            message("ERROR while reading bim file")
        }, warning=function(cond) {
            bim <- read_delim(
                paste0(bfile, '.bim'),
                col_names = c("CHR", "SNP", "CM", "BP", "A1", "A2"),
                col_types = "icdicc",
                delim = ' '
            )
            return(bim)
        }
    )
    return(bim)
}

read_fam <-function(bfile){
    fam <- tryCatch(
        {
            read_delim(
                paste0(bfile, '.fam'),
                col_names = c("FID", "IID", "PID", "MID", "SEX", "PHENO"),
                col_types = "ccccii",
                delim = '\t'
            )
        }, error=function(cond) {
            message("ERROR while reading fam file")
        }, warning=function(cond) {
            fam <- read_delim(
                paste0(bfile, '.fam'),
                col_names = c("FID", "IID", "PID", "MID", "SEX", "PHENO"),
                col_types = "ccccii",
                delim = ' '
            )
            return(fam)
        }
    )
    return(fam)
}


prepare_bigSNP <- function(bfile, inds = NULL, snps = NULL, ncores = 1){
    
    if(is.null(inds) & is.null(snps)){
        rds_path <- snp_readBed2( bedfile = paste0(bfile, '.bed'), 
                                 backingfile = tempfile(), ncores = ncores)
    }else if(is.null(inds) & !is.null(snps)){
        rds_path <- snp_readBed2( bedfile = paste0(bfile, '.bed'), 
                                 backingfile = tempfile(), ncores = ncores, ind.col = snps)
    }else if(!is.null(inds) & is.null(snps)){
        rds_path <- snp_readBed2( bedfile = paste0(bfile, '.bed'), 
                                 backingfile = tempfile(), ncores = ncores, ind.row = inds)
    }else{
        rds_path <- snp_readBed2( bedfile = paste0(bfile, '.bed'), 
                                 backingfile = tempfile(), ncores = ncores, ind.row = inds, ind.col = snps)
    }
        
    bigSNP <- snp_attach(rds_path)
    G <- snp_fastImputeSimple(
        bigSNP$genotypes,
        method = "mean2", 
        ncores = ncores
    )
    bigSNP$genotypes <- G
    rm(G)

    return(bigSNP)
}

