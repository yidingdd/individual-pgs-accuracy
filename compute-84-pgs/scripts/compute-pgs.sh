#!/bin/bash
# Step 0: Setup
# Specify the path to the directory where you have stored the downloaded data from figshare
DOWNLOADED_DATA_DIR=/opt/genomics/workingdir/yiding/projects/github-replicate
# Specify the directory where you have store the processed genotypes
GENO_DIR=/opt/genomics/workingdir/yiding/projects/github-replicate/geno
# Specify the path to the directory where you want to store pgs and its accuracy
PGS_DIR=/opt/genomics/workingdir/yiding/projects/github-replicate/pgs
mkdir -p $PGS_DIR
# Specify the trait, e.g. height
trait=height

# Step 1: Compute pgs posterior samples
dapgen score \
        --plink $GENO_DIR/all_bed.txt \
        --freq-suffix .uk-train.afreq \
        --weights $DOWNLOADED_DATA_DIR/weight-mcmc/$trait.auto.weight.tsv.gz \
        --weight-col-prefix SAMPLE_ \
        --out $PGS_DIR/${trait}.prs.tsv.gz \
        --chrom-col CHR --pos-col POS_GRCh38 \
        --alt-col EFFECT_ALLELE --ref-col OTHER_ALLELE \
        --center True \
        --threads 30 \
        --memory 60000 >> $PGS_DIR/${trait}.log
 
# Step 2: Compute pgs estimate, pgs variance, accuracy and 90% credible interval
Rscript summary.r \
    --pgs $PGS_DIR/${trait}.prs.tsv.gz \
    --param $DOWNLOADED_DATA_DIR/param/${trait}.param.txt \
    --output $PGS_DIR/${trait}.smry.tsv \
    

