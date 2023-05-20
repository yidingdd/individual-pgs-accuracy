#!/bin/bash

# Step 0: Setup
# Specify the path to the directory where you have stored the downloaded data from figshare
DOWNLOADED_DATA_DIR=
# Specify the directory where you have stored the processed genotypes
GENO_DIR=
# Specify the directory where you want to store your pca results
PCA_DIR=


# Step 1: Extract snps, align reference alleles and merge bed files
for CHROM in {1..22}
do
    plink2  --bfile $GENO_DIR/chr$CHROM \
        --extract $DOWNLOADED_DATA_DIR/pca/meansd.txt \
        --ref-allele $DOWNLOADED_DATA_DIR/snp-info.txt 6 2 \
        --make-bed --out $PCA_DIR/chr$CHROM
done
 
find $PCA_DIR -name "chr*.bim"  > $PCA_DIR/mergeList.txt
sed -i 's/.bim//g' $PCA_DIR/merge_list.txt $PCA_DIR/merge_list.txt
plink2 --pmerge-list $PCA_DIR/merge_list.txt bfile --make-bed --out $PCA_DIR/all


# Step 2: Project your genome to the pc space
flashpca --bfile $PCA_DIR/all --project \
    --inmeansd $DOWNLOADED_DATA_DIR/pca/meansd.txt \
    --inload $DOWNLOADED_DATA_DIR/pca/loadings.txt \
    --outproj $PCA_DIR/projections.txt -v
    
    
# Step 3: Compute genetic distance    
Rscript compute-genetic-distance.r \
    --pca $PCA_DIR/projections.txt \
    --output $PCA_DIR/genetic-distance.tsv
   
