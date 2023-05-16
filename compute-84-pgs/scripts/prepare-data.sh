#!/bin/bash

# Step 0: Setup
# Specify the target genome build (e.g., 38 for GRCh38 or 37 for GRCh37)
TARGET_GENOME_BUILD=
# Specify the template for your per chromosome bed file path, e.g. DIR/atlas.chr{CHROM},
# where {CHROM} will be replaced with the actual chromosome number by the script
INPUT_PER_CHROM_BED=
# Specify the path to the directory where you have stored the downloaded data from figshare
DOWNLOADED_DATA_DIR=
# Specify the directory where you want to store the processed genotypes
GENO_DIR=



# Step 1: Match SNPs in the weight file and your genotype file
# This step will output a txt file matched_snps.tsv,
# where the first column is the ID in your genotype file
# and the second column is the ID in the weight file
Rscript match.r \
    --snp_info $DOWNLOADED_DATA_DIR/snp-info.tsv \
    --target_genome_build $TARGET_GENOME_BUILD \
    --input_bed_temp $INPUT_PER_CHROM_BED \
    --GENO_DIR $GENO_DIR

# Step 2: Extract SNPs and update SNP ID
rm -f $GENO_DIR/all_bed.txt
for CHROM in {1..22}
do
  # Replace {CHROM} with the value of the CHROM variable
  bfile=${INPUT_PER_CHROM_BED/\{CHROM\}/$CHROM}

  plink2 --bfile $bfile \
      --extract $GENO_DIR/matched_snps.tsv \
      --update-name $GENO_DIR/matched_snps.tsv 2 1 \
      --make-bed --out $GENO_DIR/chr$CHROM

  echo $GENO_DIR/chr$CHROM.bed >> $GENO_DIR/all_bed.txt
done

# Step 3: Copy training allele frequency to the geno folder
cp $DOWNLOADED_DATA_DIR/maf-train/* $GENO_DIR