Here, we provide a step-by-step demonstration of how to estimate individual PGS accuracy for the 84 traits studied in our research using our pre-trained weights.

# Part 1. Prepartion

## Step 1. Install required packages
- `plink2`: https://www.cog-genomics.org/plink/2.0/
- `dask-pgen`: https://github.com/KangchengHou/dask-pgen

The `dask-pgen` package serves as a helpful wrapper built on top of plink2. It allows for efficient computation of the posterior PGS distribution using the MCMC samples stored in the weight file.

## Step 2. Download data from figshare 

```{bash}
# Make a directory to store the data
DOWNLOADED_DATA_DIR=<put your path here>
mkdir -p $DOWNLOADED_DATA_DIR
cd $DOWNLOADED_DATA_DIR

# Download data from Figshare
# The file is 11 GB, so the downloading process may take a while
https://figshare.com/ndownloader/articles/22413970

# Unzip the downloaded file
unzip 22413970
```
The downloaded data contains essential information required for computing individual PGS accuracy and genetic distance:  
- `weight-mcmc`: This folder contains 84 PGS weight files. Unlike conventional weight files that typically contain a single column of weights, the weight files we provide may contain 500-1000 (the exact number varies among traits due to QC of MCMC chains) columns of weights. Each column represents a separate MCMC sampling from the posterior distribution of genetic effects, which will be used for computing the posterior PGS distribution in subsequent steps. 
- `maf-train`: This folder contains the minor allele frequencies of PGS SNPs among the 371K UK Biobank British training individuals. 
- `param`: This folder contains the heritability and phenotye variance of the 84 traits, which are essential parameters for computing individual PGS accuracy. 
- `snp-info`: This file contains the information of SNPs used to train PGS.
- `pca`: This folder contains the `meansd.txt` and `loadings.txt` output obtained from running flashpca on training data. You can project your target genotype data to the PC space of training data and compute the genetic distance. 

## Step 3. Prepare your genotype data 
To ensure compatibility between the weight file and the genotype file in LDpred2, we strongly recommend renaming the variant IDs in your genotype file to match the IDs provided in our weight file. For your convenience, we have provided example scripts `prepare-data.sh` that can assist you in preparing your data for analysis.


# Part 2. Compute individual PGS accuracy 

## Step 1. Compute the posterior distribution of PGS
```{bash}
trait=height

dapgen score \
    --plink <path of your plink file list> \
    --freq-suffix .uk-train.afreq \
    --weights $DOWNLOADED_DATA_DIR/weight-mcmc/${trait}.auto.weight.tsv.gz \
    --weight-col-prefix SAMPLE_ \
    --out $YOUR_OUTPUT_DIR/${trait}.prs.gz \
    --chrom-col CHR --pos-col POS_GRCh37 <or 38 depending on your own genome build> \
    --alt-col EFFECT_ALLELE --ref-col OTHER_ALLELE \
    --center True \
    --threads 30
    --memory 60000
```

## Step 2. Compute individual PGS accuracy
```{bash}
Rscript summary.r \
    --pgs $YOUR_OUTPUT_DIR/${trait}.prs.gz \
    --param $DOWNLOADED_DATA_DIR/${trait}.param.txt \
    --out $YOUR_OUTPUT_DIR/${trait}.smry.txt
```

Please read scripts for details.

        
