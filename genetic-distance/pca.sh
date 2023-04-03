#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -l h_data=6G,h_rt=5:00:00 -pe shared 10
#$ -o logs/pca.$JOB_ID

plink2=/u/project/pasaniuc/yiding/downloads/plink2
flashpca=/u/project/pasaniuc/yiding/projects/subcontinental_prs_uncertainty/prs-uncertainty-experiment/pca/flashpca_x86-64
# wget https://github.com/gabraham/flashpca/blob/master/exclusion_regions_hg19.txt
long_ld_region=/u/project/pasaniuc/yiding/projects/subcontinental_prs_uncertainty/prs-uncertainty-experiment/pca/exclusion_regions_hg19.txt 

pca_dir=/u/scratch/y/yiding/prs/ukbb-atlas/pca
rm -r $pca_dir
mkdir -p $pca_dir/geno
cd $pca_dir

train_list=/u/scratch/y/yiding/prs/ukbb-atlas/meta/train.list
test_list=/u/scratch/y/yiding/prs/ukbb-atlas/meta/test.list

for CHROM in {1..22}
do
    train_bfile_per_chrom=/u/scratch/y/yiding/prs/ukbb-atlas/geno/train.chr${CHROM}
    test_bfile_per_chrom=/u/scratch/y/yiding/prs/ukbb-atlas/geno/test.chr${CHROM}

    $plink2 --bfile $train_bfile_per_chrom \
        --geno 0.1 \
        --indep-pairwise 1000 50 0.05 \
        --exclude range $long_ld_region \
        --out $pca_dir/geno/chr$CHROM

    $plink2 --bfile $train_bfile_per_chrom \
        --extract $pca_dir/geno/chr$CHROM.prune.in \
        --make-bed \
        --out $pca_dir/geno/train.chr${CHROM}

    $plink2 --bfile $test_bfile_per_chrom\
        --extract $pca_dir/geno/chr$CHROM.prune.in \
        --make-bed \
        --out $pca_dir/geno/test.chr${CHROM}
done


# step merge list and merge 
find $pca_dir/geno -name "train.chr*.bim"  > $pca_dir/geno/train.merge.list
sed -i 's/.bim//g' $pca_dir/geno/train.merge.list $pca_dir/geno/train.merge.list   
$plink2 --pmerge-list $pca_dir/geno/train.merge.list bfile --make-bed --out $pca_dir/geno/train

find $pca_dir/geno -name "test.chr*.bim"  > $pca_dir/geno/test.merge.list
sed -i 's/.bim//g' $pca_dir/geno/test.merge.list $pca_dir/geno/test.merge.list   
$plink2 --pmerge-list $pca_dir/geno/test.merge.list bfile --make-bed --out $pca_dir/geno/test


$flashpca  --bfile $pca_dir/geno/train -d 20 --outload $pca_dir/loadings.txt --outmeansd $pca_dir/meansd.txt
$flashpca --bfile $pca_dir/geno/test --project --inmeansd $pca_dir/meansd.txt \
   --outproj $pca_dir/projections.txt --inload $pca_dir/loadings.txt -v
