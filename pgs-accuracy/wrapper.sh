#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -l h_data=2G,highp,h_rt=1:00:00 -pe shared 2
#$ -o logs/run-ldpred2.$JOB_ID


export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

. /u/local/Modules/default/init/modules.sh
module load R/4.1.0-DS
module load python/3.9.6


CONFIG_FILE=$1
echo $CONFIG_FILE

###############################
######## plink ################
###############################
# get parameters
plink=$(cat $CONFIG_FILE | shyaml get-value script.plink)
out_dir=$(cat $CONFIG_FILE | shyaml get-value output.out)
pheno_name=$(cat $CONFIG_FILE | shyaml get-value data.pheno_name)
pheno_file=$(cat $CONFIG_FILE | shyaml get-value data.train.pheno)
covar_file=$(cat $CONFIG_FILE | shyaml get-value data.train.covar)
covar_name=$(cat $CONFIG_FILE | shyaml get-value data.train.covar_name)
train_list=$(cat $CONFIG_FILE | shyaml get-value data.train.list)
mkdir -p $out_dir

# plink
for CHROM in {1..22}
do
    train_bfile_per_chrom=$( cat $CONFIG_FILE  | shyaml get-value data.train.bfile_per_chrom_tmp | sed "s/CHROM/${CHROM}/" )
    $plink \
        --1 \
        --bfile ${train_bfile_per_chrom}  \
        --covar $covar_file  --covar-name $covar_name --covar-variance-standardize \
        --glm hide-covar omit-ref allow-no-covars \
        --out $out_dir/chr${CHROM} \
        --pheno $pheno_file \
        --pheno-name $pheno_name \
        --keep $train_list
done


# merge plink
rm $out_dir/chrall.$pheno_name.glm.linear
cat $out_dir/chr1.$pheno_name.glm.linear >>  $out_dir/chrall.$pheno_name.glm.linear
for CHROM in {2..22}
do 
    tail -n +2 $out_dir/chr${CHROM}.$pheno_name.glm.linear >>  $out_dir/chrall.$pheno_name.glm.linear
done
# munge plink to LDpred2
Rscript munge-plink-to-ldpred2.r --sumstats $out_dir/chrall.$pheno_name.glm.linear

# clean up
for CHROM in {1..22}
do 
    rm $out_dir/chr${CHROM}.$pheno_name.glm.linear 
    rm $out_dir/chr${CHROM}.log
done


################################
######### ldpred2 ##############
################################
# get parameters
out_dir=$(cat $CONFIG_FILE | shyaml get-value output.out)
ld_dir=$(cat $CONFIG_FILE | shyaml get-value output.ld)

pheno_name=$(cat $CONFIG_FILE | shyaml get-value data.pheno_name)
pheno_file=$(cat $CONFIG_FILE | shyaml get-value data.train.pheno)
covar_file=$(cat $CONFIG_FILE | shyaml get-value data.train.covar)
covar_name=$(cat $CONFIG_FILE | shyaml get-value data.train.covar_name)
train_list=$(cat $CONFIG_FILE | shyaml get-value data.train.list)
train_bfile=$( cat $CONFIG_FILE  | shyaml get-value data.train.bfile_per_chrom_tmp | sed "s/CHROM/1/" )

model=$(cat $CONFIG_FILE | shyaml get-value ldpred2.model)
n_core=$(cat $CONFIG_FILE | shyaml get-value ldpred2.n_core)
n_burn_in=$(cat $CONFIG_FILE | shyaml get-value ldpred2.n_burn_in)
n_iter=$(cat $CONFIG_FILE | shyaml get-value ldpred2.n_iter)
report_step=$(cat $CONFIG_FILE | shyaml get-value ldpred2.report_step)
n_chain=$(cat $CONFIG_FILE | shyaml get-value ldpred2.n_chain)
# ldpred2
Rscript ldpred2-auto.r \
    --model=$model \
    --ld_dir=$ld_dir \
    --train_bfile=$train_bfile \
    --train_sumstats=$out_dir/chrall.$pheno_name.glm.linear.ldpred2 \
    --n_core=$n_core \
    --n_burn_in=$n_burn_in \
    --n_iter=$n_iter \
    --report_step=$report_step \
    --n_chain=$n_chain \
    --out_prefix=${out_dir}/${pheno_name}.${model}



################################
######### predict ##############
################################ 
module load R/4.1.0-DS

module load python/3.9.6

dapgen=/u/project/pasaniuc/yiding/downloads/dask-pgen/bin/dapgen

out_dir=$(cat $CONFIG_FILE | shyaml get-value output.out)
all_test_bfile=$(cat $CONFIG_FILE | shyaml get-value data.test.all_bed)
test_list=$(cat $CONFIG_FILE | shyaml get-value data.test.list)
pheno_name=$(cat $CONFIG_FILE | shyaml get-value data.pheno_name)
set_name=$(cat $CONFIG_FILE | shyaml get-value data.test.set_name)

model=$(cat $CONFIG_FILE | shyaml get-value ldpred2.model)
n_core=$(cat $CONFIG_FILE | shyaml get-value ldpred2.n_core)
freq_suffix=$(cat $CONFIG_FILE | shyaml get-value data.train.freq_suffix)


$dapgen score \
    --plink $all_test_bfile \
    --freq-suffix $freq_suffix \
    --weights ${out_dir}/${pheno_name}.${model}.post.weight.tsv.gz \
    --weight-col-prefix SAMPLE_ \
    --out ${out_dir}/${pheno_name}.${model}.${set_name}.dapgen.prs.gz \
    --chrom-col CHR --pos-col POS --alt-col A1 --ref-col A2 \
    --keep-fam $test_list\
    --center True \
    --threads $n_core \
    --memory 60000


################################
######### summary ##############
################################
out_dir=$(cat $CONFIG_FILE | shyaml get-value output.out)   
pheno_name=$(cat $CONFIG_FILE | shyaml get-value data.pheno_name)
pheno_file=$(cat $CONFIG_FILE | shyaml get-value data.test.pheno)
model=$(cat $CONFIG_FILE | shyaml get-value ldpred2.model)
set_name=$(cat $CONFIG_FILE | shyaml get-value data.test.set_name)
ancestry=$(cat $CONFIG_FILE | shyaml get-value data.test.ancestry)
train_list=$(cat $CONFIG_FILE | shyaml get-value data.train.list)
test_list=$(cat $CONFIG_FILE | shyaml get-value data.test.list)
covar=$(cat $CONFIG_FILE | shyaml get-value data.train.covar)

echo $ancestry
Rscript summary.r \
    --pheno_name=$pheno_name \
    --pheno_file=$pheno_file \
    --test_prs=${out_dir}/${pheno_name}.${model}.${set_name}.dapgen.prs.gz \
    --out_prefix=${out_dir}/${pheno_name}.${model} \
    --ancestry=$ancestry \
    --train_list=$train_list \
    --test_list=$test_list \
    --model=${out_dir}/${pheno_name}.${model}.rds \
    --covariates=$covar







