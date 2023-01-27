#!/bin/bash -l
#$ -cwd
#$ -l h_data=8G,h_rt=23:00:00 -pe shared 8
#$ -j y
#$ -o ./logs/ld.$JOB_ID.$TASK_ID
#$ -t 1-22

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
. /u/local/Modules/default/init/modules.sh
module load R/4.1.0-DS

CONFIG_FILE=$1
CHROM=${SGE_TASK_ID}
outdir=$(cat $CONFIG_FILE | shyaml get-value output.ld)
train_list=$(cat $CONFIG_FILE  | shyaml get-value data.train.list)
train_bfile=$( cat $CONFIG_FILE  | shyaml get-value data.train.bfile_per_chrom_tmp | sed "s/CHROM/${CHROM}/" )

echo $outdir
echo $train_list
echo $train_bfile
mkdir -p ${outdir}

Rscript ld.r \
    --train_bfile $train_bfile \
    --train_list $train_list \
    --chrom ${CHROM} \
    --ld_dir ${outdir} \
    --n_core 8 
    
    
    
