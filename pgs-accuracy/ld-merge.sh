#!/bin/bash -l
#$ -cwd
#$ -l h_data=8G,h_rt=03:00:00 -pe shared 4
#$ -j y
#$ -o ./logs/merge-ld.$JOB_ID


export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
. /u/local/Modules/default/init/modules.sh
module load R/4.1.0-DS

CONFIG_FILE=$1

Rscript ld-merge.r \
    --ld_dir $(cat $CONFIG_FILE | shyaml get-value output.ld) 