#!/bin/bash

module load R/4.1.0-foss-2021a
module load Anaconda3/2020.11
eval "$(conda shell.bash hook)"
conda activate betalactamases

echo "Running on:"
echo $MODULE_CPU_TYPE
echo "at "`date`

export PYTHONPATH=/well/shaw/users/amu125/miniconda/$MODULE_CPU_TYPE/envs/betalactamases/bin/:${PYTHONPATH}
export PATH=/well/shaw/users/amu125/miniconda/$MODULE_CPU_TYPE/envs/betalactamases/bin/:${PATH}



#snakemake -R --until download_CARD_DB --cores 1 --rerun-incomplete
#snakemake -R --until download_CARD_DB --cores 1 --rerun-incomplete
snakemake --cores 2 --rerun-incomplete
