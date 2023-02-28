#!/bin/bash

#module load R/4.2.0-foss-2021b 
echo "aat start of script"
which python

module load Anaconda3/2020.11
eval "$(conda shell.bash hook)"
conda deactivate
echo "after conda deactivate"
which python

conda activate betalactamases
echo "after conda activate betalactamases"
which python
echo "Running on:"

echo $MODULE_CPU_TYPE
echo "at "`date`

if [ $MODULE_CPU_TYPE == "skylake" ]; then
        module load R/4.2.0-foss-2021b
fi

if [ $MODULE_CPU_TYPE == "ivybridge" ]; then
        module load R/4.2.1-foss-2022a
fi

#module load R/4.2.0-foss-2021b 
#export PYTHONPATH=/well/shaw/users/amu125/miniconda/$MODULE_CPU_TYPE/envs/betalactamases/bin/:${PYTHONPATH}
#export PATH=/well/shaw/users/amu125/miniconda/$MODULE_CPU_TYPE/envs/betalactamases/bin/:${PATH}
#which python
#module load Anaconda3/2020.11
#eval "$(conda shell.bash hook)"

echo $PYTHONPATH
echo $PATH
which python

# 
#export PATH=/well/shaw/users/amu125/programs/ncbi-blast-2.13.0+/bin:$PATH
#echo "after add blast to path"
#which python

# Julia path
#export PATH=$HOME/bin:/usr/local/bin:/well/shaw/users/amu125/programs/julia-1.7.2/bin:/well/shaw/users/amu125/programs/pangraph/bin:$PATH
#export JULIA_DEPOT_PATH="/well/shaw/users/amu125/programs/julia-1.7.2/depot:$JULIA_DEPOT_PATH"

#echo "after julia paths"
#which python

#snakemake -R --until download_CARD_DB --cores 1 --rerun-incomplete
#snakemake -R --until download_CARD_DB --cores 1 --rerun-incomplete
snakemake --cores 2 --rerun-incomplete
