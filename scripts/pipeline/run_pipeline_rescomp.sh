#!/bin/bash

#module load R/4.2.0-foss-2021b 
echo "at start of script, python is"
which python


echo "Running on:"
echo $MODULE_CPU_TYPE
echo "at "`date`

# Load different R module 
if [ $MODULE_CPU_TYPE == "skylake" ]; then
        module load R/4.2.0-foss-2021b
fi

if [ $MODULE_CPU_TYPE == "ivybridge" ]; then
        module load R/4.2.1-foss-2022a
fi

module load Anaconda3/2020.11
eval "$(conda shell.bash hook)"
conda deactivate
echo "after conda deactivate"
which python

conda activate betalactamases
echo "after conda activate betalactamases"
which python

echo $PYTHONPATH
echo $PATH
which python

# Crucial: add the correct python path
# with blast location as well
export PATH=/apps/well/ncbi-blast+/2.7.1/bin/:/well/shaw/users/amu125/miniconda/$MODULE_CPU_TYPE/envs/betalactamases/bin/:${PATH}
which python
# Crucial: add the location of the 'betalacatamases' site packages to the PYTHONPATH
export PYTHONPATH=/well/shaw/users/amu125/miniconda/$MODULE_CPU_TYPE/envs/betalactamases/lib/python3.11/site-packages/:$PATH  


#snakemake -R --until download_CARD_DB --cores 1 --rerun-incomplete
#snakemake -R --until download_CARD_DB --cores 1 --rerun-incomplete
snakemake --cores 2 --rerun-incomplete --configfile cluster_config.yaml
