#!/bin/bash

module load R/4.1.0-foss-2021a
module load Anaconda3/2020.11
eval "$(conda shell.bash hook)"
conda activate betalactamases

#snakemake -R --until download_CARD_DB --cores 1 --rerun-incomplete
#snakemake -R --until download_CARD_DB --cores 1 --rerun-incomplete
snakemake --cores 2 --rerun-incomplete
