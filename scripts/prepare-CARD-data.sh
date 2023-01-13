#!/bin/bash

Rscript prepare-01_summarise-CARD-data.R 

python prepare-02_make-matrix-CARD-data.py 

source prepare-03a_align-families.sh 

python prepare-03b_make-heatmap.py 
