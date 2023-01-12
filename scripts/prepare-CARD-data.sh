#!/bin/bash

Rscript 01_summarise-CARD-data.R 

python 02_make-matrix-card-data.py 

source 03a_align-families.sh 

python 03b_make-heatmap.py 
