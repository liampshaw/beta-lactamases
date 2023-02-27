# Beta-lactamase gene evolution

Repository of scripts and data for analysis of the genomic contexts of clinically important beta-lactamases genes from twelve major enzyme families. 

## Preparation

*Warning.* This repo is a work-in-progress and not currently set-up with wider use in mind.   

Using CARD prevalence data, we obtain sequences containing beta-lactamases from twelve clinically important beta-lactamase families (`beta-lactamases.txt`) and analyse the surrounding regions of genes using [pangraph](https://github.com/neherlab/pangraph/).  

This means installing pangraph following the instructions on the `pangraph` wiki. 

To create a conda environment with the necessary packages for other scripts excluding the pangraph analysis, run:

```
conda env create -f betalactamases.yml 
#conda create -n betalactamases
#conda activate betalactamases 
#conda install scipy pandas snp-dists numpy seaborn biopython ncbi-acc-download ncbi-genome-download fasttree seqkit untangle entrez-direct
# Optional to use mmseqs version of pipeline:
#conda install -c conda-forge -c bioconda mmseqs2 
```

To run scripts in `scripts` directory will also require downloading prevalence data files from the helpful Comprehensive Antimicrobial Resistance Database [CARD](https://card.mcmaster.ca) and extracting in the `data` directory. Only three CARD files are needed: see `README.md` in `data` directory for more information.  

