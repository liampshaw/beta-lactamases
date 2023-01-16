# Beta-lacatamase gene evolution

Analysing the evolution of beta-lactamases over time using pangraph.

Using CARD prevalence data, we obtain sequences containing beta-lactamases from twelve clinically important beta-lactamase families (`beta-lactamases.txt`) and analyse the surrounding regions of genes using [pangraph](https://github.com/neherlab/pangraph/). 

To create a conda environment with the necessary packages for other scripts excluding the pangraph analysis, run:

```
conda create -n betalactamases
conda activate betalactamases 
conda install scipy pandas snp-dists numpy seaborn biopython ncbi-acc-download ncbi-genome-download fasttree seqkit untangle entrez-direct
# Optional: conda install -c conda-forge -c bioconda mmseqs2 
```

To run scripts in `scripts` directory will also require downloading CARD data from [here](https://card.mcmaster.ca/download/6/prevalence-v3.1.0.tar.bz2) and extracting in the `data` directory. Only three files are needed: see README in `data` directory for more information. 
