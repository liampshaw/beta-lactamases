# Beta-lacatamase gene evolution

Analysing the evolution of beta-lactamases over time using pangraph.

Using CARD prevalence data, we obtain sequences containing beta-lactamases from twelve clinically important beta-lactamase families (`beta-lactamases.txt`) and analyse the surrounding regions of genes using [pangraph](https://github.com/neherlab/pangraph/). 

To create a conda environment with the necessary packages for other scripts excluding the pangraph analysis, run:

```
conda create -n betalactamases
conda activate betalactamases 
conda install scipy pandas snp-dists numpy seaborn biopython ncbi-acc-download ncbi-genome-download fasttree seqkit untangle entrez-direct
```
