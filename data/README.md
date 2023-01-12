# Data sources

Description of data files:

* `accessions.txt` - accessions of all the NCBI chromosomes and plasmids (excluding WGS-only contigs) that contain a hit for one of the analysed beta-lactamase families (CARD prevalence data)
* `metadata.csv` - cleaned metadata for all the sequences listed in the above. See `metadata-processing` for details of how to do this
* `nucleotide_fasta_protein_homolog_model.fasta` - CARD v3.2.2 nucleotide fasta of AMR genes which confer resistance through their presence (includes all beta-lactamases) 

Sequence data can be downloaded from NCBI directly using the accessions. 

Two large data files from CARD Prevalence Data are required. Specifically, after downloading the prevalence data, you will need to extract the archive in this folder and uncompress the following two files (the other files are not needed and can be deleted): 

* `card-genomes.txt`  
* `index-for-model-sequences.txt`

In the current analysis (Jan 2023) we are using the v3.1.0 release which can be downloaded from [here](https://card.mcmaster.ca/download/6/prevalence-v3.1.0.tar.bz2). See CARD's [download page](https://card.mcmaster.ca/download) for other releases.

The CARD nucleotide fasta for gene sequences is from CARD Data v3.2.2, which can be downloaded [here](https://card.mcmaster.ca/download/0/broadstreet-v3.2.2.tar.bz2). 

CARD citation:

Alcock et al. "CARD 2020: antibiotic resistome surveillance with the Comprehensive Antibiotic Resistance Database" Nucleic Acids Research, 48: D517-D525. PMID: [31665441](https://www.ncbi.nlm.nih.gov/pubmed/31665441)

