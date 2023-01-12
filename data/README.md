# Data sources

In the current analysis (Jan 2023) we are using the v3.1.0 release which can be downloaded from [here](https://card.mcmaster.ca/download/6/prevalence-v3.1.0.tar.bz2). See CARD's [download page](https://card.mcmaster.ca/download) for other releases.

* `accessions.txt` - accessions of all the NCBI chromosomes and plasmids (excluding WGS-only contigs) that contain a hit for one of the analysed beta-lactamase families (CARD prevalence data)
* `metadata.csv` - cleaned metadata for all the sequences listed in the above. See `metadata-processing` for details of how to do this

Sequence data can be downloaded from NCBI directly using the accessions. 

Further data from CARD is required. Specifically, after downloading the prevalence data, you will need to extract the archive in this folder and uncompress the following three files (the other files are not needed and can be deleted): 

* `card-genomes.txt`  
* `index-for-model-sequences.txt`
* `nucleotide_fasta_protein_homolog_model.fasta` 

CARD citation:

Alcock et al. "CARD 2020: antibiotic resistome surveillance with the Comprehensive Antibiotic Resistance Database" Nucleic Acids Research, 48: D517-D525. PMID: [31665441](https://www.ncbi.nlm.nih.gov/pubmed/31665441)

