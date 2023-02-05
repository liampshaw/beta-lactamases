# Metadata summary for a gene
library(ape)

metadata = read.csv('../data/metadata.csv',
                    header=T,
                    stringsAsFactors = F,
                    row.names = 1)

gene = "CMY-2"

fasta = read.dna(paste0('../../data/2023-01-18-', gene, '-mmseqs2-polish-u5000-d5000/', gene, '-mmseqs2-polish.all_u5000_d5000.fa'),
                 format = 'fasta')
gene.metadata = metadata[gsub(" .*", "", names(fasta)),]
