# Split metadata
# Find all isolates with a gene, and split it 

args = commandArgs(trailingOnly=TRUE)

metadata = read.csv('../data/metadata.csv',
                    header=T,
                    stringsAsFactors = F,
                    row.names = 1)

gene = args[1]
genus = args[2]

chromosome = rownames(metadata)[grepl(gene, metadata$Gene.hits) & 
                     metadata$TaxGenus==genus &
                       metadata$Contig=="chromosome"]
plasmid = rownames(metadata)[grepl(gene, metadata$Gene.hits) & 
                               metadata$TaxGenus==genus &
                               metadata$Contig=="plasmid"]

cat(chromosome, sep="\n", file=paste0('tmp.', gene, '.chromosomes.txt'))
cat(plasmid, sep="\n", file=paste0('tmp.', gene, '.plasmids.txt'))