dna <- read.dna('SHV/all_u5k_d5k_focal_gene.dedup.aln', format='fasta')
tree <- nj(ape::dist.dna(dna, model = "N"))

metadata <- read.csv('SHV/all_u5k_d5k_focal_gene.dedup.txt', 
                     sep='\t', 
                     header=F, 
                     col.names = c("n", "isolates"))
rownames(metadata) <- sapply(metadata$isolates, 
                             function(x) gsub(",.*", "", x))
metadata$isolates <- NULL
metadata$isolate <- rownames(metadata)

metadata.single <- data.frame(isolate=tree$tip.label[which(!tree$tip.label %in% metadata$isolate)],
                              n=1)

metadata <- rbind(metadata, metadata.single)
metadata <- metadata[,c("isolate", "n")]
p <- ggtree(tree, layout = 'fan')
p %<+% metadata+
  geom_tippoint(aes(size=n))+
  geom_tiplab(aes(label=n), offset=0.3)
