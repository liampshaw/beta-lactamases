library(ape)
library(ggtree)
library(ggplot2)

dna <- read.dna('all_u5k_d5k_focal_gene.dedup.aln', format='fasta')
dist.obj <- ape::dist.dna(dna, model = "N") # distance in SNPs (diffs, actual number)
tree <- ape::bionj(dist.obj)

metadata <- read.csv('all_u5k_d5k_focal_gene.dedup.txt', 
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
rownames(metadata) <- metadata$isolate

metadata$focal <- as.factor(c(1, rep(0, nrow(metadata)-1)))# just store whether focal gene or not
# Add SNPs compared to focal gene
dist.mat <- as.matrix(dist.obj)
metadata$SNPs <- dist.mat[metadata$isolate[1],metadata$isolate]
metadata$SNPS.categorical <- sapply(metadata$SNPs, function(x)
  ifelse(x<8, x, ">7"))
metadata$categoricalSNPs <- ordered(metadata$SNPS.categorical,
                              levels=c(seq(0,7), ">7"))
snps.categorical.colour.palette <- RColorBrewer::brewer.pal(name="RdYlBu", n=8)

metadata$n.label <- sapply(metadata$n, 
                           function(x) ifelse(x==1, "", x))

p <- ggtree(tree, layout = 'equal_angle')
pdf('NJ-tree-gene.pdf')
p %<+% metadata+
  geom_tippoint(aes(size=n, fill=categoricalSNPs), shape=21, colour="black")+
  #geom_tiplab2(aes(label=n.label), geom="label",size=4, hjust = 0)+
  scale_fill_manual(values=snps.categorical.colour.palette)+
  scale_size_continuous(range=c(2, 10), breaks=c(1, 10, 50, 100), name="No. sequences")+
  geom_treescale(width = 1, offset = 0.1, linesize = 1)+
  labs(fill="SNPs")
dev.off()


treePlot <- function(deduplicated_prefix){
  dna <- read.dna(paste0(deduplicated_prefix,'.aln'), format='fasta')
  dist.obj <- ape::dist.dna(dna, model = "N") # distance in SNPs (diffs, actual number)
  tree <- ape::bionj(dist.obj)
  
  metadata <- read.csv(paste0(deduplicated_prefix,'.txt'), 
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
  rownames(metadata) <- metadata$isolate
  
  metadata$focal <- as.factor(c(1, rep(0, nrow(metadata)-1)))# just store whether focal gene or not
  # Add SNPs compared to focal gene
  dist.mat <- as.matrix(dist.obj)
  metadata$SNPs <- dist.mat[metadata$isolate[1],metadata$isolate]
  metadata$SNPS.categorical <- sapply(metadata$SNPs, function(x)
    ifelse(x<8, x, ">7"))
  metadata$categoricalSNPs <- ordered(metadata$SNPS.categorical,
                                      levels=c(seq(0,7), ">7"))
  snps.categorical.colour.palette <- RColorBrewer::brewer.pal(name="RdYlBu", n=8)
  
  metadata$n.label <- sapply(metadata$n, 
                             function(x) ifelse(x==1, "", x))
  
  p <- ggtree(tree, layout = 'equal_angle')
  p <- p %<+% metadata+
    geom_tippoint(aes(size=n, fill=categoricalSNPs), shape=21, colour="black")+
    #geom_tiplab2(aes(label=n.label), geom="label",size=4, hjust = 0)+
    scale_fill_manual(values=snps.categorical.colour.palette)+
    scale_size_continuous(range=c(2, 10), breaks=c(1, 10, 50, 100), name="No. sequences")+
    geom_treescale(width = 1, offset = 0.1, linesize = 1)+
    labs(fill="SNPs")
  return(p)
}
