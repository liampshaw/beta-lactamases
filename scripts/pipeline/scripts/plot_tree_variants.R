library(ape)
library(ggplot2)
library(ggtree)


args = commandArgs(trailingOnly=TRUE)


variant.names = read.csv(args[3], header=F, stringsAsFactors=F)
print(variant.names)

dedup.names = read.csv(args[4], header=F, stringsAsFactors=F, sep='\t')
print(dedup.names)
dedup.names$rep = gsub(",.*", "", dedup.names$V2)
rownames(dedup.names) = dedup.names$rep


variant_alignment_file = args[1]
tree_file = args[2]

root=NULL

variants = read.csv(args[3], header=F, stringsAsFactors = F, row.names = 1)

  if (!is.null(root)){
    tree = root(tree, outgroup=root)
  }
  
  #labels(dna)
   dna <- read.dna(args[], format='fasta')
  isolates = gsub(" .*", "",labels(dna))

  variant.names = variants[isolates, "V2"]
  n.values = dedup.names[isolates,"V1"]
  n.values[is.na(n.values)] = 1 # if they only have one rep, not in dedup.txt file
    metadata = data.frame(isolate=isolates, 
                        variant=variant.names,
                        n=n.values)

  #tree$tip.label = gsub(" \\|.*", "", tree$tip.label)
  rownames(metadata) <- metadata$isolate
  
  #metadata$focal <- as.factor(c(1, rep(0, nrow(metadata)-1)))# just store whether focal gene or not
  # Add SNPs compared to focal gene
  #dist.mat <- as.matrix(dist.obj)
  #rownames(dist.mat) = gsub(" \\|.*", "", rownames(dist.mat))
  #colnames(dist.mat) = rownames(dist.mat)
  #metadata$SNPs <- dist.mat[metadata$isolate[1],metadata$isolate]
  #metadata$SNPS.categorical <- sapply(metadata$SNPs, function(x)
  #  ifelse(x<8, x, ">7"))
  #metadata$categoricalSNPs <- ordered(metadata$SNPS.categorical,
  #                                    levels=c(seq(0,7), ">7"))
  #snps.categorical.colour.palette <- RColorBrewer::brewer.pal(name="RdYlBu", n=9)
  # 
  #metadata$n.label <- sapply(metadata$n, 
  #                           function(x) ifelse(x==1, "", x))
  metadata$variant[which(metadata$variant=="unnamed")] = ""
  metadata$variant[which(metadata$variant=="truncated")] = "*"
  
  metadata$isolate2 = metadata$isolate
  metadata$n = as.numeric(as.character(metadata$n))

  
  p <- ggtree(tree, layout = 'rectangular')
  
  p <- p %<+% metadata+
    geom_tippoint(aes(size=n), shape=21, colour="black", fill="white")+
    geom_tiplab(aes(label=variant), size=2, offset = 0.2)+
    scale_size_continuous(range=c(2, 10), breaks=c(1, 10, 50, 100), name="No. sequences")+
    scale_fill_manual(values=snps.categorical.colour.palette)+
    theme(legend.title= element_text(size=8), 
          legend.text = element_text(size=6))+
    scale_size_continuous(range=c(2, 10), breaks=c(1, 10, 50, 100), name="No. sequences")+
    geom_treescale(width = 1, offset = 0.1, linesize = 1)+
    guides(fill="none")




pdf(args[5], width=8, height=8)
p
dev.off()

