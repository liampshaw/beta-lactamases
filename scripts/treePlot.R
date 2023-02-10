treePlot <- function(variant_alignment_file, tree_file="", root=""){
  dna <- read.dna(variant_alignment_file, format='fasta')
  dist.obj <- ape::dist.dna(dna, model = "N") # distance in SNPs (diffs, actual number)
  if (tree_file!=""){
    tree = read.tree(tree_file)
  }
  if (tree_file==""){
    tree <- ape::bionj(dist.obj)
    
  }
  
  if (root!=""){
    tree = root(tree, outgroup=root)
  }
  
  #labels(dna)

    metadata = data.frame(isolate=gsub(" \\|.*", "",labels(dna)),
                        variant=gsub(" .*", "", gsub(".*\\| ", "",labels(dna))),
                        n=gsub(".* ", "", labels(dna)))

  #tree$tip.label = gsub(" \\|.*", "", tree$tip.label)
  rownames(metadata) <- metadata$isolate
  
  metadata$focal <- as.factor(c(1, rep(0, nrow(metadata)-1)))# just store whether focal gene or not
  # Add SNPs compared to focal gene
  dist.mat <- as.matrix(dist.obj)
  rownames(dist.mat) = gsub(" \\|.*", "", rownames(dist.mat))
  colnames(dist.mat) = rownames(dist.mat)
  metadata$SNPs <- dist.mat[metadata$isolate[1],metadata$isolate]
  metadata$SNPS.categorical <- sapply(metadata$SNPs, function(x)
    ifelse(x<8, x, ">7"))
  metadata$categoricalSNPs <- ordered(metadata$SNPS.categorical,
                                      levels=c(seq(0,7), ">7"))
  snps.categorical.colour.palette <- RColorBrewer::brewer.pal(name="RdYlBu", n=9)
  
  #metadata$n.label <- sapply(metadata$n, 
  #                           function(x) ifelse(x==1, "", x))
  metadata$variant[which(metadata$variant=="unnamed")] = ""
  metadata$variant[which(metadata$variant=="truncated")] = "*"
  
  metadata$isolate2 = metadata$isolate
  metadata$n = as.numeric(as.character(metadata$n))

  
  p <- ggtree(tree, layout = 'rectangular')
  
  p <- p %<+% metadata+
    geom_tippoint(aes(size=n,   fill=categoricalSNPs), shape=21, colour="black")+
    geom_tiplab(aes(label=variant), size=2, offset = 0.2)+
    theme_tree2()+
    scale_size_continuous(range=c(2, 10), breaks=c(1, 10, 50, 100), name="No. sequences")+
    scale_fill_manual(values=snps.categorical.colour.palette)+
    theme(legend.title= element_text(size=8), 
          legend.text = element_text(size=6))+
    scale_size_continuous(range=c(2, 10), breaks=c(1, 10, 50, 100), name="No. sequences")+
    geom_treescale(width = 1, offset = 0.1, linesize = 1)+
    guides(fill="none")
  return(p)
}