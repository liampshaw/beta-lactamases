
library(ape)
library(ggplot2)
library(RColorBrewer)
suppressMessages(library(ggtree))
options(ignore.negative.edge=TRUE) # for ggtree
options(warn = -1)

args = commandArgs(trailingOnly=TRUE)

# Two arguments
# 1 - input alignment
# 2 - output pdf
# 3 width
# 4 height

treePlot <- function(variant_alignment_file){
  #gene = 'CMY-2'
  #deduplicated_prefix = paste0('../../data/2023-02-05-', gene, '-mmseqs2-polish-u5000-d5000/', gene, '-mmseqs2-polish.all_u5000_d5000_focal_gene.dedup')
  #dna <- read.dna(paste0(deduplicated_prefix,'.aln'), format='fasta')
  dna <- read.dna(variant_alignment_file, format='fasta')
  dist.obj <- ape::dist.dna(dna, model = "N") # distance in SNPs (diffs, actual number)
  tree <- ape::bionj(dist.obj)
  
  #metadata <- read.csv(paste0(deduplicated_prefix,'.txt'), 
  #                     sep='\t', 
  #                     header=F, 
  #                     col.names = c("n", "isolates"))
  #rownames(metadata) <- sapply(metadata$isolates, 
  #                             function(x) gsub(",.*", "", x))
  #metadata$isolates <- NULL
  #metadata$isolate <- rownames(metadata)
  
  #if (length(which(!tree$tip.label %in% metadata$isolate))>0){

  #metadata.single <- data.frame(isolate=tree$tip.label[which(!tree$tip.label %in% metadata$isolate)],
#                                n=1)
  
  #metadata <- rbind(metadata, metadata.single)
#}
  metadata = data.frame(isolate=gsub(" \\|.*", "", tree$tip.label),
                        variant=gsub(" .*", "", gsub(".*\\| ", "", tree$tip.label)),
                        n=gsub(".* ", "", tree$tip.label))
  #metadata <- metadata[,c("isolate", "n")]
  
  tree$tip.label = gsub(" \\|.*", "", tree$tip.label)
  rownames(metadata) <- tree$tip.label
  
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
  #metadata$isolate = as.character(rownames(metadata))
  #metadata = metadata[,c("isolate", "variant", "n")]
  p <- ggtree(tree, layout = 'rectangular')
  
  p <- p %<+% metadata+
    geom_tippoint(aes(size=n,   fill=categoricalSNPs), shape=21, colour="black")+
    geom_tiplab(aes(label=variant), size=2)+
    theme_tree2()+
    scale_size_continuous(range=c(2, 10), breaks=c(1, 10, 50, 100), name="No. sequences")+
    scale_fill_manual(values=snps.categorical.colour.palette)+
 theme(legend.title= element_text(size=8), 
          legend.text = element_text(size=6))+
    scale_size_continuous(range=c(2, 10), breaks=c(1, 10, 50, 100), name="No. sequences")+
    geom_treescale(width = 1, offset = 0.1, linesize = 1)+
    guides(fill="none")+
    ggtitle("NJ tree of central gene")
  return(p)
}


pdf(file=args[2], width=args[3], height=args[4])
p.tree = treePlot(args[1])
p.tree
dev.off()
