
require(ggplot2)
require(ggridges)
require(cowplot)
require(RColorBrewer)
require(ape)
require(ggtree)
#require(magick)

args = commandArgs(trailingOnly=TRUE)
# Args:
# [1] pangraph_all_u5k_d5k.gfa.output_dists.csv 
# [2] representative comparison sequence: pangraph_all_u5k_d5k.gfa.most_frequent_path_representative.txt 
# [3] deduplicated gene sequences, list: all_u5k_d5k_focal_gene.dedup.txt
# pangraph_all_u5k_d5k.gfa.unique_paths.txt
# This version of script uses metadata to include only 'different' comparisons
# i.e. exclude same country & same year comparisons, which may be result of outbreaks and repeated sequencing
metadata <- read.csv('../data/metadata.csv', header=T, stringsAsFactors = F, row.names = 1)

d <- read.csv(args[1],
              header=T,
              stringsAsFactors = F)

d$snps.categorical <- sapply(d$snps, function(x)
  ifelse(x<8, x, ">7"))
d$snps.categorical <- ordered(d$snps.categorical,
                              levels=c(seq(0,7), ">7"))
snps.categorical.colour.palette <- RColorBrewer::brewer.pal(name="RdYlBu", n=9)

# Add metadata information
d$country1 = metadata[d$seq1, "Country"]
d$country2 = metadata[d$seq2, "Country"]
d$year1 = metadata[d$seq1, "Year"]
d$year2 = metadata[d$seq2, "Year"]
d$genus1 = metadata[d$seq1, "TaxGenus"]
d$genus2 = metadata[d$seq2, "TaxGenus"]




# get representative of the most common path
representative = read.csv(args[2], header=F, stringsAsFactors=F)$V1

print("before subsetting")
print(nrow(d))
print("rows")

d.subset = d[which(d$seq1==representative | d$seq2==representative),]
print("After subsetting")
print(nrow(d.subset))
print("rows")


# Also get most common gene sequence - assuming it exists
seq.metadata <- read.csv(args[3], 
                     sep='\t', 
                     header=F, 
                     col.names = c("n", "isolates"))
focal.gene.isolates = strsplit(seq.metadata$isolates[1], split=", ")[[1]]

d.subset.2 = d[which(d$seq1 %in% focal.gene.isolates | d$seq2 %in% focal.gene.isolates),]


dist.max <- max(c(d$dist.up, d$dist.down))


makePlots <- function(df){
  p.upstream <-  ggplot(df, aes( -dist.up, group=snps.categorical, colour=snps.categorical, y = 1 - ..y..))+
    stat_ecdf()+
    theme_bw()+
    theme(legend.position = "none")+
    labs(colour="SNPs")+
    ylab("cdf")+
    xlab("distance from gene (bp)")+
    theme(panel.grid = element_blank())+
    scale_color_manual(values=snps.categorical.colour.palette)+
    ggtitle("Upstream")+
    theme(axis.line.x = element_line(colour = "black"),
          axis.line.y=element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    xlim(c(-dist.max, 0))+
    theme(plot.title=element_text(hjust=0.5))+
    scale_y_reverse()
  p.downstream <- ggplot(df, aes( dist.down,group=snps.categorical, colour=snps.categorical))+
    stat_ecdf()+
    theme_bw()+
    labs(colour="SNPs")+
    ylab("cdf")+
    xlab("distance from gene (bp)")+
    theme(panel.grid = element_blank())+
    scale_color_manual(values=snps.categorical.colour.palette)+
    xlim(c(0,dist.max))+
    theme(axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank())+
    ggtitle("Downstream")+
    theme(axis.line.x = element_line(colour = "black"),
          axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    xlim(c(0, dist.max))+
    theme(plot.title=element_text(hjust=0.5))+
    scale_y_reverse()
  p = cowplot::plot_grid(p.upstream, p.downstream, rel_widths = c(0.8, 1))
  return(p)
}

# For tree of SNPs
treePlot <- function(deduplicated_prefix){
  dna <- read.dna(paste0(deduplicated_prefix,'.aln'), format='fasta')
  dist.obj <- ape::dist.dna(dna, model = "N") # distance in SNPs (diffs, actual number)
  tree <- ape::bionj(dist.obj)
  
  seq.metadata <- read.csv(paste0(deduplicated_prefix,'.txt'), 
                       sep='\t', 
                       header=F, 
                       col.names = c("n", "isolates"))
  rownames(seq.metadata) <- sapply(seq.metadata$isolates, 
                               function(x) gsub(",.*", "", x))
  seq.metadata$isolates <- NULL
  seq.metadata$isolate <- rownames(seq.metadata)
  
  if (length(which(!tree$tip.label %in% seq.metadata$isolate))>0){
    
    seq.metadata.single <- data.frame(isolate=tree$tip.label[which(!tree$tip.label %in% seq.metadata$isolate)],
                                  n=1)
    
    seq.metadata <- rbind(seq.metadata, seq.metadata.single)
  }
  seq.metadata <- seq.metadata[,c("isolate", "n")]
  rownames(seq.metadata) <- seq.metadata$isolate
  
  seq.metadata$focal <- as.factor(c(1, rep(0, nrow(seq.metadata)-1)))# just store whether focal gene or not
  # Add SNPs compared to focal gene
  dist.mat <- as.matrix(dist.obj)
  seq.metadata$SNPs <- dist.mat[seq.metadata$isolate[1],seq.metadata$isolate]
  seq.metadata$SNPS.categorical <- sapply(seq.metadata$SNPs, function(x)
    ifelse(x<8, x, ">7"))
  seq.metadata$categoricalSNPs <- ordered(seq.metadata$SNPS.categorical,
                                      levels=c(seq(0,7), ">7"))
  snps.categorical.colour.palette <- RColorBrewer::brewer.pal(name="RdYlBu", n=9)
  
  seq.metadata$n.label <- sapply(seq.metadata$n, 
                             function(x) ifelse(x==1, "", x))
  
  p <- ggtree(tree, layout = 'rectangular')
  p <- p %<+% seq.metadata+
    geom_tippoint(aes(size=n, fill=categoricalSNPs), shape=21, colour="black")+
    #geom_tiplab2(aes(label=n.label), geom="label",size=4, hjust = 0)+
    scale_fill_manual(values=snps.categorical.colour.palette)+
    theme(legend.title= element_text(size=8), 
          legend.text = element_text(size=6))+
    scale_size_continuous(range=c(2, 10), breaks=c(1, 10, 50, 100), name="No. sequences")+
    geom_treescale(width = 1, offset = 0.1, linesize = 1)+
    guides(fill="none")+
    ggtitle("NJ tree of central gene")
  return(p)
}



pdf(paste0(args[1], '.flanking-plot-output-', 'all-not-same-country-year.pdf'), width=8, height=4)
d.omit.same.country.year = d[which(d$country1!=d$country2 & d$year1!=d$year2),]
makePlots(d.omit.same.country.year)
dev.off()
