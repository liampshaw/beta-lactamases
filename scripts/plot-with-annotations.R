require(ggplot2)
require(cowplot)
require(gggenes)
require(ggdendro)
require(reshape2)
require(vegan) 

# Example data
#args= c('../output/GES-24/plotting/GES.all_u5000_d5000_pangraph.json.blocks.csv',
#        'ESELPCQPHB',
#        '../output/GES-24/plotting/GES-test-abricate-hits.tsv',
#        '../output/GES-24/plotting/GES-test-annotation-output.pdf',
#         10,
#         6)

args = commandArgs(trailingOnly=TRUE)
# 1 - gfa csv
# 2 - name of block with gene in (for centering)
# 3 - name of file with annotations (from abricate)
# 4 - output pdf
# 5 - output pdf width (inches)
# 6 - output pdf height (inches)


genome.blocks = read.csv(args[1], header=T, stringsAsFactors = F)
genome.blocks$forward <- ifelse(genome.blocks$strand=="+", TRUE, FALSE)
block.counts <- table(genome.blocks$block)
blocks.which.need.colours <- names(block.counts)[which(block.counts>1)]

genome.blocks$block.coloured <- sapply(genome.blocks$block,
                                       function(x) ifelse(x %in% blocks.which.need.colours,
                                                          x,
                                                          "_other"))
block.colours <- unique(genome.blocks$colour)
names(block.colours) <- unique(genome.blocks$block.coloured)


gene.block.locations <- genome.blocks[which(genome.blocks$block==args[2]),c("start", "end")]
rownames(gene.block.locations) <- genome.blocks[which(genome.blocks$block==args[2]),"genome"]
transformAnnotationsBlocks <- function(genome, default_offset=9926){
  offset <- default_offset-gene.block.locations[genome, c("start")]
  return(offset)
}

genome.blocks$new.start <- apply(genome.blocks, MARGIN=1, function(x) as.numeric(x["start"])+transformAnnotationsBlocks(x["genome"]))-10001
genome.blocks$new.end <- apply(genome.blocks, MARGIN=1, function(x) as.numeric(x["end"])+transformAnnotationsBlocks(x["genome"]))-10001

# Get all the genome paths (in terms of blocks)
genome.paths <- sapply(unique(genome.blocks$genome), function(x) paste(genome.blocks[which(genome.blocks$genome==x), "block"], collapse=","))
genome.paths <- sort(genome.paths, decreasing = TRUE)

# Subset to unique paths
# Keep one representative genome for each
genome.path.reps <- names(genome.paths)[!duplicated(genome.paths)]
genome.blocks.unique <- genome.blocks[which(genome.blocks$genome %in% genome.path.reps),]
# And store the number of examples of each
genome.blocks.unique$genome.path <- genome.paths[genome.blocks.unique$genome]
genome.blocks.unique$n.reps <- sapply(genome.blocks.unique$genome.path,
                                      function(x) table(genome.paths)[x])
genome.blocks.unique$genome.path.name <- paste0("Type", as.numeric(as.factor(genome.blocks.unique$genome.path )))
genome.blocks.unique$genome.n <- sapply(genome.blocks.unique$n.reps, function(x) ifelse(x==1,
                                                                                        "", paste0("n=", x)))

# Now use a dendrogram to order the genomes
m <- acast(genome ~ block, data=genome.blocks.unique, fill=0, fun.aggregate=length)
m.dist <- vegdist(m, method="jaccard") # jaccard distances based on block presence/absence
dendro <- as.dendrogram(hclust(m.dist))
dendro_order <- order.dendrogram(dendro)
genome_labels <- dendro_data(dendro)$labels$label
genome.blocks.unique$genome.ordered <- ordered(genome.blocks.unique$genome,
                                               levels=genome_labels)

# Make the linear block plot
p.blocks <- ggplot(genome.blocks.unique, aes(xmin = new.start, xmax = new.end, forward = forward, y = genome.ordered, fill = block.coloured)) +
  geom_gene_arrow(arrow_body_height = unit(2, "mm"),
                  arrowhead_height = unit(2, "mm"),
                  arrowhead_width = unit(1, "mm"),
                  alpha=0.7) +
  theme_genes()+
  scale_fill_manual(values=block.colours)+
  ylab("")+
  theme(legend.position = "none")+
  scale_y_discrete(breaks=genome.blocks.unique$genome.ordered, labels=genome.blocks.unique$genome.n)+
  ggtitle("Linearized blocks")+
  theme(plot.title=element_text(hjust=0.5))+
  xlab("Position (bp)")

# Add IS/AMR hits
annotation.hits = read.csv(args[3], sep='\t', header=T, stringsAsFactors = F)
# Need to readjust the positions
annotation.hits$new.start <- apply(annotation.hits, MARGIN=1, function(x) as.numeric(x["START"])+transformAnnotationsBlocks(x["SEQUENCE"]))-10001
annotation.hits$new.end <- apply(annotation.hits, MARGIN=1, function(x) as.numeric(x["END"])+transformAnnotationsBlocks(x["SEQUENCE"]))-10001
annotation.hits$genome.ordered = ordered(annotation.hits$SEQUENCE,
                                 levels=genome_labels)
#annotation.hits.threshold = annotation.hits[which(annotation.hits$coverage>0.8),]
#annotation.hits.threshold = annotation.hits[which(annotation.hits$X.IDENTITY>90),]
annotation.hits.threshold = annotation.hits # Currently no thresholding, could using X.COVERAGE or X.IDENTITY (%)

# We have genome.blocks which has the start/end position of all blocks
# we have the start/end within each block of the IS hits
# Now for each genome we need to add in the IS hits
# Just pick one
annotation.hits.threshold.plot = annotation.hits.threshold[which(annotation.hits.threshold$genome %in% annotation.hits.threshold$genome.ordered),]

annotation.hits.threshold.plot$genome.ordered.numeric = as.numeric(annotation.hits.threshold.plot$genome.ordered)

# Use numeric genome
p.blocks.numeric = ggplot(genome.blocks.unique, aes(xmin = new.start, xmax = new.end, forward = forward, y = as.numeric(genome.ordered), fill = block.coloured)) +
  geom_gene_arrow(arrow_body_height = unit(2, "mm"),
                  arrowhead_height = unit(2, "mm"),
                  arrowhead_width = unit(1, "mm"),
                  alpha=0.7) +
  theme_genes()+
  scale_fill_manual(values=block.colours)+
  ylab("")+
  theme(legend.position = "none")+
  scale_y_discrete(breaks=as.numeric(genome.blocks.unique$genome.ordered), labels=genome.blocks.unique$genome.n)+
  ggtitle("Linearized blocks")+
  theme(plot.title=element_text(hjust=0.5))+
  xlab("Position (bp)")

# Make type into an offset so lines don't overplot
#annotation.hits.threshold.plot$type.offset = as.numeric(as.factor(annotation.hits.threshold.plot$type))-1
#annotation.hits.threshold.plot$type.offset[annotation.hits.threshold.plot$type.offset==0] = -1
#annotation.hits.threshold.plot$genome.ordered.numeric = as.numeric(annotation.hits.threshold.plot$genome.ordered)+
#  (annotation.hits.threshold.plot$type.offset)/8
pdf(args[4], height=as.numeric(args[5]), width=as.numeric(args[6]))
p.blocks.numeric+geom_segment(aes( x=new.start, xend=new.end, y=genome.ordered.numeric, yend=genome.ordered.numeric), 
                      inherit.aes = FALSE, 
                      size=0.5,
                      data=annotation.hits.threshold.plot)+
  geom_text(aes(label=GENE, x=(new.start+new.end)/2, y=genome.ordered.numeric), nudge_y=0.4, size=2,inherit.aes = FALSE, data=annotation.hits.threshold.plot)+
  scale_colour_manual(values=c("red", "black"))
dev.off()


