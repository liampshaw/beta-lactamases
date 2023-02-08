args = c("../../data/2023-02-05-CMY-2-mmseqs2-polish-u5000-d5000/CMY-2-mmseqs2-polish.all_u5000_d5000_pangraph.json.blocks.csv","SWUPAMVWBM", "none" ,"test_pangraph_blocks_plot.pdf", "20", "20")

# Genome blocks approach
genome.blocks <- read.csv(args[1], header=T, stringsAsFactors = F)
#colnames(genome.blocks) <- c("genome", "block", "strand", "start", "end", "colour")
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

# Make
p = ggplot(genome.blocks.unique, aes(xmin = new.start, xmax = new.end, forward = forward, ymin = as.numeric(genome.ordered), ymax=as.numeric(genome.ordered)+0.5, fill = block.coloured)) +
  geom_rect() +scale_fill_manual(values=block.colours)+theme_bw()+
  theme(legend.position = "none")  +
  theme(panel.border = element_blank())+
  theme(axis.line.y=element_blank())
plotly(p)
