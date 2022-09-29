
require(ggplot2)
require(ggridges)
require(cowplot)

args = commandArgs(trailingOnly=TRUE)
# Args:
# [1] pangraph_all_u5k_d5k.gfa.output_dists.csv 
# [2] representative comparison sequence: pangraph_all_u5k_d5k.gfa.most_frequent_path_representative.txt 
# [3] deduplicated gene sequences, list: all_u5k_d5k_focal_gene.dedup.txt
# pangraph_all_u5k_d5k.gfa.unique_paths.txt

d <- read.csv(args[1],
              header=T,
              stringsAsFactors = F)

d$snps.categorical <- sapply(d$snps, function(x)
  ifelse(x<8, x, ">7"))
d$snps.categorical <- ordered(d$snps.categorical,
                      levels=c(seq(0,7), ">7"))

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
metadata <- read.csv(args[3], 
                     sep='\t', 
                     header=F, 
                     col.names = c("n", "isolates"))
focal.gene.isolates = strsplit(metadata$isolates[1], split=", ")[[1]]

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
      scale_color_brewer(palette="RdYlBu")+
      ggtitle("Upstream")+
      theme(axis.line.x = element_line(colour = "black"),
        axis.line.y=element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
      xlim(c(-dist.max, 0))+
      theme(plot.title=element_text(hjust=0.5))
    p.downstream <- ggplot(df, aes( dist.down,group=snps.categorical, colour=snps.categorical))+
      stat_ecdf()+
      theme_bw()+
      labs(colour="SNPs")+
      ylab("cdf")+
      xlab("distance from gene (bp)")+
      theme(panel.grid = element_blank())+
      scale_color_brewer(palette="RdYlBu")+
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
        theme(plot.title=element_text(hjust=0.5))
    p = cowplot::plot_grid(p.upstream, p.downstream, rel_widths = c(0.8, 1))
    return(p)
}




pdf(paste0(args[1], '.flanking-plot-output-', 'all.pdf'), width=8, height=4)
makePlots(d)
dev.off()

pdf(paste0(args[1],'.flanking-plot-output-', representative, '.pdf'), width=8, height=4)
makePlots(d.subset)
dev.off()

pdf(paste0(args[1],'.flanking-plot-output-focal-gene-seq.pdf'), width=8, height=4)
makePlots(d.subset.2)
dev.off()
