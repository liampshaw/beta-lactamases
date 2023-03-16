# Summarise metadata for each gene variant

library(ggplot2)
library(dplyr)
library(RColorBrewer)
source('plot-theme.R')

# Read in metadata
metadata.df = read.csv('../data/metadata.csv',
                       header=T,
                       stringsAsFactors = F,
                       row.names = 1)

genus.counts = sort(table(metadata.df$TaxGenus))
low.abundance.genera = names(genus.counts[which(genus.counts<10)])
# assign 'other' for these genera
metadata.df$TaxGenus.simple = sapply(metadata.df$TaxGenus,
                                     function(x) ifelse(x %in% low.abundance.genera, "Other", x))

plotSubset <- function(isolates){
  # subset to isolates requested
  metadata.subset = metadata.df[isolates,]
  # By year (ignore missing data)
  #metadata.subset$missing.year = ifelse(metadata.subset$Year=="missing", "missing", "present")
  year.counts = metadata.subset %>% filter(Year!="missing") %>% group_by(Year) %>%
    summarise(n=length(Year))
  year.counts$Year = as.numeric(as.character(year.counts$Year))
  p.year = ggplot(year.counts, aes(Year, n))+
    geom_bar(stat="identity", fill="black")+
    theme_basic()+
    xlab("")+
    theme(legend.position = "none")
  # By genus and chrom/plasmid
  genus.df = metadata.subset %>% group_by(TaxGenus.simple, Contig) %>%
    summarise(n=length(TaxGenus.simple))
  totals = genus.df %>% group_by(TaxGenus.simple) %>% summarise(n=sum(n))
  genus.df$TaxGenus.simple = ordered(genus.df$TaxGenus.simple,
                                         levels=totals$TaxGenus.simple[order(totals$n)])
  p.genus = ggplot(genus.df, aes(TaxGenus.simple, n, fill=Contig))+
    scale_x_discrete()+
    geom_bar(stat="identity", position="stack")+
    theme_basic()+
    xlab("")+
    scale_fill_manual(values=c("black", "red"),
                        breaks=c("chromosome", "plasmid"),
                      drop=FALSE)+
    coord_flip()+
    theme(axis.text.y=element_text(face="italic"))+
    theme(legend.position = c(0.6, 0.2),
          legend.title = element_blank(),
          legend.text = element_text(size=6))
  # Regions
  region.df = metadata.subset %>% group_by(Region) %>%
    summarise(n=length(Region))
  region.df$Region = ordered(region.df$Region,
                             levels=REGIONS)
  p.region = ggplot(region.df, aes(Region, n, fill=Region))+
    scale_x_discrete(drop=FALSE)+
    geom_bar(stat="identity", position="stack")+
    theme_basic()+
    xlab("")+
    scale_fill_manual(values=c("grey", RColorBrewer::brewer.pal(7, "Set3")), drop=FALSE)+
    coord_flip()+
    theme(legend.position = "none")
  
  
  p.top = cowplot::plot_grid(p.genus+ggtitle("(a)"), p.region+ggtitle("(b)"),
                         nrow=1,
                         align='h')
  p.bottom = p.year+ggtitle("(c)")
  p = cowplot::plot_grid(p.top, p.bottom, nrow=2, rel_widths = c(1, 0.8))
  p = cowplot::plot_grid(p.genus+ggtitle("(a) Genus"), p.region+ggtitle("(b) Region"),
                         p.year+ggtitle("(c) Year"), nrow=1, rel_widths = c(1, 1.2, 0.7))
  
  return(p)
  
}

genes = read.csv('../data/genes.txt', header=F)$V1
for (gene in genes){
  p = plotSubset(read.csv(paste0('../data-processing/analysis-clusters/', gene, '.txt'), header=F)$V1)
  ggsave(p, file=paste0('../output/fig-metadata-', gene, '.pdf'), width=8, height=3)
  print(gene)
  
}
