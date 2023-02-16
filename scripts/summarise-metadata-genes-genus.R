# Summarise metadata with plots
library(ggplot2)
library(dplyr)
source('plot-theme.R')

# Read in metadata
metadata.df = read.csv('../data/metadata.csv',
                       header=T,
                       stringsAsFactors = F,
                       row.names = 1)
# >4 for no. of beta-lactamases
metadata.df$N.betalactamase.hits.simple = sapply(metadata.df$N.betalactamase.hits, 
                                                 function(x) ifelse(x>4, ">4", x))

# Genus table
genus.counts = sort(table(metadata.df$TaxGenus))
low.abundance.genera = names(genus.counts[which(genus.counts<10)])
# assign 'other' for these genera
metadata.df$TaxGenus.simple = sapply(metadata.df$TaxGenus,
                                     function(x) ifelse(x %in% low.abundance.genera, "Other", x))

# I want an overall overview of the dataset? 
# What things can I do?
# Summarise by genus and

# 
# Summarise number of genes by genus
gene.hits.by.genus = metadata.df %>% group_by(TaxGenus.simple, N.betalactamase.hits.simple) %>%
  summarise(n=length(TaxGenus.simple))

# Sort by total number in genus
gene.hits.by.genus$TaxGenus.simple = ordered(gene.hits.by.genus$TaxGenus.simple,
                                      levels=c("Other", names(genus.counts[which(genus.counts>10)])))

# Summarise number of genes by genus and chromosome/plasmid
gene.hits.by.genus.chrom = metadata.df %>% group_by(TaxGenus.simple, Contig, N.betalactamase.hits.simple) %>%
  summarise(n=length(TaxGenus.simple))
# Order by total amounts
gene.hits.by.genus.chrom$TaxGenus.simple = ordered(gene.hits.by.genus.chrom$TaxGenus.simple,
                                                                              levels=c("Other", names(genus.counts[which(genus.counts>10)])))
gene.hits.by.genus.chrom$N.betalactamase.hits.simple = ordered(gene.hits.by.genus.chrom$N.betalactamase.hits.simple,
                                                               levels=c("1", "2", "3", "4", ">4"))

gene.hits.by.genus.chrom$pseudo.n = gene.hits.by.genus.chrom$n+0.001
p.genes.genus = ggplot(gene.hits.by.genus.chrom, aes(TaxGenus.simple, pseudo.n, fill=N.betalactamase.hits.simple))+
  geom_bar(stat="identity", position="stack", colour="black", size=0.1)+
  theme_basic()+
  coord_flip()+
  scale_y_continuous(breaks=seq(0,2500, 500),
                     minor_breaks=seq(0,2700,100),
                     position = "right")+
  theme(panel.grid.major.x = element_line(colour="grey"),
        panel.grid.minor.x = element_line(colour="grey"))+
  ylab("")+
  xlab("")+
  labs(fill="No. of\nbeta-lactamases\ncarried")+
  theme(axis.text.y=element_text(face="italic"))+
  scale_fill_brewer(palette="Reds")+
  theme(panel.border  = element_blank(),
        axis.line=element_line(colour="black"),
        axis.text.y=element_text(size=8))+
  facet_wrap(~Contig)+
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size=10))

ggsave(file="../output/fig-metadata-genes-by-genus.pdf",
       p.genes.genus,
       width=8, height=4)

# Different types of enzyme family carried...
# But probably do this at the level of a subset





