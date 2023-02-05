# Genome blocks and year
library(dplyr)
library(ggplot2)
GENE = "NDM-1"

genome.blocks.file = paste0("../../data/2023-01-18-", GENE, "-mmseqs2-polish-u5000-d5000/", GENE, "-mmseqs2-polish.all_u5000_d5000_pangraph.json.blocks.csv")
genome.blocks =  read.csv(genome.blocks.file, header=T, stringsAsFactors = F)
genome.blocks$length = abs(genome.blocks$end-genome.blocks$start)

metadata = read.csv('../data/metadata.csv', 
                    stringsAsFactors = F,
                    header = T,
                    row.names = 1)
genome.blocks$year = metadata[genome.blocks$genome, "Year"]
plot.df = genome.blocks %>% filter(year!="missing") %>% group_by(genome) %>% 
  summarise(n=length(genome), 
            year=as.numeric(unique(year)))

ggplot(plot.df, aes(n, year))+
  geom_point()+
  stat_smooth(method="lm")
))



ggplot(genome.blocks, aes(length, group=genome))+
  geom_density()+
  facet_wrap(~year)
