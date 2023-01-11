library(ggplot2)
library(dplyr)


df <- read.csv('SHV/pangraph_all_u5k_d5k.gfa.output_dists.csv', header=T)

metadata <- read.csv('SHV/all_u5k_d5k_focal_gene.dedup.txt', 
                     sep='\t', 
                     header=F, 
                     col.names = c("n", "isolates"))
focal.gene.isolates = strsplit(metadata$isolates[1], split=", ")[[1]]

d.subset = df[which(df$seq1 %in% focal.gene.isolates | df$seq2 %in% focal.gene.isolates),]

d.small.subset <- df[which(df$seq1 %in% focal.gene.isolates[1] | df$seq2 %in% focal.gene.isolates[1]),]

d.subset %>% group_by(snps) %>%
  summarise(median=median(dist.up), 
            upper=quantile(dist.down,p=0.75),
            lower=quantile(dist.down,p=0.25))

ggplot(d.small.subset, aes(snps, dist.up))+
  geom_jitter(aes(group=snps), width=0.25, shape=21, colour='black', fill='white')+
  theme_bw()+
  theme(panel.grid = element_blank())


#  geom_smooth(method="lm", formula= (y ~ exp(-x)), linetype = 1) 