# 
library(ggplot2)
library(ggbeeswarm)
dist.df = read.csv('NDM-1-dists.csv', header=T)
ggplot(dist.df, aes(snps.gene, breakpoints))+
  geom_boxplot(aes(group=snps.gene))

ggplot(dist.df, aes(as.numeric(snps.gene), shared.seq))+
  ggbeeswarm::geom_quasirandom(aes(group=snps.gene), 
                               alpha=0.2)+
  stat_summary(fun=median, aes(group=snps.gene), colour="red")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ylim(c(0,10000))


