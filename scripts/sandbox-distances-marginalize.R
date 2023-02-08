# 
library(ggplot2)
library(ggbeeswarm)
dist.df = read.csv('TEM-1-dists.csv', header=T)
ggplot(dist.df, aes(snps.gene, breakpoints))+
  geom_boxplot(aes(group=snps.gene))

ggplot(dist.df, aes(as.numeric(snps.gene), shared.seq))+
  ggbeeswarm::geom_quasirandom(aes(group=snps.gene), 
                               alpha=0.2)+
  stat_summary(fun=median, aes(group=snps.gene), colour="red")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ylim(c(0,10000))


# 
for (gene in c("TEM-1", "CTX-M-15", "CTX-M-65", "CMY-2", "OXA-1")){
  print(gene)
  dist.df = read.csv(paste0('salmonella-', gene, '-dists.csv'), header=T)
  print(c(median(dist.df$shared.seq[dist.df$snps.gene==0]), 
          median(dist.df$shared.seq[dist.df$snps.gene>0])))
}


