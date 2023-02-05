library(ggplot2)
library(cowplot)
library(dplyr)


#entropy.df = read.csv('entropy.txt', sep=' ', header=F, stringsAsFactors = F)
entropy.df = read.csv('entropy-upstream-downstream.txt', sep=' ', header=F, stringsAsFactors = F)

colnames(entropy.df) = c("gene", "direction", "position", "entropy")

entropy.max = entropy.df %>% group_by(gene) %>%
  summarise(max=max(entropy))
new.order.genes = entropy.max$gene[order(entropy.max$max)]

entropy.df <- within(entropy.df, gene <- factor(gene, levels =new.order.genes))


p = ggplot(entropy.df, aes(position, entropy,group=direction, colour=direction))+
    geom_line()+
    theme_bw()+
    theme(panel.grid = element_blank())+
    ylab("Positional entropy (block diversity)")+
    xlab("Position (kb)\nrelative to central gene")+
    theme(panel.border = element_blank())+
    theme(axis.line = element_line())+
    facet_wrap(~gene, ncol=6)+
  scale_x_continuous(breaks=c(0,1000,2000,3000,4000,5000),
                     labels=c("0", "1", "2", "3", "4", "5"))+
  ylim(c(0,1))


colnames(entropy.df) = c("gene", "direction", "position", "entropy")
# Plot upstream-downstream at each point
library(dplyr)
entropy.df.diff = entropy.df %>% group_by(gene, position) %>%
  summarise(entropy.diff=entropy[direction=="upstream"]-entropy[direction=="downstream"])

order.genes.by.diff = entropy.df.diff %>% group_by(gene) %>%
  summarise(n=sum(entropy.diff))
new.order = order.genes.by.diff$gene[order(order.genes.by.diff$n)]
entropy.df.diff$gene = ordered(entropy.df.diff$gene,
                               levels=new.order)
entropy.df.diff <- within(entropy.df.diff, gene <- factor(gene, levels =new.order))

p.diff = ggplot(entropy.df.diff, aes(position, entropy.diff,group=gene))+
  geom_hline(yintercept = 0, colour='black')+
  geom_line()+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ylab("Positional entropy difference (upstream-downstream)")+
  xlab("Position (kb)\nrelative to central gene")+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line())+
  facet_wrap(~gene, ncol=6)+
  scale_x_continuous(breaks=c(0,1000,2000,3000,4000,5000),
                     labels=c("0", "1", "2", "3", "4", "5"))+
  theme(legend.position = "none")

ggsave(file='../output/plot-positional-entropy-diff.pdf', p.diff)

# and in close proximity?
p.diff.small = ggplot(entropy.df.diff, aes(position, entropy.diff,group=gene, colour=gene))+
  geom_hline(yintercept = 0, colour='black')+
  geom_line()+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ylab("Positional entropy difference (upstream-downstream)")+
  xlab("Position (kb)\nrelative to central gene")+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line())+
  facet_wrap(~gene, ncol=6)+
  xlim(c(0,1000))+
  theme(legend.position = "none")



p.short = ggplot(entropy.df, aes(V2, V3, group=V1, colour=V1))+
  geom_line()+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ylab("Block entropy")+
  xlab("Position (kb)\ncentral gene starts @ 5kb")+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line())+
  facet_wrap(~V1)+
  theme(legend.position = "none")+
  xlim(c(3000,7000))+
  scale_x_continuous(breaks=c(3000,4000,5000,6000,7000),
                     labels=c("3", "4", "5", "6", "7"),
                    limits = c(3000,7000))




cowplot::plot_grid(`p.CMY-2`, `p.VIM-1`)
get("p.CMY-2")
