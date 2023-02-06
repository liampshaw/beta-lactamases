library(ggplot2)
library(cowplot)
library(dplyr)


plotForGenus <- function(GENUS){
  #entropy.df = read.csv('entropy.txt', sep=' ', header=F, stringsAsFactors = F)
  entropy.df = read.csv(paste0('../output/positional-entropies-', GENUS, '.csv'), sep=',', header=F, stringsAsFactors = F)
  
  colnames(entropy.df) = c("gene", "contig", "direction","N", "position", "entropy")
  
  # entropy.max = entropy.df %>% group_by(gene) %>%
  #   summarise(max=max(entropy))
  # new.order.genes = entropy.max$gene[order(entropy.max$max)]
  # 
  # entropy.df <- within(entropy.df, gene <- factor(gene, levels =new.order.genes))
  
  entropy.df$gene.N = paste(entropy.df$gene, entropy.df$N, sep="\n")
  
  N.df = entropy.df %>% group_by(gene, contig) %>%
    summarise(n=paste0("n=",unique(N)))
  
  
  p = ggplot(entropy.df, aes(position, entropy,group=direction, colour=direction))+
    geom_line()+
    theme_bw()+
    theme(panel.grid = element_blank())+
    ylab("Positional entropy (block diversity)")+
    xlab("Position (kb)\nrelative to central gene")+
    theme(panel.border = element_blank())+
    theme(axis.line = element_line())+
    facet_grid(contig~gene)+
    scale_x_continuous(breaks=c(0,1000,2000,3000,4000,5000),
                       labels=c("0", "1", "2", "3", "4", "5"))+
    ylim(c(0,1))+
    geom_text(x=1000, y=0.8, inherit.aes = FALSE,data=N.df, aes(label=n),
              size=2)
  
  ggsave(p, file=paste0('../output/positional-entropies-', GENUS, '.pdf'),
      width=15, height=4)
  #p
  #dev.off()

}

for (genus in c("escherichia", "klebsiella", "acinetobacter", 
                "pseudomonas", "salmonella", "enterobacter",
                "campylobacter", "citrobacter")){
  plotForGenus(genus)
}

# For klebsiella
entropy.df = read.csv('../output/positional-entropies-klebsiella.csv', sep=',', header=F, stringsAsFactors = F)

colnames(entropy.df) = c("gene", "contig", "direction","N", "position", "entropy")

# entropy.max = entropy.df %>% group_by(gene) %>%
#   summarise(max=max(entropy))
# new.order.genes = entropy.max$gene[order(entropy.max$max)]
# 
# entropy.df <- within(entropy.df, gene <- factor(gene, levels =new.order.genes))

entropy.df$gene.N = paste(entropy.df$gene, entropy.df$N, sep="\n")

N.df = entropy.df %>% group_by(gene, contig) %>%
  summarise(n=paste0("n=",unique(N)))


p = ggplot(entropy.df, aes(position, entropy,group=direction, colour=direction))+
  geom_line()+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ylab("Positional entropy (block diversity)")+
  xlab("Position (kb)\nrelative to central gene")+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line())+
  facet_grid(contig~gene)+
  scale_x_continuous(breaks=c(0,1000,2000,3000,4000,5000),
                     labels=c("0", "1", "2", "3", "4", "5"))+
  ylim(c(0,1))+
  geom_text(x=1000, y=0.8, inherit.aes = FALSE,data=N.df, aes(label=n),
            size=2)

ggsave(p, file='../output/positional-entropies-klebsiella.pdf',
       width=15, height=4)
