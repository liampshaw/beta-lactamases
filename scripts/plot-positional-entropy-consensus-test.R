entropy.df.consensus = read.csv('test-consensus.csv', sep=',', header=F, stringsAsFactors = F)
colnames(entropy.df.consensus) = c("gene", "contig", "direction","N", "position", "entropy")

entropy.df = read.csv('test.csv', sep=',', header=F, stringsAsFactors = F)
colnames(entropy.df) = c("gene", "contig", "direction","N", "position", "entropy")


p1 = ggplot(entropy.df, aes(position, entropy,group=direction, colour=direction))+
  geom_line()+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ylab("Positional entropy (block diversity)")+
  xlab("Position (kb)\nrelative to central gene")+
  ggtitle("Actual lengths")+
  theme(legend.position = "none")

p2 =ggplot(entropy.df.consensus, aes(position, entropy,group=direction, colour=direction))+
  geom_line()+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ylab("Positional entropy (block diversity)")+
  xlab("Position (kb)\nrelative to central gene")+
  ggtitle("Consensus")+
  theme(legend.position = "none")


cowplot::plot_grid(p1, p2)
