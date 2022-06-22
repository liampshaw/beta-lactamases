d <- read.csv('output_dists.csv', 
              header=F, 
              stringsAsFactors = F)

library(ggplot2)
library(ggridges)
library(cowplot)

p.upstream <- ggplot(d, aes(x=as.factor(V5), y=V3))+
  geom_violin()+
  xlab("SNPs in CTX-M")+
  ylab("Distance to breakpoint")+
  ggtitle("Upstream")+
  theme_bw()+
  theme(panel.grid = element_blank())



p.downstream <- ggplot(d, aes(x=as.factor(V5), y=V4))+
  geom_violin()+
  xlab("SNPs in CTX-M")+
  ylab("Distance to breakpoint")+
  ggtitle("Downstream")+
  theme_bw()+
  theme(panel.grid = element_blank())



png('CTX-M-snps.png', res=300, width=6*300, height=4*300, units="cm")
cowplot::plot_grid(p.upstream, p.downstream)
dev.off()

ggplot(d, aes( V3, colour=as.factor(V5)))+
  stat_ecdf()
p.upstream <-  ggplot(d, aes( -V3, colour=as.factor(V5), y = 1 - ..y..))+
  stat_ecdf()+
  theme_bw()+
  theme(legend.position = "none")+
  labs(colour="SNPs")+
  ylab("cdf")+
  xlab("distance from CTX-M")+
  theme(panel.grid = element_blank())

p.downstream <- ggplot(d, aes( V4, colour=as.factor(V5)))+
  stat_ecdf()+
  theme_bw()+
  theme(legend.position = c(0.8,0.2))+
  labs(colour="SNPs")+
  ylab("cdf")+
  xlab("distance from CTX-M")+
  theme(panel.grid = element_blank())
cowplot::plot_grid(p.upstream, p.downstream)
