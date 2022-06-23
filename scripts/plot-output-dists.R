
require(ggplot2)
require(ggridges)
require(cowplot)

args = commandArgs(trailingOnly=TRUE)

d <- read.csv(args[1],
              header=T, 
              stringsAsFactors = F)

d$snps.categorical <- sapply(d$snps, function(x) 
  ifelse(x<8, x, "8 or more"))

dist.max <- max(c(d$dist.up, d$dist.down))

p.upstream <-  ggplot(d, aes( -dist.up, group=as.factor(snps.categorical), colour=as.factor(snps.categorical), y = 1 - ..y..))+
  stat_ecdf()+
  theme_bw()+
  theme(legend.position = "none")+
  labs(colour="SNPs")+
  ylab("cdf")+
  xlab("distance from gene")+
  theme(panel.grid = element_blank())+
  scale_color_brewer(palette="RdYlBu")+
  ggtitle("Upstream")
p.downstream <- ggplot(d, aes( dist.down,group=as.factor(snps.categorical), colour=as.factor(snps.categorical)))+
  stat_ecdf()+
  theme_bw()+
  labs(colour="SNPs")+
  ylab("cdf")+
  xlab("distance from gene")+
  theme(panel.grid = element_blank())+
  scale_color_brewer(palette="RdYlBu")+
  xlim(c(0,dist.max))+
  theme(axis.line.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.title.y=element_blank())+
  ggtitle("Downstream")

pdf(args[2], width=8, height=4)
p = cowplot::plot_grid(p.upstream, p.downstream, rel_widths = c(0.7, 1))
p
dev.off()
