
require(ggplot2)
require(ggridges)
require(cowplot)

args = commandArgs(trailingOnly=TRUE)

d <- read.csv(args[1],
              header=T,
              stringsAsFactors = F)

d$snps.categorical <- sapply(d$snps, function(x)
  ifelse(x<8, x, ">7"))
d$snps.categorical <- ordered(d$snps.categorical,
                      levels=c(seq(0,7), ">7"))

dist.max <- max(c(d$dist.up, d$dist.down))

p.upstream <-  ggplot(d, aes( -dist.up, group=snps.categorical, colour=snps.categorical, y = 1 - ..y..))+
  stat_ecdf()+
  theme_bw()+
  theme(legend.position = "none")+
  labs(colour="SNPs")+
  ylab("cdf")+
  xlab("distance from gene")+
  theme(panel.grid = element_blank())+
  scale_color_brewer(palette="RdYlBu")+
  ggtitle("Upstream")+
  theme(axis.line.x = element_line(colour = "black"),
    axis.line.y=element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())+
  xlim(c(-dist.max, 0))+
  theme(plot.title=element_text(hjust=0.5))
p.downstream <- ggplot(d, aes( dist.down,group=snps.categorical, colour=snps.categorical))+
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
  ggtitle("Downstream")+
  theme(axis.line.x = element_line(colour = "black"),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())+
    xlim(c(0, dist.max))+
    theme(plot.title=element_text(hjust=0.5))



pdf(args[2], width=8, height=4)
p = cowplot::plot_grid(p.upstream, p.downstream, rel_widths = c(0.8, 1))
p
dev.off()
