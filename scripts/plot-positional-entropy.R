library(ggplot2)
library(cowplot)


entropy.df = read.csv('entropy.txt', sep=' ', header=F, stringsAsFactors = F)

p = ggplot(entropy.df, aes(V2, V3, group=V1, colour=V1))+
    geom_line()+
    theme_bw()+
    theme(panel.grid = element_blank())+
    ylab("Block entropy")+
    xlab("Position (kb)\ncentral gene starts @ 5kb")+
    theme(panel.border = element_blank())+
    theme(axis.line = element_line())+
    facet_wrap(~V1)+
  theme(legend.position = "none")+
  scale_x_continuous(breaks=c(0,2500,5000,7500, 10000),
                     labels=c("0", "2.5", "5", "7.5", "10"))


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
