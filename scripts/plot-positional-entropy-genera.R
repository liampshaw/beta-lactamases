library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(ggrepel)


plotForGenus <- function(GENUS, output_pdf, GENES="all", CONTIGS="all"){
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
  
  if (genes=="ALL"){
    entropy.df.plot = entropy.df
    N.df.plot = N.df
  }
  else{
    entropy.df.plot = entropy.df[which(entropy.df$gene %in% GENES),]
    N.df.plot = N.df[which(N.df$gene %in% GENES),]
  }
  if (CONTIGS!="ALL"){
    entropy.df.plot = entropy.df.plot[which(entropy.df.plot$contig %in% CONTIGS),]
    N.df.plot = N.df.plot[which(N.df.plot$contig %in% CONTIGS),]
    
  }
  p = ggplot(entropy.df.plot, aes(position, entropy,group=direction, colour=direction))+
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
    geom_text(x=1000, y=0.8, inherit.aes = FALSE,data=N.df.plot,
              aes(label=n),
              size=2)
  
  ggsave(p, file=paste0('../output/', output_pdf),
      width=15, height=4)
  #p
  #dev.off()
  return(p)

}

for (genus in c("escherichia", "klebsiella", "acinetobacter", 
                "pseudomonas", "salmonella", "enterobacter",
                "campylobacter", "citrobacter")){
  plotForGenus(genus, paste0('plot-entropies-', genus, '.pdf'))
}


p = plotForGenus("escherichia", 
             output_pdf='fig-entropies-escherichia.pdf',
             GENES=c("CMY-2", "CTX-M-15", 
                     "CTX-M-65", "KPC-2", "NDM-1",
                     "OXA-1",  "TEM-1"),
             CONTIGS="plasmid")
pdf('../manuscript/figs/fig-entropies-escherichia-gene-subset-plasmid.pdf', 
    width=7.5, height=3)
p + scale_color_manual(values=c("#f1a340", "#998ec3"))+
  theme(axis.text=element_text(colour="black"))
dev.off()

# Any correlation between rough date of emergence and structural diversity?
GENUS = "escherichia"
GENES.OF.INTEREST = c("CMY-2", "CTX-M-15", 
                      "CTX-M-65", "KPC-2", "NDM-1",
                      "OXA-1",  "TEM-1")
entropy.df = read.csv(paste0('../output/positional-entropies-', GENUS, '.csv'), sep=',', header=F, stringsAsFactors = F)

colnames(entropy.df) = c("gene", "contig", "direction","N", "position", "entropy")

entropy.df$gene.N = paste(entropy.df$gene, entropy.df$N, sep="\n")

N.df = entropy.df %>% group_by(gene, contig) %>%
  summarise(n=paste0("n=",unique(N)))

metadata.df = read.csv('../data/metadata.csv',
                       header=T,
                       row.names = 1,
                       stringsAsFactors = F)

# Earliest date for a gene
earliestObservation <- function(data.df, gene){
  min(data.df$Year[grepl(gene, data.df$Gene.hits)])
}
# Earliest observations for Escherichia
earliest.observations.dataset = as.numeric(sapply(GENES.OF.INTEREST, 
       function(x) earliestObservation(metadata.df[which(tolower(metadata.df$TaxGenus)==GENUS),], 
                   x)))
# From PubMed search
earliest.observations.pubmed = c(1990, 2001, 2009, 2003, 2009, 1983, 1962)
# These are correlated ~ 0.9
cor.test(earliest.observations.dataset, 
         earliest.observations.pubmed)

# Does it correlate with entropy at all?
sum.of.entropies = entropy.df %>% filter(gene %in% GENES.OF.INTEREST,
                          contig %in% "plasmid") %>% 
      group_by(direction, gene) %>%
  summarise(entropy=sum(entropy))
entropies.by.direction = pivot_wider(sum.of.entropies, names_from = "direction", values_from="entropy" )
entropies.by.direction$earliest.observation = earliest.observations.dataset
entropies.by.direction$earliest.observation.pubmed = earliest.observations.pubmed
entropies.by.direction$discrepancy = entropies.by.direction$downstream-entropies.by.direction$upstream
plot(entropies.by.direction$earliest.observation, entropies.by.direction$downstream-entropies.by.direction$upstream )
plot(entropies.by.direction$earliest.observation.pubmed, entropies.by.direction$downstream-entropies.by.direction$upstream )


cor.test(entropies.by.direction$earliest.observation.pubmed, 
         entropies.by.direction$downstream-entropies.by.direction$upstream ,
         method="spearman")



p.entropy.discrepancy.and.date = ggplot(entropies.by.direction, aes(earliest.observation.pubmed, discrepancy))+
  geom_point(shape=21, fill="black")+
  ggrepel::geom_text_repel(aes(label=gene), size=2)+
  theme_basic()+
  xlab("Earliest occurrence (pubmed)")+
  ylab("Area under entropy curve difference\n(downstream - upstream)")+
  annotate(geom="text", x=2000, y=25,
           label=paste0("r=",formatC( cor.test(entropies.by.direction$earliest.observation.pubmed,
                                               entropies.by.direction$discrepancy,
                                               method="spearman")$estimate, digits=2 )))


p.upstream = ggplot(entropies.by.direction, aes(earliest.observation.pubmed, upstream))+
  geom_point(shape=21, fill="black")+
  ggrepel::geom_text_repel(aes(label=gene), size=2)+
  theme_basic()+
  xlab("Earliest occurrence (pubmed)")+
  ylab("Area under entropy curve \n(upstream)")+
  annotate(geom="text", x=2000, y=90,
           label=paste0("r=",formatC( cor.test(entropies.by.direction$earliest.observation.pubmed,
                          entropies.by.direction$upstream,
                          method="spearman")$estimate, digits=2 )))
p.downstream = ggplot(entropies.by.direction, aes(earliest.observation.pubmed, downstream))+
  geom_point(shape=21, fill="black")+
  ggrepel::geom_text_repel(aes(label=gene), size=2)+
  theme_basic()+
  xlab("Earliest occurrence (pubmed)")+
  ylab("Area under entropy curve \n(downstream)")+
  annotate(geom="text", x=2000, y=70,
           label=paste0("r=",formatC( cor.test(entropies.by.direction$earliest.observation.pubmed,
                                               entropies.by.direction$downstream,
                                               method="spearman")$estimate, digits=2 )))

pdf('../manuscript/figs/fig-entropies-escherichia-gene-subset-plasmid-dates.pdf', width=3, height=3)
p.entropy.discrepancy.and.date
dev.off()


pdf('../manuscript/figs/fig-entropies-escherichia-gene-subset-plasmid-dates.pdf',
    width=8, height=3)
cowplot::plot_grid(p.upstream+ggtitle("(a)"),
                   p.downstream+ggtitle("(b)"),
                     p.entropy.discrepancy.and.date+ggtitle("(c)"),
                   nrow=1)
dev.off()
