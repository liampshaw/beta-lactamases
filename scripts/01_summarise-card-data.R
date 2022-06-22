require(dplyr)
require(ggplot2)
require(cowplot)

DATA_DIR = '../data-processing/'

d <- read.csv('../CARD/index-for-model-sequences.txt',
              header=T,
              stringsAsFactors = F,
              sep='\t')
d.strict <- d[which(d$data_type=="ncbi_contig"),] # allow both rgi_criteria - Perfect or Stricti
d.strict$genus <- gsub(" .*", "", d.strict$species_name)
d.gene.genus <- d.strict %>% group_by(card_short_name, genus) %>%
  summarise(n=length(genus))
d.gene.genus$card_gene_group <- gsub("-[0-9].*", "", d.gene.genus$card_short_name)

# Clinically interesting beta-lactamases
beta.lactamases <- as.character(read.csv("../beta-lactamases.txt", header=F, stringsAsFactors = F)$V1 )
gene.groups <- beta.lactamases

bl.df <- d.gene.genus[which(d.gene.genus$card_gene_group %in% beta.lactamases),]


bl.df.summary <- bl.df %>% group_by(card_gene_group) %>%
  summarise(n=sum(n))
# order by number of genomes
bl.df.summary <- bl.df.summary[order(bl.df.summary$n, decreasing = T),]

bl.df$card_gene_group <- ordered(bl.df$card_gene_group,
                                                levels=bl.df.summary$card_gene_group)

genera.sums <- bl.df %>% group_by(genus) %>%
  summarise(n=sum(n))
top.genera <- genera.sums$genus[which(genera.sums$n>1000)]
bl.df$genus.simplified <- sapply(bl.df$genus,
                                            function(x) ifelse(x %in% top.genera, x, "_Other"))
genus.colours <- RColorBrewer::brewer.pal(n=length(top.genera), name="Set3")
genus.colours <- c("grey50", genus.colours)

pdf(paste0(DATA_DIR,'clinical-betalactamase-groups.pdf'), width=12, height=4)
p1 <- ggplot(bl.df, aes(card_gene_group, n, fill=genus.simplified))+
  geom_bar(stat="identity", position="stack")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text=element_text(colour="black"))+
  xlab("")+
  ylab("NCBI genomes (n)")+
  scale_fill_manual(values=genus.colours)+
  coord_flip()+
  labs(fill="Genus")+
  theme(legend.position = c(0.9,0.55))
small.groups <- bl.df.summary$card_gene_group[which(bl.df.summary$n<1000)]
p2.zoom <-  ggplot(bl.df[which(bl.df$card_gene_group %in% small.groups),], aes(card_gene_group, n, fill=genus.simplified))+
  geom_bar(stat="identity", position="stack")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text=element_text(colour="black"))+
  xlab("")+
  ylab("NCBI genomes (n)")+
  scale_fill_manual(values=genus.colours)+
  coord_flip()+
  theme(legend.position = "none")+
  ylim(c(0,1000))
cowplot::plot_grid(p1, p2.zoom, rel_widths = c(2,1) )
dev.off()

png(paste0(DATA_DIR,'clinical-betalactamase-groups.png'), width=12*300, height=4*300, units = "px", res=300)
cowplot::plot_grid(p1, p2.zoom, rel_widths = c(2,1) )
dev.off()

# Summarise the prevalences
bl.df.per.gene <- bl.df %>% group_by(card_short_name) %>% 
  summarise(n=sum(n))
bl.df.per.gene$group <- gsub("-[0-9].*", "", bl.df.per.gene$card_short_name)
write.csv(bl.df.per.gene, file=paste0(DATA_DIR, 'gene-prevalence-ncbi_wgs.csv'), row.names = F)

d.plasmid <- d[which(d$data_type=="ncbi_plasmid"),] # allow both rgi_criteria - Perfect or Stricti
d.plasmid$genus <- gsub(" .*", "", d.plasmid$species_name)
d.gene <- d.plasmid %>% group_by(card_short_name) %>%
  summarise(n=length(genus))
d.gene$group <- gsub("-[0-9].*", "", d.gene$card_short_name)
write.csv(d.gene, file=paste0(DATA_DIR,'gene-prevalence-ncbi_plasmid.csv'), row.names = F)

# Chromosome
d.chromosome <- d[which(d$data_type=="ncbi_chromosome"),] # allow both rgi_criteria - Perfect or Stricti
d.chromosome$genus <- gsub(" .*", "", d.chromosome$species_name)
d.gene <- d.chromosome %>% group_by(card_short_name) %>%
  summarise(n=length(genus))
d.gene$group <- gsub("-[0-9].*", "", d.gene$card_short_name)
write.csv(d.gene, file=paste0(DATA_DIR,'gene-prevalence-ncbi_chromosome.csv'), row.names = F)

