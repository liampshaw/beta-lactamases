correlationForGene <- function(GENE){
  dir = paste0('/Users/Liam/Dropbox/_Projects/pangraph/beta-lactamase/scripts/', GENE)
  distances = read.csv(paste0(dir, '/pangraph_all_u5k_d5k.gfa.output_dists.csv'),
                       header=T,
                       stringsAsFactors = F)
  #return(c(cor(distances$dist.up, distances$snps)^2,
  #       cor(distances$dist.down, distances$snps)^2))
  return(cor(distances$dist.up, distances$dist.down)^2)
}

plotCorrelationForGene <- function(GENE){
  dir = paste0('/Users/Liam/Dropbox/_Projects/pangraph/beta-lactamase/scripts/', GENE)
  distances = read.csv(paste0(dir, '/pangraph_all_u5k_d5k.gfa.output_dists.csv'),
                       header=T,
                       stringsAsFactors = F)
  #return(c(cor(distances$dist.up, distances$snps)^2,
  #       cor(distances$dist.down, distances$snps)^2))
  p = ggplot(distances, aes(dist.up, dist.down))+
    geom_point()
  return(p)
}

meanBreakpointForGene <- function(GENE){
  dir = paste0('/Users/Liam/Dropbox/_Projects/pangraph/beta-lactamase/scripts/', GENE)
  distances = read.csv(paste0(dir, '/pangraph_all_u5k_d5k.gfa.output_dists.csv'),
                       header=T,
                       stringsAsFactors = F)
   return(mean(na.omit(distances[which(distances$snps==0),]$dist.up-distances[which(distances$snps==0),]$dist.down)))
  #return(c(cor(distances$dist.up, distances$snps)^2,
  #       cor(distances$dist.down, distances$snps)^2))
  
}


gene.list = read.csv('../scripts/gene-lists.csv', 
                     header=F, 
                     stringsAsFactors = F,
                     col.names = c("family", "gene"))
genes = c("CMY-2"  ,  "CTX-M-15" , "CTX-M-65", "GES-24"  , "IMP-4",
          "NDM-1"  ,  "OXA-48" ,"OXA-10" ,"OXA-58"  , "OXA-23"  , "PER-1"  ,
          "SHV-134" , "TEM-1"   , "VEB-1", "VIM-1" ,"KPC-2"  )
genes = sort(genes)
results = matrix(nrow=0, ncol=2)
for (gene in genes){
  results = rbind(results, c(gene, meanBreakpointForGene(gene)))
}
results = data.frame(results)
results$X2 = as.numeric(results$X2)
gene.types = read.csv('gene-classes.csv',
                      header=T, stringsAsFactors = F, row.names = 1)
results$family = gsub("-[0-9]", "", gsub("-[0-9][0-9]$", "", results$X1))
results$family[results$family=="SHV34"] = "SHV"
results$class = gene.types[results$family,"class"]
ggplot(results, aes(X1, X2, group=family))+
  geom_bar(stat="identity", fill="black", width=0.9)+
  facet_grid(. ~ class, scales = "free", space = "free")+
  theme_bw()+
  xlab("Gene")+
  ylab("Mean upstream - downstream")+
  theme(panel.grid = element_blank())+
  theme(axis.text.x=element_text(angle=45, hjust=1))




# Proportion chromosome/plasmid
for (gene in genes){
  counts = table(read.csv(paste0('../scripts/', gene, "_within_25_diffs_accs.txt"), 
            header=F, 
            stringsAsFactors = F)$V2)
  print(c(gene, counts))
}

