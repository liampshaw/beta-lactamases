genes = read.csv('../data/genes.txt', header=F, stringsAsFactors = F)$V1


abricate = read.csv('abricate-IS-hits.tsv', header=F, stringsAsFactors = F, sep='\t')

abricate$gene=gsub("-mmseqs2.*", "", gsub(".*/", "", abricate$V1))



ggplot(abricate[which(abricate$V10>80),], aes(x=V4, group=gene))+
  geom_histogram()+
  facet_wrap(~gene, scales="free_y")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  xlab("position (bp)")+
  ylab("IS count")
