# I am comparing to NZ_CP023845.1 (2009, plasmid) since it looks like it has a nice signal between length
# of longest subsequence and SNP density

GENE = "ndm"
plotGeneResults = function(
    GENE){
  dir = paste0('/Users/Liam/Dropbox/_Projects/pangraph/beta-lactamase/scripts/', GENE)
  test.dists = read.csv(#€'/Users/Liam/github/pangenome_evolution/sandbox/comparisons.csv', 
    paste0(dir, '/comparisons.csv'),
    header=T,
    stringsAsFactors = F)
  # break = shared sequence before first breakpoint in both directions
  # lcs = longest common subsequence, any number of insertions/deletions allowed
  # colnames(test.dists) = c("id1", "id2", "n.blocks.lcs", "length.lcs", "snps.lcs", "snp.density.lcs",
  #                          "n.blocks.break", "length.break", "snps.break","snp.density.break", "diffs.break", "diff.density.break",
  #                          "n.blocks.upstream", "length.upstream", "snps.upstream","snp.density.upstream", "diffs.upstream", "diff.density.upstream",
  #                          "n.blocks.downstream", "length.downstream", "snps.downstream","snp.density.downstream", "diffs.downstream", "diff.density.downstream",
  #                          "diffs.focal.block")
  #                          
  
  metadata = read.csv("2022-10-26-metadata-country-and-year.csv",
                      header=T,
                      stringsAsFactors = F,
                      row.names = 1)
  
  
  
  
  
  
  p.upstream = ggplot(test.dists, aes(L_upstream, y=as.factor(N_single_base_focal)), fill = factor(stat(quantile))) +
    geom_vline(xintercept = 0, linetype='dashed', colour='grey')+
    geom_density_ridges(
      aes(point_color = as.factor(N_single_base_focal)), 
      alpha = .2, point_alpha = 1, jittered_points = TRUE, scale = 0.6
    ) +
    #scale_point_color_hue(l = 40)+ 
    theme_bw()+
    theme(panel.grid=element_blank())+
    ylab("Single base differences in focal block\n(SNPs, small indels)")+
    xlab("Bases until first breakpoint (upstream)")+
    theme(axis.text=element_text(colour="black"))+
    theme(legend.position = "none")+
    scale_x_reverse()
  
  p.downstream = ggplot(test.dists, aes(L_downstream, y=as.factor(N_single_base_focal)), fill = factor(stat(quantile))) +
    geom_vline(xintercept = 0, linetype='dashed', colour='grey')+
    geom_density_ridges(
      aes(point_color = as.factor(N_single_base_focal)), 
      alpha = .2, point_alpha = 1, jittered_points = TRUE, scale = 0.6, size=0.2,
    ) +
    theme_bw()+
    theme(panel.grid=element_blank())+
    ylab("")+
    xlab("Bases until first breakpoint (downstream)")+
    theme(axis.text=element_text(colour="black"))+
    theme(legend.position = "none")
  #pdf(paste0(Sys.Date(), "-CTX-M-15-focal-block-first-breakpoints.pdf"), width=8, height=4)
  p.lengths = cowplot::plot_grid(p.upstream+ggtitle(GENE), p.downstream+ggtitle(""), align='v')
  #dev.off()
  
  # Other differences - larger scale?
  
  
  #cor.test(test.dists$L_upstream, test.dists$L_downstream, method="spearman")
  
  p.correlation = ggplot(test.dists,aes(L_upstream, L_downstream))+
    #geom_point(aes(colour=as.factor(N_single_base_focal)))+
    geom_point()+
    stat_smooth(method="lm", colour="black")+
    coord_fixed()+
    theme_bw()+
    annotate(geom="text", x=4000, y=1000,
             label=paste0("Pearson's corr=", round(cor(test.dists$L_upstream, test.dists$L_downstream), 
                                      3)))+
    theme(panel.grid = element_blank())
    #facet_wrap(~as.factor(N_single_base_focal), nrow=1)
  p = cowplot::plot_grid(p.lengths, p.correlation, nrow=2, rel_widths = c(2,1))
  return(p)
}

p.ndm = plotGeneResults("ndm")
ggsave(plot=p.ndm, file='../scripts/ndm/plot-length-correlations.pdf')
p.ctx = plotGeneResults("CTX-M-15")
ggsave(plot=p.ctx, file='../scripts/CTX-M-15//plot-length-correlations.pdf')

