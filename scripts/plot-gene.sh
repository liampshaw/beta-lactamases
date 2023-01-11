GENE=$1
cd $GENE
geneBlock=$(blastn -query gene.fa -db pangraph_all_u5k_d5k.fa -outfmt 6 | cut -f 2)

/Applications/Bandage.app/Contents/MacOS/Bandage image pangraph_all_u5k_d5k.gfa.coloured.gfa \
	pangraph_all_u5k_d5k.gfa.png \
	--height 4000 --width 7000 \
	--colour custom

Rscript ../plot-blocks.R \
	pangraph_all_u5k_d5k.gfa.blocks.csv \
	$geneBlock pangraph_all_u5k_d5k.gfa.png \
	pangraph_graph_plot.pdf

cd ../
