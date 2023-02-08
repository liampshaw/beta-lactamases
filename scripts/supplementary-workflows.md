Notes on supplementary workflows:

# Positional entropy

We calculate the positional entropy for data files for each gene (downloaded from cluster run)
```
while read f;
do
# get gene start/end positions
makeblastdb -in ../../data/2023-01-18-$f-mmseqs2-polish-u5000-d5000/$f-mmseqs2-polish.all_u5000_d5000.fa -dbtype 'nucl'
grep -A 1 $f" " ../data/genes.fasta > tmp.fa
blastn -max_target_seqs 10000 -query tmp.fa -db ../../data/2023-01-18-$f-mmseqs2-polish-u5000-d5000/$f-mmseqs2-polish.all_u5000_d5000.fa -outfmt "6 sseqid sstart send" > ../../data/2023-01-18-$f-mmseqs2-polish-u5000-d5000/$f-mmseqs2-polish.all_u5000_d5000.blast_hits_gene.txt
done < ../data/genes.txt 

while read f;  
do
geneblock=$(cat ../../data/2023-01-18-$f-mmseqs2-polish-u5000-d5000/$f-mmseqs2-polish.all_u5000_d5000.gene_block.txt); python positionalEntropy.py --json ../../data/2023-01-18-$f-mmseqs2-polish-u5000-d5000/$f-mmseqs2-polish.all_u5000_d5000_pangraph.json --genelocations ../../data/2023-01-18-$f-mmseqs2-polish-u5000-d5000/$f-mmseqs2-polish.all_u5000_d5000.blast_hits_gene.txt --name $f --normalise
done < ../data/genes.txt > entropy-upstream-downstream.txt
```


# Gene variants

A script `assignGeneVariants.py` computes the enzyme variant for a given enzyme. It uses the translated protein to assign whether it is a known variant (with a name), an unnamed variant, or a truncated protein (i.e. non-functional - either genuinely or because of misassembly). This is useful for plotting a NJ tree of the gene and seeing the variation. 

```
while read f;
do
python assignGeneVariants.py --gene $f --outputprefix ../../data/2023-02-05-$f-mmseqs2-polish-u5000-d5000/$f-variants --inputprefix ../../data/2023-02-05-$f-mmseqs2-polish-u5000-d5000/$f-mmseqs2-polish.all_u5000_d5000
Rscript plot-NJ-tree.R ../../data/2023-02-05-$f-mmseqs2-polish-u5000-d5000/$f-variants.aln ../../data/2023-02-05-$f-mmseqs2-polish-u5000-d5000/$f-variants-NJ-tree.pdf 6 8
echo $f
done < ../data/genes.txt
```

# Marginalize statistics

Compute marginalized statistics for pairwise comparisons of unique paths. For example, for Salmonella plasmids with TEM-1:

```
 python subsetMetadata.py --genefamily TEM --TaxGenus Salmonella --Contig plasmid > tmp.list.txt
python marginalizeStatistics.py --pangraph ../../data/2023-02-05-TEM-1-mmseqs2-polish-u5000-d5000/TEM-1-mmseqs2-polish.all_u5000_d5000_pangraph.json  --subset tmp.list.txt --gene TEM-1
```

E.g. to look at Escherichia plasmids for all genes

```
python subsetMetadata.py --TaxGenus Salmonella --Contig plasmid > tmp.salmonella.plasmids.txt
while read f;
do
	python marginalizeStatistics.py --pangraph ../../data/2023-02-05-$f-mmseqs2-polish-u5000-d5000/$f-mmseqs2-polish.all_u5000_d5000_pangraph.json  --subset tmp.salmonella.plasmids.txt --gene $f --outputprefix salmonella-$f
done < ../data/genes.txt
