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
