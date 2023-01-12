#!/bin/bash

while read gene;
do
	echo $gene
	grep -A 1 $gene ../data/nucleotide_fasta_protein_homolog_model.fasta  | awk '$1!="--"' | sed -e 's/.*|/>/g' | sed -e 's/ .*//g' > ../data-processing/"$gene".fasta
	#clustalo -i ../data-processing/"$gene".fasta -o ../data-processing/"$gene".aln
	mafft --auto ../data-processing/"$gene".fasta > ../data-processing/"$gene".aln
	snp-dists -a -b ../data-processing/"$gene".aln > ../data-processing/"$gene"_diffs_aln.tsv # differences including ins/del (not just SNPs)
done < ../beta-lactamases.txt
