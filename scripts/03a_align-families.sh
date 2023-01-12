#!/bin/bash

mkdir -p ../data-processing/fasta
mkdir -p ../data-processing/fasta_aln
mkdir -p ../data-processing/fasta_diffs
mkdir -p ../data-processing/plots

while read gene;
do
	echo $gene
	grep -A 1 $gene ../data/nucleotide_fasta_protein_homolog_model.fasta  | awk '$1!="--"' | sed -e 's/.*|/>/g' | sed -e 's/ .*//g' > ../data-processing/fasta/"$gene".fasta
	#clustalo -i ../data-processing/"$gene".fasta -o ../data-processing/"$gene".aln
	mafft --auto ../data-processing/"$gene".fasta > ../data-processing/fasta_aln/"$gene".aln
	snp-dists -a -b ../data-processing/"$gene".aln > ../data-processing/fasta_diffs/"$gene"_diffs_aln.tsv # differences including ins/del (not just SNPs)
done < ../beta-lactamases.txt
