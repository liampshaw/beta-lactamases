while read f;
do
	makeblastdb -in ../../data/2023-02-05-$f-mmseqs2-polish-u5000-d5000/$f-mmseqs2-polish.all_u5000_d5000.fa -dbtype 'nucl'
	grep -A 1 $f" " ../data/genes.fasta > tmp.fa
	blastn -query tmp.fa -db ../../data/2023-02-05-$f-mmseqs2-polish-u5000-d5000/$f-mmseqs2-polish.all_u5000_d5000.fa -outfmt '6 sseqid sstart send' > ../../data/2023-02-05-$f-mmseqs2-polish-u5000-d5000/$f-mmseqs2-polish.all_u5000_d5000.blast_hits_gene.txt 
done < ../data/genes.txt
