from Bio import SeqIO
import sys
seqs = SeqIO.to_dict(SeqIO.parse(snakemake.input[0], "fasta"))
db = snakemake.params.db
gene = snakemake.params.gene_name
print(gene)
with open(snakemake.output[0], 'w') as f:
	for seqid in seqs:
		if gene in seqs[seqid].description:
			if db=="CARD":
				gene_id = seqs[seqid].description.split("|")[-1].split(" ")[0]
			elif db=="NCBI":
				gene_id = seqs[seqid].description.split(",")[0].split(" ")[-1]
			if gene in gene_id:
				f.write(">"+gene_id+"\n")
				f.write(str(seqs[seqid].seq)+"\n")
