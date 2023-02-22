import pandas as pd

metadata = pd.read_csv(snakemake.input.metadata, index_col=0)

GENES = ["CMY",\
			"CTX-M",\
			"GES",\
			"IMP",\
			"KPC",\
			"NDM",\
			"PER",\
			"SHV",\
			"TEM",\
			"VEB",\
			"VIM"]


for gene in GENES:
	with open("data/"+gene+".txt", 'w') as f:
		hits = metadata.loc[[(str(gene) in x) for x in list(metadata["Gene.hits"])]].index
		for h in hits:
			f.write("%s\n" % h)