# Make a matrix of the CARD data

import pandas as pd
import numpy as np

DATA_DIR = '../data-processing/'

df = pd.read_csv('../CARD/card-genomes.txt', sep='\t', index_col=0)
# replace NaN with empty string
df = df.fillna("")

# Make a list of perfect hits
df = df.assign(perfect_hits_list=[x.split(', ') for x in df['perfect_hits']],
                strict_hits_list=[x.split(', ') for x in df['strict_hits']])
df = df.assign(all_hits_list = [list(set(a).union(set(b))) for a, b in zip(df['strict_hits_list'], df['perfect_hits_list'])])
# get unique gene names
unique_genes = sorted(list(set([x for subl in df['all_hits_list'] for x in subl])))
unique_genes.pop(0) # remove empty string

#
df_chrom = df[df['data_source']=='ncbi_chromosome' ]
df_plas = df[df['data_source']=='ncbi_plasmid' ]

df_chrom_plas = pd.concat([df_chrom, df_plas])
# make the p/a matrix
presence_absence = np.zeros((len(df_chrom_plas), len(unique_genes)))

for i in range(0, len(df_chrom_plas)):
    hits = df_chrom_plas.loc[df_chrom_plas.index[i]]['all_hits_list']
    for h in hits:
        if h!='':
            presence_absence[i, unique_genes.index(h)] = 1

# Turn into a dataframe
pa_df = pd.DataFrame(presence_absence,
    index=df_chrom_plas.index,
    columns=unique_genes)

# Write to file
pa_df.to_csv(DATA_DIR+'CARD-all-hits-presence-absence-chrom-plasmid.csv')
