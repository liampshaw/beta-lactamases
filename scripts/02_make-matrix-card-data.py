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

def makePA(a_df):
    """makes presence absence matrix from df"""
    presence_absence = np.zeros((len(a_df), len(unique_genes)))
    for i in range(0, len(a_df)):
        hits = a_df.loc[a_df.index[i]]['all_hits_list']
        for h in hits:
            if h!='':
                presence_absence[i, unique_genes.index(h)] = 1
    # Turn into a dataframe
    converted_pa_df = pd.DataFrame(presence_absence,
        index=a_df.index,
        columns=unique_genes)
    return converted_pa_df



# Write to file
pa_df = makePA(df_chrom_plas)
pa_df.to_csv(DATA_DIR+'CARD-all-hits-presence-absence-chrom-plasmid.csv')

# Make for plasmids and chromosomes
pa_plas_df = makePA(df_plas)
pa_plas_df.to_csv(DATA_DIR+'CARD-all-hits-presence-absence-ncbi_plasmid.csv')
pa_chrom_df = makePA(df_chrom)
pa_chrom_df.to_csv(DATA_DIR+'CARD-all-hits-presence-absence-ncbi_chromosome.csv')
# Below not working atm
#df_contig = df[df['data_source']=='ncbi_contig' ]
#pa_contig_df = makePA(df_contig)
#pa_contig_df.to_csv(DATA_DIR+'CARD-all-hits-presence-absence-ncbi_contig.csv')
