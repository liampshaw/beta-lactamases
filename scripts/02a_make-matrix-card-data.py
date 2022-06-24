# Make a matrix of the CARD data

import pandas as pd
import numpy as np

DATA_DIR = '../data-processing/'

df = pd.read_csv('../CARD/index-for-model-sequences.txt', sep='\t', index_col=0)
# replace NaN with empty string
df = df.fillna("")

df_perfect = df[df["rgi_criteria"]=="Perfect"]
df_strict = df[df["rgi_criteria"]=="Strict"]

# Make a list of perfect hits
#df = df.assign(perfect_hits_list=[x.split(', ') for x in df['perfect_hits']],
#                strict_hits_list=[x.split(', ') for x in df['strict_hits']])
#df = df.assign(all_hits_list = [list(set(a).union(set(b))) for a, b in zip(df['strict_hits_list'], df['perfect_hits_list'])])
# get unique gene names
unique_genes = sorted(list(set(df["card_short_name"])))

#
df_chrom = df[df['data_type']=='ncbi_chromosome' ]
df_chrom.index = range(len(df_chrom))
df_plas = df[df['data_type']=='ncbi_plasmid' ]
df_plas.index = range(len(df_plas))
df_chrom_plas = pd.concat([df_chrom, df_plas])
df_chrom_plas.index = range(len(df_chrom_plas))

def makePA(a_df, pa=True):
    """makes presence absence matrix from long form df"""
    accessions = sorted(list(set(a_df["ncbi_accession"])))
    presence_absence = np.zeros((len(accessions), len(unique_genes)))
    for i in range(0, len(a_df)):
        accession = a_df["ncbi_accession"][i]
        hit = a_df["card_short_name"][i]
        if pa==True:
            presence_absence[accessions.index(accession), unique_genes.index(hit)] = 1 # doesn't count >1 hit on sequence
        elif pa==False:
            presence_absence[accessions.index(accession), unique_genes.index(hit)] += 1 # counts all hits
    # Turn into a dataframe
    converted_pa_df = pd.DataFrame(presence_absence,
        index=accessions,
        columns=unique_genes)
    return converted_pa_df



# Write to file
print("Make chromosome-plasmid matrices...")
pa_df = makePA(df_chrom_plas)
pa_df.to_csv(DATA_DIR+'CARD-all-hits-presence-absence-chrom-plasmid.csv')
hits_df = makePA(df_chrom_plas, False)
hits_df.to_csv(DATA_DIR+'CARD-all-hits-ncbi_chromosome_and_plasmid.csv')

# Make for plasmids and chromosomes
print("Make plasmid matrices...")
pa_plas_df = makePA(df_plas)
pa_plas_df.to_csv(DATA_DIR+'CARD-all-hits-presence-absence-ncbi_plasmid.csv')
hits_df_plas = makePA(df_plas, False)
hits_df_plas.to_csv(DATA_DIR+'CARD-all-hits-ncbi_plasmid.csv')
print("Make chromosome matrices...")
pa_chrom_df = makePA(df_chrom)
pa_chrom_df.to_csv(DATA_DIR+'CARD-all-hits-presence-absence-ncbi_chromosome.csv')
hits_df_chrom = makePA(df_chrom, False)
hits_df_chrom.to_csv(DATA_DIR+'CARD-all-hits-ncbi_chromosome.csv')


# Below not working atm
#df_contig = df[df['data_source']=='ncbi_contig' ]
#pa_contig_df = makePA(df_contig)
#pa_contig_df.to_csv(DATA_DIR+'CARD-all-hits-presence-absence-ncbi_contig.csv')
