# downloads all chromosomes/plasmids which have a gene present in them
# going off a cluster of genomes
# conda env: scipy

import pandas as pd
import scipy.cluster.hierarchy as hc


pa_df = pd.read_csv('CARD-all-hits-presence-absence-chrom-plasmid.csv',index_col=0)

genome_df = pd.read_csv('card-genomes.txt', sep='\t', index_col=0)

def getNeighbours(family_name, gene_name, distance=10, dataset='ncbi_wgs'):
    """get the neighbours using clustering within a certain distance threshold of the gene"""
    snps = pd.read_csv('../card-data/'+family_name+'_diffs_aln.tsv', sep='\t', index_col=0)
    #plt.close()
    df = snps
    # and cluster
    linkage_mat = hc.linkage(df, method='average')
    # Get clusters
    maximum_distance = distance
    clusters = hc.fcluster(linkage_mat, maximum_distance, criterion='distance')
    cluster_of_gene = clusters[list(df.index).index(gene_name)]
    names = [list(df.index)[i] for i in [x for x in range(len(clusters)) if clusters[x]==cluster_of_gene]]
    return(names)

# Gene clusters
gene_cluster_plas = getNeighbours("CTX-M", 'CTX-M-55', 25, dataset='ncbi_plasmid')
gene_cluster_chrom = getNeighbours("CTX-M", 'CTX-M-55', 25, dataset='ncbi_chromosome')
# combine
gene_cluster = list(set(gene_cluster_plas+gene_cluster_chrom))

#accessions = []:xv
#for g in gene_cluster:
#    if g in pa_df.columns:
#        results =  list(pa_df.index[pa_df[g]!=0])
#        accessions += results
plasmids = list(genome_df[genome_df['data_source']=='ncbi_plasmid'].index)

plasmid_df = pa_df.loc[plasmids]
accessions = {}
for g in gene_cluster_plas:
    if g in plasmid_df.columns:
        results = list(plasmid_df.index[plasmid_df[g]!=0])
        for r in results:
            if r in accessions.keys():
                accessions[r].append(g)
            else:
                accessions[r] = [g]

with open('CTX-M-55/accs_plasmid.txt', 'w') as f:
    for k, v in accessions.items():
        for value in v:
            _ = f.write('%s,%s\n' % (value, k))

chroms = list(genome_df[genome_df['data_source']=='ncbi_chromosome'].index)

chrom_df = pa_df.loc[chroms]
accessions = {}
for g in gene_cluster_plas:
    if g in chrom_df.columns:
        results = list(chrom_df.index[chrom_df[g]!=0])
        for r in results:
            if r in accessions.keys():
                accessions[r].append(g)
            else:
                accessions[r] = [g]

with open('CTX-M-55/accs_chromosome.txt', 'w') as f:
    for k, v in accessions.items():
        for value in v:
            _ = f.write('%s,%s\n' % (value, k))
