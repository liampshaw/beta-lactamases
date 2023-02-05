# downloads all chromosomes/plasmids which have a gene present in them
# going off a cluster of genomes
# conda env: scipy


import pandas as pd
import scipy.cluster.hierarchy as hc
import argparse
import warnings
warnings.filterwarnings("ignore") # because seaborn gives warning as we are passing a distance matrix, even though we are doing this deliberately

DATA_DIR = '../data-processing/'
def get_options():
    parser = argparse.ArgumentParser(description='download all sequences containing a sequence variant similar to a named beta-lactamase')
    parser.add_argument('family', help='Name of enzyme family (e.g. CTX-M)', type=str)
    parser.add_argument('gene', help='Name of sequence (e.g. CTX-M-55)', type=str)
    parser.add_argument('threshold', help='Distance threshold (diffs in clustalo MSA)', type=int)
    parser.add_argument('--singlehits', help='Only return accessions with a single hit', action='store_true', required=False, default=False)
    return parser.parse_args()


def getNeighbours(family_name, gene_name, distance=10):
    """get the neighbours using clustering within a certain distance threshold of the gene"""
    snps = pd.read_csv(DATA_DIR+'/fasta_diffs/'+family_name+'_diffs_aln.tsv', sep='\t', index_col=0)
    #plt.close()
    df = snps
    # 
    variant_names = list(snps.index[snps[gene_name]<10])
    # and cluster
    #linkage_mat = hc.linkage(df, method='average')
    # Get clusters
    #maximum_distance = distance
    #clusters = hc.fcluster(linkage_mat, maximum_distance, criterion='distance')
    #cluster_of_gene = clusters[list(df.index).index(gene_name)]
    #names = [list(df.index)[i] for i in [x for x in range(len(clusters)) if clusters[x]==cluster_of_gene]]
    #print(names)
    return(variant_names)

def main():
    args = get_options()

    # Gene clusters
    gene_cluster = getNeighbours(args.family, args.gene, args.threshold)

    # Plasmid accessions
    plasmid_df = pd.read_csv(DATA_DIR+'/CARD-all-hits-ncbi_plasmid.csv',index_col=0)

    # Dict of accessions and whether chromosome/plasmid
    accessions = {}
    for g in gene_cluster:
        if g in plasmid_df.columns:
            if args.singlehits==True:
                results = list(plasmid_df.index[plasmid_df[g]==1])
            elif args.singlehits==False:
                results = list(plasmid_df.index[plasmid_df[g]!=0])
            for r in results:
                if r not in accessions.keys():
                    accessions[r] = 'plasmid'


    chrom_df = pd.read_csv(DATA_DIR+'/CARD-all-hits-ncbi_chromosome.csv',index_col=0)
    for g in gene_cluster:
        if g in chrom_df.columns:
            if args.singlehits==True:
                results = list(chrom_df.index[chrom_df[g]==1])
            elif args.singlehits==False:
                results = list(chrom_df.index[chrom_df[g]!=0])
            for r in results:
                if r not in accessions.keys():
                    accessions[r] = 'chromosome'

    # Write to stdout
    #with open(args.gene+'_within_'+str(args.threshold)+'_diffs_accs.txt', 'w') as f:
    #    for k, v in accessions.items():
    for k, v in accessions.items():
        _ = print('%s,%s' % (k, v))

if __name__ == "__main__":
    main()
