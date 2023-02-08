# Subset metadata

import pandas as pd
import argparse
# Use metadata to subset to strains with:
# - a particular gene
# - a particular genus
# - a particular species
# - chromosome/plasmid

# Functions should be callable from other scripts i.e.
# don't write all outputs
def get_options():
    parser = argparse.ArgumentParser(description='subset metadata')
    parser.add_argument('--metadata', help='input metadata csv', type=str, default='../data/metadata.csv') 
    parser.add_argument('--genefamily', help='gene', type=str, required=False, default='')
    parser.add_argument('--TaxGenus', help='genus', type=str, required=False, default='')
    parser.add_argument('--Contig', help='contig type (chromosome/plasmid)', type=str, required=False, default='')
    return parser.parse_args()

def subset(metadata_file, values_filter, gene_family=''):#gene_family='', genus='', species='', contig=''):
	# We've passed a dictionary of filter values
	metadata_df = pd.read_csv(metadata_file, index_col=0) 
	if values_filter!={}:
		metadata_df = metadata_df.loc[(metadata_df[list(values_filter)] == pd.Series(values_filter)).all(axis=1)]
	if gene_family!='':
		metadata_df = metadata_df.loc[ [gene_family in x for x in metadata_df['Gene.hits']]]
	return(list(metadata_df.index))


def main():
	args = get_options()
	values_filter = {}
	if args.TaxGenus!='':
		values_filter['TaxGenus'] = args.TaxGenus
	if args.Contig!='':
		values_filter['Contig'] = args.Contig


	subset_list = subset(args.metadata, values_filter, args.genefamily)
	for s in subset_list:
		print(s)

if __name__=="__main__":
	main()