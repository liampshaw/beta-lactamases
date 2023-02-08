import subprocess
import itertools as iter
import json
import convertPangraphToBlockList as cp
import mostFrequentPathGFA as mf
import pandas as pd
import glob
import re
import argparse
import os
import subsetMetadata as sm

def get_options():
    parser = argparse.ArgumentParser(description='Marginalize statistics on unique paths')
    parser.add_argument('--pangraph', help='input pangraph json', type=str) 
    parser.add_argument('--gene', help='gene', type=str) 
    parser.add_argument('--subset', help='file with subset of strains', type=str, required=False, default='') 
    parser.add_argument('--outputdir', help='outputdir', type=str, required=False, default='./') 
    return parser.parse_args()



def marginalizePangraph(pangraph_json, output_prefix, strain_list=''): # output_prefix is made as dir if no strain list
	if strain_list=='':
		marginalize_command = 'pangraph marginalize '+pangraph_json+' -o '+output_prefix
	else:
		strain_string = ','.join(strain_list)
		marginalize_command = 'pangraph marginalize --strains '+strain_string+' '+pangraph_json+' > '+output_prefix+'.json'
	print("Command:", marginalize_command)
	marginalize_output = subprocess.call(marginalize_command, shell=True)
	print("Return code:", marginalize_output)
	return


def exportPangraph(pangraph_json, prefix, output_dir='./'):
	export_command = 'pangraph export '+pangraph_json+' -p '+prefix+' -o '+output_dir
	print("Command:", export_command)
	export_output = subprocess.call(export_command, shell=True)
	print("Return code:", export_output)
	return

def blockStatistics(pair_json):
	'''Block statistics for a pairwise pangraph'''
	with open(str(pair_json), 'r') as f:
		pangraph = json.load(f)
	paths = cp.extractPaths(pangraph['paths'])
	# Number of shared blocks +1 = number of breakpoints
	block_dict = {k:[b[0] for b in paths[k]] for k in paths.keys() }
	names = list(paths.keys())
	pair_member_0 = names[0]# pick random member of pair
	pair_member_1 = names[1]
	shared_blocks = [b for b in block_dict[pair_member_0] if b in block_dict[pair_member_1]]
	breakpoints = len(shared_blocks)+1
	shared_seq = sum([abs(b[3]-b[2]) for b in paths[pair_member_0] if b[0] in shared_blocks])
	return([breakpoints, shared_seq])
	# Other statistics that could be added: 
	# - location of the breakpoints
	# - alignment distance within the shared sequence


def main():
	args = get_options()
	GENE  = args.gene
	pangraph_json = args.pangraph

	# If we have a subset pre-defined, then we first marginalize the pangraph 
	if args.subset!='':
		subset_ids = [line.strip('\n') for line in open(args.subset, 'r').readlines()]
		marginalizePangraph(pangraph_json, 'tmp', subset_ids)
		pangraph_json = 'tmp.json'

	# Get unique paths 
	# a) first export pangraph to gfa
	exportPangraph(pangraph_json, prefix='tmp')
	# b) then get unique paths using mostFrequentPathGFA.py (bad name for script...)
	mf.computePaths('tmp.gfa')
	# this give output:
	unique_paths_file = 'tmp.gfa.unique_paths.txt'

	# Should compute these output distances directly from the fasta file...
	# we get a matrix of SNP distances between all the genes
	output_dists = '../../data/2023-02-05-'+GENE+'-mmseqs2-polish-u5000-d5000/'+GENE+'-mmseqs2-polish.all_u5000_d5000_pangraph.json.output_dists.csv' # contains SNPs
	dist_df = pd.read_csv(output_dists)
	# convert long format to a matrix we can access with strain1,strain2 indices
	idx = sorted(set(dist_df['seq1']).union(dist_df['seq2']))
	# reshape into matrix
	dist_mat = (dist_df.pivot(index='seq1', columns='seq2', values='snps')
	   .reindex(index=idx, columns=idx)
	   .fillna(0, downcast='infer')
	   .pipe(lambda x: x+x.values.T)
	 )

	# Subset further, just to representatives of unique paths
	isolates = [line.strip('\n').split('\t')[0].split(',')[0] for line in open(unique_paths_file, 'r').readlines()]
	# Marginalize to just these unique-path isolates
	marginalizePangraph(pangraph_json, 'test', isolates)
	# Then marginalize to get all pairwise combinations
	marginalize_output_dir = 	args.outputdir+'/marginalize'
	marginalizePangraph('test.json', marginalize_output_dir)

	# Compute the pairwise block statistics for these pairs of isolates
	with open(args.outputdir+'/'+GENE+'-dists.csv', 'w') as f:
		f.write('id1,id2,breakpoints,shared.seq,snps.gene\n')
		total_files = len(glob.glob(marginalize_output_dir+'/*json'))
		i = 1
		for json in glob.glob(marginalize_output_dir+'/*json'):
			strains = re.sub('.json', '', re.sub('.*/', '', json)).split('-')
			f.write(','.join(strains)+','+','.join([str(x) for x in blockStatistics(json)])+','+str(dist_mat[strains[0]].loc[strains[1]])+'\n')
			print(str(i)+'/', total_files)
			i += 1

	# Remove json files
	for g in glob.glob(marginalize_output_dir+'/*json'):
		os.remove(g)



	

if __name__=="__main__":
	main()