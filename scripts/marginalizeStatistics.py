import subprocess
import itertools as iter
import json
import convertPangraphToBlockList as cp
import mostFrequentPathGFA as mf
import pandas as pd
import glob
import re
import argparse

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

	#for k, v in paths.items():
	#	if v in block_dict[k]
	#	shared_seq += v[3]-v[2] if 
	#print(total_blocks)
	# Length of shared sequence	


def main():
	args = get_options()
#	example_strains = ['NZ_CP081345.1', 'NZ_CP074187.1'] # get from unique paths list
	GENE  = args.gene
	#unique_paths_file = '../../data/2023-02-05-'+GENE+'-mmseqs2-polish-u5000-d5000/'+GENE+'-mmseqs2-polish.all_u5000_d5000_pangraph.gfa.unique_paths.txt'
	pangraph_json = args.pangraph

	if args.subset!='':
		subset_ids = [line.strip('\n') for line in open(args.subset, 'r').readlines()]
		marginalizePangraph(pangraph_json, 'tmp', subset_ids)
		pangraph_json = 'tmp.json'

	# Get unique paths 
	# first export pangraph to gfa
	exportPangraph(pangraph_json, prefix='tmp')
	# then get unique paths
	mf.computePaths('tmp.gfa')

	unique_paths_file = 'tmp.gfa.unique_paths.txt'

	# Should compute these output distances directly from the fasta file...
	output_dists = '../../data/2023-02-05-'+GENE+'-mmseqs2-polish-u5000-d5000/'+GENE+'-mmseqs2-polish.all_u5000_d5000_pangraph.json.output_dists.csv' # contains SNPs
	dist_df = pd.read_csv(output_dists)
	#print(dist_df)
	idx = sorted(set(dist_df['seq1']).union(dist_df['seq2']))
	# reshape into matrix
	dist_mat = (dist_df.pivot(index='seq1', columns='seq2', values='snps')
	   .reindex(index=idx, columns=idx)
	   .fillna(0, downcast='infer')
	   .pipe(lambda x: x+x.values.T)
	 )
	#print(dist_mat)

	isolates = [line.strip('\n').split('\t')[0].split(',')[0] for line in open(unique_paths_file, 'r').readlines()]
	# Marginalize to just these unique-path isolates
	marginalizePangraph(pangraph_json, 'test', isolates)
	# Then marginalize again to get all of the pairwise combinations
	#marginalize_output_dir = '../../data/2023-02-05-'+GENE+'-mmseqs2-polish-u5000-d5000/marginalize-output'
	marginalize_output_dir = 	args.outputdir+'/marginalize'
	marginalizePangraph('test.json', marginalize_output_dir)

	with open(args.outputdir+'/'+GENE+'-dists.csv', 'w') as f:
		f.write('id1,id2,breakpoints,shared.seq,snps.gene\n')
		total_files = len(glob.glob(marginalize_output_dir+'/*json'))
		i = 1
		for json in glob.glob(marginalize_output_dir+'/*json'):
			strains = re.sub('.json', '', re.sub('.*/', '', json)).split('-')
			f.write(','.join(strains)+','+','.join([str(x) for x in blockStatistics(json)])+','+str(dist_mat[strains[0]].loc[strains[1]])+'\n')
			print(str(i)+'/', total_files)
			i += 1

#	print(blockStatistics('output/NZ_AP021952.1-NZ_CP074187.1.json'))

#	for a, b in iter.combinations(isolates, 2):
#		print(a, b)
	#	marginalizePangraph(pangraph_json, strains, 'test.json')
	#	exportPangraph('test.json', 'test')



	#marginalizePangraph(pangraph_json, strains, 'test.json')
	

if __name__=="__main__":
	main()