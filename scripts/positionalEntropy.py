# Calculate positional entropy along a pangraph away from a central gene
import json
import argparse
from random import seed
from random import randint
import convertPangraphToBlockList as cp
import numpy as np
from math import e

# --json ../output/GES-24/GES-24-mmseqs2-polish.all_u5000_d5000_pangraph.json
# --geneblock LXDNJCUWIP


def get_options():
    parser = argparse.ArgumentParser(description='Calculate positional entropy along a pangraph away from a central gene')
    parser.add_argument('--json', help='Input json file', type=str)
    parser.add_argument('--name', help='Prefix (if concatenating multiple files together)', type=str, required=False, default='')
    #parser.add_argument('--geneblock', help='Block to centre the analysis on', type=str)
    return parser.parse_args()

def shannonEntropy(labels, base=None):
  value,counts = np.unique(labels, return_counts=True)
  norm_counts = counts / counts.sum()
  base = e if base is None else base
  return -(norm_counts * np.log(norm_counts)/np.log(base)).sum()

def main():
    args = get_options()
    json_file = str(args.json)
    with open(json_file, 'r') as f:
      pangraph_json = json.load(f)

    genome_dict = cp.extractPaths(pangraph_json['paths'])
    # Don't care about exact lengths of blocks, will just use the consensus (average?) length
    block_dict = {x['id']:len(x['sequence']) for x in pangraph_json['blocks']}
    block_nums = {x['id']:i for i, x in enumerate(pangraph_json['blocks'])}
    block_dict_num = {block_nums[x['id']]:len(x['sequence']) for x in pangraph_json['blocks']}

    #print(block_dict)
    #print(block_nums)
    #print(block_dict_num)
    #print(genome_dict)

    #gene_block_num = block_nums[args.geneblock]
    #print('gene block is:', gene_block_num)

    possible_paths = [[block_nums[x[0]] for x in genome_dict[y]] for y in genome_dict.keys()]
    #print(possible_paths)

    # Upstream paths
    #downstream_paths = [x[x.index(gene_block_num)+1:] for x in possible_paths]
    #print(downstream_paths)
    #upstream_paths = [x[:x.index(gene_block_num)] for x in possible_paths]
    #print(upstream_paths) # could reverse so going away from gene?

    # Now use the consensus lengths of blocks to create vectors (say at 100bp scale) for each genome path
    #for d in downstream_paths:
#        print("new genome")
    #    counter = 0
    #    for x in d:
    #        counter += block_dict_num[x]
    #        print(str(x), counter)

    # Use actual positions in genomes
    genomes_as_int = [None]*len(genome_dict.keys())
    genome_i = 0
    for g, path in genome_dict.items():
        #print(g)
        genome_as_int = [None] * max([p[3] for p in path])

        genome_as_int = [item for sublist in [[block_nums[p[0]]]*(p[3]-p[2]) for p in path ] for item in sublist]
        #for p in path:
            #print(genome_as_int[p[2]:p[3]])
            #genome_as_int[p[2]:p[3]] = block_nums[p[0]]*(p[3]-p[2])
            #print(block_nums[p[0]], p[2], p[3])
        #print(genome_as_int)
        genomes_as_int[genome_i] = genome_as_int
        genome_i += 1


    #print(genomes_as_int)
    # Table of each position
    for i in range(0, 10000, 100):
        block_vector = [g[i] for g in genomes_as_int]# this is the vector we want
        if args.name=='':
            print(i, shannonEntropy(block_vector))
        else:
            print(args.name, i, shannonEntropy(block_vector))

    # Get upstream and downstream paths
    # Can use start/end points of blocks as breakpoints

    # Can we convert the blocks to ASCII characters? Perhaps too many...
    # Or integers is probably best




    # Convert to a table of counts for Shannon entropy
    # a, b, c, d, e










if __name__=="__main__":
    main()
