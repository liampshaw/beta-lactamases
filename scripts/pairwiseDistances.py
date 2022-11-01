import sys
from pkg_resources import to_filename
import json
import matplotlib.pyplot as plt
import operator
import pandas as pd
from itertools import combinations
from Bio import SeqIO
from random import sample
import argparse 
import os

#sys.path.append("..") # to use Marco's pre-existing scripts
#import src.py_utils as put
import blockFunctions as bf

def get_options():
    parser = argparse.ArgumentParser(description='make a coloured version of pangraph gfa')
    parser.add_argument('--json', help='Input pangraph json file', type=str)
    parser.add_argument('--focalblock', help='name of focal block')
    parser.add_argument('--panXdir', help='panX output directory')
    return parser.parse_args()


def sharedBlocksBeforeBreakpoint(a, b, starting_block, upstream=True, includeStartingBlock=False):
    """a, b are lists of blocks"""
    a_start = a.index(starting_block)
    b_start = b.index(starting_block)
    shared_blocks = []
    if includeStartingBlock:
        shared_blocks.append(starting_block)
    if upstream==False:
        i = 0
        while i<len(a)-a_start-1 and i<len(b)-b_start-1:
            i += 1
            if b_start<len(b):
                if a[i+a_start]==b[b_start+i]:
                    shared_blocks.append(b[b_start+i])
                else:
                    #if a[i+a_start] in shared_blocks: # append the block if it's a duplicated pre-existing one
                    #    shared_blocks.append(a[i+a_start]) # this means duplications 'don't count'
                    break
    elif upstream==True:
        i = 0
        while a_start-i > 0 and b_start-i > 0:
            i = i+1
            if a[a_start-i]==b[b_start-i]:
                shared_blocks.append(b[b_start-i])
            else:
                #if a[a_start+i] in shared_blocks:
                #    shared_blocks.append(a[a_start+i])
                break
    return(shared_blocks)

def editPositions(block_id, contig_id, block_dict):
    """Returns the positions that need to be edited for a given block in a given contig for a dictionary
    of blocks.
    Returns: dict of edits"""
    return {y:[x[1] for x in block_dict[block_id][y] if x[0]['name']==contig_id] for y in ['mutate', 'insert', 'delete']}

# Read in a block
def readBlockAln(gene_cluster_id, gene_cluster_dir="panX/vis/geneCluster/"):
    return SeqIO.to_dict(SeqIO.parse(gene_cluster_dir+gene_cluster_id+"_na_aln.fa",\
                                    format="fasta"))


def islandsOfDifference(seq_1, seq_2):
    """Returns number of regions of difference between two sequences."""
    bin_num_string = ''.join([str(int(seq_1[x]!=seq_2[x])) for x in range(len(seq_1))])
    return(len([x for x in bin_num_string.split('0') if x!='']))    # e.g. 11100001100000000100000000 has 3 diffs

def main():
    args = get_options()
    print(args)
    cwd = os.getcwd()
    json_file = str(args.json)
    #pangraph_dir = "/Users/Liam/Dropbox/_Projects/pangraph/beta-lactamase/scripts/"
    focal_block = args.focalblock #'OLNBEGPCML' #CTX-M-15
    panXdir = args.panXdir

    with open(json_file) as fh:
        G = json.load(fh)

    # The names are 
    contig_ids = [x['name'] for x in G['paths']]
    # the blocks are
    blocks = [{'id':x['id'], 'length':len(x['sequence']), 'depth':len(x['positions'])} for x in G['blocks']]

    # The blocks dictionary
    block_dict = {x['id']: {'gaps':x['gaps'], \
                        'mutate':x['mutate'],\
                        'insert':x['insert'],\
                        'delete':x['delete']} for x in G['blocks']}

    # We have panX output
    panX_gene_cluster_file = panXdir+"/vis/geneCluster.json"
    with open(panX_gene_cluster_file) as fh:
        panx = json.load(fh)

    # We need to convert between the alignment names and the block ids
    block_names = {x['ann']:x['msa'] for x in panx}


    # Read in all blocks
    aln_dict = {block: readBlockAln(block_names[block], panXdir+"/vis/geneCluster/") for block in block_names.keys()}


    # We can compute pairwise comparisons for all a, b
    N = len(contig_ids)
    total_combinations = N*(N-1)/2
    i = 0
    N_SUBSAMPLE = 5000
    comparisons = [(a, b) for a, b in combinations(range(N), 2)]
    subset_comparisons = sample(comparisons, N_SUBSAMPLE)
    #NAME = "NZ_CP023853.1"
    # random subset of all possible comparisons for speed
    with open('comparisons.csv', 'w') as f:
        # for name in ["NZ_CP027395.1", "NZ_CP013657.1", "NZ_CP014316.1" ,"NZ_CP013835.1" ,"NZ_CP025573.1", "NZ_CP025578.1", "NZ_CP041339.1", "NZ_CP018884.1"]:
        # a = contig_ids.index(name)
            #for b in range(N):
        for x in subset_comparisons:
            a, b = x[0], x[1]
            print(a, b)
            print(float(i)/N_SUBSAMPLE * 100)
            #print(i/N * 100)
            a_edits = {x: editPositions(x, contig_ids[a], block_dict) for x in block_dict.keys()}
            b_edits = {x: editPositions(x, contig_ids[b], block_dict) for x in block_dict.keys()}
            a_blocks = [x['id'] for x in G['paths'][a]['blocks']]
            b_blocks = [x['id'] for x in G['paths'][b]['blocks']]
            LCS_blocks = bf.lcs(a_blocks, b_blocks) 
            # this assumes the first block encountered, so might not handle duplications properly?
            a_lcs_block_seqs = [str(aln_dict[x][contig_ids[a]+'#1'].seq) for x in LCS_blocks]
            b_lcs_block_seqs = [str(aln_dict[x][contig_ids[b]+'#1'].seq) for x in LCS_blocks]
            lcs_block_total_length = sum([len(b) for b in a_lcs_block_seqs])
            lcs_block_diffs = sum([[a_lcs_block_seqs[y][x]!=b_lcs_block_seqs[y][x] for x in range(len(b_lcs_block_seqs[y]))].count(True) \
                    for y in range(len(a_lcs_block_seqs))])
            lcs_snp_density = float(lcs_block_diffs)/lcs_block_total_length
            # Flanking
            flanking_blocks = sharedBlocksBeforeBreakpoint(a_blocks, b_blocks, focal_block, upstream=True, includeStartingBlock=True)+\
                sharedBlocksBeforeBreakpoint(a_blocks, b_blocks, focal_block, upstream=False)
            # this assumes the first block encountered, so might not handle duplications properly?
            a_flanking_blocks_seqs = [str(aln_dict[x][contig_ids[a]+'#1'].seq) for x in flanking_blocks] 
            b_flanking_blocks_seqs = [str(aln_dict[x][contig_ids[b]+'#1'].seq) for x in flanking_blocks]
            flanking_L_total = sum([len(b) for b in a_flanking_blocks_seqs]) # Includes empty sites...
            flanking_N_diff_sites = sum([[a_flanking_blocks_seqs[y][x]!=b_flanking_blocks_seqs[y][x] for x in range(len(b_flanking_blocks_seqs[y]))].count(True) \
                    for y in range(len(b_flanking_blocks_seqs))])         
            flanking_seq_N_diff_regions = islandsOfDifference(''.join(a_flanking_blocks_seqs),\
                                                    ''.join(b_flanking_blocks_seqs))
            flanking_proportion_sites_different = float(flanking_N_diff_sites)/flanking_L_total
            flanking_seq_N_diff_region_density = float(flanking_seq_N_diff_regions)/flanking_L_total
                    # Also calculate the SNP density in the region up until the first breakpoint, rather than the whole thing
            # Sort the contig names so that the lexicographically first one comes first
            contig_first = sorted([contig_ids[a], contig_ids[b]])[0]
            contig_second = sorted([contig_ids[a], contig_ids[b]])[1]
            output_string = contig_first+','+contig_second+','+\
                str(len(LCS_blocks))+','+str(lcs_block_total_length)+','+str(lcs_block_diffs)+','+str(lcs_snp_density)+','+\
                    str(len(flanking_blocks))+','+str(flanking_L_total)+','+str(flanking_N_diff_sites)+','+str(flanking_proportion_sites_different)+','+str(flanking_seq_N_diff_regions)+','+str(flanking_seq_N_diff_region_density)
            print(output_string)
            _ = f.write(output_string+'\n')
            i += 1

if __name__=="__main__":
    main()
