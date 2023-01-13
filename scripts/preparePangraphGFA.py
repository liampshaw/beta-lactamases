# Converts a gfa to a visualization of linear blocks
import re
from random import randint
import sys
from random import seed
import argparse


def get_options():
    parser = argparse.ArgumentParser(description='make a coloured version of pangraph gfa')
    parser.add_argument('input_file', help='Input gfa', type=str)
    parser.add_argument('--all', '-a', help='colour all blocks (even unique ones)',
                    action='store_true')
    return parser.parse_args()

def readGFA(gfa_file, colour_all=True, colour_seed=1789):
    """Reads in a gfa and returns useful format."""
    output_blocks=gfa_file+".blocks.csv"
    block_colour_file=gfa_file+".colours.csv"

    node_dict = {}
    genome_dict = {}
    block_count_dict = {}
    with open(gfa_file, 'r') as gfa:
        i = 1
        for line in gfa.readlines():
            entries = line.split()
            if entries[0] == "S":
                block_count_dict[entries[1]] = 0
                node_dict[entries[1]] = int(re.sub(".*:", "", entries[3]))
                i += 1
            if entries[0] == "P":
                blocks_in_path = entries[2].split(',')
                new_path = []
                i = 0
                for b in blocks_in_path:
                    block_count_dict[b[:-1]] += 1
                    new_path.append([b[:-1], b[-1], i, node_dict[b[:-1]]+i-1])
                    i += node_dict[b[:-1]]
                genome_dict[entries[1].strip('\n')] = new_path
    # Create random colours
    block_colour_dict = {}
    if colour_all==False:
        block_num = len([x for x in block_count_dict.values() if x>1])
    if colour_all==True:
        block_num = len([x for x in block_count_dict.values()])
    block_colors = []
    seed(colour_seed) # for random colours
    for i in range(block_num):
        block_colors.append("#%06X" % randint(0, 0xFFFFFF))
    i = 0
    for block, count in block_count_dict.items():
        if count==1 and colour_all==False:
            block_colour_dict[block] = "#BDBABA" # colour grey
        else:
            block_colour_dict[block] = block_colors[i]
            i += 1
    with open(output_blocks, 'w') as output_f:
        output_f.write("genome,block,strand,start,end,colour\n")
        for g, values in genome_dict.items():
            for v in values:
                block = v[0]
                output_f.write("%s\n" % (g+","+','.join([str(x) for x in v])+","+block_colour_dict[block]))
    with open(block_colour_file, 'w') as output_f_colours:
        output_f_colours.write("block,colour,length\n")
        for block, colour in block_colour_dict.items():
            output_f_colours.write("%s,%s,%s\n" % (block, colour, str(node_dict[block])))

def mostFrequentPathGFA(gfa_file): 
    """Finds the most frequent (full-length) path in the gfa file
    and returns a random ID of sequence that has that path.
    Also saves a file with the unique paths and their IDs just in case."""
    path_dict = {}
    with open(gfa_file, 'r') as gfa:
        i = 1
        for line in gfa.readlines():
            entries = line.split()
            if entries[0] == "P":
                path_dict[entries[1]] = entries[2] # store paths as strings
    # get unique paths (including direction/strand)
    unique_paths = list(set(path_dict.values()))
    # get how many times they occur
    unique_path_counts = [list(path_dict.values()).count(x) for x in unique_paths]
    # find most frequent path (can be tied - but .index will select just one)
    most_frequent_path = unique_paths[unique_path_counts.index(max(unique_path_counts))]
    scores = {k:v for k, v in zip(unique_paths, unique_path_counts)}    
    scores_sorted = {k: v for k, v in sorted(scores.items(), key=lambda item: item[1], reverse=True)}
    # Summarise these 
    output_path_file = gfa_file+'.unique_paths.txt'
    unique_path_ids = dict.fromkeys(scores_sorted.keys()) 
    # Pick a random sequence that has that path
    for seq_id, path in path_dict.items():
        if unique_path_ids[path]==None:
            unique_path_ids[path] = [seq_id]
        else:
            unique_path_ids[path].append(seq_id)
    with open(output_path_file, 'w') as output:
        for path, ids in unique_path_ids.items():
            output.write('%s\t%s\n' % (','.join(ids), path)) 
    most_frequent_path = list(scores_sorted.keys())[0]
    representative_most_frequent = unique_path_ids[most_frequent_path][0]
    with open(gfa_file+'.most_frequent_path_representative.txt', 'w') as f:
       f.write(representative_most_frequent)








def rewriteGFA(gfa_file, colour_file, new_gfa_file):
    """Rewrites a GFA with colours."""
    block_colour_dict = {}
    with open(colour_file, 'r') as colours:
        for line in colours.readlines():
            block_colour_dict[line.split(',')[0]] = line.split(',')[1].strip('\n')
    with open(new_gfa_file, "w") as new_gfa:
        with open(gfa_file, 'r') as gfa:
            for line in gfa.readlines():
                entries = line.split()
                if entries[0] == "S":
                    entries.append("CL:z:"+block_colour_dict[entries[1]])
                new_gfa.write("%s\n" % "\t".join(entries))


def main():
    args = get_options()
    input_gfa = str(args.input_file)
    colour_all_blocks = args.all
    readGFA(input_gfa, colour_all_blocks)
    # get representative of most frequent full-length path in GFA
    most_frequent_path_id = mostFrequentPathGFA(input_gfa)
    rewriteGFA(input_gfa,
        input_gfa+".colours.csv",
        input_gfa+".coloured.gfa")


if __name__ == "__main__":
    main()
