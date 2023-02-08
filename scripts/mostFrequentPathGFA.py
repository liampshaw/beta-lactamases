# Finds the most frequent (full-length) path in the gfa file
# and returns a random ID of sequence that has that path.
# sAlso saves a file with the unique paths and their IDs just in case.
import re
import sys
import argparse

def get_options():
        parser = argparse.ArgumentParser(description='Finds most frequent path in GFA')
        parser.add_argument('--gfa', help='Input gfa', type=str)
        return parser.parse_args()

def computePaths(gfa_file):
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

def main():
    """Finds the most frequent (full-length) path in the gfa file
    and returns a random ID of sequence that has that path.
    Also saves a file with the unique paths and their IDs just in case."""
    args = get_options()
    computePaths(args.gfa)
    

if __name__=="__main__":
    main()
