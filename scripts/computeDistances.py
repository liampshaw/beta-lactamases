
#
import itertools as iter
import argparse

def get_options():
    parser = argparse.ArgumentParser(description='computes a variety of distances between sequences')
    parser.add_argument('input_file', help='Input file (pangraph json)', type=str)
    parser.add_argument('gene_block', help='Block containing gene of interest (anchor)', type=str)
    return parser.parse_args()

def distFirstBreakpoint(a, b, starting_block, upstream=True):
    """a, b are lists of blocks"""
    a_start = a.index(starting_block)
    b_start = b.index(starting_block)
    shared_blocks = []
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


def jaccard(a, b):
    """a, b are lists of blocks"""
    # Needs writing
    return

# TO ADD
# Function to call blast to get the starting block (contains gene) from the pangraph

def main():
    #block_lengths = {}
    args = get_options()
    #with open(args.input_file+'.colours.csv', 'r') as f:
    #    for i, line in enumerate(f.readlines()):
    #        if i>0:
    #            line = line.split(',')
    #            block_lengths[line[0]] = int(line[2].strip())
    #block_lengths = {}
    block_lengths = {}
    path_dict = {}
    with open(args.input_file+'.blocks.csv', 'r') as f:
        for i, line in enumerate(f.readlines()):
            if i>0:
                line = line.split(',')
                if line[0] in path_dict.keys():
                    path_dict[line[0]].append(line[1])
                else:
                    path_dict[line[0]] = [line[1]]
                if line[1] in block_lengths.keys():
                    block_lengths[line[1]][line[0]] = abs(int(line[4])-int(line[3]))
                else:
                    block_lengths[line[1]] = {line[0]: abs(int(line[4])-int(line[3]))}
    # Generate snp-dists?
    starting_block = args.gene_block

    with open(args.input_file+'.output_dists.csv', 'w') as output_file:
        output_file.write('seq1,seq2,dist.up,dist.down,snps\n')
        with open(args.input_file+'.gene.snps.tsv', 'r') as f: # generated with snp-dists -p from gene seqs
            for line in f.readlines():
                line = line.strip().split()
                a, b, snps = line[0], line[1], int(line[2])
                upstream_blocks = distFirstBreakpoint(path_dict[a], path_dict[b], starting_block)
                upstream_dist = sum([block_lengths[x][a] for x in upstream_blocks])
                downstream_blocks = distFirstBreakpoint(path_dict[a], path_dict[b], starting_block, upstream=False)
                downstream_dist = sum([block_lengths[x][a] for x in downstream_blocks])
                output_file.write('%s,%s,%d,%d,%d\n' % (a, b, upstream_dist, downstream_dist, snps))

if __name__ == "__main__":
    main()
