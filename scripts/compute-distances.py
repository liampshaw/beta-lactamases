
#
block_lengths = {}
with open('pangraph_all_u5k_d5k.gfa.colours.csv', 'r') as f:
    for i, line in enumerate(f.readlines()):
        if i>0:
            line = line.split(',')
            block_lengths[line[0]] = int(line[2].strip())

path_dict = {}
with open('pangraph_all_u5k_d5k.gfa.blocks.csv', 'r') as f:
    for i, line in enumerate(f.readlines()):
        if i>0:
            line = line.split(',')
            if line[0] in path_dict.keys():
                path_dict[line[0]].append(line[1])
            else:
                path_dict[line[0]] = [line[1]]

import itertools as iter

#for a, b in iter.combinations(path_dict.keys(), 2):
#    shared_blocks = set(path_dict[a]).intersection(set(path_dict[b]))
#    print(a,b, 1-sum([block_lengths[x] for x in shared_blocks])/sum([block_lengths[x] for x in path_dict[a]]))

# Distance until first breakpoint of pairs
gene_block = "ZZIMKFDWIU"
def distFirstBreakpoint(a, b, start="ZZIMKFDWIU", upstream=True):
    """a, b are lists of blocks"""
    a_start = a.index(start)
    b_start = b.index(start)
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
            print(i)
            if a[a_start-i]==b[b_start-i]:
                shared_blocks.append(b[b_start-i])
            else:
                #if a[a_start+i] in shared_blocks:
                #    shared_blocks.append(a[a_start+i])
                break
    return(shared_blocks)




def jaccard(a, b):
    """a, b are lists of blocks"""
    return (1-float(len(set(a).intersection(set(b)))/len(set(a).union(set(b)))))

with open('output_dists.csv', 'w') as output_file:
    with open('all_u5k_d5k.CTX-M.snpdists.tsv', 'r') as f: # generated with snp-dists -p from gene seqs
        for line in f.readlines():
            line = line.strip().split()
            a, b, snps = line[0], line[1], int(line[2])
            upstream_blocks = distFirstBreakpoint(path_dict[a], path_dict[b])
            upstream_dist = sum([block_lengths[x] for x in upstream_blocks])
            downstream_blocks = distFirstBreakpoint(path_dict[a], path_dict[b], upstream=False)
            downstream_dist = sum([block_lengths[x] for x in downstream_blocks])
            output_file.write('%s,%s,%d,%d,%d\n' % (a, b, upstream_dist, downstream_dist, snps))
