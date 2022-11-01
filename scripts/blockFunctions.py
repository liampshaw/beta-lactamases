import sys
import operator
import pandas as pd

class CircularList(list):
    """Circular list class. 
    From https://stackoverflow.com/questions/8951020/pythonic-circular-list
    Allows indexing over the boundary of list, which is useful for circular genomes."""
    def __getitem__(self, x):
        if isinstance(x, slice):
            return [self[x] for x in self._rangeify(x)]
        index = operator.index(x)
        try:
            return super().__getitem__(index % len(self))
        except ZeroDivisionError:
            raise IndexError('list index out of range')
    def _rangeify(self, slice):
        start, stop, step = slice.start, slice.stop, slice.step
        if start is None:
            start = 0
        if stop is None:
            stop = len(self)
        if step is None:
            step = 1
        return range(start, stop, step)


def listToString(block_list):
    """Converts a list (signed permutation) to a string. Removes the sign.
    Args:
        block_list (list)
            Signed permutation. 
    Returns:
        String version, with integers mapped to character
    There are at least 10,000 characters (including emojis) so this is ok. 
    # To do: add in option to have sign treated differently too. 
    """
    return ''.join([chr(ord('@')+abs(x)) for x in block_list])

def stringToUnsignedPermutation(a_string):
    """Converts a string back into an unsigned permutation"""
    return list([ord(x)-ord('@') for x in a_string])


def lcs(s1, s2):
    """Returns the longest common subsequence (LCS) between two strings
    Implementation from: https://stackoverflow.com/questions/48651891/longest-common-subsequence-in-python
    Args:
        s1, s2 (str)
            Strings
    Returns:
        lcs (str)
            Longest common subsequence
    """
    matrix = [[[] for x in range(len(s2))] for x in range(len(s1))]
    for i in range(len(s1)):
        for j in range(len(s2)):
            if s1[i] == s2[j]:
                if i == 0 or j == 0:
                    matrix[i][j] = [s1[i]]
                else:
                    matrix[i][j] = matrix[i-1][j-1] + [s1[i]]
            else:
                matrix[i][j] = max(matrix[i-1][j], matrix[i][j-1], key=len)
    lcs = matrix[-1][-1]
    return lcs

def circularLCS(blocks_1, blocks_2):
    """Returns the longest common subsequence between two plasmids,
    taking into account all possible rotations but not taking into account
    the sign of the blocks.
    Args:
        blocks_1, blocks_2 (list)
            Two plasmids as signed permutations.
    Returns:
        max_lcs (str)
            String version of LCS (in character form)
    """
    circular_blocks_1 = CircularList(blocks_1) # make one circular
    max_lcs = ''
    for i in range(len(blocks_1)):
        circular_blocks_1_instance = circular_blocks_1[i:(i+len(blocks_1))]
        s1 = listToString(circular_blocks_1_instance)
        s2 = listToString(blocks_2)
        lcs_instance = lcs(s1, s2)
        if len(lcs_instance) > len(max_lcs):
            max_lcs = lcs_instance
    # TODO: convert this lcs back to a list of signed permutations
    return stringToUnsignedPermutation(max_lcs)

def lcsIntervals(blocks_1, blocks_2):
    """Returns a list of the intervals between the LCS elements in both sequences."""
    lcs = circularLCS(blocks_1, blocks_2)
    # Make a circular list
    # Go through until find the first element of LCS
    # Then keep going through the list, keeping track of where in the LCS you are too
    # And count each position
    # Choose the larger sequence
    lcs = CircularList(lcs)
    lcs_insertions = {(lcs[i], lcs[i+1]): [] for i in range(len(lcs)) }
    for blocks in [blocks_1, blocks_2]:
        blocks = [abs(x) for x in blocks]
        blocks = CircularList(blocks)
        block_start = blocks.index(lcs[0]) if lcs[0] in blocks else blocks.index(-lcs[0]) # Start from first element
        within_lcs_coordinate = False
        for i in range(block_start, block_start+len(max_seq)):
            b = abs(max_seq[i])
            lcs_element = lcs[j]
            if b==lcs_element:
                j += 1
            else:
                lcs_insertions[(lcs[j-1], lcs[j])].append(b)
    return lcs_insertions
   
    # #max_seq = blocks_1 if len(blocks_1)>len(blocks_2) else blocks_2
    # # Remove signs
    # max_seq = [abs(x) for x in max_seq]
    # max_seq = CircularList(max_seq)
    # block_start = max_seq.index(lcs[0]) if lcs[0] in max_seq else max_seq.index(-lcs[0]) # Start from first element
    # j = 0
    
        #print(lcs[j-1], lcs[j], b)

def breakpoints(blocks_1, blocks_2):
    insertion_dict = lcsIntervals(blocks_1, blocks_2)
    return len([x for x in insertion_dict.values() if x!=[]])

def distanceMatrixFromPlasmids(plasmid_dict, normalize=False, breaks=False):
    """Returns a distance matrix using the LCS distance from a dictionary of plasmids
    stored as signed permutations.
    Args:
        plasmid_dict (dict)
            plasmid_id : signed permutation
        normalize (Bool)
            Optional. Whether to normalize the distance by length of larger plasmid or not.
    Returns:
        dist_matrix (list)
            Matrix of all pairwise distances
    """
    plasmid_ids = list(plasmid_dict.keys())
    dist_matrix = [[0 for x in range(len(plasmid_ids))] for x in range(len(plasmid_ids))]
    for i in range(len(plasmid_ids)): # Todo: use iter.combinations instead since distance is symmetric
        for j in range(len(plasmid_ids)):
            plasmid_i = plasmid_dict[plasmid_ids[i]]
            plasmid_j = plasmid_dict[plasmid_ids[j]]
            if breaks==False:
                if normalize==True: # IMPORTANT: this definition is stupid - it's dataset-dependent. Need to marginalize shared blocks. Use length?
                    dist_matrix[i][j] = 1-len(circularLCS(plasmid_i, plasmid_j))/max(len(plasmid_i), len(plasmid_j))            
                else:
                    dist_matrix[i][j] = max(len(plasmid_i), len(plasmid_j))-len(circularLCS(plasmid_i, plasmid_j))
            elif breaks==True:
                dist_matrix[i][j] = breakpoints(plasmid_i, plasmid_j) # This is OK I think. But need to take into account second plasmid
    return dist_matrix

def plasmidSystemToSignedPermutation(pangraph_json_fname):
    """Takes a pangraph json and returns the plasmids as signed permutations.
    WARNING: assumes that there all plasmids share at least one block with the minimal backbone plasmid. 
    # Can easily be fixed to try/except
    Args:
        pangraph_json_fname (str)
            Pangraph file name
    Returns:
        plasmid_S_representations (dict)
            plasmid_id : signed_permutation
    """
    system_pangraph = put.Pangraph.load_json(pangraph_json_fname)
    system_blockstatsdf = system_pangraph.to_blockstats_df()
    # Find backbone plasmids
    backbone_ids = [x for x in list(system_pangraph.strains()) if "backbone" in x]
    # Find the smallest backbone - in terms of number of blocks? Or actual size? Needs to be actual size. Use consensus length rather than reading in fasta
    backbone_plasmid_lengths = [sum(system_blockstatsdf.loc[system_pangraph.paths[x].block_ids]["len"]) for x in backbone_ids]
    # Choose the smallest (if more than one, just choose the first - doesn't matter)
    minimal_backbone_plasmid = backbone_ids[backbone_plasmid_lengths.index(min(backbone_plasmid_lengths))]
    # Condition against this plasmid: prune plasmids that do not contain these backbone blocks
    backbone_blocks = CircularList(list(system_pangraph.paths[minimal_backbone_plasmid].block_ids)) # To check: are duplications included here?
    # Convert the backbone block IDs to integers (these can subsequently be signed)
    block_id_dict = {b:i+1 for i, b in enumerate(backbone_blocks)}
    # Add strands
    strands = {True: 1, False: -1} # dict to map strand bool to +/-
    backbone_strands = list(system_pangraph.paths[minimal_backbone_plasmid].block_strands)
    # Block ids
    backbone_blocks_strands = {} # define default strand based on first occurrence of block in plasmid
    for i in range(len(backbone_blocks)):
        if backbone_blocks[i] not in backbone_blocks_strands.keys():
            backbone_blocks_strands[backbone_blocks[i]] = strands[backbone_strands[i]] 
    block_id_dict = {b:(i+1) for i, b in enumerate(backbone_blocks_strands.keys())}
    # Block basic order
    # This horrible expression allows dupliation of the blocks with potentially different strand 
    backbone_plasmid_S = [backbone_blocks_strands[backbone_blocks[i]]*block_id_dict[backbone_blocks[i]]*strands[backbone_strands[i]] for i in range(len(backbone_blocks))]
    if -1 in backbone_plasmid_S:
        backbone_plasmid_S = [-x for x in list(reversed(backbone_plasmid_S))] # reverse and flip sign
    # Rotate so first block is first
    backbone_plasmid_S = CircularList(backbone_plasmid_S)
    backbone_plasmid_S = backbone_plasmid_S[backbone_plasmid_S.index(1):backbone_plasmid_S.index(1)+len(backbone_plasmid_S)]
    # Give other (non-backbone) blocks arbitrary + strand as default
    other_blocks =  list(set(list(system_blockstatsdf.index)) - set(list(block_id_dict.keys())))
    other_block_id_dict = {b:(i+1+len(block_id_dict)) for i, b in enumerate(other_blocks)}
    # Add them to the block ID dictionary
    block_id_dict.update(other_block_id_dict)
    # For each plasmid
    plasmid_ids = list(system_pangraph.strains())
    # Check for the conservation of core block orders and strands
    plasmid_S_representations = {k: [] for k in plasmid_ids}
    for plasmid in plasmid_ids:
        plasmid_blocks = list(system_pangraph.paths[plasmid].block_ids)
        plasmid_strands = list(system_pangraph.paths[plasmid].block_strands)
        plasmid_blocks_strands = {plasmid_blocks[i]:strands[plasmid_strands[i]] for i in range(len(plasmid_blocks))}
        new_plasmid_S = [plasmid_blocks_strands[p]*block_id_dict[p] for p in plasmid_blocks_strands.keys()]
        # If the 'anchor block' 1 is on a different strand to how we have defined in backbone, then reverse and flip strands 
        if -1 in new_plasmid_S:
            new_plasmid_S = [-x for x in list(reversed(new_plasmid_S))] # reverse and flip sign
        new_plasmid_S = CircularList(new_plasmid_S)
        # Rotate plasmid so that anchor block (1) is first and has same orientation as in backbone plasmids  
        plasmid_S_representations[plasmid] = new_plasmid_S[new_plasmid_S.index(1):new_plasmid_S.index(1)+len(new_plasmid_S)]
        print(plasmid_S_representations[plasmid])
    return [plasmid_S_representations, block_id_dict]

def writeDistanceMatrices(pangraph_json_fname, output_prefix="test"):
    plasmids, block_dict = plasmidSystemToSignedPermutation(pangraph_json_fname)
    plasmid_ids = list(plasmids.keys())
    dist_matrix_norm = distanceMatrixFromPlasmids(plasmids, normalize=True)
    dist_matrix = distanceMatrixFromPlasmids(plasmids)
    dist_matrix_breakpoints = distanceMatrixFromPlasmids(plasmids, breaks=True)
    # write to csv
    pd.DataFrame(dist_matrix).to_csv(output_prefix+"_dist_matrix.csv")
    pd.DataFrame(dist_matrix_norm).to_csv(output_prefix+"_dist_matrix_norm.csv")
    pd.DataFrame(dist_matrix_breakpoints).to_csv(output_prefix+"_dist_matrix_breakpoints.csv")
    # Write plasmid sequences as metadata
    pd.DataFrame(zip(range(len(plasmid_ids)), plasmid_ids, [[' '.join([str(x) for x in y])][0] for y in list(plasmids.values())]), columns=["id", "name", "sequence"]).to_csv(output_prefix+"_dist_metadata.csv")
    # Write block stats
    system_pangraph = put.Pangraph.load_json(pangraph_json_fname)
    system_blockstatsdf = system_pangraph.to_blockstats_df()
    block_lengths = [(block_dict[x], x, system_blockstatsdf.loc[x]["len"]) for x in list(block_dict.keys())]
    for k, v in block_dict.items():
        print(k, v)
    with open(output_prefix+"_block_lengths.txt", "w") as f:
        for b,id, l in block_lengths:
            f.write("%d,%s,%d\n" % (b, id,l))


def writeBlockDistances(pangraph_json_fname, output_prefix="test"):
    system_pangraph = put.Pangraph.load_json(pangraph_json_fname)
    system_blockstatsdf = system_pangraph.to_blockstats_df()
    block_lengths = system_blockstatsdf["lengths"]