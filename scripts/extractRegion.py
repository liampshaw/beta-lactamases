# Extracts a region around gene of interest from a set of contigs.
# A useful precursor to further analysis.

from Bio import SeqIO
import os
import sys
import subprocess
import pandas as pd
import re
import logging
import numpy as np
import argparse
import glob


def get_options():
    parser = argparse.ArgumentParser(description='Extract region around a gene.',
                                     prog='extract_region')
    parser.add_argument('--gene', help='gene to search for', required=True)
    parser.add_argument('--input', help='input fasta file', required=True)
    parser.add_argument('--output', help='output prefix', required=True)
    parser.add_argument('--upstream', help='upstream bases', default=2500, required=False)
    parser.add_argument('--downstream', help='downstream bases', default=2500, required=False)
    parser.add_argument('--complete', help='only keep contigs with all requested region', action='store_true')
    parser.add_argument('--circular', help='whether input contigs are circular', action='store_true')
    parser.add_argument('--smh', help='extract contigs with multiple hits into separate file', type=str, required=False, default='')
    return parser.parse_args()

def reverse_complement(seq):
    '''Returns the reverse complement of a DNA sequence. Assumes no insertions.
    Args:
        seq (str)
            String of bases (e.g. ATCCG)
    Returns:
        bases (str)
            Reverse complemented string (e.g. CGGAT)
            (N.B. if seq contains non-standard characters, returns None)
    '''
    # Mapping of bases to complement (note that N->N)
    complement_map = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    bases = list(seq.upper()) # Make sure upper-case
    # Check if non-standard bases are present, and error if so
    if any([x not in list(complement_map.keys()) for x in bases]):
        non_standard_bases = [x for x in bases if x not in list(complement_map.keys())]
        print('Refusing to reverse complement!')
        print('Your sequence contains non-standard characters:', ''.join(set(non_standard_bases)))
        return
    else:
        bases = reversed([complement_map.get(base,base) for base in bases])
        bases = ''.join(bases)
        return(bases)

def read_fasta(fasta_file):
    '''Reads in fasta file as dict using SeqIO.
    Args:
        fasta_file (str)
            Filename of fasta
    Returns:
        fasta_dict (dict)
            Dictionary of sequences in fasta
    '''
    path_to_file = os.path.abspath(fasta_file) # get absolute path
    fasta_parsed = SeqIO.parse(path_to_file, 'fasta') # read in fasta
    fasta_dict = SeqIO.to_dict(fasta_parsed, lambda rec:rec.id) # convert to dict
    return(fasta_dict)



# Makes more sense to do a single blast once for the whole input file
# Then only keep those contigs that contain one hit only (for now - can have more sophisticated handling of multi-hit later)
# And simply take all of the sequence upstream/downstream as required
def blast_search(query_fasta, db_fasta):
    '''runs a blast search'''

    logging.debug('\nMaking blast database...')
    subprocess.check_call(['makeblastdb', '-in', db_fasta, '-dbtype', 'nucl'],\
        stderr=subprocess.DEVNULL,\
        stdout=open(os.devnull, 'w'))
    logging.debug('\nSearching for gene...')
    blast_process = subprocess.Popen(['blastn', '-db', db_fasta, \
                            '-query', query_fasta, \
                            '-outfmt', '6', \
                            '-max_target_seqs', '10000000'],
                            stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    blast_out, _ = blast_process.communicate() # Read the output from stdout
    blast_output = re.split('\n|\t',blast_out.decode()) # Decode
    #print(blast_output)
    logging.debug('\nRemoving temporary blast databases...')
    [os.remove(x) for x in glob.glob(db_fasta+'.n*')]
    # Looking through to cut out upstream region
    if blast_output == ['']:
        logging.debug('\nNo blast hit for gene!')
        return
    else:
        blast_output = blast_output[:-1]
        blast_results = np.reshape(blast_output, (-1, 12))
        contigs_with_hits = list(blast_results[:,1])
        contigs_one_hit = [x for x in set(contigs_with_hits) if contigs_with_hits.count(x)==1]
        contigs_multiple_hits = [x for x in set(contigs_with_hits) if contigs_with_hits.count(x)>1]
        results_dict_one_hit = {result[1]: [int(result[8]), int(result[9])] for result in blast_results if result[1] in contigs_one_hit} # only keep contigs with exactly one hit
        results_dict_multiple_hit = {x: 'multiple_hits' for x in contigs_multiple_hits}
        return({**results_dict_one_hit, **results_dict_multiple_hit})

def extract_regions(sequences, blast_hits, gene_length, is_circular=False, upstream_bases=1000, downstream_bases=1000):
    '''extracts the regions from a fasta file'''
    # Check for sequences with >1 hit in blast hits
    extracted_seqs = {}

    for seq_id, coords in blast_hits.items():
        if coords=="multiple_hits":
            print('WARNING: Contig '+seq_id+' has multiple blast hits. Not including in output.')
        else:
            contig_seq = str(sequences[seq_id].seq)
            # Check the length - if shorter, use whole contig sequence
            if len(contig_seq)<(upstream_bases+downstream_bases+gene_length):
                print('WARNING: Contig '+seq_id+' is shorter than flanking region requested!')
                if coords[0]<coords[1]:
                    extracted_seqs[seq_id] = contig_seq
                elif coords[0]>coords[1]:
                    extracted_seqs[seq_id] = reverse_complement(contig_seq)
            #print(region_seq)
            else:
                #print(seq_id)
                triplicate_seq = contig_seq+contig_seq+contig_seq
                if coords[0]<coords[1]: # positive strand
                    if not is_circular:
                        limits = [max(0, coords[0]-1-upstream_bases), min(coords[1]+downstream_bases, len(contig_seq))]
                        extracted_seqs[seq_id] = contig_seq[limits[0]:limits[1]]
                        #print(limits, limits[1]-limits[0])
                    elif is_circular:
                        limits = [max(coords[1], coords[0]-1-upstream_bases+len(contig_seq)), min(coords[0]+2*len(contig_seq), coords[1]+downstream_bases+len(contig_seq))]
                        extracted_seqs[seq_id] = triplicate_seq[limits[0]:limits[1]]
                        #print(limits, limits[1]-limits[0])
                elif coords[0]>coords[1]: # negative strand
                    if not is_circular: # also put in a max thing here
                        limits = [max(0, coords[1]-1-downstream_bases), min(coords[0]+upstream_bases, len(contig_seq))]
                        extracted_seqs[seq_id] = reverse_complement(contig_seq[limits[0]:limits[1]])
                        #print(limits, limits[1]-limits[0])
                    elif is_circular:
                        limits = [max(coords[0], coords[1]-1-downstream_bases+len(contig_seq)), min(coords[1]+2*len(contig_seq), coords[0]+upstream_bases+len(contig_seq))]
                        extracted_seqs[seq_id] = reverse_complement(triplicate_seq[limits[0]:limits[1]])
                        #print(limits, limits[1]-limits[0])
    return extracted_seqs

def store_multiple_hits(sequences, blast_hits, smh_file):
    """stores multiple hits to file"""
    N_multiple_hits = 0
    with open(smh_file, 'w') as output_f:
        for seq_id, coords in blast_hits.items():
            if coords=="multiple_hits":
                N_multiple_hits += 1
                output_f.write('>%s\n%s\n' % (seq_id, str(sequences[seq_id].seq)))
    return(N_multiple_hits)


def main():
    args = get_options()
    input_fasta = args.input
    input_sequences = read_fasta(input_fasta)

    output_fasta = args.output+'.fa'
    output_gene_fasta = args.output+'_focal_gene.fa' # output for just gene sequences
    gene_fasta = args.gene
    # Get length of gene
    gene_seq = read_fasta(gene_fasta)
    if len(gene_seq)>1:
        print('The fasta containing your gene of interest has more than sequence! Check the fasta.')
        exit()
    else:
        gene_length = [len(str(gene_seq[v].seq)) for v in gene_seq][0]
    #print('gene length is: ', str(gene_length))

    results = blast_search(gene_fasta, input_fasta)
    #print(results)
    extractions = extract_regions(input_sequences, results,
                                gene_length=gene_length,
                                is_circular=args.circular,
                                upstream_bases=int(args.upstream),
                                downstream_bases=int(args.downstream))
    N_seqs_written = 0
    seqs_written = []
    with open(output_fasta, 'w') as output_file:
        for k, v in extractions.items():
            #print(k+',', len(v)-gene_length, 'bases extracted around gene', '('+str(len(v))+' total)')
            if args.complete==False:
                output_file.write('>%s %sbp\n%s\n' % (k, str(len(v)), v))
                seqs_written.append(k)
                N_seqs_written += 1
                #print("...writing to file.")
            elif v is not None:
                if len(v)>(int(args.upstream)+gene_length+int(args.downstream)-1):
                    if "n" not in v and "N" not in v: # check for ambiguous characters
                        #print("...writing to file.")
                        output_file.write('>%s %sbp\n%s\n' % (k, str(len(v)), v))
                        seqs_written.append(k)
                        N_seqs_written += 1
    # Write multiple hits if requested
    if args.smh!='':
        N_multiple_hits = store_multiple_hits(input_sequences, results, args.smh)
        print("... Wrote "+str(N_multiple_hits)+" contigs (whole) with multiple blast hits to: "+args.smh)

    # Write just the gene sequences
    print(results)
    gene_extractions = extract_regions(input_sequences, results, gene_length=gene_length, is_circular=False, upstream_bases=0, downstream_bases=0)
    with open(output_gene_fasta, 'w') as output_file:
        for k, v in gene_extractions.items():
            if k in seqs_written: # only write those that have been written into regions
                output_file.write('>%s\n%s\n' % (k, v))

    print("\nSUMMARY:\n"+str(len(input_sequences))+" contigs in input fasta\n"+str(N_seqs_written)+" regions extracted (from contigs with one blast hit for gene)")

if __name__ == "__main__":
    main()
