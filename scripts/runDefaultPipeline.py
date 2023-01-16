# Runs the default pipeline on a multi-fasta file with a specified gene
import argparse
import subprocess
import os
import re
from datetime import datetime


def get_options():
    parser = argparse.ArgumentParser(description='Run pangraph pipeline on the flanking regions of a gene.',
                                     prog='runDefaultPipeline')
    parser.add_argument('--fasta', help='Multi-fasta file of input sequences containing the gene.', required=True)
    parser.add_argument('--gene', help='Fasta of focal gene to centre analysis on.', required=True)
    parser.add_argument('--outputdir', help='Output directory', required=False, default='./')
    parser.add_argument('--upstream', help='Upstream bases', required=False, default=5000)
    parser.add_argument('--downstream', help='Upstream bases', required=False, default=5000)
    parser.add_argument('--polish', help='Whether to use pangraph polish (time-intensive).', required=False, action='store_true', default=False)
    parser.add_argument('--panx', help='Whether to export panX output (time-intensive).', required=False, action='store_true', default=False)
    parser.add_argument('--bandage', help='Whether to run Bandage to get Bandage-visualized graph.', required=False, action='store_true', default=False)
    parser.add_argument('--prefix', help='Output prefix.', required=False, default=False)
    return parser.parse_args()

def run_command(command_list, output_file=False):
    if output_file==False:
        command_process = subprocess.Popen(command_list,
                        stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    else:
        command_process = subprocess.Popen(command_list,
                        stdout = output_file, stderr = subprocess.PIPE)
    command_out, _ = command_process.communicate() # Read the output from stdout
    #command_output = re.split('\n|\t',command_out.decode()) # Decode
    if command_out is not None:
        return(command_out.decode())
    else:
        return('')

def main():
    args = get_options()
    if not os.path.exists(args.outputdir):
        os.mkdir(args.outputdir)
    # Output prefixes
    if args.prefix==False:
        prefix_string = datetime.now().strftime("%Y_%m_%d.%H_%M_%S")
    else:
        prefix_string = str(args.prefix)
    output_prefix_no_dir = prefix_string+'.'+'all_u'+str(args.upstream)+'_d'+str(args.downstream)
    output_prefix = args.outputdir+'/'+output_prefix_no_dir

    # Log the parameters used and write to file
    with open(args.outputdir+'/'+'params.log', 'w') as f:
        f.write('run: '+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+'\n\n')
        for k, v in vars(args).items():
            f.write('%s: %s\n' % (str(k), str(v)))

    # Extract region around gene
    extract_region_command = ['python', 'extractRegion.py',
                                        '--gene', args.gene,
                                        '--input', args.fasta,
                                        '--upstream', str(args.upstream),
                                        '--downstream', str(args.downstream),
                                        '--complete',
                                        '--output', output_prefix]
    print(run_command(extract_region_command))

    # Align central gene
    mafft_command = ['mafft', '--quiet', output_prefix+'_focal_gene.fa']
    print(run_command(mafft_command, open(output_prefix+'_focal_gene.aln', 'w')))

    # Remove duplicates (how many variants of gene)
    rmdup_command = 'seqkit rmdup -s < '+output_prefix+'_focal_gene.aln > '+output_prefix+'_focal_gene.dedup.aln -D '+output_prefix+'_focal_gene.dedup.txt'
    rmdup_output = subprocess.call(rmdup_command, shell=True)
    print(rmdup_output)
    # SNP distances between gene variants (for NJ tree)
    snpdists_command = 'snp-dists -q -m '+output_prefix+'_focal_gene.aln > '+output_prefix+'_pangraph.gfa.gene.snps.tsv' # this is a bad filename - for plot-output-dists.R
    snpdists_output = subprocess.call(snpdists_command, shell=True)
    print(snpdists_output)

    # PANGRAPH
    # Build the initial pangraph
    if args.polish==True:
        pangraph_build_polish = 'pangraph build '+output_prefix+'.fa | pangraph polish > '+output_prefix+'_pangraph.json'
        pangraph_polish_output = subprocess.call(pangraph_build_polish, shell=True)
        print(pangraph_polish_output)
    else:
        pangraph_build = 'pangraph build '+output_prefix+'.fa > '+output_prefix+'_pangraph.json'
        pangraph_build_output = subprocess.call(pangraph_build, shell=True)
        print(pangraph_build_output)

    # Export to gfa
    pangraph_export = 'pangraph export --edge-minimum-length 0 '+output_prefix+'_pangraph.json -p '+output_prefix+'_pangraph -o ./'
    pangraph_export_output = subprocess.call(pangraph_export, shell=True)
    print(pangraph_export_output)

    if args.panx==True: # currently failing on GES-24 test case
        pangraph_panx_export = 'pangraph export '+output_prefix+'_pangraph.json -p '+output_prefix+'_pangraph -o '+args.outputdir+' --export-panX --no-export-gfa'
        pangraph_panx_export_output = subprocess.call(pangraph_panx_export, shell=True)
        print(pangraph_panx_export_output)

    # Prepare gfa for visualisation
    prepare_gfa = 'python preparePangraphGFA.py '+output_prefix+'_pangraph.gfa'
    prepare_gfa_output = subprocess.call(prepare_gfa, shell=True)
    print(prepare_gfa_output)

    # Find block in pancontigs which contains focal gene
    makedb = 'makeblastdb -in '+output_prefix+'_pangraph.fa -dbtype nucl'
    makedb_output = subprocess.call(makedb, shell=True)
    print(makedb_output)

    blast_gene = 'blastn -query '+args.gene+' -db '+output_prefix+'_pangraph.fa -outfmt 6 | cut -f 2'
    blast_gene_output = subprocess.check_output(blast_gene, shell=True)
    gene_block = blast_gene_output.decode().strip('\n')

    # Compute the distances
    compute_distance = 'python computeDistances.py '+output_prefix+'_pangraph.gfa '+gene_block
    compute_distance_output = subprocess.call(compute_distance, shell=True)
    print(compute_distance_output)

    # Plot the distances
    plot_dists = 'Rscript plot-output-dists.R '+output_prefix+'_pangraph.gfa.output_dists.csv '+output_prefix+'_pangraph.gfa.most_frequent_path_representative.txt '+output_prefix+'_focal_gene.dedup.txt'
    plot_dists_output = subprocess.call(plot_dists, shell=True)
    print(plot_dists_output)

    # If Bandage
    if args.bandage==True:
        #Make Bandage plot
        bandage = 'Bandage image '+output_prefix+'_pangraph.gfa.coloured.gfa '+\
                                output_prefix+'_pangraph.gfa.png '+\
                                '--height 4000 --width 7000 --colour custom'
        print(bandage)
        subprocess.call(bandage, shell=True)
        # Plot blocks
        plot_blocks = 'Rscript plot-blocks.R '+output_prefix+'_pangraph.gfa.blocks.csv '+\
                                            gene_block+' '+\
                                            output_prefix+'_pangraph.gfa.png '+\
                                            output_prefix+'_pangraph_blocks_plot.pdf'
        print(plot_blocks)
        subprocess.call(plot_blocks, shell=True)
    else:
        # Plot blocks without bandage
        plot_blocks = 'Rscript plot-blocks.R '+output_prefix+'_pangraph.gfa.blocks.csv '+\
                                            gene_block+' '+\
                                            output_prefix+'_pangraph.gfa.png '+\
                                            output_prefix+'_pangraph_blocks_plot.pdf'



if __name__ == "__main__":
    main()
