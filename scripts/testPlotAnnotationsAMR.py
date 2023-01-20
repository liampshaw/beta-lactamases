# Plotting linearized blocks with AMR annotations
# (from abricate)
# Offers a template from which to do other types of
# annotation e.g. ISFinder
# Currently uses abricate output
# i.e. assumes tsv with following headers
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	%COVERAGE	%IDENTITY	DATABASE	ACCESSION
import argparse
import subprocess



def get_options():
    parser = argparse.ArgumentParser(description='Plots blocks with annotations: currently AMR annotations.',
                                     prog='bandagePlot')
    parser.add_argument('--prefix', help='Prefix of files to use', required=True)
    parser.add_argument('--annotations', help='File of annotations (currently assumed generated with abricate)', required=True)
    parser.add_argument('--width', help='Width of final plot (inches)', required=False, default=6)
    parser.add_argument('--height', help='Height of final plot (inches)', required=False, default=10)
    parser.add_argument('--output', help='Output pdf', required=True)
    return parser.parse_args()

def main():
    args = get_options()
    blocks_file = str(args.prefix)+'_pangraph.json.blocks.csv'
    with open(args.prefix+'.gene_block.txt', 'r') as f:
        gene_block = f.readline().strip('')
    plot_command = 'Rscript plot-with-annotations.R '+blocks_file+' '+\
                        gene_block+' '+\
                        args.annotations+' '+\
                        args.output+' '+\
                        str(args.width)+' '+\
                        str(args.height)
    print("Command:", plot_command)
    plot_command_out = subprocess.call(plot_command, shell=True)
    print("Return code:", plot_command_out)


    #args= c('../output/GES-24/plotting/GES.all_u5000_d5000_pangraph.json.blocks.csv',
#        'ESELPCQPHB',
#        '../output/GES-24/plotting/GES-test-abricate-hits.tsv',
#        '../output/GES-24/plotting/GES-test-annotation-output.pdf',
#         10,
#         6)
if __name__=="__main__":
    main()
