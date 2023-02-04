# Convenience script to plot Bandage and final plot
# if all files already exist in folder (handy when wanting
# to make final versions of plots
import subprocess
import argparse

def get_options():
    parser = argparse.ArgumentParser(description='Plots bandage and block plot for given input files.',
                                     prog='bandagePlot')
    parser.add_argument('--prefix', help='Prefix of files to use', required=True)
    parser.add_argument('--nobandage', help='whether to not do bandage', required=False, action='store_true', default=False)
    parser.add_argument('--bandagewidth', help='Width of Bandage image (px)', required=False, default=7000)
    parser.add_argument('--bandageheight', help='Height of Bandage image (px)', required=False, default=4000)
    parser.add_argument('--width', help='Width of final plot (inches)', required=False, default=7000)
    parser.add_argument('--height', help='Height of final plot (inches)', required=False, default=7000)
    parser.add_argument('--output', help='Output prefix', required=False, default='')
    parser.add_argument('--metadata', help='whether to use metadata', required=False, action='store_true', default=False)
    return parser.parse_args()


def main():
    args = get_options()
    prefix = args.prefix
    if args.output=='':
        output_prefix = prefix
    else:
        output_prefix = args.output
    bandage = 'Bandage image '+prefix+'_pangraph.gfa.coloured.gfa '+\
                            output_prefix+'_pangraph.gfa.png '+\
                            '--height '+str(args.bandageheight)+' --width '+str(args.bandagewidth)+' --colour custom'
    print("Command:", bandage)
    bandage_output = subprocess.call(bandage, shell=True)

    print("Return code:", bandage_output)
    # Plot blocks
    bandage_string = output_prefix+'_pangraph.gfa.png '
    if args.nobandage==True:
        bandage_string = 'none'
        
    print("## LINEAR PLOT ##")
    with open(prefix+'.gene_block.txt', 'r') as f:
        gene_block = f.readline().strip('')
    if args.metadata==False:
        plot_blocks = 'Rscript plot-blocks.R '+prefix+'_pangraph.json.blocks.csv '+\
                                        gene_block+' '+\
                                        bandage_string+' '+\
                                        output_prefix+'_pangraph_blocks_plot.pdf '+\
                                        str(args.width)+' '+str(args.height)
    elif args.metadata==True:
        plot_blocks = 'Rscript plot-blocks-chromosome-plasmid.R '+prefix+'_pangraph.json.blocks.csv '+\
                                        gene_block+' '+\
                                        bandage_string+' '+\
                                        output_prefix+'_pangraph_blocks_plot.pdf '+\
                                        str(args.width)+' '+str(args.height)
    print("Command:", plot_blocks)
    plot_blocks_out = subprocess.call(plot_blocks, shell=True)
    print("Return code:", plot_blocks_out)

if __name__=="__main__":
    main()
