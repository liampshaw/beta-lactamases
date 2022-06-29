# Make a heatmap from a mash matrix
# conda activate scipy

# Libraries
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib
import re
import sys
import numpy as np
from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as hc
import random
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import warnings
warnings.filterwarnings("ignore") # because seaborn gives warning as we are passing a distance matrix, even though we are doing this deliberately

DATA_DIR = '../data-processing'

def makePlot(family_name, dataset='ncbi_wgs', omit_zeros=False):
    """makes the plot of prevalence in a given dataset.
    omit_zeros: whether to plot those with zero prevalence"""
    if omit_zeros==True:
        snps = pd.read_csv(DATA_DIR+'/'+family_name+'_diffs_aln.tsv', sep='\t', index_col=0)
        #plt.close()
        # Read in prevalence data
        #prevalence_df = pd.read_csv(DATA_DIR+'/gene-prevalence-'+dataset+'.csv')
        presence_absence_df = pd.read_csv(DATA_DIR+'/CARD-all-hits-presence-absence-'+dataset+'.csv', index_col=0)
        family_df_pa = presence_absence_df[[x for x in presence_absence_df.columns if family_name in x]]
        family_df = pd.DataFrame(list(zip(list(family_df_pa.columns), family_df_pa.sum(axis=0))) ,
                                columns=['card_short_name', 'n'])
        # plot bar plot in ax
        #family_df = prevalence_df[prevalence_df['group']==family_name]
        # Those sequence variants that are present
        present_names = list(family_df[family_df['n']!=0]['card_short_name'])
        #if omit_zeros==True:
        df = snps.loc[present_names, present_names]
        if len(df)==0:
            return
        family_df = family_df[family_df['n']!=0]
        # and cluster
        linkage_mat = hc.linkage(df, method='average')
        # Get clusters
        maximum_distance = 100
        clusters = hc.fcluster(linkage_mat, maximum_distance, criterion='distance')
        # Assign random colours
        number_of_colours = len(set(clusters))
        colours = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
                     for i in range(number_of_colours)]
        colour_dict = dict(zip(set(clusters), colours))
        row_colours = [colour_dict[x] for x in clusters]
        cluster_fig = sns.clustermap(df,
                        row_linkage=linkage_mat,
                        col_linkage=linkage_mat,
                        vmin=0,
                        xticklabels=False,
                        yticklabels=True,
                        row_colors=row_colours,
                        cbar_kws=dict(use_gridspec=False,pad=.01,shrink=.25))
        cluster_fig.ax_heatmap.set_yticklabels(cluster_fig.ax_heatmap.get_ymajorticklabels(), fontsize = 1.5)
        #enlarge figure
        cluster_fig.fig.set_size_inches(8,6)
        # make some space to the right in the figure
        cluster_fig.gs.update(right=0.95)
        # divide existing axes
        divider = make_axes_locatable(cluster_fig.ax_heatmap)
        divider2 = make_axes_locatable(cluster_fig.ax_col_dendrogram)
        # create new axes for bar plot
        ax = divider.append_axes("right", size="20%", pad=0.3)
        # create empty space of same size as bar plot axes (don't use this space)
        nax = divider2.new_horizontal(size="20%", pad=0.3)
        # Sort the values for the bar plot to have the same order as clusters
        target = [t.get_text() for t in np.array(cluster_fig.ax_heatmap.get_yticklabels())]
        ind= np.array([list(df.index.values).index(t) for t in target])
        family_df_new = family_df
        family_df_new.index = family_df['card_short_name']
        # Plot
        ax.barh(np.arange(len(target)), [family_df_new['n'][x] for x in target])
        ax.set_yticklabels([])
        ax.set_ylim(-0.5,len(df.index)-.5)
        ax.set_xscale('log')
        ax.set_xticks([10, 100, 1000])
        ax.tick_params(axis='x', which='major', labelsize=6)
        ax.set_xlabel('Prevalence (n)', fontsize=8)
        ax.invert_yaxis()
        cluster_fig.fig.suptitle(family_name,fontsize=60, y=1)
        plt.title(dataset+'\n'+'n='+str(sum(family_df_new['n'])), fontsize=20)
        cluster_fig.ax_col_dendrogram.set_visible(False) #suppress column dendrogram
        plt.savefig(DATA_DIR+'/'+'omit_zeros_'+dataset+'_'+family_name+'_clusters.pdf', bbox_inches='tight')
    else:
        snps = pd.read_csv(DATA_DIR+'/'+family_name+'_diffs_aln.tsv', sep='\t', index_col=0)
        #plt.close()
        df = snps
        # and cluster
        linkage_mat = hc.linkage(df, method='average')
        # Get clusters
        maximum_distance = 100
        clusters = hc.fcluster(linkage_mat, maximum_distance, criterion='distance')
        # Assign random colours
        number_of_colours = len(set(clusters))
        colours = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
                     for i in range(number_of_colours)]
        colour_dict = dict(zip(set(clusters), colours))
        row_colours = [colour_dict[x] for x in clusters]
        cluster_fig = sns.clustermap(df,
                        row_linkage=linkage_mat,
                        col_linkage=linkage_mat,
                        vmin=0,
                        xticklabels=False,
                        yticklabels=True,
                        row_colors=row_colours,
                        cbar_kws=dict(use_gridspec=False,pad=.01,shrink=.25))
        cluster_fig.ax_heatmap.set_yticklabels(cluster_fig.ax_heatmap.get_ymajorticklabels(), fontsize = 1.5)
        #enlarge figure
        cluster_fig.fig.set_size_inches(8,6)
        # make some space to the right in the figure
        cluster_fig.gs.update(right=0.95)
        # divide existing axes
        divider = make_axes_locatable(cluster_fig.ax_heatmap)
        divider2 = make_axes_locatable(cluster_fig.ax_col_dendrogram)
        # create new axes for bar plot
        ax = divider.append_axes("right", size="20%", pad=0.3)
        # create empty space of same size as bar plot axes (don't use this space)
        nax = divider2.new_horizontal(size="20%", pad=0.3)
        # Sort the values for the bar plot to have the same order as clusters
        target = [t.get_text() for t in np.array(cluster_fig.ax_heatmap.get_yticklabels())]
        ind= np.array([list(df.index.values).index(t) for t in target])
        # Read in prevalence data
        #prevalence_df = pd.read_csv(DATA_DIR+'/gene-prevalence-'+dataset+'.csv')
        presence_absence_df = pd.read_csv(DATA_DIR+'/CARD-all-hits-presence-absence-'+dataset+'.csv', index_col=0)
        family_df_pa = presence_absence_df[[x for x in presence_absence_df.columns if family_name in x]]
        family_df = pd.DataFrame(list(zip(list(family_df_pa.columns), family_df_pa.sum(axis=0))) ,
                                columns=['card_short_name', 'n'])
        # plot bar plot in ax
        #family_df = prevalence_df[prevalence_df['group']==family_name]
        # Zero entries to add
        zero_hits = [x for x in target if x not in list(family_df['card_short_name'])]
        zero_hits_df = pd.DataFrame([[x, 0] for x in zero_hits],
        columns=family_df.columns)
        family_df_new = pd.concat([family_df, zero_hits_df])
        family_df_new.index = family_df_new['card_short_name']
        ax.barh(np.arange(len(target)), [family_df_new['n'][x] for x in target])
        ax.set_yticklabels([])
        ax.set_ylim(-0.5,len(df.index)-.5)
        ax.set_xscale('log')
        ax.set_xticks([10, 100, 1000])
        ax.tick_params(axis='x', which='major', labelsize=6)
        ax.set_xlabel('Prevalence (n)', fontsize=8)
        ax.invert_yaxis()
        cluster_fig.fig.suptitle(family_name,fontsize=60, y=1)
        plt.title(dataset+'\n'+'n='+str(sum(family_df['n'])), fontsize=20)
        cluster_fig.ax_col_dendrogram.set_visible(False) #suppress column dendrogram
        plt.savefig(DATA_DIR+'/'+dataset+'_'+family_name+'_clusters.pdf', bbox_inches='tight')
    return

def main():
    beta_lactamases = [x.strip() for x in open('../beta-lactamases.txt', 'r').readlines()]
    for b in beta_lactamases:
        for d in ['ncbi_plasmid', 'ncbi_chromosome']: # ignoring contig/wgs for now
            print(b, d)
            makePlot(b, dataset=d, omit_zeros=True)


if __name__ == "__main__":
    main()
