#!/usr/bin/python3 -u

######################################################################################################################################################################
# Portions of this file were adapted from https://github.com/nloyfer/meth_atlas/blob/master/deconvolve.py and any reuse of this code should include that attribution #
######################################################################################################################################################################
import numpy as np
import pandas as pd
from scipy import optimize
import argparse
import os.path as op
import sys
from multiprocessing import Pool
import math
import matplotlib.pylab as plt
import matplotlib.cm
import matplotlib.colors

ATLAS_FILE = './reference_atlas.csv'
OUT_PATH = '.'

# Plotting parameters:
NR_CHRS_XTICKS = 30         # number of characters to be printed of the xticks
FIG_SIZE = (15, 7)          # figure size
COLOR_MAP = 'tab10'         # color map. See https://matplotlib.org/users/colormaps.html
#COLOR_MAP = 'Vega10'
# tissues with less than OTHERS_THRESH contribution will be clustered to 'other' (black):
OTHERS_THRESH = 0.01


####################################
#       Plotting methods           #
####################################

def hide_small_tissues(df):
    """
    tissues with very small contribution are grouped to the 'other' category.
    :return: The DataFrame with the new category ('other'),
             where the low-contribution tissues are set to 0.
    """
    others = df[ df < OTHERS_THRESH].sum()
    df[df < OTHERS_THRESH] = 0.0
    df = df.append(others.rename('other'))
    return df


def gen_bars_colors_hatches(nr_tissues):
    """
    Generate combinations of colors and hatches for the tissues bars
    Every tissue will get a tuple of (color, hatch)
    the last tuple is for the 'other' category, and is always black with no hatch.
    :return: a list of tuples, with length == nr_tissues
    """
    matplotlib.rcParams['hatch.linewidth'] = 0.3
    hatches = [None, 'xxx', '...', 'O', '++'][:nr_tissues // 7]

    nr_colors = int(math.ceil(nr_tissues / len(hatches)) + 1)

    # generate bars colors:
    cmap = matplotlib.cm.get_cmap(COLOR_MAP)
    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=float(nr_colors))
    colors = [cmap(norm(k)) for k in range(nr_colors)]

    def get_i_bar_tuple(i):
        color_ind = i % nr_colors
        hatch_ind = int(i // math.ceil(nr_tissues / len(hatches)))
        return colors[color_ind], hatches[hatch_ind]

    colors_hatches_list = [get_i_bar_tuple(i) for i in range(nr_tissues - 1)]
    return colors_hatches_list + [((0, 0, 0, 1), None)]


def plot_res(df, outpath, show=False):

    df = hide_small_tissues(df)
    nr_tissues, nr_samples = df.shape
      # generate bars colors and hatches:
    colors_hatches = gen_bars_colors_hatches(nr_tissues)

    plt.figure(figsize=FIG_SIZE)
    r = [i for i in range(nr_samples)]
    bottom = np.zeros(nr_samples)
    for i in range(nr_tissues):
        plt.bar(r, list(df.iloc[i, :]),
                edgecolor='white',
                width=0.85,
                label=df.index[i],
                bottom=bottom,
                color=colors_hatches[i][0],
                hatch=colors_hatches[i][1])
        bottom += np.array(df.iloc[i, :])

    # Custom x axis
    plt.xticks(r, [w[:NR_CHRS_XTICKS] for w in df.columns], rotation='vertical', fontsize=9)
    plt.xlabel("sample")
    plt.xlim(-.6, nr_samples - .4)

    # Add a legend and a title
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
    plt.title('Deconvolution Results\n' + op.basename(outpath))

    # adjust layout, save and show
    plt.tight_layout(rect=[0, 0, .83, 1])
    plt.savefig(outpath + '_deconv_plot_megacaryocytes.png')
    print("file was saved "+ outpath + '_deconv_plot_megacaryocytes.png')
    if show:
        plt.show()


# def load_sample(self, samp_path):
    # """
    # Read samples csv file. Reduce it to the atlas sites, and save data in self.samples
    # Note: samples file must contain a header line.
    # """

    # # validate path:
    # Deconvolve._validate_csv_file(samp_path)

    # samples = pd.read_csv(samp_path)
    # samples.rename(columns={list(samples)[0]: 'acc'}, inplace=True)
    # samples = samples.sort_values(by='acc').drop_duplicates(subset='acc').reset_index(drop=True)
    # samples = samples.merge(self.atlas['acc'].to_frame(), how='inner', on='acc')
    # return samples

def run(samp_path,samp_order):

    df = pd.read_csv(samp_path,index_col=0)
    
    #rename Erythrocyte progenitors to Megakaryocytes
    df.rename(index={'Erythrocyte progenitors':'Megakaryocytes'},inplace=True)
    
    #reindex list so lung will be the first row
    index_list=df.index.tolist()
    index_of_rows=[10,0,1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
    index_list=[index_list[i] for i in index_of_rows ]
    df=df.reindex(index_list)
    
    #reorder columns to be in wanted order, ay "sample_order"
    if samp_order is not None: #there is some wanted order
        df=df[samp_order.split(",")]
    
    #print(df)
    # Dump results
    
    #print(samp_path[:-len("_deconv_output.csv")])
   
    # Plot pie charts
    plot_res(df, samp_path[:-len("_deconv_output.csv")])


####################################
#            main                  #
####################################


def main():
    parser = argparse.ArgumentParser()
  
    parser.add_argument('samples_path',
                        help='Path to samples csv file. It must have a header line.\n'
                             'The first column must be Illumina IDs (e.g cg00000029)')

    #parser.add_argument('--out_dir', '-o', default=OUT_PATH, help='Output directory')
    parser.add_argument('--samp_order', '-s',default=None, help='sample order in dataframe,comma separated')

    args = parser.parse_args()
    run(args.samples_path, args.samp_order)


if __name__ == '__main__':
    main()
