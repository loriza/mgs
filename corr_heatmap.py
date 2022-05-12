import os
import sys
import warnings
import numpy as np
import palettable
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import getpass
import argparse


class DESCRIPTION:
    Program = "Program: %(prog)s (Draw corplots via Seaborn)"
    Version = "Version: 2022.05.07-15.40"
    Description = "%s\n%s\n%s" % (Program, Version, Contact)


USER = getpass.getuser()
MPLCONFIGDIR = f"/work/workspace/{USER}/matplotlib"
os.environ['MPLCONFIGDIR'] = MPLCONFIGDIR
mpl.use('Agg')

mpl.rcParams['font.family'] = 'Microsoft YaHei'
mpl.rcParams['font.sans-serif'] = ['Microsoft YaHei']

#mpl.rc("font",family='Microsoft YaHei')

# 设置文件夹的路径
curpath = os.path.dirname(os.path.realpath(__file__))
path1 = os.path.join(curpath, 'heatmap_plot.py')
# 改变路径和打印路径的信息
warnings.filterwarnings("ignore")

delimiter = '\t'


def df_reader(path):
    try:
        return pd.read_table(path, sep=delimiter, comment='#', index_col=0, header=0)
    except Exception:
        sys.exit(f"{path} does not exist")


def sub_df(df, min_r2, df2, filter):
    if filter is True:
        #去除所有r都小于cutoff的行/列
        drop_row = df[(df.abs() < min_r2).all(axis=1)].index.values.tolist()
        drop_col = df[df.columns[(df.abs() < min_r2).all(axis=0)]].columns.values.tolist()

        df = df.drop(drop_row, inplace=False).drop(drop_col, inplace=False, axis=1)
        df2 = df2.drop(drop_row, inplace=False).drop(drop_col, inplace=False, axis=1)

    #去除含有na的行/列
    new_df = df.dropna(axis=1, how='all', inplace=False).dropna(axis=0, how='all', inplace=False)
    new_df2 = df2.dropna(axis=1, how='all', inplace=False).dropna(axis=0, how='all', inplace=False)
    
    return new_df, new_df2

def flatten(pval):
    mydata = pval.copy()

    mask1 = np.array(mydata)
    mask2 = mask1.flatten()
    mask3 = [str(x) for x in mask2]
    for i in range(len(mask3)):
        if mask3[i] > str(0.05):
            mask3[i] = " "
        elif (mask3[i] <= str(0.05)) & ((mask3[i]) > str(0.01)):
            mask3[i] = "*"
        elif mask3[i] <= str(0.01):
            mask3[i] = "**"
        elif mask3[i] <= str(0.001):
            mask3[i] = "***"
    mask4 = mask3
    # print(mask4)
    return mask4

def sumn(data, pvalue):
    mask = flatten(pvalue)
    print(len(mask))
    if (data.shape[0]+data.shape[1]) < 100:
        label_size = 12
        tick_size = 12
        dpi = 100
        annot = np.asarray(list(mask)).reshape(data.shape[0], data.shape[1])

        lw = 0.01
    elif ((data.shape[0]+data.shape[1]) > 100) & ((data.shape[0]+data.shape[1]) < 400):
        label_size = 6
        tick_size = 6
        dpi = 200
        annot = np.asarray(list(mask)).reshape(data.shape[0], data.shape[1])

        lw = 0.01
    elif ((data.shape[0] + data.shape[1]) > 400) & ((data.shape[0] + data.shape[1]) < 800):
        label_size = 0.1
        tick_size = 0.1
        dpi = 300
        annot = np.asarray(list(mask)).reshape(data.shape[0], data.shape[1])

        lw = 0.01
    elif (data.shape[0] + data.shape[1]) > 800:
        label_size = 0.1
        tick_size = 0.1
        dpi = 500
        annot = False
        lw = 0.001

    annot_font = {'size': label_size}
    ticks_font = {'fontsize': tick_size}

    return label_size, tick_size, dpi, annot, lw

def legend_title(method):
    if method == 'pearson':
        method = "Pearson's r"
    elif method == 'spearman':
        method = "Spearman's ρ"
    elif method == 'kendall':
        method = "Kendall's τ"
    return method
def corr_heatmap(corr, prefix, height, width, fmt, methods):

    fig, ax = plt.subplots()
    fig = plt.figure(figsize=(width, height), dpi=dpi)

    sns.set(font_scale=1, style="dark")
    sns.set_context("paper")

    ax = sns.heatmap(corr,
                     vmax=1,
                     vmin=-1,
                     # mask=(corr.abs()<0.25),
                     yticklabels=1,  # plot all labels
                     xticklabels=1,
                     
                     square=True,
                     linewidths=lw,
                     cbar=True,
                     cbar_kws={"shrink": 0.2,
                               # "orientation":'horizontal',
                               "label": methods,
                               'ticks': np.linspace(-1, 1, 11, endpoint=True, dtype=float)},
                     cmap=palettable.cmocean.diverging.Curl_10.mpl_colors,
                     annot=annot,
                     fmt='',
                     annot_kws={'size': lab_size,
                                'weight': 'normal',
                                'color': '#253D24'})

    """setting colorbar"""
    cax = plt.gcf().axes[-1]
    cax.tick_params(labelsize=11)

    """setting colorbar title"""
    cbar = ax.collections[0].colorbar
    font1 = {'size': 20, 'color': '#000000'}
    cbar.set_label(label=methods,
                   fontdict=font1,
                   rotation='vertical',
                   va='top',
                   ha='center')

    """setting axis label"""
    plt.yticks(fontsize=ticks_size)
    plt.xticks(rotation=45, fontsize=ticks_size, ha='right')
    plt.tight_layout()
    outfile = prefix + '.corr.heatmap.' + fmt
    plt.savefig(outfile, bbox_inches='tight')
    print(f"Heatmap: {outfile}")
    plt.close()


def cluster_heatmap(corr, prefix, width, height, fmt, methods, clusterlist):

    fig, ax2 = plt.subplots()

    sns.set(font_scale=1, style="dark")
    sns.set_context("paper")
    ax2 = sns.clustermap(data=corr,
                         figsize=(width, height),
                         vmin=-1,
                         vmax=1,
                         row_cluster=clusterlist[0],
                         col_cluster=clusterlist[1],
                         #row_colors=row_colors,
                         yticklabels=1,  # plot all labels
                         xticklabels=1,
                         dendrogram_ratio=(.2, .05),
                         annot=annot,
                         annot_kws={'size': lab_size,
                                    'weight': 'normal',
                                    'color': '#253D24'},
                         fmt='',
                         linewidths=lw,
                         cbar=True,
                         cbar_kws={"shrink": 0.7,
                                   'aspect': 20*0.7,
                                   'ticks': np.linspace(-1, 1, 11,
                                                        endpoint=True,
                                                        dtype=float)},
                         cmap=palettable.cmocean.diverging.Curl_10.mpl_colors
                         )
    """setting colorbar"""
    ax2.fig.subplots_adjust(right=0.8)
    ax2.ax_cbar.set_position((1,.8, .02, .1))

    ax2.ax_cbar.set_title(methods, fontsize=ticks_size*2,ha='center')
    
    """setting tick labels"""
    plt.setp(ax2.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=ticks_size)
    plt.setp(ax2.ax_heatmap.get_yticklabels(), fontsize=ticks_size)

    plt.tight_layout()
    outfile = prefix + '.cluster.heatmap.' + fmt
    plt.savefig(outfile, bbox_inches='tight',dpi=dpi)
    plt.close()
    print(f"Cluster heatmap: {outfile}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=DESCRIPTION.Description)
    parser.add_argument('-c', dest="c", required=True,
                        help='input correlation data file')
    parser.add_argument('-p', dest="p",
                        help='input pvalues data file')
    parser.add_argument("-o", dest="output",
                        help="directory and prefix of out images")
    parser.add_argument('-f', dest="format", default='pdf', nargs='*',
                        choices=['png', 'pdf', 'svg', 'jpg'],
                        help='output heatmap image format: default jpg')
    parser.add_argument('-t', dest="methods", default='spearman',
                        choices=['spearman', 'pearson', 'kendall'],
                        help='methods for correlation: default pearson')
    parser.add_argument("-y", "--list1", dest="list1", action='append', type=int,
                        nargs=2, choices=[0, 1], default=[1, 1],
                        help="If 1 is provided, cluster the row or columns")
    parser.add_argument("-r",  dest="rc", type=float, default=0.3,
                        help="minimum coefficient: default 0.3")
    parser.add_argument('-e', dest="height", default=25, type=float,
                        help='image width, default 25 inches')
    parser.add_argument('-w', dest="width", default=22, type=float,
                        help='image height, default 22 inches')
    parser.add_argument("-m", "--verbose", action="store_false",
                        help=" if given, filtering based on the r cutoff will not be performed")

    args = parser.parse_args()

    if os.path.exists(os.path.dirname(args.output)) is False:
        os.makedirs(os.path.dirname(args.output))


    coef = df_reader(args.c)
    pvalues = df_reader(args.p)

    coef, pvalues = sub_df(coef, args.rc, pvalues, args.verbose)
    method = legend_title(args.methods)
    lab_size, ticks_size, dpi, annot, lw = sumn(coef, pvalues)
    
    cluster_heatmap(coef, args.output, args.width, args.height, args.format, method, args.list1)
    corr_heatmap(coef, args.output, args.height, args.width, args.format, method)
