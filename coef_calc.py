# -*- coding: UTF-8 -*-
import argparse
import os
import sys
import warnings
import numpy as np
import pandas as pd
import csv

sys.path.append(r'./heatmap_plot.py')
curpath = os.path.dirname(os.path.realpath(__file__))
print("Base Directory : ", curpath)
warnings.filterwarnings("ignore")


class DESCRIPTION:
    Program = "Program: %(prog)s (Correlation Calculation)"
    Version = "Version: 2022.05.07-18.20"
    Contact = "Contact: Liao Yanwen, liaoyw_09@hotmail.com"
    Description = "%s\n%s\n%s" % (Program, Version, Contact)


# -------------------------------------------------------------------------------
# global constants
# -------------------------------------------------------------------------------
delimiter = '\t'
o_format = '.tsv'

# -------------------------------------------------------------------------------
# utilities
# -------------------------------------------------------------------------------


def df_reader(path):
    if path:
        return pd.read_table(path, sep=delimiter, comment='#', index_col=0, header=0)


def validate_metadata_format(metadata):
    """check if an object exists in the metadata"""
    if metadata is not None:
        objc = []
        for i in metadata:
            if pd.api.types.is_object_dtype(metadata[i]):
                j = metadata[i].name
                objc.append(j)
        if len(objc):
            sys.exit('Unordered categorical variables/unsupported variables：', objc)


def transpose_feature(feature, group):
    ft = np.isin(group._stat_axis.values.tolist(),
                 feature.columns.values.tolist(),
                 invert=False)
    for i in ft:
        if i == True:
            feature = feature.T
            break
    return feature


def corr_calc(matrix, method):
    cols = matrix.shape[1]
    r = np.ones(shape=(cols, cols))
    p = np.ones(shape=(cols, cols))

    for i in range(cols):
        for j in range(i + 1, cols):
            if method == 'pearson':
                from scipy.stats import pearsonr
                r_, p_ = pearsonr(
                    matrix[:, i], matrix[:, j])
            elif method == 'spearman':
                from scipy.stats import spearmanr
                r_, p_ = spearmanr(
                    matrix[:, i], matrix[:, j], nan_policy='omit')
            elif method == 'kendall':
                from scipy.stats import kendalltau
                r_, p_ = kendalltau(
                    matrix[:, i], matrix[:, j], nan_policy='omit')
            r[i, j] = r[j, i] = r_
            p[i, j] = p[j, i] = p_
    return r, p


def rename_index_colname(df, data):
    return pd.DataFrame(df, index=data.columns, columns=data.columns)


def split_matrix(sub, matrix):
    """splitting the matrix into chunks"""
    nfeatures = sub.shape[1]
    UpperLeft = matrix.reindex(
        columns=matrix.columns[:nfeatures], index=matrix.columns[:nfeatures])
    LowerRight = matrix.reindex(columns=matrix.columns[nfeatures:],
                                index=matrix.columns[nfeatures:])
    UpperRight = matrix.reindex(
        index=matrix.columns[:nfeatures], columns=matrix.columns[nfeatures:])
    return UpperLeft, LowerRight, UpperRight


def save_file(df, outfile, prefix):
    filename = prefix + o_format
    outpath = os.path.join(outfile, filename)
    df.to_csv(outpath, sep='\t', index=True)


def act(feature, metadata, method, output):
    if feature is not None and metadata is not None:
        feature = transpose_feature(feature, metadata)
        dfs = [feature, metadata]
        merged = pd.concat(dfs, axis=1, join="inner")
        merged = merged.apply(pd.to_numeric, errors='raise')
        calc_input = np.array(merged)

        r1, p1 = corr_calc(calc_input, method=method)
        cor = rename_index_colname(r1, merged)
        pvalue = rename_index_colname(p1, merged)

        # splitting the output matrix
        CorFeature, CorMetadata, Coef = split_matrix(feature, cor)
        PvalFeature, PvalMetadata, Pval = split_matrix(feature, pvalue)

        save_file(CorFeature, output, 'feature_correlation')
        save_file(CorMetadata, output, 'metadata_correlation')
        save_file(PvalFeature, output, 'feature_pvalue')
        save_file(PvalMetadata, output, 'metadata_pvalue')
        save_file(Coef, output, 'correlation')
        save_file(Pval, output, 'pvalue')

    else:
        dfs = [feature, metadata]
        for ele in dfs:
            if ele is not None:
                calc_input = np.array(ele)
                r1, p1 = corr_calc(calc_input, method=method)
                cor = rename_index_colname(r1, ele)
                pvalue = rename_index_colname(p1, ele)

                var_name = [x for x in globals() if globals()[x] is ele][0]
                c_prefix = var_name+'_correlation'
                p_prefix = var_name+'_pvalue'
                save_file(cor, output, c_prefix)
                save_file(pvalue, output, p_prefix)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=DESCRIPTION.Description)
    parser.add_argument('-v', '--version',
                        action='version',
                        version=DESCRIPTION.Version)
    parser.add_argument('-i', dest="feature",
                        help='input feature table, expected sampleID as colnames')
    parser.add_argument('-g', dest="group",
                        help='input mappingfile file with header')
    parser.add_argument('-d', dest="metadata",
                        help='input clinical/demographical information')
    parser.add_argument('-j', dest="method", default='spearman',
                        choices=['pearson', 'spearman', 'kendall'],
                        help=' correlation coefficient, default spearman')
    parser.add_argument("-o",  dest="output",
                        help="directory of the output files")
    parser.add_argument("-t", dest='transpose',
                        help="transpose the feature table",
                        action='store_true')
    args = parser.parse_args()

    if args.feature is None and args.metadata is None:
        sys.exit('Neither feature nor metadata is specified')

    if os.path.exists(args.output) is False:
        os.makedirs(args.output)

    feature = df_reader(args.feature)
    group = df_reader(args.group)
    metadata = df_reader(args.metadata)

    if metadata is not None:
        validate_metadata_format(metadata)

    if feature is not None and args.transpose:
        feature = feature.T
    act(feature, metadata, args.method, args.output)
