#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
from collections import deque
import numpy as np
import pandas as pd


def aggregate():
    pass


def process_rdtest(rdtest):
    """Standardize rdtest column names"""

    # Drop metadata columns (available from VCF) and rename CNVID
    skip_cols = 'chr Start End SampleIDs Type'.split()
    rdtest = rdtest.drop(skip_cols, axis=1).rename(columns={'CNVID': 'name'})

    numeric_cols = 'Median_Power P 2ndMaxP Median_Rank Median_Separation'
    numeric_cols = numeric_cols.split()

    # Replace strings with NA
    for col in numeric_cols:
        repl = 'All_samples_called_CNV_no_analysis'
        rdtest[col] = rdtest[col].replace(repl, np.nan).astype(np.float)

    rdtest['log_pval'] = -np.log10(rdtest.P)
    rdtest['log_2ndMaxP'] = -np.log10(rdtest['2ndMaxP'])

    maxp = rdtest.loc[rdtest.log_pval != np.inf, 'log_pval'].max()
    max2p = rdtest.loc[rdtest.log_2ndMaxP != np.inf, 'log_2ndMaxP'].max()

    rdtest.loc[rdtest.log_pval == np.inf, 'log_pval'] = maxp + 5
    rdtest.loc[rdtest.log_2ndMaxP == np.inf, 'log_2ndMaxP'] = max2p + 5

    return rdtest


def process_srtest(srtest):
    metrics = 'log_pval called_median bg_median'.split()
    srtest = srtest.pivot_table(index='name', values=metrics, columns='coord')
    srtest.columns = ['_'.join(col[::-1]).strip()
                      for col in srtest.columns.values]
    srtest = srtest.reset_index()

    return srtest


def preprocess(df, dtype):
    if dtype == 'RD':
        return process_rdtest(df)
    elif dtype == 'SR':
        return process_srtest(df)
    else:
        return df


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-r', '--RDtest', required=True)
    parser.add_argument('-s', '--SRtest', required=True)
    parser.add_argument('-p', '--PEtest', required=True)
    parser.add_argument('fout')
    args = parser.parse_args()

    evidence = deque()

    for dtype in 'PE SR RD'.split():
        df = pd.read_table(getattr(args, dtype + 'test'))
        df = preprocess(df, dtype)
        df = df.rename(columns=lambda c: dtype + '_' + c if c != 'name' else c)
        df = df.set_index('name')
        evidence.append(df)

    evidence = list(evidence)
    evidence = evidence[0].join(evidence[1:], how='outer', sort=True)
    evidence = evidence.reset_index()

    evidence.to_csv(args.fout, index=False, sep='\t', na_rep='NA')


if __name__ == '__main__':
    main()
