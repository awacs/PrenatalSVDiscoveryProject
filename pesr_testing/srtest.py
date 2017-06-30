#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import numpy as np
import scipy.stats as ss
import pandas as pd


def summary_stats(series):
    return series.mean(), series.std(), series.shape[0]


def calc_ttest(df):
    called = df.loc[df.call_status == 'called', 'count']
    called_mu, called_sigma, called_n = summary_stats(called)

    background = df.loc[df.call_status == 'background', 'count']
    bg_mu, bg_sigma, bg_n = summary_stats(background)

    t, _p = ss.ttest_ind_from_stats(called_mu, called_sigma, called_n,
                                    bg_mu, bg_sigma, bg_n)

    if t > 0:
        p = _p / 2
    else:
        p = 1

    return pd.Series([called_mu, called_sigma, bg_mu, bg_sigma, 
                      called_n, bg_n, -np.log10(p)])


def srtest(counts):
    pvals = counts.groupby('coord dist'.split()).apply(calc_ttest)

    cols = {
        0: 'called_mean',
        1: 'called_std',
        2: 'bg_mean',
        3: 'bg_std',
        4: 'called_n',
        5: 'bg_n',
        6: 'log_pval'
    }

    pvals = pvals.rename(columns=cols).reset_index()

    return pvals


def main():
    counts = pd.read_table(snakemake.input[0])
    pvals = srtest(counts)
    pvals['name'] = snakemake.wildcards.name
    cols = ('name coord dist log_pval called_mean called_std bg_mean bg_std '
            'called_n bg_n').split()
    pvals[cols].to_csv(snakemake.output[0], sep='\t', index=False)


if __name__ == '__main__':
    main()
