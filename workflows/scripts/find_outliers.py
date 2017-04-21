#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
find_outliers.py

Identify per-algorithm, per-svtype outliers
"""


from pysam import VariantFile
from collections import defaultdict
import pandas as pd


def find_outliers(vcflist, svtype):
    """
    Locate samples which have abnormally high count of a specific variant type.

    Parameters
    ----------
    vcflist : list of str
        Filepaths to VCFs
    svtype : str
        SV class to count [DEL,DUP,INV,BND]

    Returns
    -------
    outliers : list of str
        List of outlier samples
    """

    counts = defaultdict(int)
    null_GTs = [(0, 0), (0, ), (None, None), (None, )]

    # Count calls per sample of specified svtype
    for f in vcflist:
        vcf = VariantFile(f)
        for record in vcf:
            if record.info['SVTYPE'] != svtype:
                continue
            for sample in record.samples:
                gt = record.samples[sample]['GT']
                if gt not in null_GTs:
                    counts[sample] += 1

    counts = pd.DataFrame.from_dict({'var_count': counts})

    # Write empty list if no calls of SVTYPE present
    if counts.shape[0] == 0:
        return []

    # Identify outliers
    Q1 = counts.var_count.quantile(0.25)
    Q3 = counts.var_count.quantile(0.75)
    IQR = Q3 - Q1
    cutoff = Q3 + 1.5 * IQR
    outliers = counts.loc[counts.var_count > cutoff]

    return list(outliers.index)


def main():
    vcflist = snakemake.input
    svtype = snakemake.wildcards.svtype
    source = snakemake.wildcards.source

    outliers = find_outliers(vcflist, svtype)

    with open(snakemake.output[0], 'w') as fout:
        for sample in outliers:
            fout.write(sample + '\n')


if __name__ == '__main__':
    main()
