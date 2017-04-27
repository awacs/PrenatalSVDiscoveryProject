#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import pandas as pd


def make_depth_rdtest_bed(svof):
    bed = svof['#chrA start end name svtype'.split()].drop_duplicates()

    # Add samples
    def agg_samples(samples):
        return ','.join(sorted(set(samples)))
    samples = svof.groupby('name')['sample'].agg(agg_samples)
    samples = samples.rename('samples').reset_index()
    bed = pd.merge(bed, samples, on='name', how='left')

    # Format
    rename = {'#chrA': '#chrom'}
    bed = bed.rename(columns=rename)
    bed['svtype'] = bed.svtype.str.upper()

    cols = '#chrom start end name samples svtype'.split()
    return bed[cols]

def main():
    dels = pd.read_table(snakemake.input.dels)
    dups = pd.read_table(snakemake.input.dups)
    svof = pd.concat([dels, dups]).sort_values('start')

    bed = make_depth_rdtest_bed(svof)

    bed.to_csv(snakemake.output[0], sep='\t', index=False)

if __name__ == '__main__':
    main()
