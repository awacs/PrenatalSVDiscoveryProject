#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import os
from collections import deque
import pandas as pd


def get_coord(row):
    if row.svtype == 'DEL':
        if row['clip'] == 'right':
            return 'start'
        elif row['clip'] == 'left':
            return 'end'
    elif row.svtype == 'DUP':
        if row['clip'] == 'right':
            return 'end'
        elif row['clip'] == 'left':
            return 'start'


def coord_dist(row):
    coord = row['coord']

    variant_pos = row[coord]
    return variant_pos - row.pos


def main():
    # get list of samples called in each variant
    bed = open(snakemake.input.bed)
    called_keys = deque()
    for line in bed:
        data = line.strip().split()
        name = data[3]
        samples = data[4].split(',')
        for sample in samples:
            called_keys.append('{name}__{sample}'.format(**locals()))

    # Read in each dataframe (until pysam can read from S3)
    dfs = deque()
    for fname in snakemake.input.split_counts:
        name = os.path.basename(fname).split('__')[0]
        df = pd.read_table(fname)
        df['name'] = name
        dfs.append(df)
    df = pd.concat(dfs)

    # Label called samples
    df['called_key'] = df['name'] + '__' + df['sample']
    df.loc[df.called_key.isin(called_keys), 'call_status'] = 'called'
    df.loc[~df.called_key.isin(called_keys), 'call_status'] = 'background'

    # Add call coordinates
    bed_df = pd.read_table(snakemake.input.bed)
    bed_df = bed_df.drop('samples', axis=1)
    df = pd.merge(df, bed_df, on='name', how='left')

    # Calculate distance from start/end
    df['coord'] = df.apply(get_coord, axis=1)
    df['dist'] = df.apply(coord_dist, axis=1)

    # Filter by distance to appropriate coordinate.
    # Accounts for cases where a split is within the window around one
    # coordinate, but is clipped in the opposite direction
    df = df.loc[df['dist'].abs() < snakemake.params.window].copy()

    # Output
    df = df.drop('called_key', axis=1)
    df['count'] = df['count'].astype(int)
    cols = 'name sample call_status coord dist count'.split()
    df[cols].to_csv(snakemake.output[0], index=False, sep='\t')


if __name__ == '__main__':
    main()
