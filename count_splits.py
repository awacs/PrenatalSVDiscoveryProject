#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import sys
from collections import defaultdict, deque
import pandas as pd
import pysam


def collect_splits(bam):
    """
    Extract soft-clipped, unique, primary alignments.

    Excludes unmapped reads, reads with an unmapped mate, duplicate reads, and
    secondary or supplementary alignments. Reads are considered split if their
    CIGAR string contains a soft clip operation.

    Parameters
    ----------
    bam : pysam.AlignmentFile

    Yields
    ------
    read : pysam.AlignedSegment
    """

    # Restrict to unique primary alignments with a mapped mate
    def _excluded(read):
        return (read.is_unmapped or
                read.mate_is_unmapped or
                read.is_secondary or
                read.is_duplicate or
                read.is_supplementary)

    # Soft clip indicate a candidate split read
    def _is_soft_clipped(read):
        return any([tup[0] == 4 for tup in read.cigartuples])

    for read in bam:
        if _excluded(read):
            continue
        if _is_soft_clipped(read):
            yield read


def get_split_position(read):
    """
    Calculate split coordinate based on read alignment and CIGAR operations.

    Support is only present for reads soft-clipped on one side, e.g. 100M51S,
    as the coordinate is calculated by shifting the alignment position by the
    length of the flanking match operation.

    Parameters
    ----------
    read : pysam.AlignedSegment

    Returns
    -------
    pos : int
        Adjusted split read coordinate
    side : str [RIGHT,LEFT,MIDDLE]
        Direction of soft clip
    """

    pos = read.pos

    # Right soft clip - add length of aligned sequence
    if read.cigartuples[0][0] == 0:
        for operation, length in read.cigartuples:
            # Only shift based on matches, ignore DEL/INS/clips
            if operation == 0:
                pos += length
        return pos, 'RIGHT'

    # Left soft clip - sequence is already aligned to split position
    elif read.cigartuples[-1][0] == 0:
        return pos, 'LEFT'

    # Safety check - ignore match flanked by soft clips
    else:
        return None, 'MIDDLE'


class SplitStack:
    def __init__(self, splits, window=50):
        self.splits = list(splits)

        self.right_clips = defaultdict(deque)
        self.left_clips = defaultdict(deque)

    def map_splits(self):
        for split in self.splits:
            pos, side = get_split_position(split)
            if side == 'RIGHT':
                self.right_clips[pos].append(split)
            elif side == 'LEFT':
                self.left_clips[pos].append(split)

    def count_splits(self, min_splits=2):
        def _make_count_df(split_df):
            counts = {pos: len(deq) for pos, deq in split_df.items()}
            df = pd.DataFrame.from_dict({'count': counts})
            df = df.reset_index().rename(columns={'index': 'pos'})
            df = df.loc[df['count'] >= min_splits]
            return df

        right_counts = _make_count_df(self.right_clips)
        right_counts['clip'] = 'right'

        left_counts = _make_count_df(self.left_clips)
        left_counts['clip'] = 'left'

        self.split_counts = pd.concat([right_counts, left_counts])


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('sample')
    parser.add_argument('fout', type=argparse.FileType('w'))
    parser.add_argument('--min-splits', type=int, default=2)
    parser.add_argument('bam', type=argparse.FileType('rb'),
                        nargs='?', default=sys.stdin)
    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.bam)
    splits = collect_splits(bam)

    stack = SplitStack(splits)
    stack.map_splits()
    stack.count_splits(args.min_splits)
    stack.split_counts['sample'] = args.sample
    stack.split_counts.to_csv(args.fout, sep='\t', index=False)
    args.fout.close()


if __name__ == '__main__':
    main()
