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
import os
import pandas as pd
import pysam
import boto3


class ReadLimitError(Exception):
    """Exceeded read limit"""


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
    # Equivalent to `samtools view -F 3340`
    def _excluded(read):
        return (read.is_unmapped or
                read.mate_is_unmapped or
                read.is_secondary or
                read.is_duplicate or
                read.is_supplementary)

    # Soft clip indicate a candidate split read
    def _is_soft_clipped(read):
        return any([tup[0] == 4 for tup in read.cigartuples])

    for i, read in enumerate(bam):
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

    def count_splits(self, min_splits=1):
        def _make_count_df(split_df, clip):
            # None objects are dropped by pd.concat
            if len(split_df.keys()) == 0:
                return None

            # Else count by position
            counts = {pos: len(deq) for pos, deq in split_df.items()}
            df = pd.DataFrame.from_dict({'count': counts})
            df = df.reset_index().rename(columns={'index': 'pos'})
            df = df.loc[df['count'] >= min_splits]
            df['clip'] = clip
            return df

        right_counts = _make_count_df(self.right_clips, 'right')
        left_counts = _make_count_df(self.left_clips, 'left')

        self.split_counts = pd.concat([right_counts, left_counts])


def load_s3bam(bucket, bam_path, filepath_index=None):
    s3 = boto3.client('s3')
    url = s3.generate_presigned_url(
            ClientMethod='get_object',
            Params={'Bucket': bucket, 'Key': bam_path})

    return pysam.AlignmentFile(url, filepath_index=filepath_index)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bam', help='Local or S3 path to bam')
    parser.add_argument('--min-splits', type=int, default=1)
    parser.add_argument('--index-dir', default=None,
                        help='Directory of local BAM indexes if accessing '
                        'a remote S3 bam.')
    parser.add_argument('-r', '--region',
                        help='Tabix-formatted region to parse')
    parser.add_argument('fout', type=argparse.FileType('w'),
                        nargs='?', default=sys.stdout)
    args = parser.parse_args()

    bam_path = args.bam
    if args.bam.startswith('s3://'):
        s3path = bam_path[5:]
        bucket = s3path.split('/')[0]
        bam_path = '/'.join(s3path.split('/')[1:])

        # Get index if possible
        bam_name = os.path.basename(bam_path)
        if args.index_dir is None:
            filepath_index = None
        else:
            idx1 = bam_name + '.bai'
            idx2 = os.path.splitext(bam_name)[0] + '.bai'
            if os.path.exists(idx1):
                filepath_index = idx1
            elif os.path.exists(idx2):
                filepath_index = idx2
            else:
                filepath_index = None

        bam = load_s3bam(bucket, bam_path, filepath_index)
    else:
        bam = pysam.AlignmentFile(args.bam)

    if args.region:
        chrom, pos = args.region.split(':')
        start, end = [int(x) for x in pos.split('-')]
        bam = bam.fetch(chrom, start, end)

    splits = collect_splits(bam)

    stack = SplitStack(splits)
    stack.map_splits()
    stack.count_splits(args.min_splits)
    stack.split_counts.to_csv(args.fout, sep='\t', index=False, header=False)
    args.fout.close()


if __name__ == '__main__':
    main()
