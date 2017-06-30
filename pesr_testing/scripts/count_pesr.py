#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import os
import subprocess
import sys
import gzip
from collections import defaultdict, deque
import numpy as np
import pysam
import boto3
from helpers import is_excluded, is_soft_clipped
import helpers
from s3bam import load_bam


class PESRCollection:
    def __init__(self, bam, splitfile, max_split_dist=300):
        self.bam = bam
        self.splitfile = splitfile

        # SR evidence
        self.right_split_counts = defaultdict(int)
        self.left_split_counts = defaultdict(int)
        self.prev_split_pos = None
        self.curr_chrom = None
        self.max_split_dist = max_split_dist

    def collect_pesr(self):
        """
        Collect PE and SR evidence from a BAM file.

        Excludes unmapped reads, reads with an unmapped mate, duplicate reads,
        and secondary or supplementary alignments. Reads are considered split
        if their CIGAR string contains a soft clip operation.
        """

        for read in self.bam:
            # Restrict to unique primary alignments with a mapped mate
            # Equivalent to `samtools view -F 3340`
            if is_excluded(read):
                continue

            # Soft clip indicate a candidate split read
            if is_soft_clipped(read):
                self.count_split(read)

        self.flush_split_counts()

    def count_split(self, read):
        """
        Count splits at each position.

        splits : iter of pysam.AlignedSegment
        max_dist : int
            Max distance between consecutive splits before parsing
        """

        pos, side = get_split_position(read)

        # Calculate distance to previous split and update position tracker
        # Use abs to catch contig switches
        if self.prev_split_pos is None:
            dist = 0
        else:
            dist = np.abs(pos - self.prev_split_pos)
        self.prev_split_pos = pos

        if self.curr_chrom is None:
            self.curr_chrom = read.reference_name

        # Flush aggregated split reads if we've moved beyond the max dist
        if dist > self.max_split_dist:
            self.flush_split_counts()
            self.curr_chrom = read.reference_name

        # Tally the split at its corresponding position
        if side == 'RIGHT':
            self.right_split_counts[pos] += 1
        elif side == 'LEFT':
            self.left_split_counts[pos] += 1

    def flush_split_counts(self):
        """
        Write current split counts to disk and reset dictionaries
        """

        # Compile counts collected so far
        entries = deque()
        for clip in 'left right'.split():
            df = getattr(self, '%s_split_counts' % clip)

            for pos, count in df.items():
                entries.append((self.curr_chrom, pos, clip, count))

        # Sort in chunks as we go
        entries = sorted(entries, key=lambda s: s[1])

        # Flush to disk
        fmt = '%s\t%d\t%s\t%d\n'
        for entry in entries:
            self.splitfile.write((fmt % entry).encode('utf-8'))

        # Reset split counts
        self.right_split_counts = defaultdict(int)
        self.left_split_counts = defaultdict(int)


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


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bam', help='Local or S3 path to bam')
    parser.add_argument('--index-dir', default=None,
                        help='Directory of local BAM indexes if accessing '
                        'a remote S3 bam.')
    parser.add_argument('-r', '--region',
                        help='Tabix-formatted region to parse')
    parser.add_argument('fout', help='Output file.')
    parser.add_argument('-z', '--bgzip', default=False, action='store_true',
                        help='bgzip and tabix index output')
    args = parser.parse_args()

    bam = load_bam(args.bam)
    if args.region:
        bam = bam.fetch(args.region.encode('utf-8'))

    # Open output file
    if args.fout in '- stdout'.split():
        fout = sys.stdout.buffer
    else:
        fout = open(args.fout, 'wb')

    # Pass through bgzip if requested
    if args.bgzip:
        pipe = subprocess.Popen(['bgzip', '-c'],
                                stdin=subprocess.PIPE,
                                stdout=fout)
        fout = pipe.stdin

    pesr = PESRCollection(bam, fout)
    pesr.collect_pesr()

    if args.bgzip and args.fout not in '- stdout'.split():
        stdout, stderr = pipe.communicate()
        tabix = 'tabix -f -s1 -b2 -e2 %s' % args.fout
        subprocess.call(tabix.split())


    # Get splits and pile up
    #  splits = helpers.collect_splits(bam)
    #  for counts in count_splits(splits):
        #  entry = counts + '\n'
        #  # gzip requires bytes
        #  if args.gzip:
            #  entry = entry.encode('utf-8')

        #  fout.write(entry)

    #  fout.close()


if __name__ == '__main__':
    main()
