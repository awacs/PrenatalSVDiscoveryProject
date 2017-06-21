#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2017 ec2-user <mstone5@mgh.harvard.edu>
#
# Distributed under terms of the MIT license.

"""

"""

import argparse
import sys
import os
from collections import defaultdict
import pysam
from helpers import is_excluded
from s3bam import load_bam
from count_splits import get_split_position


def collect_disc(bam):
    """
    Extract discordant read pairs.

    Excludes unmapped reads, reads with an unmapped mate, duplicate reads, and
    secondary or supplementary alignments. 

    Parameters
    ----------
    bam : pysam.AlignmentFile

    Yields
    ------
    readA, readB : pysam.AlignedSegment
        Reads in a discordant pair
    """

    pairs = {}
    name_counts = defaultdict(int)

    for read in bam:
        # Restrict to unique primary alignments with a mapped mate
        # Equivalent to `samtools view -F 3340`
        if is_excluded(read):
            continue

        # Only report discordant pairs
        if read.is_proper_pair:
            continue

        # Count observations of each read pair
        name_counts[read.query_name] += 1

        # Only two reads in a pair
        if name_counts[read.query_name] > 2:
            msg = 'Third read found for ID: {0}'.format(read.query_name)
            raise Exception(msg)

        # Return pair when both reads found, then delete the first read to
        # free space
        if name_counts[read.query_name] == 2:
            yield pairs[read.query_name], read

            # TODO: check if this breaks yield
            del pairs[read.query_name]

        # First read in pair
        else:
            pairs[read.query_name] = read


def condense_pair(readA, readB):
    """
    Convert pair of reads to BEDPE format.

    Parameters
    ----------
    readA : pysam.AlignedSegment
    readB : pysam.AlignedSegment

    Returns
    -------
    entry : str
    """
    chrA = readA.reference_name
    chrB = readB.reference_name

    posA = readA.pos
    posB = readB.pos

    strandA = '-' if readA.is_reverse else '+'
    strandB = '-' if readB.is_reverse else '+'

    #  entry = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}'
    #  entry = entry.format(chrA, posA, posA + 1,
    #                       chrB, posB, posB + 1,
    #                       readA.query_name, '.', strandA, strandB)

    entry = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}'
    entry = entry.format(chrA, posA, strandA, chrB, posB, strandB)

    return entry


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bam', help='Local or S3 path to bam')
    #  parser.add_argument('fout', help='Output file.')
    parser.add_argument('--index-dir', default=None,
                        help='Directory of local BAM indexes if accessing '
                        'a remote S3 bam.')
    parser.add_argument('-r', '--region',
                        help='Tabix-formatted region to parse')
    parser.add_argument('-s', '--sample', help='Sample ID to append to output')
    #  parser.add_argument('-z', '--gzip', default=False, action='store_true',
                        #  help='Gzip output')
    #  , type=argparse.FileType('w'),
    #                    #  nargs='?', default=sys.stdout)
    args = parser.parse_args()

    if args.index_dir:
        os.chdir(args.index_dir)

    bam = load_bam(args.bam)

    #  conn_count = 0
    #  conn_success = False
    #  MAX_CONNS = 5
    #  while conn_count < MAX_CONNS and not conn_success:
        #  try:
            #  conn_success = True
        #  except:
            #  conn_count += 1

    #  if conn_count == MAX_CONNS:
        #  raise Exception('Could not load file')

    if args.region:
        bam = bam.fetch(region=args.region.encode('utf-8'))

    for readA, readB in collect_disc(bam):
        entry = condense_pair(readA, readB)
        if args.sample: 
            entry = entry + '\t' + args.sample
        
        sys.stdout.write(entry + '\n')

if __name__ == '__main__':
    main()
