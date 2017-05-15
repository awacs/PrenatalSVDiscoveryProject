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


def add_window(chrom, start, end, window):

    startA = start - window
    startB = start + window
    endA = end - window
    endB = end + window
    
    if endA <= startB:
        region = '{chrom}:{startA}-{endB}'.format(**locals())
    else:
        region = '{chrom}:{startA}-{startB} {chrom}:{endA}-{endB}'
        region = region.format(**locals())

    return region


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bed', type=argparse.FileType('r'))
    parser.add_argument('window', type=int)
    parser.add_argument('fout', type=argparse.FileType('w'))
    args = parser.parse_args()

    for line in args.bed:
        if line.startswith('#'):
            continue

        data = line.strip().split()
        chrom = data[0]
        start, end = [int(x) for x in data[1:3]]

        region = add_window(chrom, start, end, args.window)

        name = data[3]
        svtype = data[5]
        entry = '{name}\t{svtype}\t{region}\n'
        entry = entry.format(**locals())
        args.fout.write(entry)

    args.fout.close()


if __name__ == '__main__':
    main()
