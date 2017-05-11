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
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bed', type=argparse.FileType('r'))
    parser.add_argument('quadlist', type=argparse.FileType('r'))
    parser.add_argument('fout', type=argparse.FileType('w'))
    parser.add_argument('--bg-count', type=int, default=40, help='[40]')
    args = parser.parse_args()

    quads = [q.strip() for q in args.quadlist.readlines()]

    np.random.seed(110517)
    
    for line in args.bed:
        data = line.strip().split()
        samples = data[4].split(',')
        called_quads = sorted(set([s.split('.')[0] for s in samples]))

        bg_quads = [q for q in quads if q not in called_quads]
        if len(bg_quads) >= args.bg_count:
            bg_quads = np.random.choice(quads, args.bg_count, replace=False)

        bg_samples = ['{0}.{1}'.format(quad, member) \
                      for quad in bg_quads \
                      for member in 'fa mo p1 s1'.split()]

        data[4] = ','.join(bg_samples)
        entry = '\t'.join(data)
        args.fout.write(entry + '\n')

    args.fout.close()

if __name__ == '__main__':
    main()
