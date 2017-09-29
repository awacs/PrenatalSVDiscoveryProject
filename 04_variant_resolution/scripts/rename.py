#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import sys
from collections import defaultdict
import pysam


def rename(vcf, fout):
    indexes = defaultdict(int)
    fmt = 'SSC_{svtype}_{chrom}_{idx}'

    for record in vcf:
        chrom = record.chrom
        svtype = record.info['SVTYPE']
        indexes[svtype] += 1
        idx = indexes[svtype]

        record.id = fmt.format(**locals())

        # Clean metadata
        end = record.stop
        record.ref = 'N'
        record.stop = end

        for info in 'CIPOS CIEND STRANDS RMSSTD MEMBERS'.split():
            if info in record.info.keys():
                record.info.pop(info)

        fout.write(record)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('fout')
    args = parser.parse_args()

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin) 
    else:
        vcf = pysam.VariantFile(args.vcf) 

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    rename(vcf, fout)


if __name__ == '__main__':
    main()
