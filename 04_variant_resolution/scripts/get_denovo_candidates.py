#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""


import argparse
import sys
import itertools
import pysam
import svtools.utils as svu


def get_denovo_candidates(record, max_parents=10):
    """
    Obtain list of samples which are putatively called de novo
    """
    called = svu.get_called_samples(record)
    parents = [s for s in called if s.endswith('fa') or s.endswith('mo')]

    if len(parents) > max_parents:
        return []

    denovo = []
    for quad, samples in itertools.groupby(called, lambda s: s.split('.')[0]):
        # Add putative de novo calls
        samples = list(samples)
        members = [s.split('.')[1] for s in samples]
        if 'fa' not in members and 'mo' not in members:
            denovo += samples

    return denovo


def filter_denovo_records(vcf, max_parents=10):
    for record in vcf:
        # Skip records without any Mendelian violations
        candidates = get_denovo_candidates(record, max_parents)
        if len(candidates) == 0:
            continue

        # Restrict to rare (parental VF<0.1%) variants
        c = svu.get_called_samples(record)
        parents = [s for s in c if s.endswith('fa') or s.endswith('mo')]
        if len(parents) > max_parents:
            continue

        # Skip non-stranded (wham)
        if 'STRANDS' in record.info.keys():
            if record.info['STRANDS'] not in '+- -+ ++ --'.split():
                continue

        yield record


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('fout')
    parser.add_argument('--max-parents', type=int, default=10)
    args = parser.parse_args()

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    for record in filter_denovo_records(vcf, args.max_parents):
        fout.write(record)


if __name__ == '__main__':
    main()
