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
import numpy as np
import pysam
from svtools.genomeslink import GenomeSLINK, GSNode


class _DiscPair(GSNode):
    def __init__(self, chrA, posA, strandA, chrB, posB, strandB, sample,
                 name='.'):
        self.strandA = strandA
        self.strandB = strandB
        self.sample = sample

        super().__init__(chrA, posA, chrB, posB, name)

    def clusters_with(self, other, dist):
        return (super().clusters_with(other, dist) and
                self.strandA == other.strandA and
                self.strandB == other.strandB and
                self.sample == other.sample)

    def __str__(self):
        entry = super().__str__()
        entry = entry + '\t{sample}\t{strandA}\t{strandB}'
        entry = entry.format(**self.__dict__)
        return entry


def _DiscParser(discfile, region):
    """

    Parameters
    ----------
    discfile : pysam.TabixFile
    regions : list of str
    """

    if isinstance(region, str):
        region = region.encode('utf-8')

    pairs = discfile.fetch(region=region, parser=pysam.asTuple())
    for pair in pairs:
        yield _DiscPair(*pair)


class PECluster(GenomeSLINK):
    def __init__(self, nodes, dist, record, window):
        self.record = record
        self.window = window

        super().__init__(nodes, dist)

    def filter_nodes(self):
        chrB = self.record.info['CHR2']
        posB = self.record.info['END']
        strandA, strandB = self.record.info['STRANDS']

        # Pairs were selected based on window around chrA; check chrB
        for node in self.nodes:
            if node.chrB != chrB:
                continue
            if np.abs(node.posB - posB) > self.window:
                continue
            if node.strandA != strandA or node.strandB != strandB:
                continue
            yield node


class PETest:
    def __init__(self, variants, discfile, window=1000, dist=300):
        """
        variants : pysam.VariantFile
        discfile : pysam.TabixFile
            chrA, posA, strandA, chrB, posB, strandB, sample
        window : int, optional
            Window around variant start/end to query for discordant pairs
        dist : int, optional
            Clustering distance
        """

        self.variants = variants
        self.discfile = discfile
        self.window = window
        self.dist = dist

    def run(self):
        for record in self.variants:
            # Temporary tloc filter
            if record.chrom != record.info['CHR2']:
                continue
            clusters = self.cluster(record)
            import ipdb
            ipdb.set_trace()
            yield clusters
            pass

    def cluster(self, record):
        chrA = record.chrom
        posA = record.pos
        strandA, strandB = record.info['STRANDS']

        reg = '{0}:{1}-{2}'
        region = reg.format(chrA, posA - self.window, posA + self.window)

        pairs = _DiscParser(self.discfile, region)
        slink = PECluster(pairs, self.dist, record, self.window)

        return slink.cluster()
        #  for cluster in slink.cluster():
            #  yield cluster


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('variants', help='Variants (default=VCF).')
    parser.add_argument('disc', help='Table of discordant pair coordinates.')
    parser.add_argument('-w', '--window', type=int, default=1000,
                        help='Window around breakpoint to query for '
                        'discordant pairs. [1000]')
    parser.add_argument('-d', '--dist', type=int, default=300,
                        help='Clustering distance. [300]')

    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    if args.variants in ['-', 'stdin']:
        variantfile = pysam.VariantFile(sys.stdin)
    else:
        variantfile = pysam.VariantFile(args.variants)

    discfile = pysam.Tabixfile(args.disc)

    petest = PETest(variantfile, discfile)
    for clusters in petest.run():
        pass



if __name__ == '__main__':
    main(sys.argv[1:])
