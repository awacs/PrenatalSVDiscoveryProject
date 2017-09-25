#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""


import argparse
import itertools
import numpy as np
import pysam
import svtools.utils as svu
import pandas as pd
from svtools.pesr import PESRTestRunner, SRTest, PETest


def get_denovo_candidates(record, max_parents=20):
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


class DenovoTestRunner(PESRTestRunner):
    def __init__(self, vcf, countfile, discfile, pe_fout, sr_fout,
                 n_background=160, max_parental_vf=0.01):

        super().__init__(vcf, n_background)

        self.srtest = SRTest(countfile)
        self.sr_fout = sr_fout

        self.petest = PETest(discfile)
        self.pe_fout = pe_fout

        self.max_parental_vf = max_parental_vf
        parents = [s for s in self.samples
                   if s.endswith('fa') or s.endswith('mo')]
        self.max_n_parents = max_parental_vf * len(parents)

    def run(self):
        for record in self.vcf:
            # Skip records without any Mendelian violations
            candidates = get_denovo_candidates(record, self.max_n_parents)
            if len(candidates) == 0:
                continue

            # Restrict to rare (parental VF<0.1%) variants
            c = svu.get_called_samples(record)
            parents = [s for s in c if s.endswith('fa') or s.endswith('mo')]
            if len(parents) > self.max_n_parents:
                continue

            # Skip non-stranded (wham)
            if record.info['STRANDS'] not in '+- -+ ++ --'.split():
                continue

            self.test_record(record)

    def test_record(self, record):
        candidates = get_denovo_candidates(record, self.max_n_parents)

        for child in candidates:
            called = [child]
            background = self.choose_background(record, candidates)

            self.sr_test(record, called, background)
            self.pe_test(record, called, background)

    def sr_test(self, record, called, background):
        sr_results = self.srtest.test_record(record, called, background)

        child = called[0]
        sr_results['sample'] = child

        # Check parents at predicted coordinates in child
        posA, posB = sr_results.set_index('coord')['pos'][['posA', 'posB']]
        parents = self.sr_test_parents(record, child, background,
                                       posA, posB)

        # Merge parents and children
        results = pd.concat([sr_results, parents], ignore_index=True)

        # sneak posB-posA distance into output
        dist = posB - posA
        results.loc[results.coord == 'sum', 'pos'] = dist

        cols = 'name sample coord pos log_pval called background'.split()
        results = results[cols]
        results['pos'] = results.pos.astype(int)
        results = results.rename(columns={'called': 'called_median',
                                          'background': 'bg_median'})

        results.to_csv(self.sr_fout, index=False, header=False,
                       sep='\t', na_rep='NA', float_format='%.20f')

    def pe_test(self, record, called, background):
        pe_results = self.petest.test_record(record, called, background)

        child = called[0]
        pe_results['sample'] = child

        quad = child.split('.')[0]
        p_results = []
        for member in 'fa mo'.split():
            parent = quad + '.' + member
            result = self.petest.test_record(record, [parent], background)
            result['sample'] = parent
            p_results.append(result)

        p_results = pd.concat(p_results)
        results = pd.concat([pe_results, p_results], ignore_index=True)

        cols = 'name sample log_pval called background'.split()
        results[cols].to_csv(self.pe_fout, index=False, header=False,
                             sep='\t', na_rep='NA', float_format='%.20f')

    def choose_background(self, record, candidates):
        # Exclude called samples and all candidate families from background
        quads = sorted(set([s.split('.')[0] for s in candidates]))
        members = 'fa mo p1 s1'.split()
        related = ['.'.join(s) for s in itertools.product(quads, members)]

        called = set(svu.get_called_samples(record))
        blacklist = called.union(related)

        background = [s for s in self.samples if s not in blacklist]

        if len(background) >= self.n_background:
            background = np.random.choice(background, self.n_background,
                                          replace=False).tolist()

        return background

    def sr_test_parents(self, record, child, background, posA, posB):
        quad = child.split('.')[0]
        results = []
        for member in 'fa mo'.split():
            parent = quad + '.' + member

            p_results = []
            for pos, strand in zip([posA, posB], record.info['STRANDS']):
                result = self.srtest.test(record.chrom, pos, strand,
                                          [parent], background)
                result = result.to_frame().transpose()
                result['coord'] = 'posA' if pos == posA else 'posB'
                result['pos'] = pos

                p_results.append(result)

            p_results = pd.concat(p_results, ignore_index=True)
            p_total = self.srtest._test_total(p_results)

            p_results = pd.concat([p_results, p_total], ignore_index=True)
            p_results['sample'] = parent
            results.append(p_results)

        results = pd.concat(results, ignore_index=True)
        results['name'] = record.id

        return results


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('-c', '--countfile', required=True)
    parser.add_argument('-d', '--discfile', required=True)
    parser.add_argument('--background', type=int, default=160)
    parser.add_argument('--max-parental-vf', type=float, default=0.01)
    parser.add_argument('petest', type=argparse.FileType('w'), help='fout')
    parser.add_argument('srtest', type=argparse.FileType('w'), help='fout')
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)
    countfile = pysam.TabixFile(args.countfile)
    discfile = pysam.TabixFile(args.discfile)

    header = 'name sample log_pval called_median bg_median'.split()
    args.petest.write('\t'.join(header) + '\n')

    header = 'name sample coord pos log_pval called_median bg_median'.split()
    args.srtest.write('\t'.join(header) + '\n')

    runner = DenovoTestRunner(vcf, countfile, discfile,
                              args.petest, args.srtest,
                              args.background, args.max_parental_vf)
    runner.run()


if __name__ == '__main__':
    main()
