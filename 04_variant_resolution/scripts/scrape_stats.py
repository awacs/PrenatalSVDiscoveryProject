#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import itertools
from collections import defaultdict
import pysam
import svtools.utils as svu


def get_inh(called):
    """
    Get list of children per inheritance status

    Returns
    -------
    inh_status : dict of {str: list of str}
        Maps inheritance status to list of corresponding samples
    """

    inh_status = dict(denovo=[], maternal=[], paternal=[], biparental=[])

    for q, samples in itertools.groupby(called, key=lambda s: s.split('.')[0]):
        samples = list(samples)
        members = [s.split('.')[1] for s in samples]

        if 'fa' in members and 'mo' in members:
            status = 'biparental'
        elif 'fa' in members:
            status = 'paternal'
        elif 'mo' in members:
            status = 'maternal'
        else:
            status = 'denovo'

        quad = samples[0].split('.')[0]
        if 'p1' not in members and 's1' not in members:
            continue

        if 'p1' in members:
            inh_status[status].append(quad + '.p1')
        if 's1' in members:
            inh_status[status].append(quad + '.s1')

    return inh_status


def scrape_record_stats(record, sample_keys):
    """
    record : pysam.VariantRecord
    sample_keys : dict of {str: list of str}
        {batch: [samples]}
    """
    name = record.id
    svtype = record.info['SVTYPE']

    svsize = record.stop - record.pos
    sources = ','.join(record.info['SOURCES'])

    called = svu.get_called_samples(record)
    parents = [s for s in called if s.endswith('fa') or s.endswith('mo')]
    children = [s for s in called if s.endswith('p1') or s.endswith('s1')]

    n_called = len(called)
    n_parents = len(parents)
    n_children = len(children)

    n_pilot = len([s for s in called if s in sample_keys['Pilot']])
    n_phase1 = len([s for s in called if s in sample_keys['Phase1']])

    quads = sorted(set([s.split('.')[0] for s in called]))
    n_quads = len(quads)

    inh_status = get_inh(called)
    n_denovo = len(inh_status['denovo'])
    n_maternal = len(inh_status['maternal'])
    n_paternal = len(inh_status['paternal'])
    n_biparental = len(inh_status['biparental'])

    chrom, start, end = record.chrom, record.pos, record.stop

    statline = ('{chrom}\t{start}\t{end}\t'
                '{name}\t{svtype}\t{svsize}\t{sources}\t'
                '{n_called}\t{n_quads}\t{n_parents}\t{n_children}\t'
                '{n_pilot}\t{n_phase1}\t'
                '{n_denovo}\t{n_maternal}\t{n_paternal}\t{n_biparental}\t'
                '{called}')

    statline = statline.format(**locals())
    return statline


class StatsScraper:
    def __init__(self, vcf, sample_keys, var_fout, obs_fout):
        """
        vcf : pysam.VariantRecord
        sample_keys : dict of {str: list of str}
            {batch : [samples]}
        var_fout : writable File object
            Per-variant statistics
        obs_fout : writable File object
            Per-sample statistics
        """

        self.vcf = vcf
        self.sample_keys = sample_keys
        self.samples = [s for slist in sample_keys.values() for s in slist]
        self.batches = list(self.sample_keys.keys())
        self.batch_key = {}
        for batch, samples in self.sample_keys.items():
            for sample in samples:
                self.batch_key[sample] = batch

        self.var_fout = var_fout
        self.obs_fout = obs_fout

        self.var_fout.write(self.var_header)
        self.obs_fout.write(self.obs_header)

        samples = [s for slist in sample_keys.values() for s in slist]
        samples = sorted(set(samples))
        self.sample_stats = {s: defaultdict(int) for s in samples}

    def scrape(self):
        for record in self.vcf:
            var_statline = scrape_record_stats(record, self.sample_keys)
            self.var_fout.write(var_statline + '\n')

            self.scrape_sample_stats(record)

        fmt = ('{sample}\t{batch}\t{quad}\t{member}\t'
               '{name}\t{chrom}\t{svtype}\t{source}\t{inh}\t{count}\n')

        for sample in self.samples:
            batch = self.batch_key[sample]
            quad, member = sample.split('.')

            obs_counts = self.sample_stats[sample]
            for (name, chrom, svtype, source, inh), count in obs_counts.items():
                self.obs_fout.write(fmt.format(**locals()))

    def scrape_sample_stats(self, record):
        chrom = record.chrom
        svtype = record.info['SVTYPE']
        if record.info['SOURCES'] == ('depth',):
            source = 'depth-only'
        elif 'depth' in record.info['SOURCES']:
            source = 'pesr+depth'
        else:
            source = 'pesr-only'

        called = svu.get_called_samples(record)
        inh_status = get_inh(called)
        inh_map = {}
        for status, samples in inh_status.items():
            for s in samples:
                inh_map[s] = status

        for sample in called:
            if sample.endswith('fa') or sample.endswith('mo'):
                obs = (record.id, chrom, svtype, source, 'parent')
            else:
                obs = (record.id, chrom, svtype, source, inh_map[sample])

            self.sample_stats[sample][obs] += 1

    @property
    def var_header(self):
        # TODO: add variable batches
        header = ('chrom\tstart\tend\t'
                  'name\tsvtype\tsvsize\tsources\t'
                  'n_called\tn_quads\tn_parents\tn_children\t'
                  'n_pilot\tn_phase1\t'
                  'n_denovo\tn_maternal\tn_paternal\tn_biparental\t'
                  'called\n')
        return header

    @property
    def obs_header(self):
        header = ('sample\tbatch\tquad\tmember\t'
                  'name\tchrom\tsvtype\tsource\tinh\tcount\n')
        return header


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('phase1', type=argparse.FileType('r'))
    parser.add_argument('pilot', type=argparse.FileType('r'))
    parser.add_argument('var_fout', type=argparse.FileType('w'))
    parser.add_argument('obs_fout', type=argparse.FileType('w'))
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)

    sample_keys = {}
    sample_keys['Pilot'] = [s.strip() for s in args.pilot.readlines()]
    sample_keys['Phase1'] = [s.strip() for s in args.phase1.readlines()]

    scraper = StatsScraper(vcf, sample_keys, args.var_fout, args.obs_fout)
    scraper.scrape()


if __name__ == '__main__':
    main()
