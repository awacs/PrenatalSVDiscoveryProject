#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import sys
import heapq
import itertools
from collections import deque, defaultdict
import pysam
import svtools.utils as svu


def add_samples(pesr_record, depth_record):
    # TODO: add pesr/depth FORMAT fields
    for sample in svu.get_called_samples(depth_record):
        pesr_record.samples[sample]['GT'] = (0, 1)


def merge_nested_depth_record(pesr_record, depth_record):
    """Add samples from nested depth record to PE/SR record"""

    pesr_record.info['MEMBERS'] = (pesr_record.info['MEMBERS'] +
                                   (depth_record.id, ))
    pesr_record.info['SOURCES'] = pesr_record.info['SOURCES'] + ('depth', )

    def _quad(s):
        return s.strip('.')[0]

    depth_samples = svu.get_called_samples(depth_record)
    pesr_samples = svu.get_called_samples(pesr_record)

    # If a sample is called in both the pe/sr record and the nested depth
    # record, move any relatives called in the depth record to the pe/sr record
    for quad, samples in itertools.groupby(depth_samples, _quad):
        samples = list(samples)
        if any([s in pesr_samples for s in samples]):
            for sample in samples:
                svu.set_null(depth_record, sample)
                pesr_record.samples[sample]['GT'] = (0, 1)


def merge_pesr_depth(pesr_vcf, depth_vcf, frac=0.8):
    # Memory inefficient but it's easier and shouldn't matter too much
    # now that the variants have been filtered down
    records = dict()
    records['pesr'] = {record.id: record for record in pesr_vcf}
    records['depth'] = {record.id: record for record in depth_vcf}

    # Wipe MEMBERS from prior clustering
    for source in 'pesr depth'.split():
        for ID, record in records[source].items():
            record.info['MEMBERS'] = [ID]

    # Reset for bedtool creation
    pesr_vcf.reset()
    depth_vcf.reset()
    pesr_bed = svu.vcf2bedtool(pesr_vcf, split_bnd=False,
                               include_strands=False)
    depth_bed = svu.vcf2bedtool(depth_vcf, split_bnd=False,
                                include_strands=False)

    # Merge depth records with PE/SR records if they share 80% recip overlap
    sect = pesr_bed.intersect(depth_bed, wa=True, wb=True, r=True, f=frac)

    filtered_depth_IDs = deque()
    for pair in sect.intervals:
        # Check SV types match
        if pair.fields[4] != pair.fields[9]:
            continue

        pesr_id, depth_id = pair.fields[3], pair.fields[8]

        # Add depth record's samples to PE/SR
        filtered_depth_IDs.append(depth_id)
        pesr_record = records['pesr'][pesr_id]
        depth_record = records['depth'][depth_id]

        # Update metadata and samples
        pesr_record.info['MEMBERS'] = (pesr_record.info['MEMBERS'] +
                                       (depth_record.id, ))
        pesr_record.info['SOURCES'] = pesr_record.info['SOURCES'] + ('depth', )
        add_samples(pesr_record, depth_record)

    # Remove overlapping depth records (not performed in for loop to account
    # for double overlaps
    # TODO: handle double overlap of depth calls
    for ID in set(filtered_depth_IDs):
        records['depth'].pop(ID)

    # In remaining depth-only calls, add samples to PE/SR record if the
    # record covers 90% of the depth-only call.
    sect = pesr_bed.intersect(depth_bed, wa=True, wb=True, F=0.9)

    for pair in sect.intervals:
        # Check SV types match
        if pair.fields[4] != pair.fields[9]:
            continue

        pesr_id, depth_id = pair.fields[3], pair.fields[8]

        # Skip depth records we already added with 80% reciprocal
        if depth_id in filtered_depth_IDs:
            continue

        # If sample is in both depth record and pe/sr record, remove it from
        # depth record
        depth_record = records['depth'][depth_id]
        pesr_record = records['pesr'][pesr_id]

        merge_nested_depth_record(pesr_record, depth_record)

    # Merge records together
    def _sort_key(record):
        return (record.chrom, record.pos, record.info['CHR2'], record.stop)

    pesr_records = sorted(records['pesr'].values(), key=_sort_key)
    depth_records = sorted(records['depth'].values(), key=_sort_key)
    for record in heapq.merge(pesr_records, depth_records, key=_sort_key):
        # Clean out unwanted format keys
        for key in record.format.keys():
            if key != 'GT':
                del record.format[key]

        record.info['SOURCES'] = sorted(set(record.info['SOURCES']))
        record.info['MEMBERS'] = sorted(set(record.info['MEMBERS']))

        # Skip emptied depth records
        if len(svu.get_called_samples(record)) == 0:
            continue

        yield record


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('pesr', help='VCF of PE/SR calls')
    parser.add_argument('depth', help='VCF of depth calls')
    parser.add_argument('fout')
    parser.add_argument('--prefix', default='SSC')
    args = parser.parse_args()

    pesr = pysam.VariantFile(args.pesr)
    depth = pysam.VariantFile(args.depth)

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=pesr.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=pesr.header)

    counts = defaultdict(int)

    for record in merge_pesr_depth(pesr, depth):
        svtype = record.info['SVTYPE']
        counts[svtype] += 1
        record.id = '{0}_{1}_{2}'.format(args.prefix, svtype, counts[svtype])
        fout.write(record)


if __name__ == '__main__':
    main()
