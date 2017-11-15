#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Merge variants that were linked by per-sample bedtools merge
"""

import argparse
import sys
from collections import deque
from operator import itemgetter
import numpy as np
from scipy import sparse
from scipy.sparse import csgraph
import pysam
import svtools.utils as svu


def samples_overlap(recordA, recordB, upper_thresh=0.8, lower_thresh=0.5):
    # Get lists of called samples for each record
    samplesA = set(svu.get_called_samples(recordA))
    samplesB = set(svu.get_called_samples(recordB))

    # Compute fraction of each record's samples which are shared
    shared = samplesA & samplesB
    fracA = len(shared) / len(samplesA)
    fracB = len(shared) / len(samplesB)

    min_frac, max_frac = sorted([fracA, fracB])

    return min_frac >= lower_thresh and max_frac >= upper_thresh

def slink(record_links, record_map):
    """
    Single linkage cluster IDs based on whether they were merged by bedtools

    links : list of (str, str)
    """

    # Make lists of IDs and records
    # (VariantRecords aren't hashable)
    IDs = sorted(record_map.keys())
    records = np.array([record_map[ID] for ID in IDs])
    n = len(records)

    # Permit clusters of size 1
    G = sparse.eye(n, dtype=np.uint16, format='lil')

    # Add edges between linked variants
    for r1, r2 in record_links:
        idx1, idx2 = IDs.index(r1.id), IDs.index(r2.id)
        if samples_overlap(r1, r2):
            G[idx1, idx2] = 1

    # Get indices of connected components
    n_comp, comp_list = csgraph.connected_components(G, connection='weak')
    cluster_names = np.arange(n_comp)

    # Convert indices back to lists of records
    clusters = deque()
    for cname in cluster_names:
        cluster_idx = np.where(comp_list == cname)[0]
        cluster = records[cluster_idx]
        clusters.append(sorted(cluster, key=lambda r: (r.chrom, r.pos)))

    # Then sort clusters by first pair's first read's position
    return clusters


def merge_linked_depth_calls(vcf, ID_links, flagged=[]):
    """
    vcf : pysam.VariantFile
    ID_links : list of (str, str)
    flagged : list of str
        Multiallelic sites
    """

    # Make list of linked IDs and build map to corresponding records
    linked_IDs = sorted(set([ID for link in ID_links for ID in link]))
    record_map = {}

    # If a record wasn't linked with a bedtools merge, just return it
    for record in vcf:
        if record.id not in linked_IDs:
            yield record
        else:
            record_map[record.id] = record

    # Ignore links on other chromosomes
    linked_IDs = sorted(record_map.keys())
    ID_links = [l for l in ID_links 
                if l[0] in linked_IDs and l[1] in linked_IDs]

    # Convert links from pairs of IDs to pairs of records
    record_links = np.empty([len(ID_links), 2], dtype=object)
    for i, link in enumerate(ID_links):
        record_links[i, 0] = record_map[link[0]]
        record_links[i, 1] = record_map[link[1]]

    clusters = slink(record_links, record_map)

    # Merge clusters
    for cluster in clusters:
        if len(cluster) == 1:
            yield cluster[0]
            continue

        # Take maximal region
        start = np.min([record.pos for record in cluster])
        end = np.max([record.stop for record in cluster])

        merged_record = cluster[0].copy()
        merged_record.pos = start
        merged_record.stop = end

        # members = list(record.info['MEMBERS']) + [r.id for r in cluster]
        # merged_record.info['MEMBERS'] = members

        # Take union of called samples
        for record in cluster:
            called = svu.get_called_samples(record)
            for sample in called:
                merged_record.samples[sample]['GT'] = (0, 1)

        # Flag multiallelic
        if any([record.id in flagged for record in cluster]):
            merged_record.info['MULTIALLELIC'] = True

        yield merged_record


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('links', type=argparse.FileType('r'))
    parser.add_argument('fout')
    parser.add_argument('--flag-multiallelic', type=argparse.FileType('r'))
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)
    header = vcf.header

    if args.flag_multiallelic:
        flagged = [l.strip() for l in args.flag_multiallelic.readlines()]
        header.add_line('##INFO=<ID=MULTIALLELIC,Number=0,Type=Flag,Description="Multiallelic locus">')
    else:
        flagged = []

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=header)

    links = [tuple(line.strip().split()) for line in args.links.readlines()]

    for record in merge_linked_depth_calls(vcf, links, flagged):
        fout.write(record)


if __name__ == '__main__':
    main()
