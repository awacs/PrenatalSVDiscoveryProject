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


def slink(links):
    """
    Single linkage cluster IDs based on whether they were merged by bedtools

    links : list of (str, str)
    """

    IDs = sorted(set([vID for link in links for vID in link]))
    n = len(IDs)

    # Permit clusters of size 1
    G = sparse.eye(n, dtype=np.uint16, format='lil')

    # Add edges between linked variants
    for v1, v2 in links:
        idx1, idx2 = IDs.index(v1), IDs.index(v2)
        G[idx1, idx2] = 1

    # Get indices of connected components
    n_comp, comp_list = csgraph.connected_components(G, connection='weak')
    cluster_names = np.arange(n_comp)

    # Convert indices back to lists of Nodes
    # Sort clusters internally by first read's position
    clusters = deque()
    for cname in cluster_names:
        cluster_idx = np.where(comp_list == cname)[0]
        cluster = itemgetter(*cluster_idx)(IDs)
        clusters.append(sorted(cluster))

    # Then sort clusters by first pair's first read's position
    return sorted(clusters)


def merge_linked_depth_calls(vcf, links, flagged=[]):
    """
    vcf : pysam.VariantFile
    links : list of (str, str)
    flagged : list of str
        Multiallelic sites
    """

    # Map variant IDs to index
    clustered_IDs = slink(links)
    idx_map = {}
    for i, cluster in enumerate(clustered_IDs):
        for j, vID in enumerate(cluster):
            idx_map[vID] = (i, j)

    clusters = [np.empty(len(c), dtype=object) for c in clustered_IDs]

    for record in vcf:
        if record.id not in idx_map:
            #  continue
            yield record
        else:
            i, j = idx_map[record.id]
            clusters[i][j] = record

    # Filter clusters on other chromosomes
    clusters = [c for c in clusters if c[0] is not None]

    # Merge clusters
    for cluster in clusters:
        # Take maximal region
        start = np.min([record.pos for record in cluster])
        end = np.max([record.stop for record in cluster])

        merged_record = cluster[0].copy()
        merged_record.pos = start
        merged_record.stop = end

        members = list(record.info['MEMBERS']) + [r.id for r in cluster]
        merged_record.info['MEMBERS'] = members

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
