#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import itertools
from collections import deque
import numpy as np
import scipy.sparse as sps
import pysam
import pybedtools as pbt
import natsort
import svtools.utils as svu


def cx_link(vcfpath, bkpt_window=1000):
    """
    Parameters
    ----------
    vcfpath : str
        Path to breakpoint VCF
    """

    bt = pbt.BedTool(vcfpath)
    vcf = pysam.VariantFile(vcfpath)

    # Identify breakpoints which overlap within specified window
    overlap = bt.cut(range(9)).window(bt.cut(range(9)), w=bkpt_window).saveas()

    # Exclude self-hits
    overlap = overlap.filter(lambda b: b.fields[2] != b.fields[11]).saveas()

    # Restrict to overlaps involving a BCA breakpoint
    cnvtypes = '<DEL> <DUP>'.split()
    overlap = overlap.filter(lambda b: b.fields[4] not in cnvtypes).saveas()

    # Get linked variant IDs
    links = [(b[2], b[11]) for b in overlap.intervals]
    linked_IDs = natsort.natsorted(set(itertools.chain.from_iterable(links)))
    linked_IDs = np.array(linked_IDs)

    # Map variant IDs to indices
    link_key = {ID: i for i, ID in enumerate(linked_IDs)}
    keyed_links = np.array([(link_key[a], link_key[b]) for a, b in links])

    # Extract VariantRecords from VCF (assumes VCF is sorted by variant ID)
    n_variants = len(linked_IDs)
    bkpts = np.empty(n_variants, dtype=object)
    idx = 0
    for record in vcf:
        if record.id == linked_IDs[idx]:
            bkpts[idx] = record
            idx += 1
            if idx == n_variants:
                break

    # Build sparse graph from links
    G = sps.eye(n_variants, dtype=np.uint16, format='lil')
    for i, j in keyed_links:
        if shared_samples(bkpts[i], bkpts[j], sample_overlap):
            G[i, j] = 1

    # Generate lists of clustered breakpoints
    n_comp, comp_list = sps.csgraph.connected_components(G)
    clusters = [deque() for x in range(n_comp)]
    for i, c_label in enumerate(comp_list):
        clusters[c_label].append(bkpts[i])

    return clusters
    for cluster in clusters:
        cluster = resolve_cluster(cluster)
    import ipdb
    ipdb.set_trace()


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Breakpoint VCFs.')
    args = parser.parse_args()

    cx_link(args.vcf)


if __name__ == '__main__':
    main()
