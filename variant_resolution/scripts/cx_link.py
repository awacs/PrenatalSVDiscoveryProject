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


def shared_samples(recA, recB, min_thresh=0.5, max_thresh=0.8):
    samplesA = set(svu.get_called_samples(recA))
    samplesB = set(svu.get_called_samples(recB))

    shared = samplesA & samplesB
    fracA = len(shared) / len(samplesA)
    fracB = len(shared) / len(samplesB)

    min_frac, max_frac = sorted([fracA, fracB])

    return min_frac >= min_thresh and max_frac >= max_thresh


def extract_breakpoints(vcfpath, IDs):
    """
    Extract all VCF records in list of IDs
    (Assumes VCF is sorted by variant ID)

    Parameters
    ----------
    vcfpath : str
        Path to VCF
    IDs : list of str
        Variant IDs to extract

    Returns
    -------
    bkpts : list of pysam.VariantRecord
    """

    vcf = pysam.VariantFile(vcfpath)
    n_bkpts = len(IDs)
    bkpts = np.empty(n_bkpts, dtype=object)
    idx = 0

    for record in vcf:
        if record.id in IDs:
            bkpts[idx] = record
            idx += 1
            if idx == n_bkpts:
                break

    # TODO: fix upstream VCF output so input IDs are sorted
    #  for record in vcf:
    #      if record.id == IDs[idx]:
    #          bkpts[idx] = record
    #          idx += 1
    #          if idx == n_bkpts:
    #              break

    return bkpts


def cx_link(vcfpath, bkpt_window=1000):
    """
    Parameters
    ----------
    vcfpath : str
        Path to breakpoint VCF
    """

    bt = pbt.BedTool(vcfpath)

    # Identify breakpoints which overlap within specified window
    overlap = bt.cut(range(9)).window(bt.cut(range(9)), w=bkpt_window).saveas()
    #  overlap = bt.cut(range(9)).intersect(bt.cut(range(9)), wa=True, wb=True).saveas()

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

    # Extract VariantRecords corresponding to breakpoints
    n_bkpts = len(linked_IDs)
    bkpts = extract_breakpoints(vcfpath, linked_IDs)

    # Build sparse graph from links
    G = sps.eye(n_bkpts, dtype=np.uint16, format='lil')
    for i, j in keyed_links:
        if bkpts[i] is None or bkpts[j] is None:
            print('None bkpt')
            import ipdb
            ipdb.set_trace()

        if shared_samples(bkpts[i], bkpts[j]):
            G[i, j] = 1

    # Generate lists of clustered breakpoints
    n_comp, comp_list = sps.csgraph.connected_components(G)
    clusters = [deque() for x in range(n_comp)]
    for i, c_label in enumerate(comp_list):
        clusters[c_label].append(bkpts[i])

    # Remove clusters of one variant - leftover from shared sample filtering
    clusters = [c for c in clusters if len(c) > 1]

    return clusters


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Breakpoint VCFs.')
    args = parser.parse_args()

    cx_link(args.vcf)


if __name__ == '__main__':
    main()
