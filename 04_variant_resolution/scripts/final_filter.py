#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import sys
import pandas as pd
import pysam
import svtools.utils as svu


def final_filter(vcf, bed):
    bed = bed.set_index('name')
    new_variant = (bed.operation == 'add-mosaic-variant')

    filter_names = bed.loc[~new_variant].index.values

    for record in vcf:
        if record.id in filter_names:
            op = bed.loc[record.id, 'operation']
            samples = bed.loc[record.id, 'samples'].split(',')

            if op == 'add-samples':
                for sample in samples:
                    record.samples[sample]['GT'] = (0, 1)
            elif op == 'flag-mosaic':
                for sample in samples:
                    record.samples[sample]['GT'] = (0, 1)
                record.info['MOSAIC'] = samples
            elif op == 'remove-variant':
                continue

        yield record

    for new_name in bed.loc[new_variant].index.values:
        # Stick to variants on same chromosome
        data = bed.loc[new_name]
        if data.chrom != record.chrom:
            continue

        new_record = record.copy()
        for sample in new_record.samples.keys():
            svu.set_null(new_record, sample)

        samples = data['samples'].split(',')
        for sample in samples:
            record.samples[sample]['GT'] = (0, 1)
        
        record.id = new_name
        record.chrom = data.chrom
        record.pos = data.start
        record.alts = ('<{0}>'.format(data.svtype), )
        record.info['CHR2'] = data.chrom
        record.stop = data.end

        record.info['SVTYPE'] = data.svtype
        if data.svtype == 'DEL':
            record.info['STRANDS'] = '+-'
        else:
            record.info['STRANDS'] = '-+'
        record.info['SVLEN'] = int(data.end - data.start)
        record.info['SOURCES'] = data.sources.split(',')
       
        record.info['MEMBERS'] = (new_name, )
        
        yield record


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('filter_bed')
    parser.add_argument('fout')
    args = parser.parse_args()

    bed = pd.read_table(args.filter_bed) 

    vcf = pysam.VariantFile(args.vcf)

    header = vcf.header
    header.add_line('##INFO=<ID=MOSAIC,Number=.,Type=String,Description="Samples predicted to harbor somatic or germline mosaicism">')

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=header)

    for record in final_filter(vcf, bed):
        fout.write(record)


if __name__ == '__main__':
    main()
