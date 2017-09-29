#!/bin/bash
#
# dl_cytoband.sh
#
# Download and format cytoband bed file
#
# Copyright (C) 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz

# Convert to GRCh37 naming scheme and make tabixable
zcat cytoBand.txt.gz \
  | sed -e 's/^chr//' \
  | sort -k1,1V -k2,2n \
  | bgzip -c \
  > cytobands.bed.gz

tabix -p bed cytobands.bed.gz

rm cytoBand.txt.gz
