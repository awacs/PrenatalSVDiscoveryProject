#!/bin/bash
#
# filter_blacklist.sh
#
# Print IDs which are removable based on blacklist coverage
#
# Copyright (C) 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

vcf=$1
blacklist=$2

bed=${vcf}.tmp.bed
svtools vcf2bed $vcf $bed

# 30% coverage
# svtools vcf2bed $vcf stdout \
bedtools coverage -a $bed -b $blacklist \
  | awk '($10 >= 0.3) {print $4}'

# or direct overlap with start
# svtools vcf2bed $vcf stdout \
awk -v OFS="\t" '{$3=$2+1; print;}' $bed \
  | bedtools intersect -a stdin -b $blacklist \
  | cut -f4

# or direct overlap with end
# svtools vcf2bed $vcf stdout \
awk -v OFS="\t" '{end=$3; $2=end; $3=end+1; print;}' $bed \
  | bedtools intersect -a stdin -b $blacklist \
  | cut -f4

rm $bed
