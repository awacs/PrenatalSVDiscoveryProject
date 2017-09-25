#!/bin/bash
#
# sub_pesr.sh
#
# 
#
# Copyright (C) 2016 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

set -e 

snakefile=$1
shift
rerun=$@

if [[ ! -f $snakefile ]]; then
  echo "Snakefile does not exist: $snakefile"
  exit 1
fi

snakename=$(basename $snakefile)
snakename="${snakename%.*}"

mkdir -p logs

bsub -q normal -o logs/${snakename}.out -sla miket_sc -J ${snakename}_MASTER "
snakemake \
  -s ${snakefile} \
  --printshellcmds \
  --cluster-config cluster.yaml \
  --cluster 'bsub -q {cluster.queue} -o {cluster.log} -sla miket_sc -J {cluster.jobname} {cluster.flags}' \
  --jobs 1000 \
  --stats logs/${snakename}.stats \
  --keep-going \
  ${rerun}"

