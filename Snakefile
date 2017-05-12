
configfile: 'config.yaml'

subworkflow background:
    workdir: "."
    snakefile: "background.snake"

import numpy as np
from collections import defaultdict

COORDS = {}
SAMPLES = defaultdict(list)

with open(background('sr_testing.w_background.bed')) as bedfile:
    window = config['window']
    for line in bedfile:
        if line.startswith('#'):
            continue
        data = line.strip().split()
        
        chrom = data[0]
        start = int(data[1])
        end = int(data[2])
        startA = start - window
        startB = start + window
        endA = end - window
        endB = end + window
        
        if endA >= startB:
            coords = '{chrom}:{startA}-{endB}'.format(**locals())
        else:
            coords = '{chrom}:{startA}-{startB} {chrom}:{endA}-{endB}'
            coords = coords.format(**locals())

        name = data[3]

        samples = data[4].split(',')
        for sample in samples:
            SAMPLES[name].append(sample)
            # if name != 'polymorphic_cnv_3041':
            #     continue
            # NAME_SAMPLES.append('{0}__{1}'.format(name, sample))
        
        COORDS[name] = coords

NAMES = sorted(COORDS.keys())
#NAMES = ['polymorphic_cnv_3041']

rule all:
    input:
        expand('split_counts/{name}.txt', name=NAMES)
        
rule count_splits:
    input:
        bed=background("sr_testing.w_background.bed")
    output:
        counts=touch('split_counts_samples/{name}')
    params:
        coords=lambda wildcards: COORDS[wildcards.name],
        samples=lambda wildcards: SAMPLES[wildcards.name],
        min_splits=config['min_splits'],
    shell:
        """
        count_splits=$(readlink -f count_splits.py);
        fout=$(readlink -f {output});
        cd bam_indexes;
        for sample in {params.samples}; do
          samtools view -h $(s3bam $sample) {params.coords} \
            | $count_splits --min-splits {params.min_splits} $sample ${{fout}}__${{sample}};
        done
        """

def get_split_list(wildcards):
    path = expand(rules.count_splits.output.counts, 
                  name=wildcards.name, sample=SAMPLES[wildcards.name])
    return path

rule combine_splits:
    input:
        bed=config['bed'],
        count_prefix=rules.count_splits.output.counts
    params:
        samples=lambda wildcards: SAMPLES[wildcards.name],
        window=config['window']
    output:
        'split_counts/{name}.txt'
    script:
        "combine_splits.py"
