
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
        
        if endA <= startB:
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

rule all:
    input:
        'split_counts/polymorphic_cnv_3041.txt'
        
rule count_splits:
    input:
        bed=background("sr_testing.w_background.bed")
    output:
        temp('split_counts_samples/{name}__{sample}.txt')
    params:
        coords=lambda wildcards: COORDS[wildcards.name],
        min_splits=config['min_splits'],
    shell:
        """
        count_splits=$(readlink -f count_splits.py);
        fout=$(readlink -f {output});
        cd bam_indexes;
        samtools view -h $(s3bam {wildcards.sample}) {params.coords} \
          | $count_splits --min-splits {params.min_splits} {wildcards.sample} $fout
        """

def get_split_list(wildcards):
    path = 'split_counts_samples/{name}__{{sample}}.txt'.format(name=wildcards.name)
    samples = SAMPLES[wildcards.name]
    return expand(path, sample=samples)

rule combine_splits:
    input:
        bed=config['bed'],
        split_counts=get_split_list
    params:
        window=config['window']
    output:
        'split_counts/{name}.txt'
    script:
        "combine_splits.py"
