
configfile: 'config.yaml'

import numpy as np

COORDS = {}

with open(config['bed']) as bedfile:
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
        
        COORDS[name] = coords

rule all:
    input:
        'split_counts/polymorphic_cnv_3041.txt'
        
wildcard_constraints:
    name='polymorphic_cnv_\d+',
    quad="\d{5}"


rule make_background:
    output:
        bed='sr_testing.background.bed'
    shell:
        "./choose_background.py {config[bed]} {config[quads]} {output}"

rule merge_background:
    input:
        bed=config['bed'],
        bg=rules.make_background.output.bed
    output:
        bed='sr_testing.w_background.bed'
    run:
        import pandas as pd
        bed = pd.read_table(input.bed)
        bg = pd.read_table(input.bg)
        bg_samples = bg['name samples'.split()].rename(columns={'samples': 'bg_samples'})
        bed = pd.merge(bed, bg_samples, on='name', how='left')
        bed['samples'] = bed['samples'] + ',' + bed['bg_samples']
        cols = '#chrom start end name samples svtype batch'.split()
        bed[cols].to_csv(output[0], index=False, sep='\t')


rule count_splits:
    input:
        bed=rules.merge_background.output.bed
    output:
        'split_counts/{name}.txt'
    params:
        coords=lambda wildcards: COORDS[wildcards.name],
        min_splits=config['min_splits'],
    shell:
        """
        bed=$(readlink -f {input.bed});
        fout=$(readlink -f {output});
        count_splits=$(readlink -f count_splits.py)
        cd bam_indexes;
        awk '($4=="{wildcards.name}") {{print $5}}' $bed | sed -e 's/,/\\n/g' \
          | while read sample; do
              samtools view -h $(s3bam $sample) {params.coords} \
                | $count_splits --min-splits {params.min_splits} $sample $fout;
            done;
        cd ..
        """

# rule combine_splits:
#     input:
#         bed=config['bed'],
#         split_counts=expand('split_counts/{NS}.txt', NS=NAME_SAMPLES)
#     params:
#         window=config['window']
#     output:
#         'split_counts.txt'
#     script:
#         "combine_splits.py"
