
configfile: 'config.yaml'

import numpy as np

COORDS = {}
QUAD_NAMES = [] 
NAME_SAMPLES = []

with open(config['quads']) as qlist:
    QUADS = [q.strip() for q in qlist]

with open(config['bed']) as bedfile:
    window = config['window']
    for line in bedfile:
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
        quad = data[4]
        samples = data[4].split(',')
        quads = sorted(set([s.split('.')[0] for s in samples]))
        
        COORDS[name] = coords
        
        #for quad in QUADS:
        for quad in quads:
            for member in 'fa mo p1 s1'.split():
                NAME_SAMPLES.append('{0}.{1}.{2}'.format(name, quad, member))
            QUAD_NAMES.append('{0}.{1}'.format(name, quad))

rule all:
    input:
        'split_counts.txt'
        
wildcard_constraints:
    name='polymorphic_cnv_\d+',
    quad="\d{5}"

rule make_background:
    output:
        'background.bed'
    run:
        bed = open(config['bed'])
        fout = open(output[0], 'w')
        np.random.seed(110517)
        for line in bed:
            data = line.strip().split()
            samples = data[4].split(',')
            called_quads = sorted(set([s.split('.')[0] for s in samples]))
            bg_quads = [q for q in QUADS if q not in called_quads]
            bg_quads = np.random.choice(quads, 40, replace=False)
            bg_samples = ['{0}.{1}'.format(quad, member) \
                          for quad in bg_quads \
                          for member in 'fa mo p1 s1'.split()]
            data[4] = ','.join(bg_samples)
            entry = '\t'.join(data)
            fout.write(entry + '\n')
        fout.close()

rule count_splits:
    output:
        'split_counts/{name}.{sample}.txt'
    params:
        coords=lambda wildcards: COORDS[wildcards.name],
        min_splits=config['min_splits']
    shell:
        """
        samtools view -h $(s3bam {wildcards.sample}) {params.coords} \
          | ./count_splits.py --min-splits {params.min_splits} {wildcards.sample} {output} 
        """

rule combine_splits:
    input:
        bed=config['bed'],
        split_counts=expand('split_counts/{NS}.txt', NS=NAME_SAMPLES)
    params:
        window=config['window']
    output:
        'split_counts.txt'
    script:
        "combine_splits.py"
