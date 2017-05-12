
configfile: 'config.yaml'

subworkflow background:
    workdir: "."
    snakefile: "background.snake"

from collections import defaultdict

COORDS = {}
SAMPLES = defaultdict(list)

with open(background('sr_testing.w_background.bed')) as bedfile:
    pass

REGIONS = {}
with open(background('calls/sr_windows.txt')) as regionfile:
    for line in regionfile:
        name, svtype, regions = line.strip().split('\t')
        REGIONS[name] = regions

NAMES = sorted(COORDS.keys())

rule all:
    input:
        expand('split_counts/{name}.txt', name=NAMES)
        
rule count_splits:
    input:
        sample_list=background('sample_lists/{name}.txt')
    output:
        'split_counts/{name}.txt'
    params:
        regions=lambda wildcards: REGIONS[wildcards.name],
        min_splits=config['min_splits'],
    shell:
        """
        ./count_splits.sh {input.sample_list} "{params.regions}" {output}
        """
        # """
        # count_splits=$(readlink -f count_splits.py);
        # fout=$(readlink -f {output});
        # cd bam_indexes;
        # samtools view -h $(s3bam {wildcards.sample}) {params.coords} \
        #   | $count_splits --min-splits {params.min_splits} {wildcards.sample} $fout
        # """

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
