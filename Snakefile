
configfile: 'config.yaml'

subworkflow background:
    workdir: "."
    snakefile: "background.snake"

REGIONS = {}
with open(background('calls/sr_windows.txt')) as regionfile:
    for line in regionfile:
        name, svtype, regions = line.strip().split('\t')
        REGIONS[name] = regions

NAMES = sorted(REGIONS.keys())
NAMES = ['polymorphic_cnv_3041']

rule all:
    input:
        expand('split_counts/{name}.txt', name=NAMES)
        
rule count_splits:
    input:
        sample_list=background('sample_lists/{name}.txt')
    output:
        counts='split_counts/{name}.txt',
        log='count_logs/{name}.log'
    params:
        regions=lambda wildcards: REGIONS[wildcards.name],
    shell:
        """
        ./count_splits.sh {input.sample_list} {output.counts} {output.log} "{params.regions}"
        """

# def get_split_list(wildcards):
#     path = 'split_counts_samples/{name}__{{sample}}.txt'.format(name=wildcards.name)
#     samples = SAMPLES[wildcards.name]
#     return expand(path, sample=samples)
# 
# rule combine_splits:
#     input:
#         bed=config['bed'],
#         split_counts=get_split_list
#     params:
#         window=config['window']
#     output:
#         'split_counts/{name}.txt'
#     script:
#         "combine_splits.py"
