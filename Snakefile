
configfile: 'config.yaml'

subworkflow preprocessing:
    workdir: "preprocessing"

subworkflow algorithm_integration:
    workdir: "algorithm_integration"

import pandas as pd

BATCH_KEY = pd.read_table(config['batches'], dtype=str)
BATCHES = sorted(BATCH_KEY.batch.unique())

with open(config['groups']) as glist:
    GROUPS = [g.strip() for g in glist.readlines()]

PESR_SOURCES = config['pesr_sources']
SOURCES = PESR_SOURCES + ['depth']
CNV = config['cnv_types']

with open(config['chroms']) as clist:
    CHROMS = [c.strip() for c in clist.readlines()]

wildcard_constraints:
    source='(' + '|'.join(SOURCES) + ')',
    chrom='(' + '|'.join(CHROMS) + ')'

# Skip preprocessing unless explicitly requested
if not config['preprocess']:
    touch('checkpoints/preprocessing.done')

rule all:
    input:
        'checkpoints/algorithm_integration.done'

rule preprocessing:
    input:
        vcfs=preprocessing(expand('filtered_vcfs/{source}.{group}.vcf.gz',
                                  source=PESR_SOURCES, group=GROUPS)),
        beds=preprocessing(expand('std_beds/{batch}.{svtype}.bed.gz',
                                  batch=BATCHES, svtype=CNV))
    output:
        touch('checkpoints/preprocessing.done')
    shell:
        """
        mkdir -p algorithm_integration/input_vcfs;
        mkdir -p algorithm_integration/input_beds;
        for f in {input.vcfs}; do ln -s $(readlink -f $f) algorithm_integration/input_vcfs/; done;
        for f in {input.beds}; do ln -s $(readlink -f $f) algorithm_integration/input_beds/; done;
        """

rule algorithm_integration:
    input:
        'checkpoints/preprocessing.done',
        algorithm_integration(expand('rdtest_beds/{batch}.{source}.{chrom}.bed',
                                     batch=BATCHES, source=SOURCES, chrom=CHROMS)),
    output:
        touch('checkpoints/algorithm_integration.done')
    shell:
        """
        mkdir -p rdtest/input_beds;
        for f in {input}; do ln -s $(readlink -f $f) rdtest/input_beds/; done
        """

