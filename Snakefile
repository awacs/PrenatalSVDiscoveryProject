
configfile: 'config.yaml'

PESR_SOURCES = config['pesr_sources']
DEPTH_SOURCES = config['depth_sources']
SOURCES = PESR_SOURCES + DEPTH_SOURCES

CNV = config['cnv_types']


with open(config['quads']) as qlist:
    QUADS = [q.strip() for q in qlist.readlines()]

with open(config['samples']) as slist:
    SAMPLES = [s.strip() for s in slist.readlines()]

wildcard_constraints:
    source='(' + '|'.join(SOURCES) + ')',
    sample='(' + '|'.join(SAMPLES) + ')'

include: 'rules/pesr_preprocessing.rules'
include: 'rules/depth_preprocessing.rules'

rule all:
    input:
        expand('preprocessing/filtered_vcfs/{source}.{quad}.vcf',
               source=PESR_SOURCES, quad=QUADS),
        expand('preprocessing/std_beds/merged_algs/{sample}.{svtype}.merged.bed',
               sample=SAMPLES, svtype=CNV)
