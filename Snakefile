
configfile: 'config.yaml'

SOURCES = config['pesr_sources']

with open(config['quads']) as qlist:
    QUADS = [q.strip() for q in qlist.readlines()]

wildcard_constraints:
    source='(' + '|'.join(SOURCES) + ')'

include:
    'rules/pesr_preprocessing.rules'

rule all:
    input:
        expand('preprocessing/filtered_vcfs/{source}.{quad}.vcf',
               source=SOURCES, quad=QUADS)

