
configfile: 'config.yaml'

SOURCES = config['pesr_sources']

with open(config['quads']) as qlist:
    QUADS = [q.strip() for q in qlist.readlines()]

wildcard_constraints:
    source='(' + '|'.join(SOURCES) + ')'

subworkflow preprocessing:
    workdir: "preprocessing"

rule all:
    input:
        preprocessing('all.done')

