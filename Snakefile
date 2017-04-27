
configfile: 'config.yaml'

include: 'rules/pesr_preprocessing.rules'
include: 'rules/depth_preprocessing.rules'
include: 'rules/depth_integration.rules'
include: 'rules/pesr_integration.rules'

PESR_SOURCES = config['pesr_sources']
DEPTH_SOURCES = config['depth_sources']
SOURCES = PESR_SOURCES + DEPTH_SOURCES
CNV = config['cnv_types']

with open(config['quads']) as qlist:
    QUADS = [q.strip() for q in qlist.readlines()]

with open(config['samples']) as slist:
    SAMPLES = [s.strip() for s in slist.readlines()]

with open(config['chroms']) as clist:
    CHROMS = [c.strip() for c in clist.readlines()]

wildcard_constraints:
    source='(' + '|'.join(SOURCES) + ')',
    sample='(' + '|'.join(SAMPLES) + ')'

rule all:
    input:
        expand('integration/rdtest_beds/{source}/merged.{chrom}.bed',
               source=PESR_SOURCES, chrom=CHROMS),
        expand('integration/svof/depth.{svtype}.svof', svtype=CNV),
