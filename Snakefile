
configfile: 'config.yaml'

include: 'rules/pesr_preprocessing.rules'
include: 'rules/depth_preprocessing.rules'
include: 'rules/depth_integration.rules'
include: 'rules/pesr_integration.rules'
include: 'rules/RdTest.rules'

PESR_SOURCES = config['pesr_sources']
DEPTH_SOURCES = config['depth_sources']
SOURCES = PESR_SOURCES + DEPTH_SOURCES
RDTEST_SOURCES = PESR_SOURCES + ['depth']
SVTYPES = config['svtypes']
CNV = config['cnv_types']

with open(config['quads']) as qlist:
    QUADS = [q.strip() for q in qlist.readlines()]

with open(config['samples']) as slist:
    SAMPLES = [s.strip() for s in slist.readlines()]

with open(config['chroms']) as clist:
    CHROMS = [c.strip() for c in clist.readlines()]

wildcard_constraints:
    source='(' + '|'.join(SOURCES) + '|depth' + ')',
    sample='(' + '|'.join(SAMPLES) + ')',
    svtype='(' + '|'.join(SVTYPES) + ')',
    chrom='(' + '|'.join(CHROMS) + ')'

rule all:
    input:
        'logs/cnv_depth.done', 'logs/cnv_pesr.done', 'logs/bca.done'

rule cnv_pesr_integration:
    input:
        expand('integration/rdtest_filtered/pesr/merged.{chrom}.svof', chrom=CHROMS),
    output:
        touch('logs/cnv_pesr.done')

rule cnv_depth_integration:
    input:
        expand('integration/rdtest_filtered/depth/merged.{chrom}.svof', chrom=CHROMS),
    output:
        touch('logs/cnv_depth.done')

rule bca_integration:
    input:
        expand('integration/bca/merged_algs/merged.{chrom}.vcf', chrom=CHROMS),
    output:
        touch('logs/bca.done')
# TODO: add rules per submodule
