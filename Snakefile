
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
        #expand('preprocessing/std_beds/merged.{svtype}.bed.gz', svtype=CNV),
        #expand('integration/bedcluster/merged.{chrom}.svof',
        #       chrom=CHROMS)
        #expand('integration/depth_variant_lists/merged.{svtype}.{chrom}.list',
        #       svtype=CNV, chrom=CHROMS)
        #expand('integration/depth_intersect/merged.{svtype}.{chrom}.bed.gz',
        #       svtype=CNV, chrom=CHROMS)
        #expand('integration/rdtest_beds/depth/merged.{chrom}.bed', chrom=CHROMS)
        expand('rdtest/{source}/merged.rdtest_pass.{chrom}.list',
               source=RDTEST_SOURCES, chrom=CHROMS),

# TODO: add rules per submodule
