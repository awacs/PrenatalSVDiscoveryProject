
configfile: 'config.yaml'

include: 'rules/pesr_preprocessing.rules'
include: 'rules/depth_preprocessing.rules'
include: 'rules/depth_integration.rules'
include: 'rules/pesr_alg_integration.rules'
include: 'rules/pesr_cohort_integration.rules'
include: 'rules/RdTest.rules'
include: 'rules/pesr_depth_integration.rules'

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

print(CHROMS)

rule all:
    input:
        'integration/pesr_depth/variants/cohort.22.bed.gz',
        # expand('integration/pesr/bca/cohort.{chrom}.vcf', chrom=CHROMS),

# TODO: add rules per submodule
rule clean:
    shell:
        "rm $(snakemake --summary | tail -n+2 | cut -f1)"
        
