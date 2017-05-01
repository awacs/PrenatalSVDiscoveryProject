
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
        # expand('integration/rdtest_filtered/pesr/merged.{chrom}.svof', chrom=CHROMS),
        expand('integration/rdtest_filtered/depth/merged.{chrom}.svof', chrom=CHROMS),
        # expand('integration/bca/merged_algs/merged.{chrom}.vcf', chrom=CHROMS),

# TODO: add rules per submodule
rule clean:
    shell:
        "rm $(snakemake --summary | tail -n+2 | cut -f1)"
        
