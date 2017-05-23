
configfile: 'config.yaml'

with open(config['samples']) as slist:
    SAMPLES = [s.strip() for s in slist.readlines()]

CHROMS = [str(x) for x in range(1, 23)] + 'X Y'.split()

rule all:
    input:
        expand('split_counts/{chrom}/{sample}.txt.gz',
               chrom=CHROMS, sample=SAMPLES) 
        
rule count_splits:
    output:
        'split_counts/{chrom}/{sample}.txt.gz'
    params:
        index_dir=config['bam_indexes']
    shell:
        """
        cd {params.index_dir};
        ./scripts/count_splits.py -r {wildcards.chrom} $(s3bam {wildcards.sample}) {output}
        """
