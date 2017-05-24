
configfile: 'config.yaml'

with open(config['samples']) as slist:
    SAMPLES = [s.strip() for s in slist.readlines()]

with open(config['chroms']) as clist:
    CHROMS = [c.strip() for c in clist.readlines()]

#CHROMS = [str(x) for x in range(1, 23)] + 'X Y'.split()

rule all:
    input:
        expand('split_counts/{chrom}/{sample}.txt.gz',
               chrom=CHROMS, sample=SAMPLES) 
       

# Force iterative counting over chromosomes to prevent too many simultaneous
# queries to same BAM
def recurse_counts(wildcards):
    chrom = wildcards.chrom
    sample = wildcards.sample

    path = 'split_counts/{{chrom}}/{sample}.txt.gz'.format(sample=sample)
    idx = CHROMS.index(chrom)
    chroms = CHROMS[:idx]

    return [path.format(chrom=chrom) for chrom  in chroms]
 
rule count_splits:
    input:
        recurse_counts        
    output:
        'split_counts/{chrom}/{sample}.txt.gz'
    params:
        index_dir=config['bam_indexes']
    shell:
        """
        fout=$(readlink -f {output});
        count_splits=$(readlink -f scripts/count_splits.py);
        cd {params.index_dir};
        $count_splits -z -r {wildcards.chrom} $(s3bam {wildcards.sample}) stdin \
          | bgzip -c > $fout
        """
