
configfile: 'config.yaml'

with open(config['samples']) as slist:
    SAMPLES = [s.strip() for s in slist.readlines()]

with open(config['chroms']) as clist:
    CHROMS = [c.strip() for c in clist.readlines()]

with open('failures.list') as clist:
    CS = [c.strip() for c in clist.readlines()]
#CHROMS = [str(x) for x in range(1, 23)] + 'X Y'.split()

# OUTDIRS = 'split_counts disc_counts'.split()
OUTDIRS = ['disc_counts']

rule all:
    input:
        expand('{outdir}/{chrom}/{sample}.txt.gz', outdir=OUTDIRS,
               chrom=CHROMS, sample=SAMPLES)
#        expand('split_counts/{cs}.txt.gz', cs=CS),
#        expand('disc_counts/{cs}.txt.gz', cs=CS)
       

# Force iterative counting over chromosomes to prevent too many simultaneous
# queries to same BAM
def recurse_counts(wildcards):
    chrom = wildcards.chrom

    path = '{outdir}/{{chrom}}/{sample}.txt.gz'
    path = path.format(sample=wildcards.sample, outdir=wildcards.outdir)

    idx = CHROMS.index(chrom)
    chroms = CHROMS[:idx]

    return [path.format(chrom=chrom) for chrom  in chroms]
 
rule count_splits:
    input:
        recurse_counts
    output:
        '{outdir}/{chrom}/{sample}.txt.gz'
    params:
        index_dir=config['bam_indexes']
    wildcard_constraints:
        outdir="split_counts"
    shell:
        """
        fout=$(readlink -f {output});
        count_splits=$(readlink -f scripts/count_splits.py);
        cd {params.index_dir};
        $count_splits -r {wildcards.chrom} $(s3bam {wildcards.sample}) stdout \
          | bgzip -c > $fout
        """

rule count_disc:
    input:
        recurse_counts
    output:
        '{outdir}/{chrom}/{sample}.txt.gz'
    params:
        index_dir=config['bam_indexes']
    wildcard_constraints:
        outdir="disc_counts"
    shell:
        """
        fout=$(readlink -f {output});
        count_disc=$(readlink -f scripts/count_disc.py);
        cd {params.index_dir};
        $count_disc -r {wildcards.chrom} $(s3bam {wildcards.sample}) \
          | sort -k1,1V -k2,2n -k5,5n \
          | bgzip -c > $fout
        """
