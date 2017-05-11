
configfile: 'config.yaml'

COORDS = {}
QUAD_NAMES = [] 

with open(config['quads']) as qlist:
    QUADS = [q.strip() for q in qlist]

with open(config['bed']) as bedfile:
    for line in bedfile:
        data = line.strip().split()
        coords = ' '.join(data[:3])
        name = data[3]
        samples = data[4].split(',')
        quads = sorted(set([s.split('.')[0] for s in samples]))
        
        COORDS[name] = coords
        
        #for quad in QUADS:
        for quad in quads:
            QUAD_NAMES.append('{0}.{1}'.format(name, quad))

rule all:
    input:
        expand('psl_filtered/{NQ}.psl', NQ=QUAD_NAMES),
        expand('tracks/{NQ}.bed', NQ=QUAD_NAMES),
        'igv/merged.igv',
        'blat.stopped'
        
wildcard_constraints:
    quad="\d{5}"

rule gather_splits:
    """Extract all split reads at locus in each family and reformat to fasta"""
    output:
        'fasta/{name}.{quad}.fa'
    params:
        coords=lambda wildcards: COORDS[wildcards.name]
    shell:
        """
        ./scripts/collect_reads.sh {params.coords} {wildcards.name} {wildcards.quad}
        """

rule start_server:
    """Start blat server"""
    input:
        expand('fasta/{NQ}.fa', NQ=QUAD_NAMES)
    output:
        touch('blat.started')
    shell:
        """
        ./scripts/start_server.sh \
          {config[port]} {config[seqDir]} {config[db]}
        """

rule blat:
    """Blat split reads at each locus in each family"""
    input:
        started='blat.started', 
        fasta='fasta/{name}.{quad}.fa'
    output:
        psl='psl/{name}.{quad}.psl'
    shell:
        """
        gfClient localhost {config[port]} {config[seqDir]} \
          -minScore=0 -minIdentity=0 {input.fasta} {output.psl}
        """

rule stop_server:
    """Stop blat server"""
    input:
        expand('psl/{NQ}.psl', NQ=QUAD_NAMES)
    output:
        touch('blat.stopped')
    shell:
        """
        gfServer stop localhost {config[port]}
        """

rule filter_psl:
    """Exclude blat mappings to other chromosomes"""
    input:
        'psl/{name}.{quad}.psl'
    output:
        'psl_filtered/{name}.{quad}.psl'
    params:
        chrom=lambda wildcards: COORDS[wildcards.name].split()[0]
    shell:
        """
        awk -v chrom={params.chrom} '($14==chrom)' {input} > {output}
        """

def windowed_coords(wildcards):
    """Add viewing window to SV coordinates"""
    coords = COORDS[wildcards.name]
    chrom, start, end = coords.split()
    start = int(start) - config['view_window']
    end = int(end) + config['view_window']
    return 'chr{0}:{1}-{2}'.format(chrom, start, end)

rule make_blat_track:
    """Create color-coded BED track of blat mappings"""
    input:
        'psl_filtered/{name}.{quad}.psl'
    output:
        'tracks/{name}.{quad}.bed'
    params:
        view_window=windowed_coords
    shell:
        """
        echo "browser position {params.view_window}" > {output};
        echo "track name=\"{wildcards.name} {wildcards.quad} SR\" visibility=2 itemRgb=On" >> {output};
        ./scripts/make_blat_track.py {input} | sort -k1,1V -k2,2n >> {output};
        """

def start_window(wildcards):
    """Add viewing window to SV coordinates"""
    coords = COORDS[wildcards.name]
    chrom, start, end = coords.split()
    start = int(start) - config['view_window']
    end = int(start) + config['view_window']
    return 'chr{0}:{1}-{2}'.format(chrom, start, end)

def end_window(wildcards):
    """Add viewing window to SV coordinates"""
    coords = COORDS[wildcards.name]
    chrom, start, end = coords.split()
    start = int(end) - config['view_window']
    end = int(end) + config['view_window']
    return 'chr{0}:{1}-{2}'.format(chrom, start, end)

rule make_igv_batch:
    """Create IGV batch file"""
    input: 
        'tracks/{name}.{quad}.bed'
    output:
        'igv/{name}.{quad}.igv'
    params:
        start_window=start_window,
        end_window=end_window,
    shell:
        """
        echo "new"; \
        echo "genome hg19"; \
        echo "snapshotDirectory {config[snapshots]}"; \
        echo "load $(readlink -f {input})"; \
        echo "goto {params.start_window}"; \
        echo "snapshot {wildcards.name}.{wildcards.quad}.start.png"; \
        echo "goto {params.end_window}"; \
        echo "snapshot {wildcards.name}.{wildcards.quad}.end.png" > {output}
        """

rule merge_igv_batch:
    """Merge IGV batch files"""
    input:
        expand('igv/{NQ}.igv', NQ=QUAD_NAMES)
    output:
        'igv/merged.igv'
    shell:
        """
        cat {input} > {output}
        """
