# Algorithm integration

The algorithm integration module combines variant predictions across each set
of algorithms (i.e. PE/SR and RD), then integrates across the two classes of
evidence.

## RD integration
All data based on the standardized calls provided in 
`{workdir}/preprocessing/std_beds/merged.{svtype}.bed.gz`

Per-file column details to be described later as necessary.

* `depth_intersect/merged.{svtype}.{chrom}.bed.gz`  
    Self-intersection of standardized calls, split by chromsoome. 
* `depth_variant_lists/merged.{svtype}.{chrom}.list`  
    List of variant IDs present in standardized calls, split by chromosome.
    Necessary for bedcluster's sparse graph clustering.
* `bedcluster/merged.{svtype}.{chrom}.svof`  
    Clustering of variants, linking variants based on self-intersection hits
* `svof/depth.{svtype}.svof`
    Clustered variants, merged across chromosomes and sorted.
