# Preprocessing

This workflow standardizes SV calls across a series of algorithms to a
consistent VCF format for PE/SR algorithms and a consistent BED format for
depth algorithms.

This workflow contains many elements which are specific to the SSC analyses for
which the pipeline was initially developed, and is not intended to be
generalizable to all use cases. However, it is presented here for a)
reproducibility of the SSC analyses and b) as an example of how SV calls may be
standardized to the format required by the remaining steps of the pipeline.

## Input data

### PE/SR calls
PE/SR calls should be placed in the `data/raw_vcfs/` directory and follow the
filename convention `{source}.{group}.vcf.gz`, where `source` refers to the
source algorithm which generated the calls and `group` refers to an identifier
of the group of samples on which the algorithm was run. In the SSC analyses,
each algorithm was run on a per-quad basis.

### Read depth calls
To be completed.

## Output data

### PE/SR preprocessing
The PE/SR preprocessing module standardizes VCFs and removes calls specific to
outlying samples. Outliers are determined to be samples with more than 
`(Q3 + 1.5 * IQR)` variants observed, and are calculated on a per-algorithm,
per-svtype basis (e.g. Delly inversions).

* `std_vcfs/{source}.{quad}.vcf`  
    Standardized VCFs.  
        - INFO fields for CHR2, END, STRANDS, SVTYPE, SVLEN, and SOURCE    
        - BND ALTs are converted to VCF specification  
        - Balanced events are separated into stranded breakpoints  

* `outliers/{source}.list`  
    Tables of svtype-specific outlier samples by caller.  
        - `sample`: sample ID  
        - `svtype`: type of SV for which the sample is an outlier  
        - `var_count`: number of variants observed in the sample  
        - `cutoff`: Q3 + 1.5 * IQR for that svtype/algorithm  

* `filtered_vcfs/{source}.{quad}.vcf.gz`  
    VCFs with calls specific to outlier samples or null in all samples removed.
    Sorted and tabix-indexed to prepare for clustering.


### RD preprocessing
The RD preprocessing module merges all calls made by any depth algorithm in a 
sample to create a single bed per sample, annotated with the algorithms which
contributed to each merged call.

* `std_beds/`
    * `{source}/{sample}.{svtype}.raw.bed`  
        Raw calls. Columns are standardized to chrom, start, end, svsize
    * `{source}/{sample}.{svtype}.merged.bed`  
        Raw beds subjected to bedtools merged.
    * `merged_algs/{sample}.{svtype}.raw.bed`  
        Sorted concatenation of raw calls.
    * `merged_algs/{sample}.{svtype}.merged.bed`  
        Per-sample merging of raw calls across all algorithms.
    * `merged_algs/{sample}.{svtype}.cov.bed`  
        Final per-sample standardized bed. Each merged call is queried for
        coverage by each individual algorithm's callset in order to determine 
        the algorithms which predicted the variant.  Columns are chrom, start, 
        end, sample svsize, cnmops coverage, cnvnator coverage, and genomestrip 
        coverage. 
    * `{batch}.{svtype}.bed.gz`  
        Concatenation of standardized `*.cov.bed` from all samples in a batch.