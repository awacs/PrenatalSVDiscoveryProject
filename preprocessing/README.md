# SV Preprocessing

The preprocessing module standardizes VCFs and removes calls specific to
outlying samples. Outliers are determined to be samples with more than 
`(Q3 + 1.5 * IQR)` variants observed, and are calculated on a per-algorithm,
per-svtype basis (e.g. Delly inversions).

* `raw_vcfs/`  
    Raw VCFs.
    '{source}/{source}.{quad}.vcf'

* `std_vcfs/`  
    Standardized VCFs.  
    `{source}.{quad}.vcf`  
        - INFO fields for CHR2, END, STRANDS, SVTYPE, SVLEN, and SOURCE    
        - BND ALTs are converted to VCF specification  
        - Balanced events are separated into stranded breakpoints  

* `outliers/`  
    Tables of svtype-specific outlier samples by caller.  
    `{source}.list`, with the following columns  
        - `sample`: sample ID  
        - `svtype`: type of SV for which the sample is an outlier  
        - `var_count`: number of variants observed in the sample  
        - `cutoff`: Q3 + 1.5 * IQR for that svtype/algorithm  

* `filtered_vcfs/`  
    Filtered VCFS.  
    `{source}.{quad}.vcf`  
    Calls specific to outlier samples or null in all samples are removed.  
