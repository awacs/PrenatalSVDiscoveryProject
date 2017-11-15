# 04. Variant resolution
This module merges structural variation calls across all PE/SR and depth algorithms after filtering, then resolves complex variants from the consensus inversion/translocation breakpoints and CNV intervals.
## Process through snakemake

### Module configuration
The configuration file `config.yaml` outlines the module's inputs and parameters, and should be modified accordingly to each specific project. 

* `batches` : Sample/group/batch key.

* `samples` : list of individules to be included

* `chromos` : list of chromosomes to be processed.
This file should be modified according to different reference genome. It is recommended that autosomes and allosomes are prepared differently.

* `pesr_sources` : 	Names of pair end/split read algorithms to be processed

* `depth_sources` :	Names of read depth algorithms to be processed

* `input_vcfs` : folder containing vcf files from `03_variant_filtering`
../03_variant_filtering/filtered_vcfs/

* `pe_counts` : pe_counts/{batch}.pe.sorted.txt.gz 

* `sr_counts` : sr_counts/{batch}.sr.sorted.txt.gz

* `coveragefile`: {batch}.binCov.bed.gz 

* `medianfile` : {batch}.binCov.median

* `famfile` : ../../ref/{batch}.fam

* `cutoffs`: cutoff information trained from the Random Forest

### Input
* `filtered_vcfs`: '../03_variant_filtering/filtered_vcfs/{batch}.{source}.{chrom}.vcf.gz'
The vcf files that were filtered from previous step. `{source}` includes all pesr and depth callers included in the study

## Manual process
### Input
* `filtered_vcfs`: '../03_variant_filtering/filtered_vcfs/{batch}.{source}.{chrom}.vcf.gz'
The vcf files that were filtered from previous step. `{source}` includes all pesr and depth callers included in the study

* `vcflist`: vcflists/pesr/{batch}.{chrom}.list,  vcflists/pesr_depth/{batch}.{chrom}.list
Each list file contains the filtered vcf file from different algorithm on the same chromosome. Lists under `/pesr/` contains pesr calls while `pesr_depth` contains rd calls. Here's an example:
```
../03_variant_filtering/filtered_vcfs/{batch}.delly.20.vcf.gz
../03_variant_filtering/filtered_vcfs/{batch}.lumpy.20.vcf.gz
../03_variant_filtering/filtered_vcfs/{batch}.manta.20.vcf.gz
../03_variant_filtering/filtered_vcfs/{batch}.wham.20.vcf.gz
```

### Process
Follow these steps to process through this step:
1a. Cluster VCFs across algorithms
```
svtk vcfcluster vcflists/pesr/{batch}.10.list stdout -p {batch} -d dist -f frac -x blacklist -z svsize -t svtypes | vcf-sort -c | bgzip -c > vcfcluster/pesr/{batch}.10.vcf.gz
tabix -p vcf vcfcluster/pesr/{batch}.10.vcf.gz
```

1b. Cluster VCFs from pesr and rd algorithms

```
svtk vcfcluster vcflists/pesr_depth/{batch}.10.list stdout -p {batch} -d dist -f frac -x blacklist -z svsize -t svtypes | vcf-sort -c | bgzip -c > vcfcluster/pesr_depth/{batch}.10.vcf.gz
tabix -p vcf vcfcluster/pesr_depth/{batch}.10.vcf.gz
```

2. Link complex SVs:
``` 
svtk resolve -p {batch}_CPX_{chrom} vcfcluster/pesr_depth/{batch}.{chrom}.vcf.gz complex_linking/{batch}.{chrom}.resolved.vcf -u complex_linking/{batch}.{chrom}.unresolved.vcf
```

3. Merge raw vcfs 
```
vcf-concat vcfcluster/pesr_depth/{batch}.*.vcf.gz | bgzip -c > merged_vcfs/{batch}.alg_merged.vcf.gz
tabix -p vcf merged_vcfs/{batch}.alg_merged.vcf.gz
```

### Output



