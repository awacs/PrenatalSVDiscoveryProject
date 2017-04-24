# SV detection pipeline

* `data/`  
    Raw input and final output data.
* `rules/`  
    Snakemake rules governing pipeline operation.
* `ref/`  
    Reference data and configuration files (e.g. sample lists).
* `preprocessing/`  
    Data standardization and outlier removal.
* `logs/`  
    LSF logs.
* `scripts/`  
    General purpose scripts. Scripts which invoke snakemake data can be found
    in `rules/scripts/`. 

## Input data

### PE/SR calls
PE/SR calls should be placed in the `data/raw_vcfs/` directory and follow the
filename convention `{source}.{quad}.vcf.gz`, where `source` refers to the
source algorithm which generated the calls and `quad` refers to an identifier
of the batch on which the algorithm was run. In the SSC analyses, each
algorithm was run on a per-quad basis.

### Read depth calls

## Pipeline configuration

All variables controlling pipeline operation can be modified in `config.yaml`.

* `quads` : filepath  
    Path to list of quads (or other batch identifier) on which to run the
pipeline.

* `samples` : filepath  
    Path to list of all samples represented when merging the quads.

* `pesr_sources` : list of strings  
    List of all PE/SR algorithms to merge.

* `svtypes` : list of strings  
    List of SV classes to consider when combining PE/SR variants.

## Pipeline modules

* `preprocessing/`  
    Data standardization and outlier removal. Details to follow.
