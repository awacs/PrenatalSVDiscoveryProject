# SV detection pipeline

Talkowski Lab structural variant detection pipeline. Documentation in progress.

## Dependencies

### SVtools
The pipeline requires the `svtools` Python package, which is currently
available only on github.

```
$ git clone git@github.com:talkowski-lab/svtools.git
$ cd svtools
$ pip install -e .
```

### Snakemake
The pipeline is built with `snakemake`, which can be installed through `pip` or
`conda` with one of the two following commands.

```
$ pip install snakemake
$ conda install -c bioconda snakemake
```

`snakemake` is an excellent and recommended tool for building bioinformatics
pipelines, but a comprehensive understanding of the tool is not necessary to
run this pipeline. If interested in learning more, extended `snakemake`
documentation can be found on their [Read the Docs
page](https://snakemake.readthedocs.io/en/stable/). A
[tutorial](https://snakemake.bitbucket.io/snakemake-tutorial.html) and
[demonstrative slides](http://slides.com/johanneskoester/deck-1#/) are also
available. 

## Installation and Usage
As a `snakemake` workflow, the pipeline is intended to be cloned for each
project, e.g.

```
$ git clone git@github.com:talkowski-lab/sv-pipeline.git MySVDiscoveryProject
```

After cloning the pipeline, edit `config.yaml` to update the configuration as
necessary for the project, then link or copy raw data into the `data/` or
`ref/` directories. (More detail below or to come. For current testing
purposes, symlink the `data/` and `ref/` directories in
`/data/talkowski/Samples/SFARI/deep_sv/asc_540/sv-pipeline-devel/`). 

Optionally, run a dry run of `snakemake` to test the configuration, then run
the pipeline with `snakemake`.

```
$ vim config.yaml
$ ln -s ${raw_SV_calls} data/
$ cp ${reference_data} ref/
$ snakemake -np
$ snakemake
```

The pipeline can remove all files it has generated without affecting
configuration or data files.

```
$ snakemake clean
```

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

## Directory structure
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
