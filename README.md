# ASC SV detection and analysis

This repository contains the code used in the generation and analysis of the SV
callset described in Werling et al., 2017. It includes the scripts used to run
each SV detection algorithm, any reference files used in curation, and the
customized pipeline used to combine, filter, and annotate raw SV predictions.

## Overview

* `ref/`  
    Reference files used for callset generation and curation.
   
* `environment.yaml`  
    Anaconda environment used during callset generation.

## SV pipeline

SV detection and filtering was broken into a series of self-contained modules,
each provided here as a Snakemake workflow.
 
* `00_preprocessing/`  
    Conversion of raw algorithm inputs to standard VCF and BED formats.

* `01_algorithm-integration/`  
    Per-algorithm, batch-wide clustering of algorithm predictions.

* `02_evidence-assessment/`  
    Collection of PE, SR, RD, and BAF evidence supporting each putative variant
    locus.

* `03_variant-filtering/`  
    Random-forest filtering of variants based on collected metrics.

* `04_variant-resolution/`  
    Algorithm merging, batch merging, stringent de novo filtering, per-sample
    CNV merging, and resolution of complex variant.

* `05_annotation/`  
    Annotation of resolved variants with genic and noncoding impacts.
