[![Build Status](https://travis-ci.com/maxibor/metaphlan-nf.svg?token=pwT9AgYi4qJY4LTp9WUy&branch=master)](https://travis-ci.com/maxibor/metaphlan-nf)

# metaphlan-nf

Simple [AdapterRemoval](https://github.com/MikkelSchubert/adapterremoval) - [Metaphlan2](https://bitbucket.org/biobakery/metaphlan2/src/default) (by default, version 2.7.5) Nextflow pipeline.

## Dependancies

- [conda](https://conda.io/en/latest/) 
- [Nextflow](https://www.nextflow.io/) : `conda install -c bioconda nextflow`

## Usage

```
nextflow run maxibor/metaphlan-nf --reads "/path/to/paired_end_reads_*.{1,2}.fastq.gz"
```

### Input

#### --reads

Use this to specify the location of your input FastQ files. For example:

`--reads 'path/to/data/sample_*_{1,2}.fastq'`

**Please note the following requirements:**

- The path must be enclosed in quotes
- The path must have at least one * wildcard character
- When using the pipeline with paired end data, the path must use {1,2} notation to specify read pairs.


### Output

#### metaphlan_taxon_table.csv

Taxon count table of all input samples.  
Samples in columns, Taxon in rows

## Help

```
$ nextflow run maxibor/metaphlan-nf --help
metaphlan-nf: simple metaphlan2 Nextflow pipeline
 Homepage: https://github.com/maxibor/metaphlan-nf
 Author: Maxime Borry <borry@shh.mpg.de>
=========================================
Usage:
The typical command for running the pipeline is as follows:
nextflow run maxibor/metaphlan-nf --reads '/path/to/paired_end_reads_*.{1,2}.fastq.gz'
Mandatory arguments:
  --reads                       Path to input data (must be surrounded with quotes)

Settings:
  --phred                       Specifies the fastq quality encoding (33 | 64). Defaults to 33
  --pairedEnd                   Specified if reads are paired-end (true | false). Default = true

Options:
  --results                     The output directory where the results will be saved. Defaults to ./results
  --help  --h                   Shows this help page
```
