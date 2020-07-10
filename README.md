[![Build Status](https://travis-ci.com/maxibor/metaphlan-nf.svg?token=pwT9AgYi4qJY4LTp9WUy&branch=master)](https://travis-ci.com/maxibor/metaphlan-nf)

# metaphlan-nf

[AdapterRemoval](https://github.com/MikkelSchubert/adapterremoval) - [Metaphlan3]([https://bitbucket.org/biobakery/metaphlan2/src/default) (by default, version 3.0](https://github.com/biobakery/MetaPhlAn/tree/3.0)) Nextflow pipeline.

## Dependancies

- [Nextflow](https://www.nextflow.io/)

## Usage

```bash
nextflow run maxibor/metaphlan-nf -profile [docker|conda|singularity] --reads "/path/to/paired_end_reads_*.{1,2}.fastq.gz"
```

Specify the container/environment system you would like to use with `-profile`:  

- conda: `-profile conda`
- docker: `-profile docker`
- singularity: `-profile singularity`

By default, metaphlan will use the database version `mpa_v30_CHOCOPhlAn_201901`

### Input

#### --reads

Use this to specify the location of your input FastQ files. For example:

`--reads 'path/to/data/sample_*_{1,2}.fastq'`

**Please note the following requirements:**

- The path must be enclosed in quotes
- The path must have at least one * wildcard character
- When using the pipeline with paired end data, the path must use {1,2} notation to specify read pairs.


### Output

#### mapdamage

[mapDamage](https://ginolhac.github.io/mapDamage/) output, per sample

#### metaphlan

[metaphlan3](https://github.com/biobakery/MetaPhlAn/tree/3.0) output, per sample

#### pydamage

[pydamage](https://github.com/maxibor/pydamage) output, per sample

#### metaphlan_taxon_table.csv

Taxon count table of all input samples.  
Samples in columns, Taxon in rows

## Help

```bash
$ nextflow run maxibor/metaphlan-nf --help
metaphlan-nf: simple metaphlan Nextflow pipeline
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
  --collapse                    Collapse forward and reverse read, for paired-end reads. Default = false
  --ancient                     Run MapDamage and Pydamage. Default = false
  --mpa_db_name                 Metaphlan database name. Default = mpa_v30_CHOCOPhlAn_201901
  --bt2db                       Directory to store metaphlan database files. Default = /Users/borry/Documents/GitHub/metaphlan-nf/btdb

Options:
  --results                     The output directory where the results will be saved. Defaults to ./results
  --help  --h                   Shows this help page
```
