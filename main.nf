#!/usr/bin/env nextflow




params.reads = ''
params.phred = 33
params.results = './results'
params.pairedEnd = true
params.help = false
params.mpa_db_name = "mpa_v20_m200"
params.mpa_db = "https://github.com/maxibor/metaphlan-nf/releases/download/0.1/mpa_v20_m200.fna.bz2"
params.mpa_pkl = "$baseDir/assets/mpa_v20_m200.pkl"

def helpMessage() {
    log.info"""
     metaphlan-nf: simple metaphlan2 Nextflow pipeline
     Homepage: https://github.com/maxibor/metaphlan-nf
     Author: Maxime Borry <borry@shh.mpg.de>
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run maxibor/metaphlan-nf --reads '/path/to/paired_end_reads_*.{1,2}.fastq.gz' --metaphlandb '/path/to/minimetaphlan2_v2_8GB_201904_UPDATE.tgz'
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)

    Settings:
      --phred                       Specifies the fastq quality encoding (33 | 64). Defaults to ${params.phred}
      --pairedEnd                   Specified if reads are paired-end (true | false). Default = ${params.pairedEnd}

    Options:
      --results                     The output directory where the results will be saved. Defaults to ${params.results}
      --help  --h                   Shows this help page
    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

Channel
    .fromFilePairs( params.reads, size: params.pairedEnd ? 2 : 1 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\n" }
	.set {reads_to_trim}

Channel
    .fromPath (params.mpa_db)
    .ifEmpty { exit 1, "Cannot find Metaphlan Markers fasta file: ${params.mpa_db}\n" }
    .set {mpa_db}

Channel
    .fromPath(params.mpa_pkl)
    .ifEmpty{exit 1, "Cannot find Metaphlan index pickle file: ${params.mpa_pkl}\n"}
    .set {mpa_pkl}

process build_metaphlan_db {
    tag "${params.mpa_db_name}"

    label 'intenso'

    stageInMode 'copy'

    input:
        file(fasta) from mpa_db
    output:
        file("*.bt2") into mpa_db_path
    
    script:
        decomp_fasta = fasta.toString().minus(".bz2")
        """
        bunzip2 $fasta
        bowtie2-build --threads ${task.cpus} $decomp_fasta ${params.mpa_db_name}
        """
}


process AdapterRemoval {
    tag "$name"

    label 'expresso'

    input:
        set val(name), file(reads) from reads_to_trim

    output:
        set val(name), file('*.trimmed.fastq') into trimmed_reads
        set val(name), file("*.settings") into adapter_removal_results

    script:
        out1 = name+".pair1.trimmed.fastq"
        out2 = name+".pair2.trimmed.fastq"
        se_out = name+".trimmed.fastq"
        settings = name+".settings"
        if (params.pairedEnd){
            """
            AdapterRemoval --basename $name --file1 ${reads[0]} --file2 ${reads[1]} --trimns --trimqualities --minquality 20 --minlength 30 --output1 $out1 --output2 $out2 --threads ${task.cpus} --qualitybase ${params.phred} --settings $settings
            """
        } else {
            """
            AdapterRemoval --basename $name --file1 ${reads[0]} --trimns --trimqualities --minquality 20 --minlength 30 --output1 $se_out --threads ${task.cpus} --qualitybase ${params.phred} --settings $settings
            """
        }
            
}

process get_read_count {
    tag "$name"

    label 'ristretto'

    input:
        set val(name), file(ar_settings) from adapter_removal_results
    output:
        set val(name), stdout into nb_reads_ch
    script:
        """
        grep 'Number of retained reads:' $ar_settings | cut -d : -f 2 | tr -d '\040\011\012\015'
        """
}


process metaphlan {
    tag "$name"

    label 'intenso'

    publishDir "${params.results}/metaphlan/$name", mode: 'copy'

    input:
        set val(name), file(reads) from trimmed_reads
        file(bt_index) from mpa_db_path
        file(pkl) from mpa_pkl

    output:
        set val(name), file('*.metaphlan.out') into metaphlan_out
        set val(name), file('*_bowtie.sam') into metaphlan_bowtie_out

    script:
        out = name+".metaphlan.out"
        bt_out = name+"_bowtie.sam"
        tmp_dir = baseDir+"/tmp"
        if (params.pairedEnd){
            """
            mkdir ${params.mpa_db_name}
            mv $pkl ${params.mpa_db_name}/$pkl
            mv *.bt2 ${params.mpa_db_name}
            metaphlan2.py --bowtie2db ${params.mpa_db_name} \\
                          -o $out \\
                          --input_type fastq \\
                          --bowtie2out $bt_out  \\
                          --nproc ${task.cpus} \\
                          ${reads[0]},${reads[1]}
            """    
        } else {
            """
            mkdir ${params.mpa_db_name}
            mv $pkl ${params.mpa_db_name}/$pkl
            mv *.bt2 ${params.mpa_db_name}
            metaphlan2.py --bowtie2db ${params.mpa_db_name} \\
                          -o $out \\
                          --input_type fastq \\
                          --bowtie2out $bt_out  \\
                          --nproc ${task.cpus} \\
                           $reads
            """  
        }
        
}

process metaphlan_parse {
    tag "$name"

    label 'ristretto'

    input:
        set val(name), file(metaphlan_r), val(nb_reads) from metaphlan_out.join(nb_reads_ch)

    output:
        set val(name), file('*.metaphlan_parsed.csv') into metaphlan_parsed

    script:
        out = name+".metaphlan_parsed.csv"
        """
        metaphlan_parse.py -n $nb_reads -o $out $metaphlan_r
        """    
}

process metaphlan_merge {

    label 'ristretto'

    publishDir "${params.results}", mode: 'copy'

    input:
        file(csv_count) from metaphlan_parsed.collect()

    output:
        file('metaphlan_taxon_table.csv') into metaphlan_merged

    script:
        out = "metaphlan_taxon_table.csv"
        """
        merge_metaphlan_res.py -o $out
        """    
}