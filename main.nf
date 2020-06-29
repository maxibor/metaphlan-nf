#!/usr/bin/env nextflow
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
      --collapse                    Collapse forward and reverse read, for paired-end reads. Default = ${params.collapse}
      --ancient                     Run DamageProfiler and Pydamage. Default = ${params.ancient}
      --mpa_db_name                 Metaphlan database name. Default = ${params.mpa_db_name}
      --bt2db                       Directory to store metaphlan database files. Default = ${params.bt2db}

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

fasta_db_file = "${params.bt2db}/*.fna.bz2"

// Channel
//     .fromPath (params.mpa_db_tar)
//     .ifEmpty { exit 1, "Cannot find Metaphlan database tar file: ${params.mpa_db_tar}\n" }
//     .set {mpa_db_tar}

process build_metaphlan_db {
    tag "${params.mpa_db_name}"

    label 'intenso'

    output:
        val("${params.mpa_db_name}") into mpa_db_path_wait
    
    script:
        """
        metaphlan --install -x ${params.mpa_db_name} --bowtie2db ${params.bt2db} --nproc ${task.cpus}
        """
}

mpa_db_path_wait.into{mpa_db_path_wait1;  mpa_db_path_wait2}

process AdapterRemoval {
    tag "$name"

    label 'expresso'

    input:
        set val(name), file(reads) from reads_to_trim

    output:
        set val(name), file('*.trimmed.fastq') into trimmed_reads
        set val(name), file("*.settings") into adapter_removal_results, adapter_removal_results_multiqc

    script:
        settings = name+".settings"
        if (params.pairedEnd && !params.collapse){
            out1 = name+".pair1.trimmed.fastq"
            out2 = name+".pair2.trimmed.fastq"
            """
            AdapterRemoval --basename $name --file1 ${reads[0]} --file2 ${reads[1]} --trimns --trimqualities --minquality 20 --minlength 30 --output1 $out1 --output2 $out2 --threads ${task.cpus} --qualitybase ${params.phred} --settings $settings
            """
        } else if (params.pairedEnd && params.collapse) {
            se_out = name+".trimmed.fastq"
            """
            AdapterRemoval --basename $name --file1 ${reads[0]} --file2 ${reads[1]} --trimns --trimqualities --minquality 20 --minlength 30 --collapse --outputcollapsed $se_out --threads ${task.cpus} --qualitybase ${params.phred} --settings $settings
            """
        } 
        else {
            se_out = name+".trimmed.fastq"
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
        val (mpa_db) from mpa_db_path_wait1

    output:
        set val(name), file('*.metaphlan.out') into metaphlan_out
        set val(name), path('*.sam') into mpa_aln
        path("*.bowtie2.log") into bt2_log

    script:
        out = name+".metaphlan.out"
        sam_out = name+".sam"
        bt_out = name+".bowtie2.log"
        btdb = "${params.bt2db}/${params.mpa_db_name}"
        if (params.pairedEnd && !params.collapse){
            """
            bowtie2 --no-unal --very-sensitive -S $sam_out -x $btdb -p ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} 2> $bt_out
            metaphlan $sam_out \\
                      -o $out \\
                      --bowtie2db ${params.bt2db} \\
                      -x ${params.mpa_db_name} \\
                      --input_type sam \\
                      --nproc ${task.cpus} \\
            """    
        } else {
            """
            bowtie2 --no-unal --very-sensitive -S $sam_out -x $btdb -p ${task.cpus} -U $reads --met-file 2> $bt_out
            metaphlan $sam_out \\
                      -o $out \\
                      --bowtie2db ${params.bt2db} \\
                      -x ${params.mpa_db_name} \\
                      --input_type sam \\
                      --nproc ${task.cpus} \\
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

process decompress_fasta {
    label 'expresso'

    input:
        val(wait) from mpa_db_path_wait2
    output:
        file("*.fa") into fasta_ref
    script:
        """
        bunzip2 -c $fasta_db_file > mpa_db_${params.mpa_db_name}.fa
        """
}

fasta_ref.into{fasta_ref_decomp_ch; fasta_ref_damage_profiler}

process sam2bam {
    tag "$name"

    label 'expresso'

    input:
        path(fasta) from fasta_ref_decomp_ch
        tuple val(name), path(sam) from mpa_aln
    output:
        set val(name), path('*.sorted.bam') into mpa_aln_damageprofiler, mpa_aln_pydamage
    script:
        """
        samtools view -F 4 $sam | cut -f 3 | sort | uniq > mapped_refs.txt
        samtools view -H $sam > all_refs.txt
        grep -Ff mapped_refs.txt all_refs.txt > mapped.sam
        samtools view -F 4 $sam >> mapped.sam
        samtools view -h -b -@ ${task.cpus} mapped.sam | samtools sort -@ ${task.cpus} > ${name}.sorted.bam
        """
}


if (params.ancient) {

    // process damageprofiler {
    //     tag "$name"

    //     label 'expresso'

    //     publishDir "${params.results}/damageProfiler/${name}_", mode: 'copy'

    //     input:
    //         set val(name), path(aln) from mpa_aln_damageprofiler
    //         path(fasta) from fasta_ref_damage_profiler
    //     output:
    //         path("*.dmgprof.json") into damageprofiler_result_ch
    //     script:
    //         outfile = name+".dmgprof.json"
    //         maxmem = task.memory.toGiga()
    //         """
    //         samtools index $aln
    //         samtools faidx $fasta
    //         damageprofiler -Xmx${maxmem}g -i $aln -r $fasta -o tmp
    //         mv tmp/${name}.sorted/dmgprof.json $outfile
    //         """
    // }

        process mapdamage {
        tag "$name"

        label 'long_single'

        publishDir "${params.results}/mapdamage/$name", mode: 'copy'

        input:
            set val(name), path(aln) from mpa_aln_damageprofiler
            path(fasta) from fasta_ref_damage_profiler
        output:
            path("mapdamage_out/*") into damageprofiler_result_ch
        script:
            maxmem = task.memory.toGiga()
            """
            samtools index $aln
            samtools faidx $fasta
            mapDamage --merge-reference-sequences -r $fasta -i $aln -d mapdamage_out
            """
    }
    
    process pydamage {
        tag "$name"
        
        label 'intenso'

        errorStrategy 'ignore'

        publishDir "${params.results}/pydamage/$name", mode: 'copy'

        input:
            tuple val(name), path(aln) from mpa_aln_pydamage
        output:
            tuple val(name), path("*.pydamage_results.csv") optional true into pydamage_result_ch
            path "${name}/plots", optional: true
        script:
            output = name
            if (params.pydamage_plot) {
                plot = "--plot"
            } else {
                plot = ""
            }
            """
            samtools index $aln
            pydamage --force -p ${task.cpus} -m ${params.minread} -c ${params.coverage} $plot -o $output $aln
            mv ${name}/pydamage_results.csv ${name}.pydamage_results.csv
            """
    }

} else {
    pydamage_result_ch = Channel.empty()
    damageprofiler_result_ch = Channel.empty()
}

adapter_removal_results_multiqc
    .map {it -> it[1]}
    .set {adapter_removal_results_multiqc}


process multiqc {

    label 'ristretto'
 
    publishDir "${params.results}", mode: 'copy'

    input:
        path('adapterRemoval/*') from adapter_removal_results_multiqc.collect().ifEmpty([])
        path('bowtie2/*') from bt2_log.collect().ifEmpty([])
        // path('damageProfiler/*') from damageprofiler_result_ch.collect().ifEmpty([])
    output:
        path('*multiqc_report.html')
    script:
        """
        multiqc -c ${params.multiqc_config} .
        """
}