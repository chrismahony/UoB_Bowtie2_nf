process FASTQC {

     beforeScript '''\
        module purge; module load bluebear
        module load FastQC/0.11.9-Java-11
     '''.stripIndent() 

    time '10h'
    cpus 20	
    memory { 50.GB * task.attempt }


    tag "FASTQC on $sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("*_fastqc.{zip,html}"), emit: fastqc_results

    script:
    """
    module purge; module load bluebear
    module load FastQC/0.11.9-Java-11
    fastqc -q ${reads}

    """
}

