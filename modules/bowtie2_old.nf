process BOWTIE2 {
    tag "$meta.id"
    label 'process_high'
    
    beforeScript '''\
        module purge; module load bluebear
        module load bear-apps/2022b
        module load Bowtie2/2.5.1-GCC-12.2.0
        module load SAMtools/1.17-GCC-12.2.0
    '''.stripIndent() 

    time '10h'
    cpus 20	
    memory { 50.GB * task.attempt }
    
    publishDir "${params.outdir}/bowtie2", mode: 'copy'
        
    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.bam"), emit: aligned_bam
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
    bowtie2 -p $task.cpus -x \$INDEX -1 ${reads[0]} -2 ${reads[1]} $args | \
    samtools view -bS - > ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version | sed 's/.*bowtie2-align-s version //; s/ .*//')
        samtools: \$(samtools --version | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}