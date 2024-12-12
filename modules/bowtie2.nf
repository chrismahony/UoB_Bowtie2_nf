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
    def index_name = index.name
    
    // Separate read1 and read2 files
    def read1_files = reads.findAll { it.name.contains('_1') }
    def read2_files = reads.findAll { it.name.contains('_2') }
    def read1_command = read1_files.collect { "-1 $it" }.join(' ')
    def read2_command = read2_files.collect { "-2 $it" }.join(' ')

    """
    mkdir -p bowtie2_index
    cp -L ${index}/* bowtie2_index/
    
    bowtie2 -p $task.cpus -x bowtie2_index/${index_name} \
        $read1_command $read2_command $args | \
    samtools view -bS - > ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version | sed 's/.*bowtie2-align-s version //; s/ .*//')
        samtools: \$(samtools --version | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}