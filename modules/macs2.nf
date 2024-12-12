process MACS2 {
beforeScript '''\
        module purge; module load bluebear
        module load bear-apps/2022a
        module load MACS2/2.2.9.1-foss-2022a
     '''.stripIndent() 

time '10h'
    cpus 20	
    memory { 50.GB * task.attempt }

   publishDir "${params.outdir}/macs2", mode: 'copy'


    input:
    path bams
    val macs2_gsize

    output:
    path "*.{narrowPeak,broadPeak}", emit: peak
    path "*.xls"                   , emit: xls
    path "versions.yml"            , emit: versions
    path "*.gappedPeak"            , optional:true, emit: gapped
    path "*.bed"                   , optional:true, emit: bed
    path "*.bdg"                   , optional:true, emit: bdg

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "macs2_all_samples"
    """
    macs2 callpeak \\
        $args \\
        --gsize $macs2_gsize \\
        --format BAM \\
        --name $prefix \\
        --treatment ${bams.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs2: \$(macs2 --version | sed -e "s/macs2 //g")
    END_VERSIONS
    """
}