process FEATURECOUNTS {
    beforeScript '''\
        module purge; module load bluebear
        module load bear-apps/2022a
        module load Subread/2.0.4-GCC-11.3.0
    '''.stripIndent()

    time '20h'
    cpus 40
    memory { 100.GB * task.attempt }

    tag "featureCounts on all samples"
    publishDir "${params.outdir}/featurecounts", mode: 'copy'

    input:
    path bams
    path gtf

    output:
    path "*_featurecounts.txt", emit: counts
    path "versions.yml"       , emit: versions

    script:
    """
    featureCounts -a $gtf -o all_samples_featurecounts.txt -T $task.cpus $bams

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        featureCounts: \$(featureCounts -v 2>&1 | sed -n '/\/featureCounts v/p' | sed 's/^.*featureCounts //' | sed 's/ .*\$//')
    END_VERSIONS
    """
}