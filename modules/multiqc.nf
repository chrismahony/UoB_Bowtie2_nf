process MULTIQC {
    beforeScript '''\
        module purge; module load bluebear
        module load MultiQC/1.9-foss-2019b-Python-3.7.4
    '''.stripIndent()

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path(multiqc_config)
    path('fastqc/*')

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data
    path "multiqc_report_data", optional: true

    script:
    """
    module purge; module load bluebear
    module load MultiQC/1.9-foss-2019b-Python-3.7.4

    echo "Current directory: \$(pwd)"
    echo "Contents of fastqc directory:"
    ls -R fastqc

    multiqc . \
        --outdir . \
        --filename multiqc_report.html \
        -c $multiqc_config

    # Ensure multiqc_data directory exists
    if [ ! -d "multiqc_data" ]; then
        mkdir multiqc_data
    fi

    # Move multiqc_report_data to multiqc_data if it exists
    if [ -d "multiqc_report_data" ]; then
        mv multiqc_report_data/* multiqc_data/
        rmdir multiqc_report_data
    fi

    # List contents of current directory
    echo "Contents of current directory after MultiQC:"
    ls -la
    """
}