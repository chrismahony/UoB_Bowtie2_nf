include { FASTQC } from '../modules/fastqc'
include { MULTIQC } from '../modules/multiqc'
include { BOWTIE2 } from '../modules/bowtie2'
include { MACS2 } from '../modules/macs2'

workflow RNASEQ {
    // Define channels based on params
    reads_ch = Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .map { name, files -> 
        def sampleName = name.split('_').take(2).join('_')
        return tuple([id: sampleName], files)
    }
    

    genome_dir = file(params.genome_dir)
    gtf_file = file(params.gtf_file)
    multiqc_config_file = file(params.multiqc_config)

    // Run FASTQC on all samples
    FASTQC(reads_ch)
    
    // Run Bowtie2 alignment
    BOWTIE2(reads_ch, file(params.bowtie2_index))
         
    // Debug step
    println "MultiQC config: ${params.multiqc_config ? multiqc_config_file : 'No config provided'}"

    // Collect FastQC outputs
    ch_fastqc_outputs = FASTQC.out.fastqc_results.map { meta, files -> files }.flatten().collect()

    // MultiQC config file
    ch_multiqc_config = Channel.fromPath(params.multiqc_config, checkIfExists: true)

    // Run MultiQC on FastQC outputs
    MULTIQC(ch_multiqc_config, ch_fastqc_outputs)

    // Debug: Print Bowtie2 output directory
    println "Bowtie2 output directory: ${params.outdir}/bowtie2"

    // Create BAM files channel from the results directory
    bam_files_ch = Channel.fromPath("${params.outdir}/bowtie2/*.bam").collect()
    
    // Debug: Print BAM files and genome size
    bam_files_ch.view { "BAM files for MACS2: $it" }
    println "Genome size for MACS2: ${params.genome_size}"

    // Run MACS2
    MACS2(bam_files_ch, params.genome_size)

    // Debug: Print MACS2 output
    MACS2.out.peak.view { "MACS2 peak output: $it" }
    MACS2.out.xls.view { "MACS2 XLS output: $it" }

    
   
 }


