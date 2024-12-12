process SUBMIT_SBATCH_MACS2 {
    tag "submit_macs2"
    
    input:
    val results_dir
    val genome_size
    
    output:
    val "done", emit: completed
    
    script:
    """
    job_id=\$(sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=macs2_peak_calling
#SBATCH --output=macs2_peak_calling_%j.out
#SBATCH --error=macs2_peak_calling_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Load necessary modules
module purge
module load bear-apps/2022b
module load MACS2/2.2.7.1-foss-2022b

# Change to the results directory
cd ${results_dir}

# Search for BAM files and run MACS2
for bam in \$(find . -name "*.bam"); do
    sample_name=\$(basename \$bam .bam)
    echo "Processing \$sample_name"
    
    macs2 callpeak -t \$bam \\
        -f BAM \\
        -g ${genome_size} \\
        -n \$sample_name \\
        --outdir macs2_output \\
        -q 0.05
done

EOF
)

    job_id=\${job_id##* }
    
    # Wait for the job to complete
    scontrol wait job \$job_id
    
    echo "done" > macs2_job_completed.txt
    """
}