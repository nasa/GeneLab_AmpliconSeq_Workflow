/**
 * COPY_READS
 * 
 * This process stages raw read files and renames them to standardized naming convention.
 * 
 * For URLs: Nextflow automatically downloads them during staging, then this process
 * renames the downloaded files (1.gz, 2.gz) to the expected format using configured suffixes.
 * 
 * For local files: Nextflow stages them as symlinks, then this process renames
 * them to the standardized format.
 * 
 * Input format:
 * - Paired-end: [sample_id, [R1_file, R2_file], "true"] 
 * - Single-end: [sample_id, [R1_file], "false"]
 * 
 * Output format: 
 * - Paired-end: [sample_id, [sample_R1_raw.fastq.gz, sample_R2_raw.fastq.gz], "true"]
 * - Single-end: [sample_id, [sample_R1_raw.fastq.gz], "false"]
 */
process COPY_READS {
    publishDir "${params.raw_reads_dir}",
        mode: params.publishDir_mode
    tag "${ sample_id }"

    input:
        tuple val(sample_id), path("*.gz"), val(paired)

    output:
        tuple val(sample_id), path("${sample_id}*.gz"), val(paired), emit: raw_reads

    script:
        if ( paired == 'true' ) {
        """
        cp -P 1.gz ${sample_id}${params.raw_R1_suffix}
        cp -P 2.gz ${sample_id}${params.raw_R2_suffix}
        """
        } else {
        """
        cp -P 1.gz ${sample_id}${params.raw_R1_suffix}
        """
        }
}
