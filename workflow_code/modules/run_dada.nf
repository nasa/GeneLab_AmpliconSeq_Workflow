/**************************************************************************************** 
*********************  ASV generation using Dada2 ***************************************
****************************************************************************************/

include { cleanedPrefix; assaySuffix } from './utils.nf'

process RUN_DADA2 {

    tag "Running dada2..."

    input:
        tuple path(sample_IDs_file), val(isPaired)
        path(reads)
        path(trimmed_read_counts) // Empty when primers are already trimmed.
        path(database)
    output:
        path("*_filtered.fastq.gz"),                                               emit: reads
        path({ "${cleanedPrefix()}filtered-read-counts${assaySuffix()}.tsv" }), emit: filtered_count
        path("${cleanedPrefix()}taxonomy${assaySuffix()}.tsv"),             emit: taxonomy
        path({ "${cleanedPrefix()}taxonomy-and-counts${assaySuffix()}.biom" }), emit: biom
        path("${cleanedPrefix()}ASVs${assaySuffix()}.fasta"),               emit: fasta
        path("${cleanedPrefix()}read-count-tracking${assaySuffix()}.tsv"),  emit: read_count
        path("${cleanedPrefix()}counts${assaySuffix()}.tsv"),               emit: counts
        path("${cleanedPrefix()}taxonomy-and-counts${assaySuffix()}.tsv"),  emit: taxonomy_count
        path("versions.txt"),                                                      emit: version
    script:
        def trim_primers = params.trim_primers
        def input_dir     = trim_primers ? "Trimmed_Sequence_Data/" : "raw_reads/"
        def trimmed_count = trim_primers ? "${trimmed_read_counts}" : ""
        def trimming_reverse = params.force_single_end == "R2"
        def single_end_trunc = trimming_reverse ? params.right_trunc : params.left_trunc
        def single_end_maxEE = trimming_reverse ? params.right_maxEE : params.left_maxEE
        
        def starting_reads_suffix = trim_primers ? "trimmed" : "raw"
        def R1_suffix     = "${assaySuffix()}_R1_${starting_reads_suffix}.fastq.gz"
        def R2_suffix     = "${assaySuffix()}_R2_${starting_reads_suffix}.fastq.gz"
        def single_end_suffix = isPaired == "false" ? (
            params.force_single_end == "R1" ? R1_suffix :
            params.force_single_end == "R2" ? R2_suffix :
            "${assaySuffix()}_${starting_reads_suffix}.fastq.gz"
        ) : ""
        
        def filtered_R1_suffix   = "${assaySuffix()}_R1_filtered.fastq.gz"
        def filtered_R2_suffix   = "${assaySuffix()}_R2_filtered.fastq.gz"
        def filtered_single_end_suffix = isPaired == "false" ? (
            params.force_single_end == "R1" ? filtered_R1_suffix :
            params.force_single_end == "R2" ? filtered_R2_suffix :
            "${assaySuffix()}_filtered.fastq.gz"
        ) : ""
        """
        mkdir -p ${input_dir}
        if [ ${isPaired} == true ]; then
            mv *${R1_suffix} *${R2_suffix} ${trimmed_count} ${input_dir}
        else
            mv *${single_end_suffix} ${trimmed_count} ${input_dir}
        fi

        if [ ${isPaired} == true ]; then
            Illumina-PE-R-processing.R \\
                "${params.left_trunc}" \\
                "${params.right_trunc}" \\
                "${params.left_maxEE}" \\
                "${params.right_maxEE}" \\
                "${trim_primers ? 'TRUE' : 'FALSE'}" \\
                ${sample_IDs_file} \\
                "${input_dir}" \\
                "${R1_suffix}" \\
                "${R2_suffix}" \\
                "${filtered_R1_suffix}" \\
                "${filtered_R2_suffix}" \\
                "${cleanedPrefix()}" \\
                "${params.target_region}" \\
                "${params.concatenate_reads_only ? 'TRUE' : 'FALSE'}" \\
                "${assaySuffix()}" \\
                "${database.name}"
        else
            Illumina-SE-R-processing.R \\
                "${single_end_trunc}" \\
                "${single_end_maxEE}" \\
                "${trim_primers ? 'TRUE' : 'FALSE'}" \\
                ${sample_IDs_file} \\
                "${input_dir}" \\
                "${single_end_suffix}" \\
                "${filtered_single_end_suffix}" \\
                "${cleanedPrefix()}" \\
                "${params.target_region}" \\
                "${assaySuffix()}" \\
                "${database.name}"
        fi

        (head -n 1 ${cleanedPrefix()}taxonomy-and-counts${assaySuffix()}.tsv; \\
            awk 'NR>1{print}' ${cleanedPrefix()}taxonomy-and-counts${assaySuffix()}.tsv | sort -V -k1) \\
            > temp_tax_cont.tsv && mv temp_tax_cont.tsv ${cleanedPrefix()}taxonomy-and-counts${assaySuffix()}.tsv

        R --vanilla --version | grep "R version" > versions.txt
        get_R_package_version.R
        """
}
