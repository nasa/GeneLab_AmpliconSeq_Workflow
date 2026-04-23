/**************************************************************************************** 
*********************  ASV generation using Dada2 ***************************************
****************************************************************************************/

process RUN_DADA2 {

    tag "Running dada2..."

    input:
        tuple path(sample_IDs_file), val(isPaired)
        path(reads)
        path(trimmed_read_counts) // or dummy NO_FILE if no trimming
        path(database)
    output:
        path("*_filtered.fastq.gz"),                                                    emit: reads
        path("${params.cleaned_prefix}filtered-read-counts${params.assay_suffix}.tsv"), emit: filtered_count
        path("${params.cleaned_prefix}taxonomy${params.assay_suffix}.tsv"),             emit: taxonomy
        path("${params.cleaned_prefix}taxonomy-and-counts${params.assay_suffix}.biom"), emit: biom
        path("${params.cleaned_prefix}ASVs${params.assay_suffix}.fasta"),               emit: fasta
        path("${params.cleaned_prefix}read-count-tracking${params.assay_suffix}.tsv"),  emit: read_count
        path("${params.cleaned_prefix}counts${params.assay_suffix}.tsv"),               emit: counts
        path("${params.cleaned_prefix}taxonomy-and-counts${params.assay_suffix}.tsv"),  emit: taxonomy_count
        path("versions.txt"),                                                           emit: version
    script:
        def input_dir     = params.trim_primers ? "Trimmed_Sequence_Data/" : "raw_reads/"
        def trimmed_count = params.trim_primers ? "${trimmed_read_counts}" : ""
        def R1_suffix     = params.trim_primers ? "${params.assay_suffix}_R1_trimmed.fastq.gz" : "${params.assay_suffix}_R1_raw.fastq.gz"
        def R2_suffix     = (isPaired ? (params.trim_primers ? "${params.assay_suffix}_R2_trimmed.fastq.gz" : "${params.assay_suffix}_R2_raw.fastq.gz") : "")
        def filtered_R1   = "${params.assay_suffix}_R1_filtered.fastq.gz"
        def filtered_R2   = "${params.assay_suffix}_R2_filtered.fastq.gz"
        """
        mkdir -p ${input_dir}
        if [ ${isPaired} == true ]; then
            mv *${R1_suffix} *${R2_suffix} ${trimmed_count} ${input_dir}
        else
            mv *${R1_suffix} ${trimmed_count} ${input_dir}
        fi

        if [ ${isPaired} == true ]; then
            Illumina-PE-R-processing.R \\
                "${params.left_trunc}" \\
                "${params.right_trunc}" \\
                "${params.left_maxEE}" \\
                "${params.right_maxEE}" \\
                "${params.trim_primers ? 'TRUE' : 'FALSE'}" \\
                ${sample_IDs_file} \\
                "${input_dir}" \\
                "${R1_suffix}" \\
                "${R2_suffix}" \\
                "${filtered_R1}" \\
                "${filtered_R2}" \\
                "${params.cleaned_prefix}" \\
                "${params.target_region}" \\
                "${params.concatenate_reads_only}" \\
                "${params.assay_suffix}" \\
                "${database.name}"
        else
            Illumina-SE-R-processing.R \\
                "${params.left_trunc}" \\
                "${params.left_maxEE}" \\
                "${params.trim_primers ? 'TRUE' : 'FALSE'}" \\
                ${sample_IDs_file} \\
                "${input_dir}" \\
                "${R1_suffix}" \\
                "${filtered_R1}" \\
                "${params.cleaned_prefix}" \\
                "${params.target_region}" \\
                "${params.assay_suffix}" \\
                "${database.name}"
        fi

        (head -n 1 ${params.cleaned_prefix}taxonomy-and-counts${params.assay_suffix}.tsv; \\
            awk 'NR>1{print}' ${params.cleaned_prefix}taxonomy-and-counts${params.assay_suffix}.tsv | sort -V -k1) \\
            > temp_tax_cont.tsv && mv temp_tax_cont.tsv ${params.cleaned_prefix}taxonomy-and-counts${params.assay_suffix}.tsv

        R --vanilla --version | grep "R version" > versions.txt
        get_R_package_version.R
        """
}