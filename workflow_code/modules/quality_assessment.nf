/**************************************************************************************** 
*********************  Sequence quality assessment and control processes ****************
****************************************************************************************/

// A 3-column (single-end) or 4-column (paired-end) file
//params.csv_file = "${projectDir}/PE_file.csv" 
//params.prefix = "raw"
//params.multiqc_config = "${projectDir}/config/multiqc.config"

// FastQC performed on reads
process FASTQC {

    tag "Running fastqc on ${sample_id}..."
    beforeScript "chmod +x ${projectDir}/bin/*"
    label "fastqc"

    input:
        tuple val(sample_id), path(reads), val(isPaired)
    output:
        tuple path("*.html"), path("*.zip"), emit: fastqc
        path("versions.txt"), emit: version
    script:
        // Calculate memory per thread (100MB minimum, 10000MB maximum), adapted from https://github.com/nf-core/modules/blob/master/modules/nf-core/fastqc/main.nf
        def memory_in_mb = MemoryUnit.of("${task.memory}").toUnit('MB') / task.cpus
        def fastqc_memory = memory_in_mb > 10000 ? 10000 : (memory_in_mb < 100 ? 100 : memory_in_mb)
        """
        fastqc -o . \\
        -t ${task.cpus} \\
        --memory ${fastqc_memory} \\
        ${reads}

        fastqc --version > versions.txt
        """
}

process MULTIQC {

    tag "Running multiqc on the ${prefix} files..."
    beforeScript "chmod +x ${projectDir}/bin/*"
    
    input:
        val(prefix)   
        path(multiqc_config)
        path(files)
    output:
        path("${prefix}_multiqc${params.assay_suffix}_data.zip"), emit: zipped_data
        path("${prefix}_multiqc${params.assay_suffix}.html"), emit: html
        path("${prefix}_multiqc${params.assay_suffix}_data"), emit: data
        path("versions.txt"), emit: version
    script:
        """
        multiqc \\
            --force \\
            --interactive \\
            -o . \\
            -n ${prefix}_multiqc${params.assay_suffix} \\
            --cl-config 'max_table_rows: 99999999' \\
            --config ${multiqc_config} \\
            ${files} > /dev/null 2>&1

        # Clean paths and create zip
        clean_multiqc_paths.py ${ prefix }_multiqc${ params.assay_suffix }_data .

        multiqc --version > versions.txt
        """
}
