process PACKAGE_PROCESSING_INFO {

    tag "Purging file paths and zipping processing info"

    input:
        path(processing_scripts_dir) 
    output:
        path("${params.cleaned_prefix}processing_info${params.assay_suffix}.zip"), emit: zip

    script:
        """
        cat `which clean-paths.sh` > clean-paths.sh
        chmod +x ./clean-paths.sh
        mkdir processing_info/ && \\
        cp ${processing_scripts_dir}/* processing_info

        echo "Purging file paths"
        find processing_info/ -type f -exec bash ./clean-paths.sh '{}' ${params.baseDir} \\;
        
        # Purge file paths and then zip
        zip -r ${params.cleaned_prefix}processing_info${params.assay_suffix}.zip processing_info/
        """
} 


process GENERATE_README {

    beforeScript "chmod +x ${projectDir}/bin/*"
    tag "Generating README for ${params.OSD_accession}"

    output:
        path("${params.cleaned_prefix}README${params.assay_suffix}.txt"), emit: readme

    script:
        def raw_reads = params.include_raw_data ? "--include-raw-reads" : ""
        def trimmed_reads = params.trim_primers ? "" : "--primers_already_trimmed"
        def single_end = params.single_end ? "--single-end" : ""
        """    
        GL-gen-processed-data-amplicon-readme-updated.py \\
             --output '${params.cleaned_prefix}README${params.assay_suffix}.txt' \\
             --osd-id '${params.OSD_accession}' \\
             --name '${params.name}' \\
             --email '${params.email}' \\
             --protocol-ID '${params.protocol_id}' \\
             --assay_suffix '${params.assay_suffix}' \\
             ${raw_reads} ${trimmed_reads} ${single_end}
        """

}


process VALIDATE_PROCESSING {

    tag "Running automated validation and verification...."

    input:
        val(processed_dir)
        val(target_files)
        path(runsheet)
        path(README)
        path(processing_info)

    output:
        path("${params.GLDS_accession}_${params.cleaned_prefix}amplicon-validation${params.assay_suffix}.log"), emit: log

    script:
        def raw_fastq = params.include_raw_data ? "--include_raw_fastq" : ""
        def primers_flag = params.trim_primers ? "" : "--primers-already-trimmed"
        def single_end_flag = params.single_end ? "--single-ended" : ""
        def used_R1_as_SE_flag = params.used_R1_as_SE ? "--R1-used-as-single-ended-data" : ""
        """
        GL-validate-processed-amplicon-data \\
             --output '${params.GLDS_accession}_${params.cleaned_prefix}amplicon-validation${params.assay_suffix}.log' \\
             --GLDS-ID '${params.GLDS_accession}' \\
             --runsheet '${runsheet}' \\
             --outdir '${processed_dir}' \\
             --V_V_guidelines_link '${params.V_V_guidelines_link}' \\
             --output-prefix '${params.cleaned_prefix}' \\
             --zip_targets '${target_files}' \\
             --assay_suffix '${params.assay_suffix}' \\
             --raw_suffix "${params.assay_suffix}_raw.fastq.gz" \\
             --raw_R1_suffix "${params.assay_suffix}_R1_raw.fastq.gz" \\
             --raw_R2_suffix "${params.assay_suffix}_R2_raw.fastq.gz" \\
             --primer_trimmed_suffix "${params.assay_suffix}_trimmed.fastq.gz" \\
             --primer_trimmed_R1_suffix "${params.assay_suffix}_R1_trimmed.fastq.gz" \\
             --primer_trimmed_R2_suffix "${params.assay_suffix}_R2_trimmed.fastq.gz" \\
             --filtered_suffix "${params.assay_suffix}_filtered.fastq.gz" \\
             --filtered_R1_suffix "${params.assay_suffix}_R1_filtered.fastq.gz" \\
             --filtered_R2_suffix "${params.assay_suffix}_R2_filtered.fastq.gz" \\
             --processing_zip_file '${processing_info}' \\
             --readme '${README}' \\
             ${raw_fastq} ${primers_flag} ${single_end_flag} ${used_R1_as_SE_flag}
        """
}


process GENERATE_CURATION_TABLE {

    beforeScript "chmod +x ${projectDir}/bin/*"
    tag "Generating a file association table for curation..."

    input:
        path(input_table)
        path(fastqc_outputs_dir) // passed so the process has access to raw_multiqc to get read counts from
        path(final_outputs_dir) // passed so the process has access to final outputs to get alpha/beta diversity, taxonomy plots, and DA outputs

    output:
        path("${params.GLDS_accession}_${params.cleaned_prefix}associated-file-names${params.assay_suffix}.tsv"), emit: curation_table

    script:
        def primers_flag = params.trim_primers ? "" : "--primers-already-trimmed"
        def single_end_flag = params.single_end ? "--single-ended" : ""
        def used_R1_as_SE_flag = params.used_R1_as_SE ? "--R1-used-as-single-ended-data" : ""
        def include_raw_multiqc_flag = params.include_raw_multiqc ? "--include-raw-multiqc-in-output" : ""        
        """

        GL-gen-amplicon-file-associations-table-GLfile.py --GL-file ${input_table} \\
                    --output '${params.GLDS_accession}_${params.cleaned_prefix}associated-file-names${params.assay_suffix}.tsv' \\
                    --GLDS-ID  '${params.GLDS_accession}' \\
                    --output-prefix '${params.cleaned_prefix}' \\
                    --assay_suffix '${params.assay_suffix}' \\
                    --raw_suffix '${params.assay_suffix}_raw.fastq.gz' \\
                    --raw_R1_suffix '${params.assay_suffix}_R1_raw.fastq.gz' \\
                    --raw_R2_suffix '${params.assay_suffix}_R2_raw.fastq.gz' \\
                    --primer_trimmed_suffix '${params.assay_suffix}_trimmed.fastq.gz' \\
                    --primer_trimmed_R1_suffix '${params.assay_suffix}_R1_trimmed.fastq.gz' \\
                    --primer_trimmed_R2_suffix '${params.assay_suffix}_R2_trimmed.fastq.gz' \\
                    --filtered_suffix '${params.assay_suffix}_filtered.fastq.gz' \\
                    --filtered_R1_suffix '${params.assay_suffix}_R1_filtered.fastq.gz' \\
                    --filtered_R2_suffix '${params.assay_suffix}_R2_filtered.fastq.gz' \\
                    ${primers_flag} ${single_end_flag} ${used_R1_as_SE_flag} ${include_raw_multiqc_flag}

        """
}


process GENERATE_MD5SUMS {
    
    beforeScript "chmod +x ${projectDir}/bin/*"
    tag "Generating md5sums for the files to be released on OSDR..."
 
    input:
        val(processed_dir)
        val(done_token)

    output:
        path("${params.cleaned_prefix}raw_md5sum${params.assay_suffix}.tsv"), emit: raw_md5sum, optional: true
        path("${params.cleaned_prefix}processed_md5sum${params.assay_suffix}.tsv"), emit: processed_md5sum
    script:
        def raw_md5_flag = params.include_raw_data ? "--generate_raw_md5sums" : ""
        def primers_flag = params.trim_primers ? "" : "--primers_already_trimmed"
        """
        generate_md5sums.py --outdir ${processed_dir} --assay_suffix '${params.assay_suffix}' --output_prefix '${params.cleaned_prefix}' ${raw_md5_flag} ${primers_flag}
        """
}


process GENERATE_PROTOCOL {

    beforeScript "chmod +x ${projectDir}/bin/*"
    tag "Generating your analysis protocol..."

    input:
        path(software_versions)
        val(protocol_id)
        path(rarefaction_depth)
    output:
        path("protocol.txt"), emit: protocol
    script:
        def trim_flag = params.trim_primers ? "--trim_primers" : ""
        """
        generate_protocol.py ${software_versions} ${protocol_id} \\
            --rarefaction_depth_file ${rarefaction_depth} \\
            --target_region ${params.target_region} \\
            --trunc_len "${params.left_trunc},${params.right_trunc}" \\
            --max_ee "${params.left_maxEE},${params.right_maxEE}" \\
            --workflow_version ${workflow.manifest.version} \\
            ${trim_flag}
        """
}
