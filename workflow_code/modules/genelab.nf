process CLEAN_MULTIQC_PATHS {
    tag "Purging genelab paths from MultiQC zip files in FastQC_Outputs..."

     input:
        path fastqc_outputs_dir, stageAs: "FastQC_Outputs_dir"
    output:
        path("FastQC_Outputs"), emit: clean_dir
    script:
        OUT_DIR = "FastQC_Outputs"
        """
        WORKDIR=`pwd`

        [ -d ${OUT_DIR}/ ] || mkdir  ${OUT_DIR}/ && \\
        cp -r FastQC_Outputs_dir/*  ${OUT_DIR}/
        
        [ -f ${OUT_DIR}/versions.txt ] && rm -rf ${OUT_DIR}/versions.txt

        cat `which clean-paths.sh` > \${WORKDIR}/clean-paths.sh
        chmod +x \${WORKDIR}/clean-paths.sh

        echo "Purging paths from multiqc outputs"
        cd \${WORKDIR}/${OUT_DIR}/
        echo "Cleaning raw multiqc files with path info"
        unzip ${params.cleaned_prefix}raw_multiqc${params.assay_suffix}_report.zip && rm ${params.cleaned_prefix}raw_multiqc${params.assay_suffix}_report.zip
        cd raw_multiqc_report/raw_multiqc_data/

        # No reason not to just run it on all
        echo "Purging paths in all raw QC files..."
        find . -type f -exec bash \${WORKDIR}/clean-paths.sh '{}' ${params.baseDir} \\;
        cd \${WORKDIR}/${OUT_DIR}/

        echo "Re-zipping up raw multiqc"
        zip -r ${params.cleaned_prefix}raw_multiqc${params.assay_suffix}_report.zip raw_multiqc_report/ && rm -rf raw_multiqc_report/

        echo "Cleaning filtered multiqc files with path info..."
        unzip ${params.cleaned_prefix}filtered_multiqc${params.assay_suffix}_report.zip && rm ${params.cleaned_prefix}filtered_multiqc${params.assay_suffix}_report.zip
        cd filtered_multiqc_report/filtered_multiqc_data/


        # No reason not to just run it on all
        echo "Purging paths in all filtered QC files..."
        find . -type f -exec bash \${WORKDIR}/clean-paths.sh '{}' ${params.baseDir} \\;
        cd \${WORKDIR}/${OUT_DIR}/


        echo "Re-zipping up filtered multiqc..."
        zip -r ${params.cleaned_prefix}filtered_multiqc${params.assay_suffix}_report.zip filtered_multiqc_report/ && rm -rf filtered_multiqc_report/
        cd \${WORKDIR}

        echo "Purging paths from multiqc outputs completed successfully..."

        echo "Done! Paths purged successfully."
        """

}

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
        """    
        GL-gen-processed-data-amplicon-readme-updated.py \\
             --output '${params.cleaned_prefix}README${params.assay_suffix}.txt' \\
             --osd-id '${params.OSD_accession}' \\
             --name '${params.name}' \\
             --email '${params.email}' \\
             --protocol-ID '${params.protocol_id}' \\
             --assay_suffix '${params.assay_suffix}' \\
             ${raw_reads} ${params.readme_extra}
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
        path(cleaned_multiqc_dir) // passed as done token only, script needs to check for the existence of the cleaned multiqc files in the published directory (Post_processing/FastQC_Outputs/)
    output:
        path("${params.GLDS_accession}_${params.cleaned_prefix}amplicon-validation${params.assay_suffix}.log"), emit: log

    script:
        def raw_fastq = params.include_raw_data ? "--include_raw_fastq" : ""
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
             ${raw_fastq}  ${params.validation_extra}
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
                    ${params.file_association_extra}

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
        """
        generate_md5sums.py --outdir ${processed_dir} --assay_suffix '${params.assay_suffix}' --output_prefix '${params.cleaned_prefix}' ${raw_md5_flag}
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
        """
        generate_protocol.sh ${software_versions} ${protocol_id} ${rarefaction_depth} > protocol.txt
        """
}
