#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process CLEAN_FASTQC_PATHS {
    tag "Purging genelab paths from MultiQC zip files in ${params.FastQC_Outputs}"
    input:
        path(FastQC_Outputs_dir)
    output:
        path("${OUT_DIR}"), emit: clean_dir
    script:
        OUT_DIR = "${FastQC_Outputs_dir.baseName}"
        """
        WORKDIR=`pwd`
        mv ${FastQC_Outputs_dir} FastQC_Outputs_dir

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
        val(files_and_dirs) 
    output:
        path("${params.cleaned_prefix}processing_info${params.assay_suffix}.zip"), emit: zip

    script:
        """
        cat `which clean-paths.sh` > clean-paths.sh
        chmod +x ./clean-paths.sh
        mkdir processing_info/ && \\
        cp -r ${files_and_dirs.join(" ")} processing_info/

        echo "Purging file paths"
        find processing_info/ -type f -exec bash ./clean-paths.sh '{}' ${params.baseDir} \\;
        
        # Purge file paths and then zip
        zip -r ${params.cleaned_prefix}processing_info${params.assay_suffix}.zip processing_info/
        """
} 


process GENERATE_README {

    beforeScript "chmod +x ${projectDir}/bin/*"
    tag "Generating README for ${OSD_accession}"
    input:
        tuple val(name), val(email),
              val(OSD_accession), val(protocol_id)

    output:
        path("${params.cleaned_prefix}README${params.assay_suffix}.txt"), emit: readme

    script:
        """    
        GL-gen-processed-data-amplicon-readme-updated.py \\
             --output '${params.cleaned_prefix}README${params.assay_suffix}.txt' \\
             --osd-id '${OSD_accession}' \\
             --name '${name}' \\
             --email '${email}' \\
             --protocol-ID '${protocol_id}' \\
             --assay_suffix '${params.assay_suffix}' \\
             ${params.readme_extra}
        """

}


process VALIDATE_PROCESSING {

    tag "Running automated validation and verification...."

    input:
        tuple val(GLDS_accession), val(V_V_guidelines_link), val(output_prefix),
               val(target_files), val(assay_suffix),
               val(raw_suffix), val(raw_R1_suffix), val(raw_R2_suffix),
               val(primer_trimmed_suffix), val(primer_trimmed_R1_suffix), val(primer_trimmed_R2_suffix),
               val(filtered_suffix), val(filtered_R1_suffix), val(filtered_R2_suffix)
        path(runsheet)
        path(README)
        path(processing_info) 
        path(FastQC_Outputs)
        path(Filtered_Sequence_Data)
        path(Trimmed_Sequence_Data)
        path(Final_Outputs)

    output:
        path("${GLDS_accession}_${params.cleaned_prefix}amplicon-validation${assay_suffix}.log"), emit: log

    script:
        """
        GL-validate-processed-amplicon-data \\
             --output '${GLDS_accession}_${params.cleaned_prefix}amplicon-validation${assay_suffix}.log' \\
             --GLDS-ID '${GLDS_accession}' \\
             --readme '${README}' \\
             --runsheet '${runsheet}' \\
             --V_V_guidelines_link '${V_V_guidelines_link}' \\
             --output-prefix '${params.cleaned_prefix}' \\
             --zip_targets '${target_files}' \\
             --assay_suffix '${assay_suffix}' \\
             --raw_suffix '${raw_suffix}' \\
             --raw_R1_suffix '${raw_R1_suffix}' \\
             --raw_R2_suffix '${raw_R2_suffix}' \\
             --primer_trimmed_suffix '${primer_trimmed_suffix}' \\
             --primer_trimmed_R1_suffix '${primer_trimmed_R1_suffix}' \\
             --primer_trimmed_R2_suffix '${primer_trimmed_R2_suffix}' \\
             --filtered_suffix '${filtered_suffix}' \\
             --filtered_R1_suffix '${filtered_R1_suffix}' \\
             --filtered_R2_suffix '${filtered_R2_suffix}' \\
             --processing_zip_file '${processing_info}' \\
             --fastqc_dir '${FastQC_Outputs}' \\
             --filtered_reads_dir '${Filtered_Sequence_Data}' \\
             --trimmed_reads_dir '${Trimmed_Sequence_Data}' \\
             --final_outputs_dir  '${Final_Outputs}'  ${params.validation_extra}
        """
}


process GENERATE_CURATION_TABLE {

    beforeScript "chmod +x ${projectDir}/bin/*"
    tag "Generating a file association table for curation..."

    input:
        // GeneLab accession and Suffixes
        tuple val(GLDS_accession), val(output_prefix), val(assay_suffix),
               val(raw_suffix), val(raw_R1_suffix), val(raw_R2_suffix),
               val(primer_trimmed_suffix), val(primer_trimmed_R1_suffix), val(primer_trimmed_R2_suffix),
               val(filtered_suffix), val(filtered_R1_suffix), val(filtered_R2_suffix)
        // File labels
        tuple val(processing_zip_file), val(readme)
        // Directory labels as paths
        tuple path(raw_reads_dir), path(filtered_reads_dir),
              path(trimmed_reads_dir), path(final_outputs_dir)
        path(input_table)
        path(FastQC_Outputs)

    output:
        path("${GLDS_accession}_${params.cleaned_prefix}associated-file-names${assay_suffix}.tsv"), emit: curation_table

    script:
        """

        GL-gen-amplicon-file-associations-table-GLfile.py --GL-file ${input_table} \\
                    --output '${GLDS_accession}_${params.cleaned_prefix}associated-file-names${assay_suffix}.tsv' \\
                    --GLDS-ID  '${GLDS_accession}' \\
                    --output-prefix '${params.cleaned_prefix}' \\
                    --assay_suffix '${assay_suffix}' \\
                    --raw_suffix '${raw_suffix}' \\
                    --raw_R1_suffix '${raw_R1_suffix}' \\
                    --raw_R2_suffix '${raw_R2_suffix}' \\
                    --primer_trimmed_suffix '${primer_trimmed_suffix}' \\
                    --primer_trimmed_R1_suffix '${primer_trimmed_R1_suffix}' \\
                    --primer_trimmed_R2_suffix '${primer_trimmed_R2_suffix}' \\
                    --filtered_suffix '${filtered_suffix}' \\
                    --filtered_R1_suffix '${filtered_R1_suffix}' \\
                    --filtered_R2_suffix '${filtered_R2_suffix}' \\
                    --raw_reads_dir '${raw_reads_dir}' \\
                    --fastqc_dir '${FastQC_Outputs}' \\
                    --filtered_reads_dir '${filtered_reads_dir}' \\
                    --trimmed_reads_dir '${trimmed_reads_dir}' \\
                    --final_outputs_dir '${final_outputs_dir}' \\
                    ${params.file_association_extra}

        """
}


process GENERATE_MD5SUMS {
    
    tag "Generating md5sums for the files to be released on OSDR..."
 
    input:
        path(processing_info)
        path(README)
        path(software_versions)
        val(dirs)

    output:
        path("${params.cleaned_prefix}processed_md5sum${params.assay_suffix}.tsv"), emit: md5sum
    script:
        """
        mkdir processing/
        cp -r ${dirs.join(" ")} ${processing_info} ${README} ${software_versions} processing/

        # Generate md5sums
        find -L processing/ -type f \\( \
            ! -name "versions.txt" -a \
            ! -path "*/*-trimmed-counts.tsv" -a \
            ! -path "*/*-cutadapt.log" -a \
            ! -path "*/*_fastqc.zip" -a \
            ! -path "*/*_fastqc.html" -a \
            ! -path "*/alpha_diversity/*.png" -a \
            ! -path "*/taxonomy_plots/*.png" -a \
            ! \\( -path "*/beta_diversity/*.png" -a ! -name "vsd_validation_plot*.png" \\) -a \
            ! -path "*/*_diversity/rarefaction_depth*.txt" -a \
            ! -path "*/differential_abundance/*/*volcano*.png" \
        \\) -exec md5sum '{}' \\; |
        awk -v OFS='\\t' 'BEGIN{OFS="\\t"; printf "File Path\\tFile Name\\tmd5\\n"} \\
                {N=split(\$2,a,"/"); sub(/processing\\//, "", \$2); gsub(/^Final_Outputs\\//, "", \$2); print \$2,a[N],\$1}' \\
        > ${params.cleaned_prefix}processed_md5sum${params.assay_suffix}.tsv
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
        path("protocol.txt")
    script:
        """
        generate_protocol.sh ${software_versions} ${protocol_id} ${rarefaction_depth} > protocol.txt
        """
}
