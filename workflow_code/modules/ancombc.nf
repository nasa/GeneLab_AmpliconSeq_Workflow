//params.diff_abund_method = "ancombc2"
//params.assay_suffix = "_GLAmpSeq"
//params.output_prefix = ""
//params.group  = "groups"
//params.samples_column = "Sample Name"
//params.input_file = "../GeneLab/GLDS-487_amplicon_v1_runsheet.csv"
//params.asv_table = "../workflow_output/Final_Outputs/counts_GLAmpSeq.tsv"
//params.taxonomy = "../workflow_output/Final_Outputs/taxonomy_GLAmpSeq.tsv"

process ANCOMBC {

    tag "Running ${method} for differential abundance testing..."
    label "visualization"

    input: 
        val(method)
        val(meta)
        path(feature_table)
        path(taxonomy)
        path(metadata)
        path(dummy) // dummy path to ensure dependency between this step and the step that generates this file

    output:
        path("differential_abundance/${method}/"), emit: output_dir
        path("differential_abundance/${params.cleaned_prefix}contrasts${meta.assay_suffix}.csv"), emit: contrasts_file, optional: true
        path("differential_abundance/${params.cleaned_prefix}SampleTable${meta.assay_suffix}.csv"), emit: sample_table_file, optional: true
        path("versions.txt"), emit: version

    script:      
        """
        if [  ${method} == "ancombc1" ]; then
            script_name='pairwise_ancombc1.R'
        else
            script_name='pairwise_ancombc2.R'
        fi
        
        \${script_name} \\
                  --metadata-table ${metadata} \\
                  --feature-table '${feature_table}' \\
                  --taxonomy-table '${taxonomy}' \\
                  --group '${meta.group}' \\
                  --samples-column '${meta.samples}' \\
                  --assay-suffix  '${meta.assay_suffix}' \\
                  --output-prefix  '${params.cleaned_prefix}' \\
                  --cpus ${task.cpus} \\
                  --target-region  '${meta.target_region}' \\
                  --prevalence-cutoff ${meta.prevalence_cutoff} \\
                  --library-cutoff  ${meta.library_cutoff} ${meta.struc_zero}
                    
        # Capture common packages once
        Rscript -e "VERSIONS=sprintf('ANCOMBC %s\\ntaxize %s\\nglue %s\\nphyloseq %s\\nggrepel %s\\nggplot2 %s\\ndplyr %s\\npurrr %s\\nreadr %s\\nstringr %s\\ntibble %s\\ntidyr %s\\nmia %s\\n', \\
                                    packageVersion('ANCOMBC'), \\
                                    packageVersion('taxize'), \\
                                    packageVersion('glue'), \\
                                    packageVersion('phyloseq'), \\
                                    packageVersion('ggrepel'), \\
                                    packageVersion('ggplot2'), \\
                                    packageVersion('dplyr'), \\
                                    packageVersion('purrr'), \\
                                    packageVersion('readr'), \\
                                    packageVersion('stringr'), \\
                                    packageVersion('tibble'), \\
                                    packageVersion('tidyr'), \\
                                    packageVersion('mia')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"

        # Capture method-specific package
        if [ ${method} == "ancombc1" ]; then
            Rscript -e "VERSIONS=sprintf('DescTools %s\\n', \\
                                        packageVersion('DescTools')); \\
                        write(x=VERSIONS, file='versions.txt', append=TRUE)"
        fi
        """
}