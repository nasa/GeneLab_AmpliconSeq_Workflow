//params.assay_suffix = "_GLAmpSeq"
//params.group  = "groups"
//params.samples_column = "Sample Name"
//params.output_prefix = ""
//params.input_file = "../GeneLab/GLDS-487_amplicon_v1_runsheet.csv"
//params.asv_table = "../workflow_output/Final_Outputs/counts_GLAmpSeq.tsv"
//params.taxonomy = "../workflow_output/Final_Outputs/taxonomy_GLAmpSeq.tsv"


process PLOT_TAXONOMY  {

    tag "Making taxonomy plots..."
    label "visualization"

    input:
        val(meta)
        path(feature_table)
        path(taxonomy_table)
        path(metadata)

    output:
        path("taxonomy_plots/"), emit: output_dir
        path("versions.txt"), emit: version

    script:
        """
        plot_taxonomy.R \\
                  --metadata-table ${metadata} \\
                  --feature-table '${feature_table}' \\
                  --taxonomy-table '${taxonomy_table}' \\
                  --group '${meta.group}' \\
                  --samples-column '${meta.samples}' \\
                  --assay-suffix  '${meta.assay_suffix}' \\
                  --output-prefix  '${params.cleaned_prefix}'
                 
        Rscript -e "VERSIONS=sprintf('dplyr %s\\npurrr %s\\nreadr %s\\nstringr %s\\ntibble %s\\ntidyr %s\\nglue %s\\nggplot2 %s\\n',  \\
                                    packageVersion('dplyr'), \\
                                    packageVersion('purrr'), \\
                                    packageVersion('readr'), \\
                                    packageVersion('stringr'), \\
                                    packageVersion('tibble'), \\
                                    packageVersion('tidyr'), \\
                                    packageVersion('glue'), \\
                                    packageVersion('ggplot2')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"
        """

}