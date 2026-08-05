//params.assay_suffix = "_GLAmpSeq"
//params.group  = "groups"
//params.samples_column = "Sample Name"
//params.rarefaction_depth = 500
//params.output_prefix = ""
//params.input_file = "../GeneLab/GLDS-487_amplicon_v1_runsheet.csv"
//params.asv_table = "../workflow_output/Final_Outputs/counts_GLAmpSeq.tsv"
//params.taxonomy = "../workflow_output/Final_Outputs/taxonomy_GLAmpSeq.tsv"


process ALPHA_DIVERSITY {

    tag "Running alpha diversity analysis..."
    label "visualization"

    input:
        val(meta)
        path(feature_table)
        path(taxonomy_table)
        path(metadata)

    output:
        path("alpha_diversity/"), emit: output_dir
        path("versions.txt"), emit: version

    script:
        """
        alpha_diversity.R \\
                  --metadata-table ${metadata} \\
                  --feature-table '${feature_table}' \\
                  --taxonomy-table '${taxonomy_table}' \\
                  --group '${meta.group}' \\
                  --samples-column '${meta.samples}' \\
                  --rarefaction-depth ${meta.depth}   \\
                  --assay-suffix  '${meta.assay_suffix}' \\
                  --output-prefix  '${params.cleaned_prefix}' \\
                  --prevalence-cutoff ${meta.prevalence_cutoff} \\
                  --library-cutoff  ${meta.library_cutoff} ${meta.rare}

        Rscript -e "VERSIONS=sprintf('vegan %s\\nphyloseq %s\\nglue %s\\nFSA %s\\nmultcompView %s\\nrstatix %s\\npatchwork %s\\nRColorBrewer %s\\nggplot2 %s\\ndplyr %s\\npurrr %s\\nreadr %s\\nstringr %s\\ntibble %s\\ntidyr %s\\n', \\
                            packageVersion('vegan'), \\
                            packageVersion('phyloseq'), \\
                            packageVersion('glue'), \\
                            packageVersion('FSA'), \\
                            packageVersion('multcompView'), \\
                            packageVersion('rstatix'), \\
                            packageVersion('patchwork'), \\
                            packageVersion('RColorBrewer'), \\
                            packageVersion('ggplot2'), \\
                            packageVersion('dplyr'), \\
                            packageVersion('purrr'), \\
                            packageVersion('readr'), \\
                            packageVersion('stringr'), \\
                            packageVersion('tibble'), \\
                            packageVersion('tidyr')); \\
            write(x=VERSIONS, file='versions.txt', append=TRUE)"
        """

}





process BETA_DIVERSITY {

    tag "Running beta diversity analysis..."
    label "visualization"

    input:
        val(meta)
        path(feature_table)
        path(taxonomy_table)
        path(metadata)


    output:
        path("beta_diversity/"), emit: output_dir
        path("versions.txt"), emit: version

    script:
        """
        beta_diversity.R \\
                  --metadata-table ${metadata} \\
                  --feature-table '${feature_table}' \\
                  --taxonomy-table '${taxonomy_table}' \\
                  --group '${meta.group}' \\
                  --samples-column '${meta.samples}' \\
                  --rarefaction-depth ${meta.depth}   \\
                  --assay-suffix  '${meta.assay_suffix}' \\
                  --output-prefix  '${params.cleaned_prefix}' \\
                  --prevalence-cutoff ${meta.prevalence_cutoff} \\
                  --library-cutoff  ${meta.library_cutoff} ${meta.rare}
        
        Rscript -e "VERSIONS=sprintf('vegan %s\\nphyloseq %s\\nglue %s\\nDESeq2 %s\\nggplot2 %s\\ndplyr %s\\npurrr %s\\nreadr %s\\nstringr %s\\ntibble %s\\ntidyr %s\\nggdendro %s\\nRColorBrewer %s\\nbroom %s\\nvsn %s\\nhexbin %s\\n', \\
                            packageVersion('vegan'), \\
                            packageVersion('phyloseq'), \\
                            packageVersion('glue'), \\
                            packageVersion('DESeq2'), \\
                            packageVersion('ggplot2'), \\
                            packageVersion('dplyr'), \\
                            packageVersion('purrr'), \\
                            packageVersion('readr'), \\
                            packageVersion('stringr'), \\
                            packageVersion('tibble'), \\
                            packageVersion('tidyr'), \\
                            packageVersion('ggdendro'), \\
                            packageVersion('RColorBrewer'), \\
                            packageVersion('broom'), \\
                            packageVersion('vsn'),
                            packageVersion('hexbin')); \\
            write(x=VERSIONS, file='versions.txt', append=TRUE)"
        """

}