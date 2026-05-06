process DESEQ  {

    tag "Runnning differential abundance using DESEQ2..."
    label "visualization"

    input:
        val(meta)
        path(feature_table)
        path(taxonomy_table)
        path(metadata)
        path(dummy) // dummy path to ensure dependency between this step and the step that generates this file

    output:
        path("differential_abundance/deseq2/"), emit: output_dir
        path("differential_abundance/${params.cleaned_prefix}contrasts${meta.assay_suffix}.csv"), emit: contrasts_file, optional: true
        path("differential_abundance/${params.cleaned_prefix}SampleTable${meta.assay_suffix}.csv"), emit: sample_table_file, optional: true
        path("versions.txt"), emit: version

    script:
        """
        run_deseq2.R \\
                  --metadata-table ${metadata} \\
                  --feature-table '${feature_table}' \\
                  --taxonomy-table '${taxonomy_table}' \\
                  --group '${meta.group}' \\
                  --samples-column '${meta.samples}' \\
                  --assay-suffix  '${meta.assay_suffix}' \\
                  --output-prefix  '${params.cleaned_prefix}' \\
                  --target-region  '${meta.target_region}' \\
                  --prevalence-cutoff ${meta.prevalence_cutoff} \\
                  --library-cutoff  ${meta.library_cutoff} ${meta.rare}
        
        Rscript -e "VERSIONS=sprintf('DESeq2 %s\\ntaxize %s\\nglue %s\\nphyloseq %s\\nggrepel %s\\nggplot2 %s\\ndplyr %s\\npurrr %s\\nreadr %s\\nstringr %s\\ntibble %s\\ntidyr %s\\n', \\
                            packageVersion('DESeq2'), \\
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
                            packageVersion('tidyr')); \\
            write(x=VERSIONS, file='versions.txt', append=TRUE)"   
        """
}