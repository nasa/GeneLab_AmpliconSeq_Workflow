process ZIP {
    tag "Zipping $label"
    label "zip"
    publishDir "${target_dir}", mode: params.publishDir_mode

    input:
    tuple val(label), path(files), val(target_dir)

    output:
    path("${params.cleaned_prefix}${label}${params.assay_suffix}.zip")

    script:
    """
    zip -q ${params.cleaned_prefix}${label}${params.assay_suffix}.zip ${files}

    zip -h | grep "Zip" | sed -E 's/(Zip.+\\)).+/\\1/' >> versions.txt
    """
}