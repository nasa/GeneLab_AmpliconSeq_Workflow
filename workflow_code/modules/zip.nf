process ZIP {
    tag "Zipping $label"
    label "zip"

    input:
    tuple val(label), path(files)

    output:
    path("${params.cleaned_prefix}${label}${params.assay_suffix}.zip"), emit: zip
    path("versions.txt"), emit: version

    script:
    """
    zip -q ${params.cleaned_prefix}${label}${params.assay_suffix}.zip ${files}

    zip -h | grep "Zip" | sed -E 's/(Zip.+\\)).+/\\1/' >> versions.txt
    """
}