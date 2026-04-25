process ZIP {
    tag "Zipping $label"
    label "zip"

    input:
    tuple val(label), path(files)

    output:
    path("${params.cleaned_prefix}${label}${params.assay_suffix}*.zip"), emit: zip
    path("versions.txt"), emit: version

    script:
    def is_biom = files.toString().endsWith(".biom")
    def out_name = is_biom
            ? "${params.cleaned_prefix}${label}${params.assay_suffix}.biom.zip"
            : "${params.cleaned_prefix}${label}${params.assay_suffix}.zip"
    def junk_flag = is_biom ? "-j" : ""
    """
    zip -q ${junk_flag} ${out_name} ${files}

    zip -h | grep "Zip" | sed -E 's/(Zip.+\\)).+/\\1/' >> versions.txt
    """
}