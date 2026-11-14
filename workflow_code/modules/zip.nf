include { cleanedPrefix } from './utils.nf'

process ZIP {
    tag "Zipping $label"
    label "zip"

    input:
    tuple val(label), path(files)

    output:
    path("${cleanedPrefix()}${label}${params.assay_suffix}*.zip"), emit: zip
    path("versions.txt"), emit: version

    script:
    def is_biom = files.toString().endsWith(".biom")
    def out_name = is_biom
            ? "${cleanedPrefix()}${label}${params.assay_suffix}.biom.zip"
            : "${cleanedPrefix()}${label}${params.assay_suffix}.zip"
    def junk_flag = is_biom ? "-j" : ""
    """
    zip -q ${junk_flag} ${out_name} ${files}

    zip -h | grep "Zip" | sed -E 's/(Zip.+\\)).+/\\1/' >> versions.txt
    """
}