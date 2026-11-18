def cleanedPrefix() {
    def prefix = params.output_prefix ?: ""
    return (prefix && !prefix.endsWith("_") && !prefix.endsWith("-")) ? prefix + "_" : prefix
}

def assaySuffix() {
    def suffix = params.assay_suffix ?: ""
    def regionPrefix = "_${params.target_region}"
    if (!params.prepend_target_region || suffix == regionPrefix || suffix.startsWith("${regionPrefix}_")) {
        return suffix
    }
    return "${regionPrefix}${suffix}"
}

process SOFTWARE_VERSIONS {

    tag "Writing out software versions..."
    label "fastqc"   //unix environment

    input:
        path(software_versions)

    output:
        path({ "software_versions${assaySuffix()}.txt" }), emit: software_versions_txt

    script:
        """
        # Delete white spaces and write out unique software versions 
        grep -v "^\$" ${software_versions} | sort -u > software_versions${assaySuffix()}.txt
        """
}
