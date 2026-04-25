process SOFTWARE_VERSIONS {

    tag "Writing out software versions..."
    label "fastqc"   //unix environment

    input:
        path(software_versions)

    output:
        path("software_versions${params.assay_suffix}.txt"), emit: software_versions_txt

    script:
        """
        # Delete white spaces and write out unique software versions 
        grep -v "^\$" ${software_versions} | sort -u > software_versions${params.assay_suffix}.txt 
        """
}
