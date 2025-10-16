#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process DOWNLOAD_DATABASE {

    tag "Downloading reference database for ${target_region}..."
    
    input:
        val(target_region)
        
    output:
        path("*.RData"), emit: database
        
    script:
        def db_config = [
            "16S": ["SILVA_SSU_r138_2_2024.RData", "https://figshare.com/ndownloader/files/52846199"],
            "ITS": ["UNITE_v2024_April2024.RData", "https://figshare.com/ndownloader/files/52846346"],
            "18S": ["PR2_v4_13_March2021.RData", "https://figshare.com/ndownloader/files/46241917"]
        ]
        
        def db_name = db_config[target_region][0]
        def url = db_config[target_region][1]
        
        """
        wget --user-agent="Mozilla/5.0" ${url} -O ${db_name} || exit 1
        ls -la ${db_name}
        """
}
