process DOWNLOAD_DATABASE {

    tag "Downloading reference database for ${target_region}..."
    storeDir "${params.database_store_path}"
    
    input:
        tuple val(target_region), val(db_name), val(url)
        
    output:
        path("${db_name}"), emit: database
        
    script:        
        """
        wget --user-agent="Mozilla/5.0" ${url} -O ${db_name} || exit 1
        chmod 664 ${db_name}

        ls -la ${db_name}
        """
}
