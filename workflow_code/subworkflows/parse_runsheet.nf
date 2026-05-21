def colorCodes = [
    c_back_bright_red: "\u001b[41;1m",
    c_bright_green:    "\u001b[32;1m",
    c_reset:           "\033[0m"
]

def get_runsheet_paths(LinkedHashMap row) {
    def meta = [:]
    meta.id         = row["Sample Name"]
    meta.region     = row["Parameter Value[Library Selection]"]
    meta.paired_end = row.paired_end.toBoolean()
    meta.f_primer   = row["F_Primer"]
    meta.r_primer   = meta.paired_end ? row["R_Primer"] : "null" 
    meta.groups     = row.containsKey("groups") ? row.groups : "null"

    // Extract factors
    meta.factors = row.findAll { key, value ->
        key.startsWith("Factor Value[") && key.endsWith("]")
    }.collectEntries { key, value ->
        [(key[13..-2]): value]
    }

    def raw_reads = []
    def array = []
    raw_reads.add(file(row.read1_path))
    if (meta.paired_end) {
        raw_reads.add(file(row.read2_path))
    }
    array = [meta, raw_reads]

    return array
}

def mutate_to_single_end(it, read_choice) {
    def new_meta = it[0].clone()
    new_meta.paired_end = false
    def chosen_read = read_choice == 'R1' ? it[1][0] : it[1][1]
    return [new_meta, [chosen_read]]
}

process TRUNCATE_RUNSHEET {

    input:
        path(runsheet)
        val(limit)

    output:
        path "${runsheet.baseName}_truncated.csv", emit: truncated_runsheet

    script:
    """
    head -n 1 ${runsheet} > ${runsheet.baseName}_truncated.csv
    if [ ${limit} -gt 0 ]; then
        tail -n +2 ${runsheet} | head -n ${limit} >> ${runsheet.baseName}_truncated.csv
    else
        tail -n +2 ${runsheet} >> ${runsheet.baseName}_truncated.csv
    fi
    """
}

workflow PARSE_RUNSHEET {
    take:
        runsheet_path

    main:
        sample_limit = params.limit_samples_to ? params.limit_samples_to : -1

        if (sample_limit > 0) {
            TRUNCATE_RUNSHEET(runsheet_path, sample_limit)
            ch_runsheet = TRUNCATE_RUNSHEET.out.truncated_runsheet
        } else {
            ch_runsheet = runsheet_path
        }

        // Parse runsheet into meta tuples
        ch_samples = ch_runsheet
            | splitCsv(header: true)
            | map { row -> get_runsheet_paths(row) }
            | map { it -> params.force_single_end ? mutate_to_single_end(it, params.force_single_end) : it }

        // Validate metadata consistency across samples
        ch_samples
            .map { meta, reads -> [meta.paired_end, meta.f_primer, meta.r_primer] }
            .unique()
            .count()
            .subscribe { count ->
                if (count > 1) {
                    log.error "${colorCodes.c_back_bright_red}ERROR: Inconsistent metadata across samples.${colorCodes.c_reset}"
                    exit 1
                } else {
                    println "${colorCodes.c_bright_green}Metadata consistency check passed.${colorCodes.c_reset}"
                }
            }

        // Check all read files are unique
        ch_samples
            .flatMap { meta, reads -> reads }
            .collect()
            .map { all_reads ->
                def total_count  = all_reads.size()
                def unique_count = all_reads.toSet().size()
                if (unique_count != total_count) {
                    error "${colorCodes.c_back_bright_red}ERROR: Duplicate read files detected.${colorCodes.c_reset}"
                } else {
                    println "${colorCodes.c_bright_green}All ${unique_count} read files are unique.${colorCodes.c_reset}"
                }
            }

        // Print autodetected metadata for first sample
        ch_samples.take(1) | view { meta, reads ->
            """${colorCodes.c_bright_green}Autodetected Processing Metadata:
            Target Region : ${meta.region}
            Paired End    : ${meta.paired_end}
            F Primer      : ${meta.f_primer}
            R Primer      : ${meta.r_primer}${colorCodes.c_reset}"""
        }

    emit:
        samples  = ch_samples
        runsheet = ch_runsheet
}