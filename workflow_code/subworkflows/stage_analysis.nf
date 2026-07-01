include { GET_ACCESSIONS    } from '../modules/get_accessions.nf'
include { FETCH_ISA         } from '../modules/fetch_isa.nf'
include { ISA_TO_RUNSHEET   } from '../modules/isa_to_runsheet.nf'
include { SELECT_RUNSHEET   } from '../modules/select_runsheet.nf'
include { PARSE_RUNSHEET    } from './parse_runsheet.nf'
include { COPY_READS        } from '../modules/copy_reads.nf'
include { COPY_REMOTE_READS } from '../modules/copy_reads.nf'

process WRITE_INPUT_FILE {

    input:
        val(samples)
        val(paired_end)

    output:
        path("GLfile.csv"), emit: gl_file

    exec:
    def lines = paired_end ?
        ["sample_id,forward,reverse,paired,groups"] :
        ["sample_id,forward,paired,groups"]

    samples.each { meta, reads ->
        def r1 = reads[0].toUriString()
        if (paired_end) {
            def r2 = reads[1].toUriString()
            lines << "${meta.id},${r1},${r2},true,${meta.groups}"
        } else {
            lines << "${meta.id},${r1},false,${meta.groups}"
        }
    }

    def outfile = task.workDir.resolve("GLfile.csv")
    outfile.text = lines.join('\n') + '\n'
}

/**
 * STAGE_ANALYSIS
 *
 * This subworkflow handles the initial setup of the AmpliconSeq analysis:
 *
 * Accession mode (--accession):
 *   1. Resolves OSD and GLDS accession numbers via the OSDR API
 *   2. Fetches the ISA archive (or uses a locally provided one via --isa_archive)
 *   3. Converts the ISA archive to a runsheet using dp_tools
 *   4. Selects the correct runsheet by target region (16S/18S/ITS)
 *   5. Parses the runsheet, validates metadata consistency, and generates GLfile.csv
 *   6. Stages raw reads via COPY_REMOTE_READS
 *
 * Input file mode (--input_file):
 *   1. Reads the user-provided input CSV file (same format as the GLfile.csv
 *      generated in accession mode, with columns: sample_id, forward, [reverse,]
 *      paired, groups)
 *   2. Stages raw reads via COPY_READS
 */
workflow STAGE_ANALYSIS {
    take:
        accession
        target_region
        input_file
        isa_archive_path
        api_url
        dp_tools_plugin

    main:
        channel.empty() | set { isa_archive_ch }
        channel.empty() | set { osd_accession }
        channel.empty() | set { glds_accession }
        channel.empty() | set { gl_file_ch }
        channel.empty() | set { runsheet_ch }
        channel.empty() | set { primers_ch }
        channel.empty() | set { software_versions_ch }

        if ( accession ) {
            GET_ACCESSIONS( accession, api_url )
            osd_accession  = GET_ACCESSIONS.out.accessions_txt.map { it.readLines()[0].trim() }
            glds_accession = GET_ACCESSIONS.out.accessions_txt.map { it.readLines()[1].trim() }

            // Use provided ISA archive or fetch it
            if ( isa_archive_path ) {
                isa_archive_ch = isa_archive_path
            } else {
                FETCH_ISA( osd_accession, glds_accession )
                isa_archive_ch = FETCH_ISA.out.isa_archive
            }

            ISA_TO_RUNSHEET( osd_accession, glds_accession, isa_archive_ch, dp_tools_plugin )
            software_versions_ch = ISA_TO_RUNSHEET.out.version

            SELECT_RUNSHEET( ISA_TO_RUNSHEET.out.runsheet, target_region, params.specify_runsheet ?: "" )
            runsheet_ch = SELECT_RUNSHEET.out.runsheet

            PARSE_RUNSHEET( runsheet_ch )

            // Extract primers from runsheet
            primers_ch = runsheet_ch
                .splitCsv(header: true)
                .map { row -> ["${row.F_Primer}", "${row.R_Primer}"] }
                .first()

            // Extract paired_end from first sample to pass alongside collected samples
            paired_end_ch = PARSE_RUNSHEET.out.samples
                .map { meta, reads -> meta.paired_end }
                .first()

            WRITE_INPUT_FILE(
                PARSE_RUNSHEET.out.samples.toList(),
                paired_end_ch
            )
            gl_file_ch = WRITE_INPUT_FILE.out.gl_file

            // Build reads channel from GLfile.csv
            reads_ch = WRITE_INPUT_FILE.out.gl_file
                .splitCsv(header: true)
                .map { row ->
                    row.paired.trim().toLowerCase() == 'true'
                        ? tuple("${row.sample_id}", ["${row.forward}", "${row.reverse}"], "true")
                        : tuple("${row.sample_id}", ["${row.forward}"], "false")
                }

            COPY_REMOTE_READS( reads_ch )
            staged_reads_ch = COPY_REMOTE_READS.out.raw_reads

        } else {
            // Input file mode
            runsheet_ch    = input_file
            isa_archive_ch = channel.empty()
            gl_file_ch     = channel.empty()
            primers_ch     = Channel.value([params.F_primer, params.R_primer])

            reads_ch = input_file
                .splitCsv(header: true)
                .map { row ->
                    row.paired.trim().toLowerCase() == 'true'
                        ? tuple("${row.sample_id}", [file("${row.forward}"), file("${row.reverse}")], "true")
                        : tuple("${row.sample_id}", [file("${row.forward}")], "false")
                }

            COPY_READS( reads_ch )
            staged_reads_ch = COPY_READS.out.raw_reads
        }

    emit:
        staged_reads      = staged_reads_ch
        runsheet          = runsheet_ch
        isa_archive       = isa_archive_ch
        gl_file           = gl_file_ch
        primers           = primers_ch
        software_versions = software_versions_ch
}