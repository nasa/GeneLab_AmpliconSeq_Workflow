nextflow.enable.dsl=2

include { validateParameters } from 'plugin/nf-schema'

def prefix = params.output_prefix ?: ""
params.cleaned_prefix = (prefix && !prefix.endsWith("_") && !prefix.endsWith("-")) ? prefix + "_" : prefix

validateParameters()

// Terminal text color definitions
c_back_bright_red = "\u001b[41;1m"
c_bright_green    = "\u001b[32;1m"
c_blue            = "\033[0;34m"
c_reset           = "\033[0m"

/************************************************
*********** Show pipeline parameters ************
*************************************************/
if(params.debug){

log.info """${c_blue}
         GeneLab Post Processing Pipeline: $workflow.manifest.version
         
         You have set the following parameters:

         Input/Output Options:
         GLDS Accession : ${params.GLDS_accession}
         OSD Accession : ${params.OSD_accession}
         Input Runsheet: ${params.runsheet}
         Output Directory: ${params.outdir}
         Assay Suffix: ${params.assay_suffix}
         Output Prefix: ${params.output_prefix}
         Target Files: ${params.target_files}

         Analyst Options:
         Analyst's Name : ${params.name}
         Analyst's Email : ${params.email}
         Protocol ID: ${params.protocol_id}
         V & V Link: ${params.V_V_guidelines_link}

         Processing Flags:
         Include Raw Data: ${params.include_raw_data}
         Trim Primers: ${params.trim_primers}
         Single End: ${params.single_end}
         Force Processing Single-End: ${params.force_single_end}
         Include Raw MultiQC: ${params.include_raw_multiqc}

         General Pipeline Settings:
         Profile: ${workflow.profile}
         Nextflow Directory Publishing Mode: ${params.publish_dir_mode}
         ${c_reset}"""

}


include { PACKAGE_PROCESSING_INFO; GENERATE_README; VALIDATE_PROCESSING;
           GENERATE_CURATION_TABLE; GENERATE_MD5SUMS; GENERATE_PROTOCOL} from './modules/genelab.nf'


workflow {
      main:
        // Make sure accessions numbers are set
        if(!params.GLDS_accession || !params.OSD_accession){
           error("""${c_back_bright_red}ACCESSION ERROR!. 
                    Please supply both --GLDS_accession and --OSD_accession.
                    They can be any string you choose but they must be set.
                 ${c_reset}""")
        }

        // If force_single_end is set, ignore include_raw_data and include_raw_multiqc as they should reflect paired-end data which is assumed to have been already included
        if (params.force_single_end && (params.include_raw_data || params.include_raw_multiqc)) {
            log.warn "force_single_end is set — --include_raw_data and --include_raw_multiqc will both be ignored."
        }

       // ---------------------- Input channels -------------------------------- //
       // Input files
       // Runsheet used to execute the processing workflow
       runsheet_ch = Channel.fromPath(params.runsheet, checkIfExists: true)

       def processed_dir = file(params.outdir)
       if( !processed_dir.exists() ) {
            error "Output directory '${processed_dir}' does not exist. Make sure the main workflow is run first."
       }

       processing_scripts_dir = Channel.fromPath("${processed_dir}/processing_scripts", type: 'dir')
       
       raw_multiqc_report_dir = channel.fromPath("${processed_dir}/Raw_Sequence_Data/MultiQC_Reports", type: 'dir')
       final_outputs_dir = channel.fromPath("${processed_dir}/Final_Outputs", type: 'dir')
     
       software_versions_ch   =  channel.fromPath("${processed_dir}/GeneLab/software_versions_*.txt")
       rarefaction_depth_ch =  channel
                                  .fromPath("${processed_dir}/Final_Outputs/alpha_diversity/*rarefaction_depth${params.assay_suffix}.txt")
                                  ifEmpty { file("empty_depth.txt").tap { it.text = "" } }


        // If an assay table is provided use it as the input table otherwise use the provided ISA zip file - no longer needed, using GLfile.csv for GENERATE_CURATION_TABLE
        //input_table_ch = Channel.fromPath( params.assay_table ? params.assay_table : params.isa_zip,
          //                                checkIfExists: true)

        // ---------------------- Post-processing begins -------------------------------------- //
        
        PACKAGE_PROCESSING_INFO(processing_scripts_dir)

        GENERATE_README()

        // Automatic verification and validation
        VALIDATE_PROCESSING(channel.value(processed_dir), params.target_files, runsheet_ch, 
                            GENERATE_README.out.readme,
                            PACKAGE_PROCESSING_INFO.out.zip) 
        
        // Generate md5sums  
          // Make sure md5sums are generated after the following processes are done since they are included in the md5sum generation:
          ch_ready = GENERATE_README.out.readme                 
                .combine(PACKAGE_PROCESSING_INFO.out.zip) 
                .map { true }
          
          GENERATE_MD5SUMS(channel.value(processed_dir), ch_ready)

        // Generate file association table for curation
        GENERATE_CURATION_TABLE(runsheet_ch, raw_multiqc_report_dir, final_outputs_dir)

        // Write methods
        GENERATE_PROTOCOL(software_versions_ch, params.protocol_id, rarefaction_depth_ch)

      publish:
      processing_info_ch = PACKAGE_PROCESSING_INFO.out.zip
      readme_ch = GENERATE_README.out.readme
      log_ch = VALIDATE_PROCESSING.out.log
      raw_md5sum_ch = GENERATE_MD5SUMS.out.raw_md5sum
      processed_md5sum_ch = GENERATE_MD5SUMS.out.processed_md5sum
      curation_table_ch = GENERATE_CURATION_TABLE.out.curation_table
      protocol_ch = GENERATE_PROTOCOL.out.protocol  

}

output {
    // Post-processing outputs
    processing_info_ch { path "GeneLab" }

    readme_ch { path "GeneLab" }

    log_ch { path "GeneLab" }

    raw_md5sum_ch { path "GeneLab" }

    processed_md5sum_ch { path "GeneLab" }

    curation_table_ch { path "GeneLab" }

    protocol_ch { path "GeneLab" }
}