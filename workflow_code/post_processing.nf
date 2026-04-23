nextflow.enable.dsl=2

def prefix = params.output_prefix ?: ""
params.cleaned_prefix = (prefix && !prefix.endsWith("_") && !prefix.endsWith("-")) ? prefix + "_" : prefix

// Terminal text color definitions
c_back_bright_red = "\u001b[41;1m"
c_bright_green    = "\u001b[32;1m"
c_blue            = "\033[0;34m"
c_reset           = "\033[0m"

params.help = false


/**************************************************
* HELP MENU  **************************************
**************************************************/
if(params.help){

  println()
  println("GeneLab Post Processing Pipeline: $workflow.manifest.version")
  println("USAGE:")
  println("Example: Submit and run jobs with slurm in singularity containers.")
  println("   > nextflow -C post_processing.config run post_processing.nf -resume -profile slurm,singularity")
  println()
  println("Required Parameters:")
  println("""-profile [STRING] Specifies the profile to be used to run the workflow. Options are [slurm, singularity, docker, and  conda].
	                 singularity, docker and conda will run the workflow locally using singularity, docker, and conda, respectively.
                      To combine profiles, separate two or more profiles with a comma. 
                      For example, to combine slurm and singularity profiles, pass 'slurm,singularity' as argument. """)	
  println("  --publish_dir_mode [STRING]  Specifies how nextflow handles output file publishing. Options can be found here https://www.nextflow.io/docs/latest/process.html#publishdir Default: link.")
  println("  --GLDS_accession [STRING]  A Genelab GLDS accession number. Example GLDS-487. Default: null")
  println("  --OSD_accession [STRING]  A Genelab OSD accession number. Example OSD-487. Default: null")
  println("  --name [STRING] The analyst's full name. E.g. 'FirstName A. LastName'.  Default: FirstName A. LastName")
  println("  --email [STRING] The analyst's email address. E.g. 'mail@nasa.gov'.  Default: mail@nasa.gov")
  println("  --assay_suffix [STRING]  Genelab's assay suffix. Default: _GLAmpSeq.")
  println("  --output_prefix [STRING] Unique name to tag onto output files. Adds a "_" separator to the string if it is not empty and does not end with '_' or '-'. Default: empty string.")
  println("  --V_V_guidelines_link [URL] Genelab metagenomics data validation and verification guidelines link. Default: https://github.com/nasa/GeneLab_AmpliconSeq_Workflow/blob/main/README.md.")
  println("  --target_files [STRING] A comma separated list of target files to find in processing_info.zip. Default: command.txt,nextflow_processing_info_GLAmpSeq.txt.")
  println()
  println("Extra parameters to scripts:")
  println("      --readme_extra [STRING] Extra parameters and arguments to GL-gen-processed-amplicon-data-readme command. Run 'GL-gen-processed-amplicon-readme --help' for extra parameters that can be set. Example '--raw-reads-dir  ./Raw_Sequence_Data/'. Default: empty string")
  println("      --validation_extra [STRING] Extra parameters and arguments to GL-validate-processed-amplicon-data command. Run 'GL-validate-processed-amplicon-data --help' for extra parameters that can be set. Example '--single-ended --R1-used-as-single-ended-data --skip_raw_multiqc'. Default: empty string ")
  println("      --file_association_extra [STRING] Extra parameters and arguments to GL-gen-amplicon-file-associations-table command. Run 'GL-gen-amplicon-file-associations-table --help' for extra parameters that can be set. Example '--single-ended --R1-used-as-single-ended-data'. Default: empty string ")
  println()
  println("Files:")
  println("    --runsheet  [PATH] A 3-column (single-end) or 4-column (paired-end) input file (sample_id, forward, [reverse,] paired) used to run the processing pipeline. This is the value set to the parameter --csv_file when run the processing pipeline with a csv file as input otherwise it is the GLfile.csv in the GeneLab directory if --GLDS_accession was used as input. Example './GeneLab/GLfile.csv'.  Default: null")
  println()
  println("Optional arguments:")  
  println("    --help  Print this help message and exit")
  println("    --debug [BOOLEAN] Set to true if you'd like to see the values of your set parameters printed to the terminal. Default: false.")
  println()
  println("Paths to existing conda environments to use otherwise a new one will be created using the yaml file in envs/.")
  println("      --conda_dp_tools [PATH] Path to a conda environment containing dp_tools. Default: null.")
  exit 0
}


/************************************************
*********** Show pipeline parameters ************
*************************************************/
if(params.debug){

log.info """${c_blue}
         GeneLab Post Processing Pipeline: $workflow.manifest.version
         
         You have set the following parameters:
         Profile: ${workflow.profile} 
         Analyst's Name : ${params.name}
         Analyst's Email : ${params.email}
         GLDS Accession : ${params.GLDS_accession}
         OSD Accession : ${params.OSD_accession}
         Assay Suffix: ${params.assay_suffix}
         Output Prefix: ${params.output_prefix}
         V & V Link: ${params.V_V_guidelines_link}
         Target Files: ${params.target_files}
         Nextflow Directory publishing mode: ${params.publish_dir_mode}
         
         Extra scripts parameters:
         Readme Script Extra: ${params.readme_extra}
         Validation Script Extra : ${params.validation_extra}
         File Association Script Extra: ${params.file_association_extra}

         Files:
         Input Runsheet: ${params.runsheet}
         ${c_reset}"""

}


include { CLEAN_MULTIQC_PATHS; PACKAGE_PROCESSING_INFO; GENERATE_README; VALIDATE_PROCESSING;
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

       // ---------------------- Input channels -------------------------------- //
       // Input files
       // Runsheet used to execute the processing workflow
       runsheet_ch = Channel.fromPath(params.runsheet, checkIfExists: true)

       def processed_dir = file(params.outdir)
       if( !processed_dir.exists() ) {
            error "Output directory '${processed_dir}' does not exist. Make sure the main workflow is run first."
       }

       processing_scripts_dir = Channel.fromPath("${processed_dir}/processing_scripts", type: 'dir')
       
       fastqc_outputs_dir = channel.fromPath("${processed_dir}/FastQC_Outputs", type: 'dir')
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

        CLEAN_MULTIQC_PATHS(fastqc_outputs_dir)

       // Automatic verification and validation
        VALIDATE_PROCESSING(channel.value(processed_dir), params.target_files, runsheet_ch, 
                            GENERATE_README.out.readme,
                            PACKAGE_PROCESSING_INFO.out.zip, 
                            CLEAN_MULTIQC_PATHS.out.clean_dir) 
        
        // Generate md5sums  
          // Make sure md5sums are generated after the following processes are done since they are included in the md5sum generation:
          ch_ready = GENERATE_README.out.readme                 
                .combine( CLEAN_MULTIQC_PATHS.out.clean_dir )
                .combine(PACKAGE_PROCESSING_INFO.out.zip) 
                .map { true }
          
          GENERATE_MD5SUMS(channel.value(processed_dir), ch_ready)

        // Generate file association table for curation
        GENERATE_CURATION_TABLE(runsheet_ch, fastqc_outputs_dir, final_outputs_dir)

        // Write methods
        GENERATE_PROTOCOL(software_versions_ch, params.protocol_id, rarefaction_depth_ch)

      publish:
      processing_info_ch = PACKAGE_PROCESSING_INFO.out.zip
      readme_ch = GENERATE_README.out.readme
      purged_fastqc_ch = CLEAN_MULTIQC_PATHS.out.clean_dir
      log_ch = VALIDATE_PROCESSING.out.log
      raw_md5sum_ch = GENERATE_MD5SUMS.out.raw_md5sum
      processed_md5sum_ch = GENERATE_MD5SUMS.out.processed_md5sum
      curation_table_ch = GENERATE_CURATION_TABLE.out.curation_table
      protocol_ch = GENERATE_PROTOCOL.out.protocol  

}

output {
    // Post-processing outputs
    processing_info_ch {
        path "Post_Processing"
    }

    readme_ch {
        path "Post_Processing"
    }

    purged_fastqc_ch {
        path "Post_Processing"
    }

    log_ch {
        path "Post_Processing"
    }

    raw_md5sum_ch {
        path "Post_Processing"
    }

    processed_md5sum_ch {
        path "Post_Processing"
    }

    curation_table_ch {
        path "Post_Processing"
    }

    protocol_ch {
        path "Post_Processing"
    }
}