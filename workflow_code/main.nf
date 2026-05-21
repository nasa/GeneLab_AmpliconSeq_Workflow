nextflow.enable.dsl = 2

include { validateParameters } from 'plugin/nf-schema'

def prefix = params.output_prefix ?: ""
params.cleaned_prefix = (prefix && !prefix.endsWith("_") && !prefix.endsWith("-")) ? prefix + "_" : prefix

validateParameters()

// Terminal text color definitions
c_back_bright_red = "\u001b[41;1m";
c_bright_green    = "\u001b[32;1m";
c_blue            = "\033[0;34m";
c_reset           = "\033[0m";

/************************************************
*********** Show pipeline parameters ************
*************************************************/
if(params.debug){

log.info """${c_blue}
         Nextflow AmpIllumina Consensus Pipeline: $workflow.manifest.version
         
         You have set the following parameters:

         Amplicon target region : ${params.target_region}
         GLDS or OSD accession : ${params.accession}
         Input csv file : ${params.input_file}
         Output directory: ${params.outdir}
         Database Store Directory: ${params.database_store_path}
         Genelab Assay Suffix: ${params.assay_suffix}
         Output Prefix: ${params.output_prefix}
         Trim Primers: ${params.trim_primers}

         Cutadapt Parameters:
         Forward Primer: ${params.F_primer}
         Reverse Primer: ${params.R_primer}
         Minimum Trimmed Reads length: ${params.min_cutadapt_len}
         Primers Are linked: ${params.primers_linked}
         Primers Are Anchored: ${params.anchored_primers}
         Discard Untrimmed Reads: ${params.discard_untrimmed}
 
         Dada2 Parameters:
         Truncate left: ${params.left_trunc}bp
         Truncate right: ${params.right_trunc}bp
         Max error left: ${params.left_maxEE}
         Max error right: ${params.right_maxEE}
         Concatenate Reads: ${params.concatenate_reads_only}
         
         Diversity and Differential abundance Parameters:
         Method: ${params.diff_abund_method}
         Rarefaction Depth: ${params.rarefaction_depth}
         Remove Structural Zeros: ${params.remove_struc_zeros}
         Remove Rare Taxa and Samples: ${params.remove_rare}
         Taxa Prevalence Cut Off: ${params.prevalence_cutoff}
         Sample Library Cut Off: ${params.library_cutoff}
         Groups to Comapre Column: ${params.group}
         Samples Column: ${params.samples_column}

         Debugging Options:
         Limit Samples for Testing: ${params.limit_samples_to}
         Force Processing Single-End: ${params.force_single_end}

         General Pipeline Settings:
         Nextflow Directory publishing mode: ${params.publish_dir_mode}
         MultiQC configuration file: ${params.multiqc_config}
         Nextflow Error strategy: ${params.errorStrategy}

         Conda Environments:
         dp_tools: ${params.conda_dp_tools}
         fastqc: ${params.conda_fastqc}
         multiqc: ${params.conda_multiqc}
         cutadapt: ${params.conda_cutadapt}
         R: ${params.conda_R}
         Diversity and Differential abundance : ${params.conda_diversity}
         zip: ${params.conda_zip}
         wget: ${params.conda_wget}
         ${c_reset}"""
}

// Stage analysis setup (inputs, and raw reads)
include { STAGE_ANALYSIS } from './subworkflows/stage_analysis.nf'

// Read quality check and filtering
include { FASTQC as RAW_FASTQC ; MULTIQC as RAW_MULTIQC  } from './modules/quality_assessment.nf'

// Trim primers if requested
include { CUTADAPT ; COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE } from './modules/cutadapt.nf'

// Cluster ASVs
include { DOWNLOAD_DATABASE } from './modules/download_database.nf'
include { RUN_DADA2 } from './modules/run_dada.nf'

// Filtered quality check
include { FASTQC as FILTERED_FASTQC ; MULTIQC as FILTERED_MULTIQC  } from './modules/quality_assessment.nf'

// Diversity, differential abundance and visualizations
include { ALPHA_DIVERSITY; BETA_DIVERSITY } from './modules/diversity.nf'
include { PLOT_TAXONOMY } from './modules/taxonomy_plots.nf'
include { ZIP as ZIP_BIOM; ZIP as ZIP_ALPHA; ZIP as ZIP_BETA_EUCLIDEAN; ZIP as ZIP_BETA_BRAY; ZIP as ZIP_TAXONOMY_SAMPLES; ZIP as ZIP_TAXONOMY_GROUPS } from './modules/zip.nf'
include { ANCOMBC as ANCOMBC1 } from './modules/ancombc.nf'
include { ANCOMBC as ANCOMBC2 } from './modules/ancombc.nf'
include { DESEQ } from './modules/deseq.nf'
include { ZIP as ZIP_DA; ZIP as ZIP_ANCOMBC1; ZIP as ZIP_ANCOMBC2; ZIP as ZIP_DESEQ2 } from './modules/zip.nf'
include { SOFTWARE_VERSIONS } from './modules/utils.nf'

ch_dp_tools_plugin = params.dp_tools_plugin ? channel.value(file(params.dp_tools_plugin)) : channel.value(file("$projectDir/bin/dp_tools__NF_AmpIllumina_${params.target_region}"))

ch_input_file = params.input_file ? channel.fromPath(params.input_file) : null
ch_isa_archive = params.isa_archive ? channel.fromPath(params.isa_archive) : null

// A function to delete white spaces from an input string and covert it to lower case
def deleteWS(string){

    return string.replaceAll(/\s+/, '').toLowerCase()

}


workflow {
    main:

    //  ---------------------  Sanity Checks ------------------------------------- //
    // Test input requirement
    if (!params.accession &&  !params.input_file){
       error("""${c_back_bright_red}INPUT ERROR! 
              Please supply either an accession (OSD or Genelab number) or an input CSV file
              by passing either to the --accession or --input_file parameter, respectively.
              ${c_reset}""")
    } 
    
    // Test input csv file
    if(params.input_file){
        // Test primers
        if(!params.F_primer || !params.R_primer){

            error("""${c_back_bright_red}PRIMER ERROR! 
                  When using a csv file as input (--input_file) to this workflow you must provide 
                  foward and reverse primer sequences. Please provide your forward 
                  and reverse primer sequences as arguements to the --F_primer 
                  and --R_primer parameters, respectively.
                  ${c_reset}""")
         }
    }

    // Test ISA archive and accession
    if (params.isa_archive && !params.accession) {
        error """${c_back_bright_red}INPUT ERROR!
            --isa_archive requires --accession to resolve OSD/GLDS accessions
            for the ISA-to-runsheet conversion.${c_reset}"""
    }

    software_versions_ch = channel.empty()
        
    // Stage analysis setup (inputs, and raw reads)
    STAGE_ANALYSIS(
        params.accession,
        params.target_region,
        ch_input_file,
        ch_isa_archive,
        params.api_url,
        ch_dp_tools_plugin
    )
    staged_reads_ch = STAGE_ANALYSIS.out.staged_reads
    runsheet_ch = STAGE_ANALYSIS.out.runsheet
    isa_archive_ch = STAGE_ANALYSIS.out.isa_archive
    gl_file_ch = STAGE_ANALYSIS.out.gl_file
    primers_ch = STAGE_ANALYSIS.out.primers
    
    STAGE_ANALYSIS.out.software_versions
        | mix(software_versions_ch)
        | set { software_versions_ch }


    // Read quality check and trimming
    RAW_FASTQC(staged_reads_ch)
    raw_fastqc_files = RAW_FASTQC.out.fastqc.flatten().collect()
    
    RAW_MULTIQC("raw", params.multiqc_config,raw_fastqc_files)

    RAW_FASTQC.out.version | mix(software_versions_ch) | set{software_versions_ch}
    RAW_MULTIQC.out.version | mix(software_versions_ch) | set{software_versions_ch}

    // Download reference database for taxonomic classification
    def db_config = [
            "16S": ["SILVA_SSU_r138_2_v2.RData", "https://api.figshare.com/v2/file/download/64078939"],
            "ITS": ["UNITE_v2025.RData", "https://api.figshare.com/v2/file/download/64079011"],
            "18S": ["PR2_v4_13_March2021.RData", "https://api.figshare.com/v2/file/download/46241917"]
        ]
    target_region_ch = Channel.value(params.target_region)
        .map { region -> tuple(region, db_config[region][0], db_config[region][1]) }
    
    DOWNLOAD_DATABASE(target_region_ch)

    trimmed_reads_ch = channel.empty()
    trimmed_reads_counts = channel.empty()
    cutadapt_logs = channel.empty()
    if(params.trim_primers){

        //if(!params.accession) primers_ch = channel.value([params.F_primer, params.R_primer]) // to be removed once stage analysis workflow is implemented
        CUTADAPT(staged_reads_ch, primers_ch)
        logs = CUTADAPT.out.logs.map{ sample_id, log -> file("${log}")}.collect()
        counts = CUTADAPT.out.trim_counts.map{ sample_id, count -> file("${count}")}.collect()
        trimmed_reads_ch = CUTADAPT.out.reads.map{ 
                                              sample_id, reads, isPaired -> reads instanceof List ? reads.each{file("${it}")}: file("${reads}")
                                              }.flatten().collect()

        COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE(counts, logs, runsheet_ch)
        trimmed_reads_counts = COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE.out.counts
        cutadapt_logs = COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE.out.logs

        isPaired_ch = CUTADAPT.out.reads.map{ 
                                              sample_id, reads, isPaired -> isPaired
                                              }.first()

        samples_ch = runsheet_ch.first()
                     .concat(isPaired_ch)
                     .collate(2)
        
        
        // Run dada2
        RUN_DADA2(samples_ch, trimmed_reads_ch, trimmed_reads_counts, DOWNLOAD_DATABASE.out.database)

        CUTADAPT.out.version | mix(software_versions_ch) | set{software_versions_ch}
    }else{
        raw_reads_ch = staged_reads_ch.map{
                          sample_id, reads, isPaired -> reads instanceof List ? reads.each{file("${it}")}: file("${reads}")
                          }.flatten().collect()

        isPaired_ch = staged_reads_ch.map{sample_id, reads, isPaired -> isPaired}.first()
        samples_ch = runsheet_ch.first()
                     .concat(isPaired_ch)
                     .collate(2)
        
        // Run dada2 without primer trimming
        RUN_DADA2(samples_ch, raw_reads_ch, file("NO_FILE"), DOWNLOAD_DATABASE.out.database)
    }

    dada_counts = RUN_DADA2.out.counts
    dada_taxonomy = RUN_DADA2.out.taxonomy
    dada_biom = RUN_DADA2.out.biom
    filtered_count = RUN_DADA2.out.filtered_count

    filtered_reads_ch = RUN_DADA2.out.reads
            .flatten()
            .map { file ->
                    // derive sample_id from filename
                    def sample_id
                    if (file.name.endsWith("${params.assay_suffix}_R1_filtered.fastq.gz")) {
                            sample_id = file.name.replace("${params.assay_suffix}_R1_filtered.fastq.gz", "")
                    } else if (file.name.endsWith("${params.assay_suffix}_R2_filtered.fastq.gz")) {
                            sample_id = file.name.replace("${params.assay_suffix}_R2_filtered.fastq.gz", "")
                    }

                    tuple(sample_id, file)
            }
            .groupTuple(by:0)  // group R1/R2 by sample_id
            .map { sample_id, files ->
                    def pathFiles = files.collect { it instanceof String ? file(it) : it }  // ensure Path objects
                    def isPaired = pathFiles.size() > 1
                    tuple(sample_id, pathFiles, isPaired)
            }

    FILTERED_FASTQC(filtered_reads_ch)
    	filtered_fastqc_files = FILTERED_FASTQC.out.fastqc.flatten().collect()

    FILTERED_MULTIQC("filtered", params.multiqc_config, filtered_fastqc_files)

    RUN_DADA2.out.version | mix(software_versions_ch) | set{software_versions_ch}
    FILTERED_FASTQC.out.version | mix(software_versions_ch) | set{software_versions_ch}
    FILTERED_MULTIQC.out.version | mix(software_versions_ch) | set{software_versions_ch}

    // Zip biom file
    dada_biom
        .map { biom -> tuple("taxonomy-and-counts", biom) }
        | ZIP_BIOM

    ZIP_BIOM.out.version | mix(software_versions_ch) | set{software_versions_ch}


   
    // Diversity, differential abundance testing and their corresponding visualizations
    if(params.accession){

        values = ["samples": "Sample Name",
                "group" : "groups",
                "depth" : params.rarefaction_depth,
                "assay_suffix" : params.assay_suffix,
                "output_prefix" : params.cleaned_prefix,
                "target_region" : params.target_region,
                "library_cutoff" : params.library_cutoff,
                "prevalence_cutoff" : params.prevalence_cutoff,
                "rare" : params.remove_rare ? "--remove-rare" : "",
                "struc_zero": params.remove_struc_zeros ? "--remove-structural-zeros" : ""
                ]
        
        metadata  =  runsheet_ch

    }else{

        values = ["samples": params.samples_column,
                "group" : params.group,
                "depth" : params.rarefaction_depth,
                "assay_suffix" : params.assay_suffix,
                "output_prefix" : params.cleaned_prefix,
                "target_region" : params.target_region,
                "library_cutoff" : params.library_cutoff,
                "prevalence_cutoff" : params.prevalence_cutoff,
                "rare" :  params.remove_rare ? "--remove-rare" : "",
                "struc_zero": params.remove_struc_zeros ? "--remove-structural-zeros" : ""
                ]
        
        metadata  =  ch_input_file

    }
    meta  = channel.of(values)
    
    // Diversity analysis
    ALPHA_DIVERSITY(meta, dada_counts, dada_taxonomy, metadata)
    BETA_DIVERSITY(meta, dada_counts, dada_taxonomy, metadata)

    // Zipping diversity plots
    // Alpha diversity (if rarefaction succeeded)
    ALPHA_DIVERSITY.out.output_dir
    	.map { dir ->
		    def pngs = file(dir).listFiles()?.findAll { it.name.endsWith('.png') }
        	pngs ? tuple(
            	"alpha_diversity_plots",
            	pngs
        	) : null
    	}
    	.filter { it != null }
    	| ZIP_ALPHA

    // Beta diversity - euclidean distance
    BETA_DIVERSITY.out.output_dir
	    .map { dir ->
		    def pngs = file(dir).listFiles()?.findAll { it.name.contains('euclidean') && it.name.endsWith('.png') }
        	pngs ? tuple(
            	"euclidean_distance_plots",
            	pngs
        	) : null
    	}
	    .filter { it != null }
        | ZIP_BETA_EUCLIDEAN

    // Beta diversity - bray curtis (if rarefaction succeeded)
    BETA_DIVERSITY.out.output_dir
        .map { dir ->
            def pngs = file(dir).listFiles()?.findAll { it.name.contains('bray') && it.name.endsWith('.png') }
            pngs ? tuple(
                "bray_curtis_plots",
                pngs
            ) : null
        }
        .filter { it != null }
    	| ZIP_BETA_BRAY

    // Taxonomy plotting
    PLOT_TAXONOMY(meta, dada_counts, dada_taxonomy, metadata)

    // Zipping taxonomy plots
   // Sample plots
    PLOT_TAXONOMY.out.output_dir
        .map { dir ->
            def pngs = file(dir).listFiles()?.findAll { it.name.contains('samples') && it.name.endsWith('.png') }
            pngs ? tuple(
                "sample_taxonomy_plots",
                pngs
            ) : null
        }
        .filter { it != null }
        | ZIP_TAXONOMY_SAMPLES

   // Group taxonomy plots
    PLOT_TAXONOMY.out.output_dir
        .map { dir ->
            def pngs = file(dir).listFiles()?.findAll { it.name.contains('groups') && it.name.endsWith('.png') }
            pngs ? tuple(
                "group_taxonomy_plots",
                pngs
            ) : null
        }
        .filter { it != null }
        | ZIP_TAXONOMY_GROUPS
    
    ALPHA_DIVERSITY.out.version | mix(software_versions_ch) | set{software_versions_ch}
    BETA_DIVERSITY.out.version | mix(software_versions_ch) | set{software_versions_ch}
    PLOT_TAXONOMY.out.version | mix(software_versions_ch) | set{software_versions_ch}
    
     // Differential abundance testing
     da_contrasts_ch = channel.empty()
     da_sampleTable_ch = channel.empty()
     ancombc1_ch = channel.empty()
     zip_ancombc1_ch = channel.empty()
     ancombc2_ch = channel.empty()
     zip_ancombc2_ch = channel.empty()
     deseq2_ch = channel.empty()
     zip_deseq2_ch = channel.empty()

     method = channel.of(params.diff_abund_method)
     if (params.diff_abund_method == "deseq2"){
    
        DESEQ(meta, dada_counts, dada_taxonomy, metadata, filtered_count)
        deseq2_ch = DESEQ.out.output_dir
        da_contrasts_ch = DESEQ.out.contrasts_file
        da_sampleTable_ch = DESEQ.out.sample_table_file
        DESEQ.out.version | mix(software_versions_ch) | set{software_versions_ch}
        // Zipping DESeq2 plots
	    DESEQ.out.output_dir
		    .map { dir ->
                def pngs = file(dir).listFiles()?.findAll { it.name.contains('volcano') && it.name.endsWith('.png') }
                pngs ? tuple(
                    "deseq2_volcano_plots",
                    pngs
                ) : null
            }
            .filter { it != null }
  	        | ZIP_DESEQ2
        zip_deseq2_ch = ZIP_DESEQ2.out.zip
    
    }else if (params.diff_abund_method == "ancombc1"){
    
        ANCOMBC1(method, meta, dada_counts, dada_taxonomy, metadata, filtered_count)
        ancombc1_ch = ANCOMBC1.out.output_dir
        da_contrasts_ch = ANCOMBC1.out.contrasts_file
        da_sampleTable_ch = ANCOMBC1.out.sample_table_file
        ANCOMBC1.out.version | mix(software_versions_ch) | set{software_versions_ch}
        // Zipping ANCOMBC1 plots
	    ANCOMBC1.out.output_dir
            .map { dir ->
                def pngs = file(dir).listFiles()?.findAll { it.name.contains('volcano') && it.name.endsWith('.png') }
                pngs ? tuple(
                    "ancombc1_volcano_plots",
                    pngs
                ) : null
            }
            .filter { it != null }
            | ZIP_ANCOMBC1
        zip_ancombc1_ch = ZIP_ANCOMBC1.out.zip

    }else if (params.diff_abund_method == "ancombc2"){

        ANCOMBC2(method, meta, dada_counts, dada_taxonomy, metadata, filtered_count)
        ancombc2_ch = ANCOMBC2.out.output_dir
        da_contrasts_ch = ANCOMBC2.out.contrasts_file
        da_sampleTable_ch = ANCOMBC2.out.sample_table_file
        ANCOMBC2.out.version | mix(software_versions_ch) | set{software_versions_ch}
        // Zipping ANCOMBC2 plots
	    ANCOMBC2.out.output_dir
            .map { dir ->
                def pngs = file(dir).listFiles()?.findAll { it.name.contains('volcano') && it.name.endsWith('.png') }
                pngs ? tuple(
                    "ancombc2_volcano_plots",
                    pngs
                ) : null
            }
            .filter { it != null }
            | ZIP_ANCOMBC2
        zip_ancombc2_ch = ZIP_ANCOMBC2.out.zip

    }else{

        ANCOMBC1("ancombc1", meta, dada_counts, dada_taxonomy, metadata, filtered_count)
        ancombc1_ch = ANCOMBC1.out.output_dir
        da_contrasts_ch = ANCOMBC1.out.contrasts_file
        da_sampleTable_ch = ANCOMBC1.out.sample_table_file
        ANCOMBC1.out.version | mix(software_versions_ch) | set{software_versions_ch}

        ANCOMBC2("ancombc2", meta, dada_counts, dada_taxonomy, metadata, ANCOMBC1.out.output_dir)
        ancombc2_ch = ANCOMBC2.out.output_dir
        ANCOMBC2.out.version | mix(software_versions_ch) | set{software_versions_ch}

        DESEQ(meta, dada_counts, dada_taxonomy, metadata, ANCOMBC2.out.output_dir)
        deseq2_ch = DESEQ.out.output_dir
        DESEQ.out.version | mix(software_versions_ch) | set{software_versions_ch}

        // Zipping DA plots
	    //ANCOMBC1
	    ANCOMBC1.out.output_dir
            .map { dir ->
                def pngs = file(dir).listFiles()?.findAll { it.name.contains('volcano') && it.name.endsWith('.png') }
                pngs ? tuple(
                    "ancombc1_volcano_plots",
                    pngs
                ) : null
            }
            .filter { it != null }
            | ZIP_ANCOMBC1
        zip_ancombc1_ch = ZIP_ANCOMBC1.out.zip
	    //ANCOMBC2
	    ANCOMBC2.out.output_dir
            .map { dir ->
                def pngs = file(dir).listFiles()?.findAll { it.name.contains('volcano') && it.name.endsWith('.png') }
                pngs ? tuple(
                    "ancombc2_volcano_plots",
                    pngs
                ) : null
            }
            .filter { it != null }
            | ZIP_ANCOMBC2
        zip_ancombc2_ch = ZIP_ANCOMBC2.out.zip
	    // DESeq2
	    DESEQ.out.output_dir
            .map { dir ->
                def pngs = file(dir).listFiles()?.findAll { it.name.contains('volcano') && it.name.endsWith('.png') }
                pngs ? tuple(
                    "deseq2_volcano_plots",
                    pngs
                ) : null
            }
            .filter { it != null }
            | ZIP_DESEQ2
        zip_deseq2_ch = ZIP_DESEQ2.out.zip
    }
    

     // Software Version Capturing - combining all captured software versions
     nf_version = "Nextflow Version ".concat("${nextflow.version}")
     nextflow_version_ch = channel.value(nf_version)

     //  Write software versions to file
     software_versions_ch | map { it.text.strip() }
                          | unique
                          | mix(nextflow_version_ch)
                          | collectFile({it -> it}, newLine: true, cache: false)
                          | SOFTWARE_VERSIONS

    publish:
    // Metadata
    runsheet = runsheet_ch
    isa_archive = isa_archive_ch
    gl_file = gl_file_ch

    // Raw reads
    raw_reads = staged_reads_ch

    // Trimmed reads
    trimmed_reads = trimmed_reads_ch
    trimmed_count = trimmed_reads_counts
    cutadapt_logs = cutadapt_logs
    
    // Filtered reads
    filtered_reads = filtered_reads_ch
    filtered_count = filtered_count

    // FastQC
    raw_fastqc = RAW_FASTQC.out.fastqc
    filtered_fastqc = FILTERED_FASTQC.out.fastqc

    // MultiQC
    zip_multiqc_raw = RAW_MULTIQC.out.zipped_data
    html_multiqc_raw = RAW_MULTIQC.out.html
    zip_multiqc_filtered = FILTERED_MULTIQC.out.zipped_data
    html_multiqc_filtered = FILTERED_MULTIQC.out.html

    // Dada2 outputs
    asv = RUN_DADA2.out.fasta
    counts = RUN_DADA2.out.counts
    taxonomy = RUN_DADA2.out.taxonomy
    taxonomy_counts = RUN_DADA2.out.taxonomy_count
    biom_zip = ZIP_BIOM.out.zip
    read_count_tracking = RUN_DADA2.out.read_count

    // Alpha and beta diversity outputs
    alpha_diversity = ALPHA_DIVERSITY.out.output_dir
    zip_alpha_plots = ZIP_ALPHA.out.zip
    beta_diversity = BETA_DIVERSITY.out.output_dir
    zip_beta_euclidean_plots = ZIP_BETA_EUCLIDEAN.out.zip
    zip_beta_bray_plots = ZIP_BETA_BRAY.out.zip

    // Taxonomy plots
    taxonomy_plots = PLOT_TAXONOMY.out.output_dir
    zip_taxonomy_samples = ZIP_TAXONOMY_SAMPLES.out.zip
    zip_taxonomy_groups = ZIP_TAXONOMY_GROUPS.out.zip

    // Differential abundance outputs
    da_contrasts = da_contrasts_ch
    da_sampleTable = da_sampleTable_ch
    ancombc1 = ancombc1_ch
    zip_ancombc1 = zip_ancombc1_ch
    ancombc2 = ancombc2_ch
    zip_ancombc2 = zip_ancombc2_ch
    deseq2 = deseq2_ch
    zip_deseq2 = zip_deseq2_ch

    // GeneLab
    software_versions = SOFTWARE_VERSIONS.out.software_versions_txt

}

output {
    // Metadata
    runsheet { path "Metadata" }
    isa_archive { path "Metadata" }
    gl_file { path "Metadata" }

    // Raw reads
    raw_reads { path "Raw_Sequence_Data" }

    // Trimmed reads
    trimmed_reads { path "Trimmed_Sequence_Data" }
    trimmed_count { path "Trimmed_Sequence_Data" }
    cutadapt_logs { path "Trimmed_Sequence_Data" }
    
    // Filtered reads
    filtered_reads { path "Filtered_Sequence_Data" }
    filtered_count { path "Filtered_Sequence_Data" }

    // FastQC
    raw_fastqc { path {html, zip -> "Raw_Sequence_Data/FastQC_Outputs" } }
    filtered_fastqc { path {html, zip -> "Filtered_Sequence_Data/FastQC_Outputs" } }

    // MultiQC
    zip_multiqc_raw { path "Raw_Sequence_Data/MultiQC_Reports" }
    html_multiqc_raw { path "Raw_Sequence_Data/MultiQC_Reports" }

    zip_multiqc_filtered { path "Filtered_Sequence_Data/MultiQC_Reports" }
    html_multiqc_filtered { path "Filtered_Sequence_Data/MultiQC_Reports" }

    // Dada2 outputs
    asv { path "Final_Outputs" }
    counts { path "Final_Outputs" }
    taxonomy { path "Final_Outputs" }
    taxonomy_counts { path "Final_Outputs" }
    biom_zip { path "Final_Outputs" }
    read_count_tracking { path "Final_Outputs" }

    // Alpha and beta diversity outputs
    alpha_diversity { path "Final_Outputs" }
    zip_alpha_plots { path "Final_Outputs/alpha_diversity" }

    beta_diversity { path "Final_Outputs" }
    zip_beta_euclidean_plots { path "Final_Outputs/beta_diversity" }
    zip_beta_bray_plots { path "Final_Outputs/beta_diversity" }

    // Taxonomy plots
    taxonomy_plots { path "Final_Outputs" }
    zip_taxonomy_samples { path "Final_Outputs/taxonomy_plots" }
    zip_taxonomy_groups { path "Final_Outputs/taxonomy_plots" }

    // Differential abundance outputs
    da_contrasts { path "Final_Outputs" }
    da_sampleTable { path "Final_Outputs" }

    ancombc1 { path "Final_Outputs" }
    zip_ancombc1 { path "Final_Outputs/differential_abundance/ancombc1" }

    ancombc2 { path "Final_Outputs" }
    zip_ancombc2 { path "Final_Outputs/differential_abundance/ancombc2" }

    deseq2 { path "Final_Outputs" }
    zip_deseq2 { path "Final_Outputs/differential_abundance/deseq2" }

    // GeneLab
    software_versions { path "GeneLab" }
}