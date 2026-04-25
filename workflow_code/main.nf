nextflow.enable.dsl = 2

def prefix = params.output_prefix ?: ""
params.cleaned_prefix = (prefix && !prefix.endsWith("_") && !prefix.endsWith("-")) ? prefix + "_" : prefix

// Terminal text color definitions
c_back_bright_red = "\u001b[41;1m";
c_bright_green    = "\u001b[32;1m";
c_blue            = "\033[0;34m";
c_reset           = "\033[0m";


params.help = false

/**************************************************
* HELP MENU  **************************************
**************************************************/
if (params.help) {
  println()
  println("Nextflow AmpIllumina Consensus Pipeline: $workflow.manifest.version")
  println("USAGE:")
  println("Example 1: Submit and run jobs with slurm in singularity containers.")
  println("   > nextflow run main.nf -resume -profile slurm,singularity --input_file PE_file.csv --target_region 16S --F_primer AGAGTTTGATCCTGGCTCAG --R_primer CTGCCTCCCGTAGGAGT")
  println()
  println("Example 2: : Submit and run jobs with slurm in conda environments.")
  println("   > nextflow run main.nf -resume -profile slurm,conda --input_file SE_file.csv --target_region 1TS --F_primer AGAGTTTGATCCTGGCTCAG --R_primer CTGCCTCCCGTAGGAGT")
  println()
  println("Example 3: Run jobs locally in conda environments, supplying a GLDS or OSD accession, and specifying the path to an existing conda environment")
  println("   > nextflow run main.nf -resume -profile mamba --accession GLDS-487 --target_region 16S --conda_multiqc <path/to/existing/conda/environment>")
  println()
  println("Required arguments:")
  println("""   -profile [STRING] What profile should be used to run the workflow. Options are [singularity, docker, conda, slurm].
	                 singularity, docker and conda will run the pipelne locally using singularity, docker, and conda, respectively.
                         To combine profiles, pass them together separated by comma. For example, to run jobs using slurm in singularity containers use 'slurm,singularity' . """)			 
  println("     --input_file  [PATH] A 4-column (single-end) or 5-column (paired-end) input file (sample_id, forward, [reverse,] paired, groups). Mandatory if a GLDS or OSD accession is not provided.")
  println("                   Please see the files: SE_file.csv and PE_file.csv for single-end and paired-end examples, respectively.")
  println("                   The sample_id column should contain unique sample ids.")
  println("                   The forward and reverse columns should contain the absolute or relative path to the sample's forward and reverse reads.")
  println("                   The paired column should be true for paired-end or anything else for single-end reads.")
  println("                   The groups column contain group levels / treatments to be compared during diversity and differential abundance testing analysis. Default: null")
  println("     --target_region [STRING] What is the amplicon target region to be analyzed. Options are one of [16S, 18S, ITS]. Default: 16S.")
  println("     --trim_primers [STRING] Should primers be trimmed? TRUE or FALSE. Default: TRUE.") 
  println()
  println("Cutadapt (trimming) parameters:")
  println("	    --F_primer [STRING] Forward primer sequence e.g. AGAGTTTGATCCTGGCTCAG. Default: null.")
  println("	    --R_primer [STRING] Reverse primer sequence e.g. CTGCCTCCCGTAGGAGT. Default: null.")
  println("	    --min_cutadapt_len [INTEGER] What should be the minimum read length after quality trimming with cutadapt. Default: 130.")
  println("	    --primers_linked [STRING] Are the primers linked?. https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-paired-end-reads. Default: TRUE. ")
  println("     --anchored_primers [STRING] Are the primers anchored? https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-paired-end-reads. Default: TRUE.")
  println("	    --discard_untrimmed [STRING] Should untrimmed reads be discarded? Any supplied string except TRUE will not discard them. Default: TRUE.")
  println()	
  println("Optional arguments:")  
  println("  --help  Print this help message and exit.")
  println("  --debug [BOOLEAN] Set to true if you'd like to see the values of your set parameters printed to the terminal. Default: false.")
  println("  --publish_dir_mode [STRING]  How should nextflow publish file outputs. Options can be found here https://www.nextflow.io/docs/latest/process.html#publishdir. Default: link.")
  println("  --errorStrategy [STRING] How should nextflow handle errors. Options can be found here https://www.nextflow.io/docs/latest/process.html#errorstrategy. Default: terminate")
  println("  --multiqc_config [PATH] Path to a custom multiqc config file. Default: config/multiqc.config.")
  println()
  println("Dada2 parameters passed to filterAndTrim() function:")
  println("	    --left_trunc [INTEGER] truncate the sequences to the left by this number of bases. Default: 0.") 
  println("	    --right_trunc [INTEGER] truncate the sequences to the right by this number of bases. Default: 0.") 
  println("	    --left_maxEE [INTEGER] Maximum allowed errors to the left. Default: 1.")
  println("	    --right_maxEE [INTEGER] Maximum allowed errors to the right. Default: 1.")
  println("	    --concatenate_reads_only [STRING] Concatenate only with dada2 instead of merging paired reads if TRUE.")
  println("      This is typically used with primers like 515-926, that captured 18S fragments that are typically too long to merge.")
  println("      Note that 16S and 18S should have been separated already prior to running this workflow. This should likely be left as FALSE for any option other than 18S above.") 	    
  println("	     Values are TRUE or FALSE Default: FALSE.")
  println()
  println("Diversity and Differential abundance testing parameters:")
  println("         --diff_abund_method [STRING] The method to use for differential abundance testing. Either ['all', 'ancombc1', 'ancombc2', or 'deseq2'] respectively. Default: 'all' ")
  println("         --rarefaction_depth [INTEGER] The Minimum desired sample rarefaction depth for diversity analysis. Default: 500.")
  println("         --group [STRING] Column in input csv file with treatments to be compared. Default: 'groups' ")
  println("         --samples_column [STRING] Column in input csv file with sample names belonging to each treatment group. Default: 'sample_id' ")
  println("         --remove_struc_zeros [BOOLEAN] Should structural zeros (a.k.a ASVs with zeros count in at least one group) be removed? default is false i.e. structural zeros will be retained. Options are true or false. Default: false.")
  println("         --remove_rare [BOOLEAN] Should rare features be filtered out prior to analysis? If true, rare features will be removed. Options are true or false. Default: false.")
  println("         --prevalence_cutoff [FLOAT] If --remove_rare is true, a numerical fraction between 0 and 1. Taxa with prevalences(the proportion of samples in which the taxon is present) less than this will be excluded from diversity and differential abundance analysis. Default is 0 , i.e. do not exclude any taxa. For example, to exclude taxa that are not present in at least 15% of the samples set it to 0.15.")
  println("         --library_cutoff [INTEGER] If --remove-rare is true, a numerical threshold for filtering samples based on library sizes. Samples with library sizes less than this number will be excluded in the analysis. Default is 0 i.e do not remove any sample. For example, if you want to discard samples with library sizes less than 100, then set to 100.")
  println()
  println("Genelab specific arguements:")
  println("      --accession [STRING]  A Genelab accession number if the --input_file parameter is not set. If this parameter is set, it will ignore the --input_file parameter.")
  println("      --assay_suffix [STRING]  Genelabs assay suffix. Default: _GLAmpSeq.")
  println("      --output_prefix [STRING] Unique name to tag onto output files. Automatically appends '_' if not empty and does not end with '_' or '-'. Default: empty string.")
  println()
  println("Paths to existing conda environments to use otherwise a new one will be created using the yaml file in envs/")
  println("      --conda_fastqc [PATH] Path to a conda environment containing fastqc. Default: null.")
  println("      --conda_multiqc [PATH] Path to a conda environment containing multiqc. Default: null.")
  println("      --conda_zip [PATH] Path to a conda environment containing zip. Default: null.")
  println("      --conda_R [PATH] Path to a conda environment containing R along with the packages decipher and biomformat installed. Default: null.")
  println("      --conda_dp_tools  [PATH] Path to a conda environment containing dp_tools. Default: null.")
  println("      --conda_cutadapt [PATH] Path to a conda environment containing cutadapt. Default: null.")
  println("      --conda_diversity [PATH] Path to a conda environment containing R packages required for diversity and differential abundance testing (ANCOMBC and DESeq2). Default: null.")
  println()
  print("Advanced users can edit the nextflow.config file for more control over default settings such container choice, number of cpus, memory per task etc.")
  exit 0
  }


if(params.debug){
log.info """${c_blue}
         Nextflow AmpIllumina Consensus Pipeline: $workflow.manifest.version
         
         You have set the following parameters:
         Input csv file : ${params.input_file}
         GLDS or OSD accession : ${params.accession}
         Amplicon target region : ${params.target_region}
         Nextflow Directory publishing mode: ${params.publish_dir_mode}
         Trim Primers: ${params.trim_primers}
         Nextflow Error strategy: ${params.errorStrategy}
         MultiQC configuration file: ${params.multiqc_config}

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
         
 
         Genelab Assay Suffix: ${params.assay_suffix}
         Output Prefix: ${params.output_prefix}

         Conda Environments:
         fastqc: ${params.conda_fastqc}
         multiqc: ${params.conda_multiqc}
         zip: ${params.conda_zip}
         R: ${params.conda_R}
         dp_tools: ${params.conda_dp_tools}
         cutadapt: ${params.conda_cutadapt}
         Diversity and Differential abundance : ${params.conda_diversity}
         ${c_reset}"""
}

// Create GLDS runsheet
include { GET_RUNSHEET } from "./modules/create_runsheet.nf"

// Stage raw reads
include { COPY_READS } from './modules/copy_reads.nf'
include { COPY_REMOTE_READS } from './modules/copy_reads.nf'

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

    
   // Capture software versions
   software_versions_ch = channel.empty()

   if(params.accession){

       values = channel.of([params.accession, params.target_region])

       GET_RUNSHEET(values)
       GET_RUNSHEET.out.input_file
           .splitCsv(header:true)
           .set{file_ch}

       target_region = GET_RUNSHEET.out.runsheet
                           .splitCsv(header:true)
                           .map{row -> "${row.'Parameter Value[Library Selection]'}"}.first()
       primers_ch = GET_RUNSHEET.out.runsheet
                           .splitCsv(header:true)
                           .map{row -> ["${row.F_Primer}", "${row.R_Primer}"]}
                           .first()                 

      GET_RUNSHEET.out.version | mix(software_versions_ch) | set{software_versions_ch}

      file_ch.map{
            row -> deleteWS(row.paired)  == 'true' ? tuple( "${row.sample_id}", ["${row.forward}","${row.reverse}"], deleteWS(row.paired)) : 
                                tuple( "${row.sample_id}", [("${row.forward}")], deleteWS(row.paired))
            }.set{reads_ch} 


   }else{

        channel.fromPath(params.input_file, checkIfExists: true)
           .splitCsv(header:true)
           .set{file_ch}
           file_ch.map{
                row -> deleteWS(row.paired)  == 'true' ? tuple( "${row.sample_id}", [file("${row.forward}"), file("${row.reverse}")], deleteWS(row.paired)) : 
                                    tuple( "${row.sample_id}", [file("${row.forward}")], deleteWS(row.paired))
                }.set{reads_ch} 
   }
    

    // Use original runsheet to preserve sample order
    if(params.accession){
        runsheet_ch = GET_RUNSHEET.out.runsheet
        isa_archive_ch = GET_RUNSHEET.out.zip
        gl_file_ch = GET_RUNSHEET.out.input_file
    }else{
        runsheet_ch = channel.fromPath(params.input_file, checkIfExists: true)
        isa_archive_ch = channel.empty()
        gl_file_ch = channel.empty()
    }

    // Stage raw reads with standard naming
    if(params.accession){
        COPY_REMOTE_READS(reads_ch)
        staged_reads_ch = COPY_REMOTE_READS.out.raw_reads
        
    }else{
        COPY_READS(reads_ch)
        staged_reads_ch = COPY_READS.out.raw_reads
    }


    // Read quality check and trimming
    RAW_FASTQC(staged_reads_ch)
    raw_fastqc_files = RAW_FASTQC.out.fastqc.flatten().collect()
    
    RAW_MULTIQC("raw", params.multiqc_config,raw_fastqc_files)

    RAW_FASTQC.out.version | mix(software_versions_ch) | set{software_versions_ch}
    RAW_MULTIQC.out.version | mix(software_versions_ch) | set{software_versions_ch}

    // Download reference database for taxonomic classification
    DOWNLOAD_DATABASE(params.target_region)

    trimmed_reads_ch = channel.empty()
    trimmed_reads_counts = channel.empty()
    cutadapt_logs = channel.empty()
    if(params.trim_primers){

        if(!params.accession) primers_ch = channel.value([params.F_primer, params.R_primer])
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

    meta  = channel.of(values)
    
    metadata  =  GET_RUNSHEET.out.runsheet

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

    meta  = channel.of(values)
    
    metadata  =  channel.fromPath(params.input_file, checkIfExists: true)

    }
    
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
    runsheet {
        path "Metadata"
    }

    isa_archive {
        path "Metadata"
    }

    gl_file {
        path "Metadata"
    }

    // Raw reads
    raw_reads {
        path "Raw_Sequence_Data"
    }

    // Trimmed reads
    trimmed_reads {
        path "Trimmed_Sequence_Data"
    }

    trimmed_count {
        path "Trimmed_Sequence_Data"
    }

    cutadapt_logs {
        path "Trimmed_Sequence_Data"
    }
    
    // Filtered reads
    filtered_reads {
        path "Filtered_Sequence_Data"
    }

    filtered_count {
        path "Filtered_Sequence_Data"
    }

    // FastQC
    raw_fastqc {
        path {html, zip -> "Raw_Sequence_Data/FastQC_Outputs" }
    }

    filtered_fastqc {
        path {html, zip -> "Filtered_Sequence_Data/FastQC_Outputs" }
    }

    // MultiQC
    zip_multiqc_raw {
        path "FastQC_Outputs"
    }

    html_multiqc_raw {
        path "FastQC_Outputs"
    }

    zip_multiqc_filtered {
        path "FastQC_Outputs"
    }

    html_multiqc_filtered {
        path "FastQC_Outputs"
    }

    // Dada2 outputs
    asv {
        path "Final_Outputs"
    }

    counts {
        path "Final_Outputs"
    }

    taxonomy {
        path "Final_Outputs"
    }

    taxonomy_counts {
        path "Final_Outputs"
    }

    biom_zip {
        path "Final_Outputs"
    }

    read_count_tracking {
        path "Final_Outputs"
    }

    // Alpha and beta diversity outputs
    alpha_diversity {
        path "Final_Outputs"
    }

    zip_alpha_plots {
        path "Final_Outputs/alpha_diversity"
    }

    beta_diversity {
        path "Final_Outputs"
    }

    zip_beta_euclidean_plots {
        path "Final_Outputs/beta_diversity"
    }

    zip_beta_bray_plots {
        path "Final_Outputs/beta_diversity"
    }

    // Taxonomy plots
    taxonomy_plots {
        path "Final_Outputs"
    }
    zip_taxonomy_samples {
        path "Final_Outputs/taxonomy_plots"
    }
    zip_taxonomy_groups {
        path "Final_Outputs/taxonomy_plots"
    }

    // Differential abundance outputs
    da_contrasts {
        path "Final_Outputs"
    }

    da_sampleTable {
        path "Final_Outputs"
    }

    ancombc1 {
        path "Final_Outputs"
    }

    zip_ancombc1 {
        path "Final_Outputs/differential_abundance/ancombc1"
    }

    ancombc2 {
        path "Final_Outputs"
    }

    zip_ancombc2 {
        path "Final_Outputs/differential_abundance/ancombc2"
    }

    deseq2 {
        path "Final_Outputs"
    }

    zip_deseq2 {
        path "Final_Outputs/differential_abundance/deseq2"
    }

    // GeneLab
    software_versions {
        path "GeneLab"
    }
}