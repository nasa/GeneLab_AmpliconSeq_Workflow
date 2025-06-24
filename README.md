# GeneLab Amplicon Illumina Sequencing Data Processing Workflow

> GeneLab, part of NASA's [Open Science Data Repository (OSDR)](https://www.nasa.gov/osdr/), has wrapped each step of the Illumina amplicon sequencing data processing pipeline ([AmpIllumina](https://github.com/nasa/GeneLab_Data_Processing/tree/master/Amplicon)), starting with pipeline version C, into a Nextflow workflow. This repository contains information about the workflow along with instructions for installation and usage. Exact workflow run info and AmpIllumina version used to process specific datasets hosted on the [OSDR data repository](https://osdr.nasa.gov/bio/repo/) are provided alongside their processed data in OSDR under 'Files' -> 'GeneLab Processed Diversity Amplicon Files' -> 'Processing Info'. 

<br>

## General Workflow Info

### Implementation Tools

The current GeneLab Illumina amplicon sequencing data processing pipeline (AmpIllumina), [GL-DPPD-7104-C.md](../../Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-C.md), is implemented as a [Nextflow](https://nextflow.io/) DSL2 workflow and utilizes [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/introduction.html) containers, [Docker](https://docs.docker.com/get-started/) containers, or [conda](https://docs.conda.io/en/latest/) environments to install/run all tools. This workflow is run using the command line interface (CLI) of any unix-based system.  While knowledge of creating workflows in Nextflow is not required to run the workflow as is, [the Nextflow documentation](https://nextflow.io/docs/latest/index.html) is a useful resource for users who want to modify and/or extend this workflow.   

### Resource Requirements <!-- omit in toc -->

The table below details the default maximum resource allocations for individual Nextflow processes.

| CPU Cores | Memory |
|-------------------|----------------|
| 2                 | 5 GB           |

> **Note:** These per-process resource allocations are defaults. They can be adjusted by modifying `cpus` and `memory` directives in the configuration file: [`nextflow.config`](workflow_code/nextflow.config).

<br>

## Utilizing the Workflow

1. [Install Nextflow, Singularity, and Conda](#1-install-nextflow-singularity-and-conda)  
   1a. [Install Nextflow and Conda](#1a-install-nextflow-and-conda)  
   1b. [Install Singularity](#1b-install-singularity)  

2. [Download the Workflow Files](#2-download-the-workflow-files)  

3. [Fetch Singularity Images](#3-fetch-singularity-images)  

4. [Run the Workflow](#4-run-the-workflow)  
   4a. [Approach 1: Start with OSD or GLDS accession as input](#4a-approach-1-start-with-an-osd-or-glds-accession-as-input)  
   4b. [Approach 2: Start with a runsheet csv file as input](#4b-approach-2-start-with-a-runsheet-csv-file-as-input)  
   4c. [Modify parameters and compute resources in the Nextflow config file](#4c-modify-parameters-and-compute-resources-in-the-nextflow-config-file)  

5. [Workflow Outputs](#5-workflow-outputs)  
   5a. [Main outputs](#5a-main-outputs)  
   5b. [Resource logs](#5b-resource-logs)  

6. [Post-processing](#6-post-processing)  

<br>

---

### 1. Install Nextflow, Singularity, and Conda

#### 1a. Install Nextflow and Conda

Nextflow can be installed either through [Anaconda](https://anaconda.org/bioconda/nextflow) or as documented on the [Nextflow documentation page](https://www.nextflow.io/docs/latest/getstarted.html).

> Note: If you want to install Anaconda, we recommend installing a Miniconda version appropriate for your system, as instructed by [Happy Belly Bioinformatics](https://astrobiomike.github.io/unix/conda-intro#getting-and-installing-conda).  
> 
> Once conda is installed on your system, you can install the latest version of Nextflow by running the following commands:
> 
> ```bash
> conda install -c bioconda nextflow
> nextflow self-update
> ```
> You may also install [mamba](https://mamba.readthedocs.io/en/latest/index.html) first which is a faster implementation of conda and can be used as a drop-in replacement:
> ```bash
> conda install -c conda-forge mamba
> ```

<br>

#### 1b. Install Singularity

Singularity is a container platform that allows usage of containerized software. This enables the GeneLab workflow to retrieve and use all software required for processing without the need to install the software directly on the user's system.

We recommend installing Singularity on a system wide level as per the associated [documentation](https://docs.sylabs.io/guides/3.10/admin-guide/admin_quickstart.html).

> Note: Singularity is also available through [Anaconda](https://anaconda.org/conda-forge/singularity).

> Note: Alternatively, Docker can be used in place of Singularity. To get started with Docker, see the [Docker CE installation documentation](https://docs.docker.com/engine/install/).

<br>

---

### 2. Download the Workflow Files

All files required for utilizing the NF_AmpIllumina GeneLab workflow for processing amplicon Illumina data are in the [workflow_code](workflow_code) directory. To get a copy of the latest *NF_AmpIllumina* version on to your system, the code can be downloaded as a zip file from the release page then unzipped after downloading by running the following commands: 

```bash
wget https://github.com/nasa/GeneLab_AmpliconSeq_Workflow/releases/download/NF_AmpIllumina_1.0.0/NF_AmpIllumina_1.0.0.zip
unzip NF_AmpIllumina_1.0.0.zip && cd NF_AmpIllumina_1.0.0
```

<br>

---

### 3. Fetch Singularity Images

Although Nextflow can fetch Singularity images from a url, doing so may cause issues as detailed [here](https://github.com/nextflow-io/nextflow/issues/1210).

To avoid this issue, run the following command to fetch the Singularity images prior to running the NF_AmpIllumina workflow:

> Note: This command should be run in the location containing the `NF_AmpIllumina_1.0.0` directory that was downloaded in [step 2](#2-download-the-workflow-files) above. Depending on your network speed, fetching the images will take ~20 minutes. Approximately 4GB of RAM is needed to download and build the Singularity images.

```bash
bash ./bin/prepull_singularity.sh nextflow.config
```

Once complete, a `singularity` folder containing the Singularity images will be created. Run the following command to export this folder as a Nextflow configuration environment variable to ensure Nextflow can locate the fetched images:

```bash
export NXF_SINGULARITY_CACHEDIR=$(pwd)/singularity
```

<br>

---

### 4. Run the Workflow

> ***Note:** All the commands in this step must be run from within the `NF_AmpIllumina_1.0.0` directory that was downloaded in [step 2](#2-download-the-workflow-files) above.*

For options and detailed help on how to run the workflow, run the following command:

```bash
nextflow run main.nf --help
```

> Note: Nextflow commands use both single hyphen arguments (e.g. -help) that denote general Nextflow arguments and double hyphen arguments (e.g. --input_file) that denote workflow specific parameters.  Take care to use the proper number of hyphens for each argument.

<br>

#### 4a. Approach 1: Start with OSD or GLDS accession as input

```bash
nextflow run main.nf \
   -resume \
   -profile singularity \
   --target_region 16S \
   --accession OSD-487 
```

<br>

#### 4b. Approach 2: Start with a runsheet csv file as input

```bash
nextflow run main.nf \
   -resume \
   -profile singularity \
   --target_region 16S \
   --input_file PE_file.csv \
   --F_primer AGAGTTTGATCCTGGCTCAG \
   --R_primer CTGCCTCCCGTAGGAGT 
```

<br>


**Required Parameters For All Approaches:**

* `main.nf` - Instructs Nextflow to run the NF_AmpIllumina workflow 
* `-resume` - Resumes  workflow execution using previously cached results
* `-profile` – Specifies the configuration profile(s) to load (multiple options can be provided as a comma-separated list)
   * Software environment profile options (choose one):
      * `singularity` - instructs Nextflow to use Singularity container environments
      * `docker` - instructs Nextflow to use Docker container environments
      * `conda` - instructs Nextflow to use conda environments via the conda package manager
        > *Note: By default, Nextflow will create environments at runtime using the yaml files in the [workflow_code/envs](workflow_code/envs/) folder. You can change this behavior by using the `--conda_*` workflow parameters or by editing the [nextflow.config](workflow_code/nextflow.config) file to specify a centralized conda environments directory via the `conda.cacheDir` parameter.*
      * `mamba` - instructs Nextflow to use conda environments via the mamba package manager 
   * Other option (can be combined with the software environment option above using a comma, e.g. `-profile slurm,singularity`):
      * `slurm` - instructs Nextflow to use the [Slurm cluster management and job scheduling system](https://slurm.schedmd.com/overview.html) to schedule and run the jobs on a Slurm HPC cluster
* `--target_region` - Specifies the amplicon target region to be analyzed (type: string, options: 16S, 18S, or ITS)


**Additional Required Parameters For Approach 1:** 

* `--accession` - The OSD or GLDS accession number specifying the [OSDR](https://osdr.nasa.gov/bio/repo/) dataset to process, e.g. OSD-487 or GLDS-487 (type: string)
  > *Note: Not all datasets have the same OSD and GLDS number, so make sure the correct OSD or GLDS number is specified*


**Additional Required Parameters For Approach 2:** 

* `--input_file` –  A single-end or paired-end runsheet csv file containing assay metadata for each sample, including sample_id, forward (path to forward read), [reverse (path to reverse read, for paired-end only),] paired (boolean, TRUE | FALSE), groups (specifies sample treatment group name). Please see the [runsheet documentation](./examples/runsheet) in this repository for examples on how to format this file. (type: string)

* `--F_primer` - Forward primer sequence (type: string)

* `--R_primer` - Reverse primer sequence (type: string)


**Additional [Optional] Parameters For All Approaches**
> *Note: See `nextflow run -h` and [Nextflow's CLI run command documentation](https://nextflow.io/docs/latest/cli.html#run) for more options and details on how to run Nextflow.*

* `--errorStrategy` - Error handling strategy for Nextflow processes. If processes fail, use "ignore" to allow the workflow to continue running (type: string, default: "terminate")
* `--trim_primers` - Whether primers should be trimmed (type: string, default: "TRUE")
* `--primers_linked` - Whether forward and reverse primers are linked (type: string, default: "TRUE")
* `--anchored_primers` - Whether primers are anchored at the start of reads (type: string, default: "TRUE")
* `--min_cutadapt_len` - Minimum length of reads to keep after Cutadapt trimming (type: integer, default: 130) 
* `--discard_untrimmed` - Whether to discard untrimmed reads (type: string, default: "TRUE")
* `--left_trunc` - Truncate forward reads after this many bases. Reads shorter than this are discarded (type: integer, default: 0)
* `--right_trunc` - Truncate reverse reads after this many bases. Reads shorter than this are discarded (type: integer, default: 0)
* `--left_maxEE` - Maximum expected errors allowed in forward reads (type: integer, default: 1)
* `--right_maxEE` - Maximum expected errors allowed in reverse reads (type: integer, default: 1)
* `--concatenate_reads_only` - Whether to concatenate paired reads end-to-end instead of merging based on overlapping regions (type: string, default: "FALSE")
* `--rarefaction_depth` - The minimum desired sample rarefaction depth for beta diversity analysis (type: integer, default: 500)
* `--diff_abund_method` - Differential abundance testing method to use (type: string, default: "all")
* `--group` - Column name in input CSV file containing groups to be compared (type: string, default: "groups")
* `--samples_column` - Column name in input CSV file containing sample names (type: string, default: "sample_id")
* `--remove_struc_zeros` - Whether to remove structural zeros when running ANCOMBC (type: boolean, default: false)
* `--remove_rare` - Whether to filter out rare features and samples with low library sizes (type: boolean, default: false)
* `--prevalence_cutoff` - Taxa with prevalence below this fraction will be excluded (type: float, default: 0)
* `--library_cutoff` - Samples with library sizes below this threshold will be excluded (type: integer, default: 0)
* `--output_prefix` - Prefix to add to output filenames, e.g. "Study1_" (type: string, default: "")
* `--assay_suffix` - Suffix to add to output filenames (type: string, default: "_GLAmpSeq")
* `--use_conda` - Whether Conda environments should be used to run Nextflow processes (type: boolean, default: false)
* `--conda_cutadapt` - Path to existing Cutadapt conda environment (type: string, default: null)
* `--conda_diversity` - Path to existing R diversity analysis conda environment (type: string, default: null)
* `--conda_dp_tools` - Path to existing dp_tools conda environment (type: string, default: null)
* `--conda_fastqc` - Path to existing FastQC conda environment (type: string, default: null)
* `--conda_multiqc` - Path to existing MultiQC conda environment (type: string, default: null)
* `--conda_R` - Path to existing R conda environment (type: string, default: null)
* `--conda_zip` - Path to existing zip conda environment (type: string, default: null)

<br>

#### 4c. Modify parameters and compute resources in the Nextflow config file

Additionally, all parameters and workflow resources can be directly specified in the [nextflow.config](./workflow_code/nextflow.config) file. For detailed instructions on how to modify and set parameters in the config file, please see the [documentation here](https://www.nextflow.io/docs/latest/config.html).

Once you've downloaded the workflow template, you can modify the parameters in the `params` scope and cpus/memory requirements in the `process` scope in your downloaded version of the [nextflow.config](workflow_code/nextflow.config) file as needed in order to match your dataset and system setup. Additionally, if necessary, you can modify each variable in the [nextflow.config](workflow_code/nextflow.config) file to be consistent with the study you want to process and the computer you're using for processing.

<br>

---

### 5. Workflow Outputs

#### 5a. Main Outputs

The outputs from this pipeline are documented in the [GL-DPPD-7104-C](../../Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-C.md) processing protocol.

#### 5b. Resource Logs

Standard Nextflow resource usage logs are also produced as follows:

**Nextflow Resource Usage Logs**
   - Resource_Usage/execution_report_{timestamp}.html (an html report that includes metrics about the workflow execution including computational resources and exact workflow process commands)
   - Resource_Usage/execution_timeline_{timestamp}.html (an html timeline for all processes executed in the workflow)
   - Resource_Usage/execution_trace_{timestamp}.txt (an execution tracing file that contains information about each process executed in the workflow, including: submission time, start time, completion time, cpu and memory used, machine-readable output)

> Further details about these logs can also found within [this Nextflow documentation page](https://www.nextflow.io/docs/latest/tracing.html#execution-report).

<br>

---

### 6. Post-processing

> Please note that to run the post-processing workflow successfully, you MUST run the processing workflow above via the [launch.sh](workflow_code/launch.sh) script first. Please see the [script](workflow_code/launch.sh) for how to run it and make sure to edit the place holders before running it.

The post-processing workflow generates a README file, a protocols file, an md5sums table, and a file association table suitable for uploading to OSDR.

For options and detailed help on how to run the post-processing workflow, run the following command:

```bash
nextflow run post_processing.nf --help
```

To generate the post-processing files after running the main processing workflow successfully, modify and set the parameters in [post_processing.config](workflow_code/post_processing.config), then run the following command:

```bash
nextflow run post_processing.nf \
   -c post_processing.config \ 
   -resume \
   -profile singularity
``` 

The outputs of the post-processing workflow are described below:
> *Note: The outputs will be in a directory called `Post_Processing` by default*

**Post processing workflow output files** 
 - Post_processing/FastQC_Outputs/filtered_multiqc_GLAmpSeq_report.zip (Filtered sequence MultiQC report with paths purged) 
 - Post_processing/FastQC_Outputs/raw_multiqc_GLAmpSeq_report.zip (Raw sequence MultiQC report with paths purged)
 - Post_processing/<GLDS_accession>_-associated-file-names.tsv (File association table for OSDR curation)
 - Post_processing/<GLDS_accession>_amplicon-validation.log (Automated verification and validation log file)
 - Post_processing/processed_md5sum_GLAmpSeq.tsv (md5sums for the files published on OSDR)
 - Post_processing/processing_info_GLAmpSeq.zip  (Zip file containing all files used to run the workflow and required logs with paths purged) 
 - Post_processing/protocol.txt  (File describing the methods used by the workflow)
 - Post_processing/README_GLAmpSeq.txt (README file listing and describing the outputs of the workflow)

<br>

---

## License

The software for the Amplicon Seq workflow is released under the [NASA Open Source Agreement (NOSA) Version 1.3](License/Amplicon_NOSA_License.pdf).


### 3rd Party Software Licenses

Licenses for the 3rd party open source software utilized in the Amplicon Seq workflow can be found in the [3rd_Party_Licenses sub-directory](License/3rd_Party_Licenses/Amplicon_3rd_Party_Software.md). 

<br>

---

## Notices

Copyright © 2021 United States Government as represented by the Administrator of the National Aeronautics and Space Administration.  All Rights Reserved.

### Disclaimers

No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."

Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.
