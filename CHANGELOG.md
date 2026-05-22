# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.10](https://github.com/nasa/GeneLab_AmpliconSeq_Workflow/tree/NF_AmpIllumina_1.0.10)

### Added

- Caching support for reference databases to avoid redundant downloads on repeated runs of the main workflow
- A Bioconductor package (microbiome 1.32.0) that is required by ANCOMBC >= 2.12
- Support for `--primers-already-trimmed` behavior in the post-processing workflow by adding `trim_primers` parameter to control inclusion of trimmed-data-related information across its processes
- Validation for alpha/beta diversity outputs, including failure files and rarefaction depth consistency across both analyses, for taxonomy plots, and for outputs of the three differential abundance methods, including per-method volcano plot counts against number of contrasts
- `post_processing_schema.json` for parameter validation and automated help menu generation


### Changed
- Migrated output handling to Nextflow `workflow output {}` block in both main and post-processing workflows; removed deprecated suffix and directory parameters from `nextflow.config` and `post_processing.config`
- Switched output files under `Metadata/` and `GeneLab/` to comply with other GeneLab workflows
- Replaced `GET_RUNSHEET` process and associated workflow logic with a new staging analysis subworkflow supporting both accession-based and input-file-based execution modes
- Consolidated `CUTADAPT_ANCHORED` and `CUTADAPT_UNANCHORED` into a single `CUTADAPT` process
- Consolidated `RUN_R_TRIM` and `RUN_R_NOTRIM` into a single `RUN_DADA2` process
- Updated `CUTADAPT` and `RUN_DADA2` to accurately support both true single-end datasets and forced single-end processing in addition to the existing paired-end processing
- Updated `GENERATE_README`, `GENERATE_MD5SUMS`, `VALIDATE_PROCESSING`, and `GENERATE_CURATION_TABLE` to account for the `force_single_end` parameter
- Consolidated `ZIP_BIOM` into the general `ZIP` process using file extension detection for biom-specific flags
- Split MultiQC zip output into separate `data.zip` and `report.html` outputs
- Moved MultiQC path cleaning into the `MULTIQC` process via `clean_multiqc_paths.py`; removed standalone `MULTIQC_ZIP` process in main workflow and `CLEAN_MULTIQC_PATHS` process in post-processing workflow
- Moved `MultiQC_Reports/` from under `FastQC_Outputs/` to under their respective folders (`Raw_Sequence_Data` or `Filtered_Sequence_Data`); removed `FastQC_Outputs/`
- Converted `generate_protocol.sh` to a Python script for automated handling of reference/filtering parameters and special cases (primers already trimmed, processing only forward or reverse reads)
- Replaced bash code that generates md5sums file with `generate_md5sums.py` that provides md5sums for files published on OSDR only and optionally generates md5sums for raw data (FASTQ and MultiQC files)
- Moved post-processing output files from under `Post_Processing/` to under `GeneLab/`; removed `Post_Processing/`
- Modified `nextflow_schema.json` to reflect updated pipeline parameters with improved parameter descriptions and group definitions 
- Wired up `nf-schema` plugin to `nextflow_schema.json` for parameter validation and automated help menu generation
- Updated `post_processing.config` to support running the post-processing workflow with `-C` instead of `-c` (i.e., disregards the default `nextflow.config` and uses only `post_processing.config` ); modified launch scripts (`launch.sh` and `slurm_submit.slurm`) accordingly 

- Updated database versions:
  - SILVA_SSU_r138.2_v2, downloaded from DECIPHER and uploaded to Figshare (https://doi.org/10.6084/m9.figshare.32118616) in April 2026 
  - UNITE_v2025, downloaded from DECIPHER and uploaded to Figshare (https://doi.org/10.6084/m9.figshare.32118652) in April 2026

- Updated versions of Cutadapt and certain R/Bioconductor packages to match GeneLab Amplicon Pipeline version [GL-DPPD-7104-D](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Amplicon/Illumina/Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-D.md)
- Replaced tidyverse package with its individually-used component packages across R scripts to ensure specific versions of these packages are being used:
  - ggplot2 4.0.2
  - dplyr 1.2.1
  - readr 2.2.0
  - stringr 1.6.0
  - purrr 1.2.1
  - tibble 3.3.1
  - tidyr 1.3.2
- Removed explicit `library(tools)` loading from R scripts and its respective version capture, as it is automatically loaded with R base
- Bumped `r-diversity` and `r-dada-decipher-biomformat` singularity images from 1.1 to 1.2 to match R package updates 
- Updated third-party licenses

### Fixed
- Fixed `assay_suffix` chained parameter evaluation, which was previously evaluated at parse time instead of runtime

### Removed
- Removed `run_workflow.py` script as it was deemed unnecessary for the Nextflow workflow setup; cleaned up related commands in `slurm_submit.slurm`

<br>

---

## [1.0.9](https://github.com/nasa/GeneLab_AmpliconSeq_Workflow/tree/DEV_1.0.9)

### Fixed

- Fix NCBI taxon lookup issue by wrapping `get_ncbi_ids()` in per-taxon tryCatch with retries and returning NA on failure

<br> 

---


## [1.0.8](https://github.com/nasa/GeneLab_AmpliconSeq_Workflow/tree/NF_AmpIllumina_1.0.8)

### Fixed

- Fixed broken figshare download links by switching from the ndownloader URL to the API endpoint

<br>

---

## [1.0.7](https://github.com/nasa/GeneLab_AmpliconSeq_Workflow/tree/NF_AmpIllumina_1.0.7)

### Changed

- Increased ncbi_sleep from 0.8 to 2 in ANCOMBC1, ANCOMBC2, DESeq2 scripts

<br>

---

## [1.0.6](https://github.com/nasa/GeneLab_AmpliconSeq_Workflow/tree/NF_AmpIllumina_1.0.6)

### Added

- `dpt-isa-to-runsheet` now uses local plugin directories
- Versions of "hexbin" and "vsn" are now in software_versions.txt

### Changed

- Moved reference taxonomy database download from DADA R script into a dedicated Nextflow process (DOWNLOAD_DATABASE)
- Updated filtered FastQC/MultiQC steps to use filtered reads instead of trimmed reads
- Removed `rarefaction_depth.txt` from README, processed_md5sum, and curation table
- Excluded files not intended for OSDR from processed_md5sum
- Revised curation table to match agreed-upon format with the curation team
- Updated protocol text to remove unsupported characters
- Updated differential abundance (ANCOMBC1, ANCOMBC2, DESEQ2) output table column headers to use hyphens instead of dots, sanitize sample names before saving output tables. 

<br>

---

## [1.0.5](https://github.com/nasa/GeneLab_AmpliconSeq_Workflow/tree/NF_AmpIllumina_1.0.5)

### Added

- Added retries to remote raw read staging: for `--accession` runs, raw reads are downloaded with `wget` instead of being staged in Nextflow as paths. 

<br>

---

## [1.0.4](https://github.com/nasa/GeneLab_AmpliconSeq_Workflow/tree/NF_AmpIllumina_1.0.4)

### Changed

- Added dynamic legend sizing in taxonomy plots to prevent long taxa names from reducing plot area
- Added retries and load validation around Figshare reference database downloads in R processing scripts

<br>

---

## [1.0.3](https://github.com/nasa/GeneLab_AmpliconSeq_Workflow/tree/NF_AmpIllumina_1.0.3)

### Added

- Zipping plots generated by subdirectories of Final_outputs

### Changed

- Added output_prefix to README, processed md5sum, and processing_info.zip
- Used assay_suffix param instead of hardcoded value in config files
- Fixed README typo and pass expected output filename
- Updated file associations table to include outputs of Final_outputs subdirectories
- Updated README file and processed md5sum file to replace individual plot files with their respective zipped files and reflect OSDR's Files tree structure


### Fixed

- Fixed regex pattern in ANCOMBC scripts so factor conditions with underscores are parsed correctly
- Fixed runsheet handling so params.input_file paths with spaces no longer break processes

<br>

---

## [1.0.2](https://github.com/nasa/GeneLab_AmpliconSeq_Workflow/tree/NF_AmpIllumina_1.0.2)


### Changed

- Beta diversity script now uses the same rarefaction depth logic as the alpha diversity script; removed separate rarefaction adjustment.
- Added "GLAmpSeq" suffixes to all output files
- Added 3 retry attempts for OSDR file downloads


<br>

---

## [1.0.1](https://github.com/nasa/GeneLab_AmpliconSeq_Workflow/tree/NF_AmpIllumina_1.0.1)

### Added

- Rarefaction depth checks to alpha diversity with warning/error based on depth levels
- File output for rarefaction depth values for post-processing protocol use
- Raw reads are now staged in main.nf and renamed with COPY_READS using sample names and raw read suffix parameters

### Changed

- Updated FastQC container source to biocontainers
- Synchronized license table format with RNAseq license table
- Fixed relative paths of output directories and Nextflow-generated reports (timeline, report, trace)
- Fixed relative paths of input/output directories for post-processing workflow
- Removed "samples" parameter ("unique-sample-IDs.txt") from post-processing workflow, now uses "GLfile.csv"
- File associations table generation now uses sample names from assay table 
- An underscore ("\_") is now appended to the output prefix if params.output_prefix is not empty and does not end with "\_" or "-"
- Updated README generation script and linked it to workflow
- Updated protocol generation script with rarefaction depth checks

### Fixed

- Fixed issue where alpha diversity fails after rarefaction plot
- Fixed issue with diversity.df mangled sample name characters caused by estimate_richness()

### Removed

- Removed `GLparams_file.csv`, primer and target region info is now read from the runsheet

<br>

---

## [1.0.0](https://github.com/nasa/GeneLab_AmpliconSeq_Workflow/tree/NF_AmpIllumina_1.0.0) 

### Added

- Group-wise and Sample-wise Taxonomic Summary Plots
- Observed and Simpson alpha diversity metrics
- Bray-Curtis beta diversity metrics
- ANCOMBC1 differential abundance analysis
- ANCOMBC2 differential abundance analysis
- Added workflow parameters for finer-grained control of diversity and differential abundance analyses
- Improved error handling and robustness in R scripts, plots
- Software version tracking
- Added software versions:
  - ANCOMBC 2.8.0
  - broom 1.0.7
  - DescTools 0.99.59
  - dp_tools 1.3.8
  - FSA 0.9.6
  - ggdendro 0.2.0
  - glue 1.8.0
  - hexbin 1.28.3
  - mia 1.14.0
  - taxize 0.10.0
  - vsn 3.74.0
- Added reference database support for PR2 v4.13
- Added persistent reference links to DECIPHER databases on Figshare

### Changed

- Workflow converted from Snakemake to Nextflow
- Updated software versions:
  - MultiQC 1.27.1
  - Cutadapt 5.0
  - R-base 4.4.2
  - DADA2 1.34.0
  - DECIPHER 3.2.0
  - biomformat 1.34.0
  - DESeq2 1.46.0
  - ggrepel 0.9.6
  - phyloseq 1.50.0
- Updated database versions:
  - SILVA SSU r138_2 2024
  - UNITE v2024 April2024

<br>

---

> ***Note:** All previous workflow changes associated with the previous version of the GeneLab Amplicon Pipeline
[GL-DPPD-7104-B](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Amplicon/Illumina/Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-B.md) and can be found in the
[change log of the Snakemake workflow (SW_AmpIllumina-B)](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Amplicon/Illumina/Workflow_Documentation/SW_AmpIllumina-B/CHANGELOG.md).*
