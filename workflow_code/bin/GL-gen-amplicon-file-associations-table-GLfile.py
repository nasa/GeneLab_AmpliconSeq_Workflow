#!/usr/bin/env python

"""
This is a program for generating the file-associations table needed by Curation for newly processed datasets.
"""

import os
import sys
import argparse
import textwrap
import pandas as pd
import zipfile
import re

parser = argparse.ArgumentParser(description = "This program generates the file-assocations table needed by Curation for \
                                               newly processed amplicon datasets. It is intended to be run after `GL-validate-processed-data` \
                                               has been run successfully.")
required = parser.add_argument_group('required arguments')
required.add_argument("-g", "--GLDS-ID", help = 'GLDS ID (e.g. "GLDS-276")', action = "store", required = True)
parser.add_argument("--GL-file", 
                    help = 'A 3-column (single-end) or 4-column (paired-end) input file (sample_id, forward, [reverse,] paired) used to run the processing pipeline. This is the GLfile.csv in the GeneLab directory if --GLDS_accession was used as input. ', 
                    action = "store", default = "")
parser.add_argument("--runsheet",
                    help = """
                           Input csv runsheet file used to run nextflow. This argument must be set when running the workflow with an OSD/GLDS accession as input as opposed to passing an input csv file.
                           This argument is used to get raw input file names that are used to retrieve raw read depths per sample. 
                           """,
                    action = "store", default = "")
parser.add_argument("--type", help = 'Specify if ASVs or OTUs (default: "ASVs"; only relevant for Amplicon)', 
                    action = "store", choices = ["ASVs", "OTUs"], default = "ASVs")
parser.add_argument("--map", help = 'Mapping file if samples come from more than one primer set (only relevant for Amplicon; tab-delimited, first column holds sample IDs, \
                                    second column holds the filename prefix of the outputs specific to that sample)', action = "store")
parser.add_argument("--output", 
                    help = 'Name of output log file (default: "<GLDS-ID>_[<output_prefix>]-associated-file-names.tsv", with appended prefix if one is provided)',
                    default = "", action = "store")
parser.add_argument("-p", "--output-prefix", help = "Output additional file prefix if there is one", action = "store", default = "")
parser.add_argument("--additional-string-to-remove-from-unique-filenames",
                    help = "If there is any additional text to remove from unqiue filenames, it can be provided here.",
                    action = "store")
parser.add_argument("--assay_suffix", help = "Genelab assay suffix", action = "store", default = "_GLAmpSeq")
parser.add_argument("--raw_file_prefix", help = "Prefix to be added to the raw data files alone (Default: <GLDS_ID>_Amplicon_)", action = "store", default ="")
parser.add_argument("--file_prefix", help = "Prefix to be added to all files except the raw files (Default: <GLDS_ID>_GAmplicon_)", action = "store", default ="")
parser.add_argument("--raw_suffix", help = "Raw reads suffix", action = "store", default ="_raw.fastq.gz")
parser.add_argument("--raw_R1_suffix", help = "Raw forward reads suffix", action = "store", default = "_R1_raw.fastq.gz")
parser.add_argument("--raw_R2_suffix", help = "Raw reverse reads suffix", action = "store", default = "_R2_raw.fastq.gz")
parser.add_argument("--primer_trimmed_suffix", help = "Trimmed reads suffix", action = "store", default = "_trimmed.fastq.gz")
parser.add_argument("--primer_trimmed_R1_suffix", help = "Trimmed forward reads suffix", action = "store", default = "_R1_trimmed.fastq.gz")
parser.add_argument("--primer_trimmed_R2_suffix", help = "Trimmed reverse reads suffix", action = "store", default = "_R2_trimmed.fastq.gz")
parser.add_argument("--filtered_suffix", help = "Filtered reads suffix", action = "store", default = "_filtered.fastq.gz")
parser.add_argument("--filtered_R1_suffix", help = "Filtered forward reads suffix", action = "store", default = "_R1_filtered.fastq.gz")
parser.add_argument("--filtered_R2_suffix", help = "Filtered reverse reads suffix", action = "store", default = "_R2_filtered.fastq.gz")
parser.add_argument("--raw_reads_dir", help = "Specifies the name of the raw reads directory if they are to be included",
                    action = "store", default = "Raw_Sequence_Data/")
parser.add_argument("--fastqc_dir", help = "Specifies the name of fastqc and multiqc reports directory", 
                    action = "store", default = "FastQC_Outputs/")
parser.add_argument("--filtered_reads_dir", help = "Specifies the name of the filtered reads directory", 
                    action = "store", default = "Filtered_Sequence_Data/")
parser.add_argument("--trimmed_reads_dir", help = "Specifies the name of the trimmed reads directory", 
                    action = "store", default = "Trimmed_Sequence_Data/")
parser.add_argument("--final_outputs_dir", help = "Specifies the name of the final outputs directory.", 
                    action = "store", default = "Final_Outputs/")
parser.add_argument("--single-ended", help = "Add this flag if data are single-end sequencing.", action = "store_true")
parser.add_argument("--primers-already-trimmed", help = "Add this flag if primers were trimmed prior to GeneLab processing, \
                    therefore there are no trimmed sequence data", action = "store_true")
parser.add_argument("--R1-used-as-single-ended-data", help = "Provide this flag if processing only R1 reads as single-end (as the expected raw \
                    filename suffixes will have 'R1' in there)", 
                    action = "store_true")
parser.add_argument("--include-raw-multiqc-in-output",
                    help = "Provide this flag if wanting to include the raw multiqc zip in the file-associations output table (may be wanted for older datasets)", 
                    action = "store_true")
parser.add_argument("--use-sample-names-from-assay-table",
                    help = "Provide this flag if the unique filename strings in the processed outputs are taken directly from the \
                           'Sample Name' column of the input assay table.", action = "store_true")


if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_args()



# Setting some colors
tty_colors = {
    'green' : '\033[0;32m%s\033[0m',
    'yellow' : '\033[0;33m%s\033[0m',
    'red' : '\033[0;31m%s\033[0m'
}


######################### Aesthetic functions #########################
def color_text(text, color='green'):
    if sys.stdout.isatty():
        return tty_colors[color] % text
    else:
        return text


def wprint(text):
    """ Print wrapper """

    print(textwrap.fill(text, width=80, initial_indent="  ", 
          subsequent_indent="  ", break_on_hyphens=False))

#################### End of Aesthetic functions #########################


def report_failure(message, color = "yellow"):
    print("")
    wprint(color_text(message, color))
    print("\nCuration file-associations table generation failed.\n")
    sys.exit(1)


def check_for_file_and_contents(file_path):
    """Checks if file exists and that it is not empty"""
    if not os.path.exists(file_path):
        report_failure("The expected file '" + str(file_path) + "' does not exist.")
    if not os.path.getsize(file_path) > 0:
        report_failure("The file '" + str(file_path) + "' is empty.")


def get_assay_table_from_ISA(isa_zip):
    """ Tries to find a single assay table in an isa object """

    zip_file = zipfile.ZipFile(isa_zip)
    isa_files = zip_file.namelist()

    # Getting wanted filename (those that start with "a_"  and contain the word "amplicon" seem to be what we want)
    wanted_file_list = [item for item in isa_files if item.startswith("a_") and item.find("amplicon") != -1]
    if len(wanted_file_list) != 1:
        report_failure("We couldn't find the correct assay table in the ISA object, consider passing it directly to the '-a' argument.")

    wanted_file = wanted_file_list[0]

    df = pd.read_csv(zip_file.open(wanted_file), sep = "\t")

    return(df)



def parse_amplicon_names(name, raw_file_prefix, raw_R1_suffix, raw_R2_suffix, raw_suffix):
    """ This removes expected prefixes and suffixes """

    # Removing expected prefix
    curr_name = name.replace(raw_file_prefix, "")

    # Removing potential suffixes (also checking R2 in case they are not 
    # in the appropriate order in the sample table, e.g. R2 before R1)
    curr_name = curr_name.replace(raw_R1_suffix, "")
    curr_name = curr_name.replace(raw_R2_suffix, "")
    curr_name = curr_name.replace(raw_suffix, "")

    return(curr_name)

def get_sample_names_and_unique_filenames(assay_table,  raw_file_prefix, raw_R1_suffix,
                                          raw_R2_suffix, raw_suffix,
                                         use_sample_names_from_assay_table,
                                         additional_string_to_remove_from_unique_filenames):
    """
    This gets the sample names ('Sample Name' column) from the assay table,
    and tries to get what would have been the unique filename prefixes generated from
    what's in the Raw Data File column of the assay table.

    Unless the --use-sample-names-from-assay-table flag was provided, then it just uses what's
    in the 'Sample Name' column.
    """
    sample_names = assay_table["Sample Name"].tolist()
    
    if use_sample_names_from_assay_table:
        unique_filename_prefixes = sample_names
        return(sample_names, unique_filename_prefixes)

    all_filenames = assay_table["Raw Data File"]

    unique_filename_prefixes = []
    
    # Attempting to split if they have multiple files (like paired-end)
    # and also removing the common prefixes and suffixes intending to create the same 
    # unique filenames used for processing

    for entry in all_filenames:

        # splitting if there are more than one (like with paired-end)
        curr_name = entry.split(",")[0]

        curr_name = parse_amplicon_names(curr_name, raw_file_prefix, raw_R1_suffix, raw_R2_suffix, raw_suffix)

        unique_filename_prefixes.append(curr_name)


    if additional_string_to_remove_from_unique_filenames:

        unique_filename_prefixes = [x.replace(additional_string_to_remove_from_unique_filenames, "") for x in unique_filename_prefixes]

    return(sample_names, unique_filename_prefixes)


def get_read_counts_from_raw_multiqc(mapping_tab, raw_multiqc_stats_file_path,
                                      fastqc_dir, output_prefix,  raw_multiqc_zip):

    # These are in multiple files if there was a mapping input table
    if isinstance(mapping_tab, pd.DataFrame):
        unique_prefixes = mapping_tab.prefix.unique()

        # Initializing a list to hold all dataframes we'll read in
        df_list = []

        # Working through each one
        for prefix in unique_prefixes:

            curr_file_path = os.path.join(fastqc_dir, output_prefix + prefix + raw_multiqc_zip)

            # Making sure there are multiqc files for each unique prefix given in the mapping table
            check_for_file_and_contents(curr_file_path)

            # Reading in
            zip_file = zipfile.ZipFile(curr_file_path)
            #curr_df = pd.read_csv(zip_file.open(raw_multiqc_stats_file_path), sep = "\t", usecols = [0,6])
            curr_df = pd.read_csv(zip_file.open(raw_multiqc_stats_file_path), sep = "\t")
            curr_df = curr_df.iloc[:,[0,-1]] # retrieve the samples column[0] and  last column[-1] which is reads counts column

            # Test to see if the sequence counts are presented as raw counts
            # or as decimal representation of a million
            max_count = curr_df.iloc[:,1].max()
            if not isinstance(max_count, int):
                curr_df.iloc[:,-1] = curr_df.iloc[:,-1] * 1000000  # convert to raw counts

            curr_df.columns = ["sample", "counts"]
            curr_df.set_index("sample", inplace = True)

            # Adding to list
            df_list.append(curr_df)

        # Combining tables
        df = pd.concat(df_list, axis = 0)

        return(df)

    else:
        input_zip = os.path.join(fastqc_dir, output_prefix + raw_multiqc_zip)
        zip_file = zipfile.ZipFile(input_zip)
        #df = pd.read_csv(zip_file.open(raw_multiqc_stats_file_path), sep = "\t", usecols = [0,6])
        df = pd.read_csv(zip_file.open(raw_multiqc_stats_file_path), sep = "\t")
        df = df.iloc[:,[0,-1]] # retrieve the samples column[0] and  last column[-1] which is reads counts column
        # Test to see if the sequence counts are presented as raw counts
        # or as decimal representation of a million
        max_count = df.iloc[:,1].max()
        if not isinstance(max_count, int):
            df.iloc[:,-1] = df.iloc[:,-1] * 1000000  # convert to raw counts
        df.columns = ["sample", "counts"]
        df.set_index("sample", inplace = True)

        return(df)


def get_read_count_from_df(sample_name, read_counts_tab, 
                           raw_suffix, raw_R1_suffix, 
                           single_ended, sample_raw_prefix_dict):

    if sample_raw_prefix_dict != "":
        return(round(read_counts_tab.at[sample_raw_prefix_dict[sample_name], "counts"]))

    if single_ended:
        return(round(read_counts_tab.at[str(sample_name) + \
                     raw_suffix.replace("_raw.fastq.gz", ""), "counts"]))
    else:
        return(round(read_counts_tab.at[str(sample_name) + \
                    raw_R1_suffix.replace("_raw.fastq.gz", ""), "counts"]))



def write_colnames(raw_reads_dir, trimmed_reads_dir,
                  filtered_reads_dir, include_raw_multiqc_in_output):

    ## Builds as if primers were trimmed by the workflow (with Trimmed column),
    #  but that is removed later if
    ## --primers-already-trimmed argument was provided
    colnames = ["Sample Name", 
                f"Parameter Value[{raw_reads_dir}]",
                "Parameter Value[Read Depth]",
                "Unit",
                "Parameter Value[MultiQC File Names]",
                f"Parameter Value[{trimmed_reads_dir}]",
                f"Parameter Value[{filtered_reads_dir}]",
                f"Parameter Value[Filtered Sequence Data/MultiQC Reports]",
                "Parameter Value[Taxonomy and ASV Counts Data]",
                "Parameter Value[Alpha Diversity Data]",
                "Parameter Value[Beta Diversity Data]",
                "Parameter Value[Beta Diversity Data/Bray-Curtis]",
                "Parameter Value[Beta Diversity Data/Euclidean_distance]",
                "Parameter Value[Taxonomy Plots]",
                "Parameter Value[Differential Abundance]",
                "Parameter Value[Differential Abundance/ANCOMBC1]",
                "Parameter Value[Differential Abundance/ANCOMBC2]",
                "Parameter Value[Differential Abundance/DESeq2]"]
    
    if not include_raw_multiqc_in_output:
        colnames.remove("Parameter Value[MultiQC File Names]")

    return colnames


def create_constants(include_raw_multiqc_in_output, raw_multiqc_zip,
                     filtered_multiqc_zip, Type, combined_prefix, assay_suffix):
    """A function to create lists of contants to be in creating a file association table"""
    if include_raw_multiqc_in_output:
        fastqc = [combined_prefix + raw_multiqc_zip, 
                  combined_prefix + filtered_multiqc_zip]
    else:
        fastqc = [combined_prefix + filtered_multiqc_zip]

    if Type == "ASVs":
        rep_seq_output = combined_prefix + f"ASVs{assay_suffix}.fasta"
    else:
        rep_seq_output = combined_prefix + f"OTUs{assay_suffix}.fasta"

    final_outputs = [rep_seq_output, 
                     combined_prefix + f"counts{assay_suffix}.tsv", 
                     combined_prefix + f"read-count-tracking{assay_suffix}.tsv", 
                     combined_prefix + f"taxonomy-and-counts{assay_suffix}.biom.zip", 
                     combined_prefix + f"taxonomy-and-counts{assay_suffix}.tsv", 
                     combined_prefix + f"taxonomy{assay_suffix}.tsv"]
 
    return fastqc, final_outputs

def collect_final_outputs_columns(final_outputs_dir, file_prefix, output_prefix, assay_suffix):
    """
    Returns a dict with each new column as key and a comma-separated string of files as value
    """

    results = {}

    # Primary files
    primary_files = [
        f"{output_prefix}ASVs{assay_suffix}.fasta",
        f"{output_prefix}counts{assay_suffix}.tsv",
        f"{output_prefix}taxonomy{assay_suffix}.tsv",
        f"{output_prefix}taxonomy-and-counts{assay_suffix}.tsv",
        f"{output_prefix}taxonomy-and-counts{assay_suffix}.biom.zip",
        f"{output_prefix}read-count-tracking{assay_suffix}.tsv"
    ]
    primary_files = [file_prefix + f for f in primary_files if os.path.exists(os.path.join(final_outputs_dir, f))]
    results["Taxonomy and ASV Counts Data"] = ",".join(primary_files)
    
    # Alpha diversity
    alpha_dir = os.path.join(final_outputs_dir, "alpha_diversity")
    all_alpha_files = [file_prefix + f for f in sorted(os.listdir(alpha_dir))] if os.path.isdir(alpha_dir) else []
    if all_alpha_files:
        alpha_files = [f for f in all_alpha_files if ".png" not in f.lower() and "rarefaction_depth" not in f.lower()]
        results["Alpha Diversity Data"] = ", ".join(alpha_files)

    # Beta diversity (vsd, bray-curtis, euclidean directly under beta_diversity)
    beta_dir = os.path.join(final_outputs_dir, "beta_diversity")
    if os.path.isdir(beta_dir):
        beta_files = os.listdir(beta_dir)

        # top-level (vsd validation)
        top_files = [file_prefix + f for f in beta_files if "vsd" in f.lower()]
        if top_files:
            results["Beta Diversity Data"] = ", ".join(top_files)

        # bray-curtis group
        bray = sorted([file_prefix + f for f in beta_files if "bray" in f.lower() and ".png" not in f.lower()])
        failure = sorted([file_prefix + f for f in beta_files if "failure" in f.lower()])
        if bray:
            results["Beta Diversity Data/Bray-Curtis"] = ", ".join(bray)
        elif failure:
            results["Beta Diversity Data/Bray-Curtis"] = ", ".join(failure)
        else:
            results["Beta Diversity Data/Bray-Curtis"] = ""
            print("Beta diversity with rarefaction failed without generating a failure text file.")

        # euclidean group
        euclidean = sorted([file_prefix + f for f in beta_files if "euclidean" in f.lower() and ".png" not in f.lower()])
        results["Beta Diversity Data/Euclidean_distance"] = ", ".join(euclidean)

    # Taxonomy plots
    tax_plot_dir = os.path.join(final_outputs_dir, "taxonomy_plots")
    all_tax_plot_files = [file_prefix + f for f in sorted(os.listdir(tax_plot_dir))] if os.path.isdir(tax_plot_dir) else []
    if all_tax_plot_files:
        tax_plot_files = [f for f in all_tax_plot_files if ".png" not in f.lower()]
        results["Taxonomy Plots"] = ", ".join(tax_plot_files)


    # Differential abundance
    da_dir = os.path.join(final_outputs_dir, "differential_abundance")
    if os.path.isdir(da_dir):
        # top-level DA files
        da_files = [f for f in os.listdir(da_dir) if os.path.isfile(os.path.join(da_dir, f))]
        if da_files:
            results["Differential Abundance"] = ", ".join(file_prefix + f for f in da_files)

        # method-specific
        for method in ["ANCOMBC1", "ANCOMBC2", "DESeq2"]:
            subdir = os.path.join(da_dir, method.lower())
            if os.path.isdir(subdir):
                method_files = sorted([file_prefix + f for f in os.listdir(subdir) if not ("volcano" in f.lower() and f.lower().endswith(".png"))])
                if method_files:
                    results[f"Differential Abundance/{method}"] = ", ".join(method_files)

    return results



def runsheet_to_dict(runsheet):
    """ Reads the input nextflow runsheet into a dataframe and converts it to 
        a dictionary with sample names as keys and raw reads forward prefix used
        by multiqc as values
    """
    def get_prefix(string):
        basename  = os.path.basename(string)
        index     = basename.rfind("_")
        return(basename[0:index])
    df                    = pd.read_csv(runsheet, usecols=["sample_id", "forward"])
    df['forward']         = df.forward.apply(lambda row : get_prefix(row))
    sample_to_prefix_dict = {k:v['forward'] for k,v in df.set_index("sample_id").T.to_dict().items()}
    return(sample_to_prefix_dict)



def create_association_table(header_colnames, fastqc,
                             unique_filename_prefixes, read_count_tab, 
                             sample_file_dict, file_prefix,  output_prefix, combined_prefix,
                             assay_suffix,  raw_file_prefix, 
                             raw_suffix, raw_R1_suffix, raw_R2_suffix,
                             primer_trimmed_suffix, primer_trimmed_R1_suffix, primer_trimmed_R2_suffix,
                             filtered_suffix, filtered_R1_suffix, filtered_R2_suffix,
                             single_ended, R1_used_as_single_ended_data, sample_raw_prefix_dict,
                             include_raw_multiqc_in_output, read_count_unit = "read"):
    """Create association table and add data rows to it"""

    # Initialize association table
    association_df = pd.DataFrame(columns = header_colnames)
    # Create row
    for sample in unique_filename_prefixes:
        # Single-end (Paired-end data where only the forward reads were analyzed)
        if single_ended and R1_used_as_single_ended_data:
            # If only forward read was used, still want to include both foward and reverse read names 
            # in the "Raw Data" column because it is tied to the hosted raw data, not just what was used here
            curr_raw_data = [raw_file_prefix + sample + raw_R1_suffix,
                             raw_file_prefix + sample + raw_R2_suffix]
            # If only forward read was used, then we still want to include the _R1 portion of the filename 
            curr_trimmed_data = [file_prefix + sample + primer_trimmed_R1_suffix]
            curr_filt_data = [file_prefix + sample + filtered_R1_suffix]
        # Single-end without reverse reads
        elif single_ended:
            curr_raw_data = [raw_file_prefix + sample + raw_suffix]
            curr_trimmed_data = [file_prefix + sample + primer_trimmed_suffix]
            curr_filt_data = [file_prefix + sample + filtered_suffix]
        # Paired-end
        else:
            curr_raw_data = [raw_file_prefix + sample + raw_R1_suffix,
                             raw_file_prefix + sample + raw_R2_suffix]
            curr_trimmed_data = [file_prefix + sample + primer_trimmed_R1_suffix, 
                                 file_prefix + sample + primer_trimmed_R2_suffix]
            curr_filt_data = [file_prefix + sample + filtered_R1_suffix, 
                              file_prefix + sample + filtered_R2_suffix]
        # Get sample raw read count
        curr_read_count = get_read_count_from_df(sample, read_count_tab, raw_suffix,
                                                 raw_R1_suffix, single_ended, sample_raw_prefix_dict)

        final_outputs_dict = collect_final_outputs_columns(args.final_outputs_dir, file_prefix, output_prefix, assay_suffix)

        curr_row_as_list = [sample_file_dict[sample],
                            ", ".join(curr_raw_data),
                            curr_read_count, 
                            read_count_unit]
        
        # Divide fastqc into raw and filtered multiqc
        if include_raw_multiqc_in_output:
            raw_fastqc, filtered_fastqc = fastqc
            curr_row_as_list.append(raw_fastqc)  # Add only if flag is True
        else:
            filtered_fastqc = fastqc

        curr_row_as_list.extend([", ".join(curr_trimmed_data),
                                 ", ".join(curr_filt_data),
                                 filtered_fastqc,
                                 final_outputs_dict["Taxonomy and ASV Counts Data"],
                                 final_outputs_dict["Alpha Diversity Data"],
                                 final_outputs_dict["Beta Diversity Data"],
                                 final_outputs_dict["Beta Diversity Data/Bray-Curtis"],
                                 final_outputs_dict["Beta Diversity Data/Euclidean_distance"],
                                 final_outputs_dict["Taxonomy Plots"],
                                 final_outputs_dict["Differential Abundance"],
                                 final_outputs_dict["Differential Abundance/ANCOMBC1"],
                                 final_outputs_dict["Differential Abundance/ANCOMBC2"],
                                 final_outputs_dict["Differential Abundance/DESeq2"]])

        # Apped row to the association dataframe
        association_df.loc[len(association_df)] = curr_row_as_list

    return association_df


def write_association_table(outfile, association_df, trimmed_reads_dir='Trimmed Sequence Data', primers_already_trimmed=False):
    """Write to csv file"""
    # Removing trimmed column if primers were already removed and no primer-removal was done
    if primers_already_trimmed:
        association_df.drop(f"Parameter Value[{trimmed_reads_dir}]", axis = 1, inplace = True)
    # Writing out
    association_df.to_csv(outfile, sep = "\t", index = False)



def main():

    ### Set variables  ###
    
    # Directories
    fastqc_dir = str(args.fastqc_dir)
    raw_reads_dir = str(args.raw_reads_dir)
    filtered_reads_dir = str(args.filtered_reads_dir)
    trimmed_reads_dir = str(args.trimmed_reads_dir)
    final_outputs_dir = str(args.final_outputs_dir)
    raw_reads_dir = raw_reads_dir.replace("_", " ").rstrip("/")
    trimmed_reads_dir = trimmed_reads_dir.replace("_", " ").rstrip("/")
    filtered_reads_dir = filtered_reads_dir.replace("_", " ").rstrip("/")
    fastqc_dir = fastqc_dir.replace("_", " ").rstrip("/")
    final_outputs_dir = final_outputs_dir.replace("_", " ").rstrip("/")

    # Suffixes
    filtered_suffix = str(args.filtered_suffix)
    filtered_R1_suffix = str(args.filtered_R1_suffix)
    filtered_R2_suffix = str(args.filtered_R2_suffix)
    raw_suffix = str(args.raw_suffix)
    raw_R1_suffix = str(args.raw_R1_suffix)
    raw_R2_suffix = str(args.raw_R2_suffix)
    if args.R1_used_as_single_ended_data:
        raw_suffix = raw_R1_suffix
        # Just in case user only specified --R1-used-as-single-ended, but didn't specify --single-ended
        args.single_ended = True 
    primer_trimmed_suffix = str(args.primer_trimmed_suffix)
    primer_trimmed_R1_suffix = str(args.primer_trimmed_R1_suffix)
    primer_trimmed_R2_suffix = str(args.primer_trimmed_R2_suffix)
    assay_suffix = str(args.assay_suffix)


    # This one is only used for the raw data files
    raw_file_prefix = f"{args.GLDS_ID}_Amplicon_" if args.raw_file_prefix == "" else  str(args.raw_file_prefix)
    file_prefix = f"{args.GLDS_ID}_GAmplicon_" if args.file_prefix == "" else  str(args.file_prefix)
    raw_multiqc_zip = f"raw_multiqc{assay_suffix}_report.zip"
    filtered_multiqc_zip = f"filtered_multiqc{assay_suffix}_report.zip"
    output_prefix = str(args.output_prefix)
    combined_prefix = file_prefix + output_prefix
    raw_multiqc_stats_file_path = "raw_multiqc_report/raw_multiqc_data/multiqc_general_stats.txt"

    if args.map:
        map_tab = pd.read_csv(args.map, sep = "\t", names = ["sample", "prefix"])
        map_tab.set_index("sample", inplace = True)
    else:
        map_tab = None

    # Set output file name
    if args.output == "" and output_prefix != "":
        outfile = f"{args.GLDS_ID}_{output_prefix}-associated-file-names.tsv"
    elif args.output == "":
        outfile = f"{args.GLDS_ID}-associated-file-names.tsv"
    else:
        outfile = args.output
    
    include_raw_multiqc_in_output = args.include_raw_multiqc_in_output
    Type = str(args.type)


    # Ensure that GLfile.csv is passed as argument
    if args.GL_file == "":
        report_failure("This program requires a runsheet (such as GLfile.csv) as a required argument")

    
    # Check if assay_table exists and that it isn't empty
    check_for_file_and_contents(args.GL_file)

    #Read GLfile into dataframe and set sample names and unique prefixes
    GL_file = pd.read_csv(args.GL_file, sep = ",")
    sample_names = unique_filename_prefixes = GL_file["sample_id"].tolist()

    #sample_names, unique_filename_prefixes = get_sample_names_and_unique_filenames(assay_table,  raw_file_prefix, raw_R1_suffix,
     #                                                                              raw_R2_suffix, raw_suffix, 
      #                                                                             args.use_sample_names_from_assay_table,
       #                                                                            args.additional_string_to_remove_from_unique_filenames)

    sample_file_dict = dict(zip(unique_filename_prefixes, sample_names))

    read_counts_df = get_read_counts_from_raw_multiqc(map_tab, raw_multiqc_stats_file_path,
                                                      args.fastqc_dir, output_prefix,  raw_multiqc_zip)
    
    ###################################  Write file association table ##########################################
    header = write_colnames(raw_reads_dir, trimmed_reads_dir,
                  filtered_reads_dir, include_raw_multiqc_in_output)
    
    fastqc, final_outputs = create_constants(include_raw_multiqc_in_output, raw_multiqc_zip,
                     filtered_multiqc_zip, Type, combined_prefix, assay_suffix)
    
    # Retrieve a dictionary with sample names as keys and raw fatqfile prefix as values 
    sample_raw_prefix_dict = runsheet_to_dict(args.runsheet) if args.runsheet != "" else ""

    association_df = create_association_table(header, fastqc,
                             unique_filename_prefixes, read_counts_df, 
                             sample_file_dict, file_prefix,  output_prefix, combined_prefix,
                             assay_suffix,  raw_file_prefix, 
                             raw_suffix, raw_R1_suffix, raw_R2_suffix,
                             primer_trimmed_suffix, primer_trimmed_R1_suffix, primer_trimmed_R2_suffix,
                             filtered_suffix, filtered_R1_suffix, filtered_R2_suffix,
                             args.single_ended, args.R1_used_as_single_ended_data, 
                             sample_raw_prefix_dict, include_raw_multiqc_in_output, read_count_unit = "read")
    

    write_association_table(outfile, association_df, trimmed_reads_dir, args.primers_already_trimmed)


if __name__ == "__main__":
    main()
