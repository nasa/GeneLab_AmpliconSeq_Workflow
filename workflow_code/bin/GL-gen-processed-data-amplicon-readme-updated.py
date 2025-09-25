#!/usr/bin/env python

"""
This is a program for generating a README.txt file for GeneLab processed amplicon datasets.
"""

import os
import sys
import argparse
import textwrap
import zipfile
import datetime

parser = argparse.ArgumentParser(
    description="This program generates the corresponding README file for GeneLab processed amplicon datasets. It is intended to \
                                             be run before running `GL-validate-processed-amplicon-data` and after processing_info.zip has been created.")

required = parser.add_argument_group('required arguments')
required.add_argument("--osd-id", help='OSDR osd_id ID (e.g. "OSD-69")', action="store", required=True)

parser.add_argument("--output", "--output-name", default="README",
                    help='Name of output file (default: "README")')
parser.add_argument("--name", default="", required=True,
                    help='Name of individual who performed the processing (default: "")')
parser.add_argument("--email", default="", required=True,
                    help='Email address of individual who performed the processing (default: "")')
parser.add_argument("--protocol-ID", default="",
                    help='Protocol document ID followed')
parser.add_argument("--assay_suffix", help = "Genelab assay suffix", action = "store", default = "_GLAmpSeq")
parser.add_argument("--primers_already_trimmed", help = "Add this flag if primers were trimmed prior to GeneLab processing, \
                    therefore there are no trimmed sequence data", action = "store_true")
parser.add_argument("--single-end", action="store_true", help="Is the data single-end?")
parser.add_argument("--date", action="store", default=datetime.date.today(), type=datetime.date.fromisoformat,
                    help="Date the processed data was generated in ISO format (YYYY-MM-DD). (default: today's date)")

parser.add_argument("--include-raw-reads", action="store_true", #???
                    help="Include the raw read data (Merged sequence data)?")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_args()


################################################################################

def main():
    with open(output_file, "w") as output:
        write_header(output, args.osd_id, args.name, args.email, args.protocol_ID, args.date)

        write_amplicon_body(output)


################################################################################

# setting some colors
tty_colors = {
    'green': '\033[0;32m%s\033[0m',
    'yellow': '\033[0;33m%s\033[0m',
    'red': '\033[0;31m%s\033[0m'
}


def color_text(text, color='green'):
    if sys.stdout.isatty():
        return tty_colors[color] % text
    else:
        return text


def wprint(text):
    """ print wrapper """

    print(textwrap.fill(text, width=80, initial_indent="  ",
                        subsequent_indent="  ", break_on_hyphens=False))


def report_failure(message, color="yellow"):
    print("")
    wprint(color_text(message, color))
    print("\nREADME-generation failed.\n")

    sys.exit(1)


def check_for_file_and_contents(file_path):
    """ used by get_processing_zip_contents function """

    if not os.path.exists(file_path):
        report_failure(f"The expected file '{file_path}' does not exist.")
    if not os.path.getsize(file_path) > 0:
        report_failure("The file '{file_path}' is empty.")


def get_processing_zip_contents():
    """ this gets the filenames that are in the processing_info.zip to add them to the readme """

    check_for_file_and_contents(processing_zip_file)

    with zipfile.ZipFile(processing_zip_file) as zip_obj:
        entries = zip_obj.namelist()
        entries.sort()

    return entries


def write_header(output, osd_id, name, email, protocol_id, date=datetime.date.today()):
    if not args.osd_id.startswith("OSD-"):
        print("OSDR osd_id must start with 'OSD-'")
        sys.exit(1)

    header = ["################################################################################\n",
              "{:<77} {:>0}".format(f"## This directory holds processed data for NASA {osd_id}", "##\n"),
              "{:<77} {:>0}".format(f"## https://osdr.nasa.gov/bio/repo/data/studies/{osd_id}/", "##\n"),
              "{:<77} {:>0}".format("##", "##\n"),
              "{:<77} {:>0}".format(f"## Processed by {name} ({email}) on {date}", "##\n"),
              "{:<77} {:>0}".format(f"## Based on {protocol_id}", "##\n"),
              "################################################################################\n\n",
              "Summary of contents:\n\n"]

    output.writelines(header)


# up_and_left: '┌',
# up_and_right: '┐',
# down_and_left: '└',
# down_and_right: '┘',
# vertical: '│',
# horizontal: '─',
# vertical_and_horizontal: '┼',
# down_and_horizontal: '┬',
# up_and_horizontal: '┴',
# top_connection: None,
# bottom_connection: None,
# HorizT='├'
# HLine='─'


def add_level(file_name, file_description, max_offset, output, continues, is_last=False):
    """
    This function prints a single row of a tree, similar to the Unix `tree` command,
    showing connectors (├──, └──) and vertical bars (│) to indicate hierarchy.

    Args:
        file_name (str): Name of the file or directory to display.
        file_description (str): Description of the file or directory. If empty, no description is shown.
        max_offset (int): Target width for the file name column to align descriptions.
        output (file-like object): File or buffer to write the formatted line to.
        continues (list of bool): A list of booleans representing ancestor nodes.
            Each element corresponds to an ancestor level (from root to parent):
                True  -> that ancestor has more siblings after this node (draw '│')
                False -> that ancestor is the last child (draw spaces)
        is_last (bool, optional): True if this node is the last child at its level.
            Determines whether the connector is '└──' (last) or '├──' (not last).
            Defaults to False.

    Returns:
        None: Writes the formatted line directly to the output.
    """

    # Build the prefix based on ancestors
    indent_parts = ["   "]
    for cont in continues[:-1]:
        indent_parts.append("│   " if cont else "    ")
    
    # Choose connector for this node
    connector = "└──" if is_last else "├──"
    indent_parts.append(connector)
    
    # Join all parts to form the full prefix for the line
    prefix = "".join(indent_parts)

    # Calculate spacing for the file name to align descriptions, ensure spacing >= 1
    spacing = max(1, max_offset - len(prefix))

    # Write the formatted line to the output buffer
    if file_description == "":
        output.write(f"{prefix} {file_name:<{spacing}}\n")
    else:
        output.write(f"{prefix} {file_name:<{spacing}}   - {file_description}\n")



def add_spacer(output):
    output.write("   │\n")


def write_amplicon_body(output):
    longest_filename = f"richness_and_diversity_estimates_by_sample{assay_suffix}.png"
    # length of padding is the length of the longest file + the indent_level + some extra
    pad = len(longest_filename) + (4 * 4) + 2

    # this file
    add_level(output_file, "this file", pad, output, continues=[True])

    add_spacer(output)

    # trimmed reads
    if not args.primers_already_trimmed:
        
        add_level(trimmed_reads_dir, "trimmed fastq files", pad, output, continues=[True])
        add_level(f"cutadapt{assay_suffix}.log", "log file of standard output and error from cutadapt", pad, output, continues=[True, True])
        add_level(f"trimmed-read-counts{assay_suffix}.tsv", "per sample read counts before and after trimming", pad, output, continues=[True, True])
        # Trimmed Sequence Data
        if args.single_end:
            add_level("*_trimmed.fastq.gz", "trimmed fastq files", pad, output, continues=[True, False], is_last=True)
        else:
            add_level("*_R1_trimmed.fastq.gz", "read1 trimmed fastq files", pad, output, continues=[True, True])
            add_level("*_R2_trimmed.fastq.gz", "read2 trimmed fastq files", pad, output, continues=[True, False], is_last=True)

    add_spacer(output)

    # quality-filtered reads
    add_level(filtered_reads_dir, "quality-filtered fastq files", pad, output, continues=[True])
    add_level(f"filtered-read-counts{assay_suffix}.tsv", "per sample read counts before and after quality-filtering", pad, output, continues=[True, True])
    # Filtered Sequence Data
    if args.single_end:
        add_level("*_filtered.fastq.gz", "filtered fastq files", pad, output, continues=[True, False], is_last=True)
    else:
        add_level("*_R1_filtered.fastq.gz", "read1 filtered fastq files", pad, output, continues=[True, True])
        add_level("*_R2_filtered.fastq.gz", "read2 filtered fastq files", pad, output, continues=[True, False], is_last=True)

    add_spacer(output)

    #fastqc outputs
    add_level(fastqc_outputs_dir, "multiQC summary reports of FastQC runs", pad, output, continues=[True])

    add_spacer(output)

    # final outputs
    add_level(tax_asv_outputs_dir, "", pad, output, continues=[True])
    add_level(f"ASVs{assay_suffix}.fasta", "fasta file of recovered sequences", pad, output, continues=[True, True])
    add_level(f"counts{assay_suffix}.tsv", "count table of sequences across samples", pad, output, continues=[True, True])
    add_level(f"taxonomy{assay_suffix}.tsv", "assigned taxonomy of recovered sequences", pad, output, continues=[True, True])
    add_level(f"taxonomy-and-counts{assay_suffix}.tsv", "combined table of counts and taxonomy", pad, output, continues=[True, True])
    add_level(f"taxonomy-and-counts{assay_suffix}.biom.zip", "biom-formatted output of counts and taxonomy", pad, output, continues=[True, True])
    add_level(f"read-count-tracking{assay_suffix}.tsv", "read counts at each processing step", pad, output, continues=[True, False], is_last=True)
    # Alpha diversity Reports
    add_level(a_diversity_dir, "alpha diversity measurements and plots using observed features estimates and Shannon diversity indices", pad, output, continues=[True])
    add_level(f"alpha_diversity_plots{assay_suffix}.zip", "", pad, output, continues=[True, True])
                            ### How to represent failure files? ###
    add_level(f"statistics_table{assay_suffix}.csv", "", pad, output, continues=[True, True])
    add_level(f"summary_table{assay_suffix}.csv", "", pad, output, continues=[True, True], is_last=True)
    # Beta diversity Reports
    add_level(b_diversity_dir, "", pad, output, continues=[True])
    add_level("vsd_validation_plot.png", "VST transformation validation diagnostic plot", pad, output, continues=[True, True])
    #Bray-Curtis results
    add_level("Bray-Curtis/", "beta diversity measurements and plots using Bray-Curtis dissimilarity", pad, output, continues=[True, True])
    add_level(f"bray_curtis_plots{assay_suffix}.zip", "", pad, output, continues=[True, True, True]) 
    add_level(f"bray_adonis_table{assay_suffix}.csv", "", pad, output, continues=[True, True, True])
    add_level(f"bray_variance_table{assay_suffix}.csv", "", pad, output, continues=[True, True, True], is_last=True)
    #Euclidean results
    add_level("Euclidean_distance/", "beta diversity measurements and plots using Euclidean distance", pad, output, continues=[True, False], is_last=True)
    add_level(f"euclidean_distance_plots{assay_suffix}.zip", "", pad, output, continues=[True, False, True]) 
    add_level(f"euclidean_adonis_table{assay_suffix}.csv", "", pad, output, continues=[True, False, True])
    add_level(f"euclidean_variance_table{assay_suffix}.csv", "", pad, output, continues=[True, False, True], is_last=True)    
    # Taxonomy plots
    add_level(taxonomy_plots_dir, "", pad, output, continues=[True])
    add_level(f"sample_taxonomy_plots{assay_suffix}.zip", "relative abundance plots of taxa for each sample", pad, output, continues=[True, True])
    add_level(f"group_taxonomy_plots{assay_suffix}.zip", "relative abundance plots of taxa for each group", pad, output, continues=[True, True], is_last=True)
    #Differential abundance
    add_level(diff_abundance_dir, "", pad, output, continues=[True])
    add_level(f"SampleTable{assay_suffix}.csv", "table containing samples and their respective groups", pad, output, continues=[True, True])
    add_level(f"contrasts{assay_suffix}.csv", "table specifying the group contrasts used for differential abundance analysis", pad, output, continues=[True, True])
    #ancombc1
    add_level(ancombc1_dir, "differential abundance analysis using ANCOMBC", pad, output, continues=[True, True])
    add_level(f"ancombc1_differential_abundance{assay_suffix}.csv", "normalized ASV counts and differential abundance statistics for pairwise group comparisons", pad, output, continues=[True, True, True])
    add_level(f"ancombc1_volcano_plots{assay_suffix}.zip", "volcano plots of pairwise group comparisons", pad, output, continues=[True, True, True], is_last=True)
    #ancombc2
    add_level(ancombc2_dir, "differential abundance analysis using ANCOM-BC2", pad, output, continues=[True, True])
    add_level(f"ancombc2_differential_abundance{assay_suffix}.csv", "normalized ASV counts and differential abundance statistics for pairwise group comparisons", pad, output, continues=[True, True, True])
    add_level(f"ancombc2_volcano_plots{assay_suffix}.zip", "volcano plots of pairwise group comparisons", pad, output, continues=[True, True, True], is_last=True)
    #deseq2
    add_level(deseq2_dir, "differential abundance analysis using DESeq2", pad, output, continues=[True, False], is_last=True)
    add_level(f"deseq2_differential_abundance{assay_suffix}.csv", "normalized ASV counts and differential abundance statistics for pairwise group comparisons", pad, output, continues=[True, False, True])
    add_level("asv_sparsity_plot.png", "diagnostic plot of ASV sparsity", pad, output, continues=[True, False, True])
    add_level(f"deseq2_volcano_plots{assay_suffix}.zip", "volcano plots of pairwise group comparisons", pad, output, continues=[True, False, True], is_last=True)

    add_spacer(output)

    # processing info
    add_level(processing_zip_file, "zip archive holding info related to processing (workflow files and metadata)", pad, output, continues=[False], is_last=True)

    output.write("\n")


# variable setup #

# universal settings
assay_suffix = args.assay_suffix

processing_zip_file = f"processing_info{assay_suffix}.zip"

merged_reads_dir = "Merged Sequence Data/"
#multiqc_dir = "MultiQC Reports/" ###Change FastQC Outputs to MultiQC Outputs?###
fastqc_outputs_dir = "FastQC Outputs/"
trimmed_reads_dir = "Trimmed Sequence Data/"
filtered_reads_dir = "Filtered Sequence Data/"
tax_asv_outputs_dir = "Taxonomy and ASV Counts Data/"
a_diversity_dir = "Alpha Diversity Data/"
b_diversity_dir = "Beta Diversity Data/"
taxonomy_plots_dir = "Taxonomy Plots/"
diff_abundance_dir = "Differential Abundance/"
ancombc1_dir = "ANCOMBC1/"
ancombc2_dir = "ANCOMBC2/"
deseq2_dir = "DESeq2/"

output_file = f"{str(args.output)}" 

if __name__ == "__main__":
    main()
