"""
Example call:

python3 whole_genome_annotation.py \
--bowtie2_path bowtie2-2.5.4-sra-linux-x86_64/ \
--unannotated_genome_FAFSA_file whole_genome_annotation_sample_files/genome_W303_Liti_2009.fsa \
--labeled_sequences_FAFSA_file whole_genome_annotation_sample_files/orf_genomic_all.fasta \
--gff_output_file_name whole_genome_annotation_sample_files/W303_Liti_2009.gff

This code creates a number of "EGA-WGA-temp-*" files, which should be automatically deleted at the end
"""


import subprocess
import os
import sys
import argparse
import re
import pysam


def preprocess_fasta(input_fasta, output_fasta):
    """
    Preprocess a FASTA file to remove HTML tags and unwanted strings like "-->".
    
    Args:
        input_fasta (str): Path to the input FASTA file with labeled sequences.
        output_fasta (str): Path to the output FASTA file with cleaned sequences.
    """
    # Compile a regular expression to match HTML tags and the "-->" string
    html_tag_pattern = re.compile(r'<[^>]*>')  # Matches anything between < and >
    unwanted_string = "-->"                   # String to remove

    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        for line in infile:
            # Remove HTML tags
            cleaned_line = re.sub(html_tag_pattern, "", line)
            # Remove the "-->" string
            cleaned_line = cleaned_line.replace(unwanted_string, "-to-")
            # Write the cleaned sequence line
            outfile.write(cleaned_line)


if __name__ == "__main__":

    # STEP 1
    # Create the argument parser
    parser = argparse.ArgumentParser(description="This script requires user-supplied variables.")

    # Add required arguments
    parser.add_argument("--bowtie2_path", required=True, help="Path to bowtie2 folder.")
    parser.add_argument("--unannotated_genome_FAFSA_file",  required=True, help="Path and name of unannotated FAFSA file.")
    parser.add_argument("--labeled_sequences_FAFSA_file", required=True, help="Path and name of FAFSA file with labeled sequences.")
    parser.add_argument("--gff_output_file_name", required=True, help="Path and name of gff output file.")

    # Parse arguments
    args = parser.parse_args()

    # Access the arguments
    bowtie2_path = args.bowtie2_path
    unannotated_genome_FAFSA_file = args.unannotated_genome_FAFSA_file
    labeled_sequences_FAFSA_file = args.labeled_sequences_FAFSA_file
    gff_output_file_name	= args.gff_output_file_name

    # If no arguments are supplied, argparse will exit with an error message.
    print(f"Expecting bowtie2-build and bowtie2 in the folder: {bowtie2_path}")
    print(f"Genome FAFSA file to annotate: {unannotated_genome_FAFSA_file}")
    print(f"FAFSA file with annotated sequences: {labeled_sequences_FAFSA_file}")
    print(f"Output file: {gff_output_file_name}")


    # STEP 2
    # Clean up labeled_sequences_FAFSA_file, which often has HTML tags and the "-->" symbol
    labeled_sequences_FAFSA_file_cleaned = "EGA-WGA-temp-labels.fsa"
    preprocess_fasta(labeled_sequences_FAFSA_file, labeled_sequences_FAFSA_file_cleaned)
    print(f"Cleaned up the FAFSA file containing labels, removing HTML tags and '-->' symbol, and written to {labeled_sequences_FAFSA_file_cleaned}.")


    # STEP 3
    # Build Bowtie2 index
    print(f"--------------------------------Building Bowtie2 index for {unannotated_genome_FAFSA_file}...")
    index_base = "EGA-WGA-temp-index"
    subprocess.run([os.path.join(bowtie2_path, "bowtie2-build"),unannotated_genome_FAFSA_file,index_base], check=True)

    # Verify that expected files have been produced, else exit
    expected_index_files = [
        f"{index_base}.1.bt2",
        f"{index_base}.2.bt2",
        f"{index_base}.3.bt2",
        f"{index_base}.4.bt2",
        f"{index_base}.rev.1.bt2",
        f"{index_base}.rev.2.bt2"
    ]
    missing_files = [f for f in expected_index_files if not os.path.exists(f)]
    if missing_files:
        print(f"Error: The following Bowtie2 index files are missing: {missing_files}")
        sys.exit(1)  # Exit the script with an error
    print("Bowtie2 index files successfully created.")


    # STEP 4
    # Align sequences with Bowtie2
    print(f"--------------------------------Aligning {labeled_sequences_FAFSA_file} to {unannotated_genome_FAFSA_file}...")
    output_sam = "EGA-WGA-temp-SAM.sam"
    subprocess.run([os.path.join(bowtie2_path, "bowtie2"),"-f","-x",index_base,"-U",labeled_sequences_FAFSA_file_cleaned,"-S",output_sam], check=True)
    if not os.path.exists(output_sam):
        print(f"Error: SAM file {output_sam} was not created.")
        sys.exit(1)


    # STEP 5
    # Read SAM file
    sam_file = pysam.AlignmentFile(output_sam, "r")

    sam_data = []
    low_score = 0
    not_found = 0
    mapped = 0

    for read in sam_file:
        gene_info = {
            "Read name": read.query_name,
            "Reference name": read.reference_name,
            "Start Position": read.reference_start + 1,  # gff files are 1-based, sam files are 0-based
            "matches": read.get_cigar_stats()[0][0],
            "insertions": read.get_cigar_stats()[0][1],
            "deletions": read.get_cigar_stats()[0][2],
            "End Position": read.reference_end,
            "Mapping Quality": read.mapping_quality,
            "Strand": "+" if not read.is_reverse else "-",
        }

        sam_data.append(gene_info)

        if read.mapping_quality < 10:
            low_score += 1

    sam_file.close()


    # STEP 6
    # Write gff file, appending the genome sequence
    with open(gff_output_file_name, "w") as gff_output_file:
        gff_output_file.write("##gff-version 3\n")
        gff_output_file.write("# This genome annotation file was created with Easy Genome Annotation by L. Hug, V. Gligorovski, S. J. Rahi.\n")
        for sam_entry in sam_data:
            seq_id = str(sam_entry["Reference name"])

            if seq_id == "None":
                not_found += 1
                continue  # Skip this entry if seq_id is None

            source = "."
            feature_type = "feature"
            start = int(sam_entry["Start Position"])
            stop = int(sam_entry["End Position"]) if sam_entry["End Position"] is not None else 0
            score = str(sam_entry["Mapping Quality"])
            strand = sam_entry["Strand"]
            phase = "."
            attributes = "ID={};Name={}".format(sam_entry["Read name"], sam_entry["Read name"])
            gff_entry = "\t".join([seq_id, source, feature_type, str(start), str(stop), source,
                                   strand, source, attributes]) + "\n"
            gff_output_file.write(gff_entry)
            mapped += 1

        gff_output_file.write("###\n")
        gff_output_file.write("##FASTA\n")  # Add FASTA section marker

        # Read and append the content of the FASTA file
        with open(unannotated_genome_FAFSA_file, "r") as fasta_file:
            for line in fasta_file: 
                if line.startswith(">"):
                    seq_id = line.strip() #.split(".")[1]  SJR: delete 2024-11-15
                    gff_output_file.write(">"+seq_id + "\n")
                else: 
                    gff_output_file.write(line)

    print("Mapped: ", mapped)
    print("Of which with low mapQ score:", low_score)
    print("Not found: ", not_found)


    # STEP 7
    # Delete all temporary files
    print("Deleting temporary files...")
    [os.remove(f) for f in os.listdir() if f.startswith("EGA-WGA-temp") and os.path.isfile(f)]
    
