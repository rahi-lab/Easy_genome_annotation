"""
Example call:

python3 annotate_genomic_modifications.py \
--orig_genome_gff annotate_genomic_modifications_sample_files/genome_W303_Liti_2009_formatted.gff \
--modifications_gff annotate_genomic_modifications_sample_files/additional_seqs.gff \
--final_genome_gff annotate_genomic_modifications_sample_files/genome_W303_Liti_2009_formatted_updated.gff

The code puts everything in the final genome gff file, including ther updated genomic sequence.
"""


import argparse


# Parameters
homo_n = 20		# amount of homology to look for
text_width = 60



def parse_gff(gff_file):
    """
    Function to parse a GFF (General Feature Format) file and extract CDS features and sequences.

    Args:
        gff_file (str): Path to the GFF file.

    Returns:
        tuple: A tuple containing a list of CDS features and a dictionary of sequences.
    """
    cds_features = []
    sequences = {}
    with open(gff_file, "r") as f:
        in_fasta = False
        for line in f:
           
            if line.strip().startswith("##FASTA"):
                in_fasta = True
                continue
            if line.strip().startswith("#"):
                continue
            if in_fasta:
                if line.startswith(">"):
                    contig = line.strip()[1:]
                    sequences[contig] = []
                    
                else:
                    sequences[contig].append(line.strip())
            else:
                fields = line.strip().split("\t")
                
                if fields[2] == "gene":
                    contig = fields[0]   
                    start = int(fields[3])
                    end = int(fields[4])
                    orientation = fields[6]
                    attributes = dict(item.split("=")
                                      for item in fields[8].split(";"))
                    gene_name = attributes.get("gene", "")
                   
                    cds_name = attributes.get("Name", "")
                   
                    cds_features.append(
                        [contig, start, end, orientation, cds_name, gene_name])
    return cds_features, sequences


def reverse_complement(sequence):
    """
    Function to find the reverse complement of a DNA sequence.

    Args:
        sequence (str): DNA sequence.

    Returns:
        str: Reverse complement of the input DNA sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    reverse_comp = "".join(complement[base.upper()] for base in sequence[::-1])
    return reverse_comp


def compare_strings_last_agreement(string1, string2):
    last_agreement = -1
    for i in range(min(len(string1),len(string2))):
        if string1[i] == string2[i]:
            last_agreement = i
        else:
            break
    return last_agreement
    
    
    
def compare_strings_first_agreement(string1, string2, last_agreement):
    first_agreement = -1
    for i in range(min(len(string1),len(string2))):
        i_translated_from_start_string1 = len(string1) - i - 1 # -1 just because string indexing starts with 0
        i_translated_from_start_string2 = len(string2) - i - 1
        if string1[-i-1] == string2[-i-1] and i_translated_from_start_string1 > last_agreement and i_translated_from_start_string2 > last_agreement: # this was necessary since first_agreement could be to the left of last_agreement by chance agreement of nucleotides. messing things up
            first_agreement = i
        else:
            break
    return first_agreement



def compare_new_old_genome_seq(new_oneline_genome_seq,oneline_genome_seq):
    last_agreement  = compare_strings_last_agreement (new_oneline_genome_seq,oneline_genome_seq)
    first_agreement = compare_strings_first_agreement(new_oneline_genome_seq,oneline_genome_seq, last_agreement) #had to add last_agreement here, so we don't go further back than last_agreement
    first_agreement_in_orig_genome = len(oneline_genome_seq) - first_agreement - 1
    first_agreement_in_new_genome = len(new_oneline_genome_seq) - first_agreement - 1
    return last_agreement, first_agreement, first_agreement_in_orig_genome, first_agreement_in_new_genome



def add_new_cds_features(add_features, add_contig, num_repeats, match, match_homopiece_1_left, len_seq_added, cds_features):
    for feature_to_add in add_features:
        if feature_to_add[0] in add_contig:
            for repeat in range(num_repeats):
                if (match[3] == "+" and feature_to_add[3] == "+") or (match[3] == "-" and feature_to_add[3] == "-"):
                    new_orient= "+"
                else:
                    new_orient= "-"
    
                if match[3] == "+":
                    new_start = match_homopiece_1_left + feature_to_add[1] + repeat*len_seq_added
                    new_end   = match_homopiece_1_left + feature_to_add[2] + repeat*len_seq_added
                else:
                    new_start = match_homopiece_1_left + (len_seq_added - feature_to_add[2] + 1) + repeat*len_seq_added
                    new_end   = match_homopiece_1_left + (len_seq_added - feature_to_add[1] + 1) + repeat*len_seq_added
                new_feature   = [match[2].split(" ")[0], new_start, new_end, new_orient, feature_to_add[4], feature_to_add[5]]
                cds_features.append(new_feature)
                print("Adding new feature: " + feature_to_add[4])



def change_existing_cds_features(cds_features, last_agreement, first_agreement_in_orig_genome, len_seq_added):
    delete_cds_features = []
    duplicate_cds_features = []
    for i in range(len(cds_features)):
        if cds_features[i][0] in match[2]:
            if   cds_features[i][1] - 1 <= last_agreement and \
                 cds_features[i][2] - 1 <= last_agreement:
                print("Case 1: To remain as is: " + cds_features[i][4])
            elif cds_features[i][1] - 1 >= first_agreement_in_orig_genome and \
                 cds_features[i][2] - 1 >= first_agreement_in_orig_genome:
                cds_features[i][1] = cds_features[i][1] + len_seq_added
                cds_features[i][2] = cds_features[i][2] + len_seq_added
                print("Case 2: To be moved: " + cds_features[i][4])
            elif cds_features[i][1] - 1 <= last_agreement and \
                 cds_features[i][2] - 1 >  last_agreement and \
                 cds_features[i][2] - 1 <  first_agreement_in_orig_genome:
                cds_features[i][2] = last_agreement + 1 # + 1 because last_agreement counts string index starting from 0
                cds_features[i][4] = cds_features[i][4] + "-truncated"
                if cds_features[i][5]:
                    cds_features[i][5] = cds_features[i][5] + "-truncated"
                print("Case 3: To be truncated: " + cds_features[i][4])
            elif cds_features[i][1] - 1 >  last_agreement and \
                 cds_features[i][2] - 1 >= first_agreement_in_orig_genome and \
                 cds_features[i][1] - 1 <  first_agreement_in_orig_genome:
                cds_features[i][1] = first_agreement_in_new_genome + 1
                cds_features[i][2] = cds_features[i][2] + len_seq_added
                cds_features[i][4] = cds_features[i][4] + "-truncated"
                if cds_features[i][5]:
                    cds_features[i][5] = cds_features[i][5] + "-truncated"
                print("Case 4: To be truncated: " + cds_features[i][4])
            elif cds_features[i][1] - 1 >  last_agreement and \
                 cds_features[i][2] - 1 <  first_agreement_in_orig_genome:
                ### delete cds
                delete_cds_features.append(i)
                print("Case 5: To be deleted: " + cds_features[i][4])
            elif cds_features[i][1] - 1 <= last_agreement and \
                 cds_features[i][2] - 1 >= first_agreement_in_orig_genome:
                duplicate_cds_features.append([i, first_agreement_in_new_genome + 1, cds_features[i][2] + len_seq_added])
                cds_features[i][2] = last_agreement + 1 # cds_features[i][2] + len_seq_added
                cds_features[i][4] = cds_features[i][4] + "-inter"
                if cds_features[i][5]:
                    cds_features[i][5] = cds_features[i][5] + "-inter"
                print("Case 6: To be interrupted: " + cds_features[i][4])
    # have to duplicate first, and append, will not interfere with delete
    for feature_to_duplicate in duplicate_cds_features:
        i = feature_to_duplicate[0]
        cds_features_orig = cds_features[i][:]
        cds_features[i][4]   = cds_features[i][4]   + "-1"
        cds_features_orig[4] = cds_features_orig[4] + "-2"
        if cds_features[i][5]:
            cds_features[i][5]   = cds_features[i][5]   + "-1"
            cds_features_orig[5] = cds_features_orig[5] + "-2"
        cds_features_orig[1] = feature_to_duplicate[1]
        cds_features_orig[2] = feature_to_duplicate[2]
        cds_features.append(cds_features_orig)
    for i in sorted(delete_cds_features,reverse=True):
        del cds_features[i]
        
#                # change existing cds features
#                for i in range(len(cds_features)):
#                    if cds_features[i][0] in match[2]:
#                        if   cds_features[i][1] - 1 <  match_homopiece_1_left and \
#                             cds_features[i][2] - 1 <  match_homopiece_1_left: # I think < is the correct comparison in the last condition
#                            print("Case 1")
#                            continue
#                        elif cds_features[i][1] - 1 >= match_homopiece_1_left and \
#                             cds_features[i][2] - 1 >= match_homopiece_1_left:
#                            cds_features[i][1] = cds_features[i][1] + len_seq_added_repeated
#                            cds_features[i][2] = cds_features[i][2] + len_seq_added_repeated
#                            print("Case 2")
#                        elif cds_features[i][1] - 1 <  match_homopiece_1_left and \
#                             cds_features[i][2] - 1 >= match_homopiece_1_left:
#                            if cds_features[i][3] == "+":
#                                cds_features[i][2] = min(last_agreement + 1,cds_features[i][2]) # + 1 because last_agreement counts string index starting from 0
#                                print("Case 3")
#                            else:
#                                cds_features[i][1] = max(len(new_oneline_genome_seq) - first_agreement, cds_features[i][1] + len_seq_added_repeated)
#                                cds_features[i][2] = cds_features[i][2] + len_seq_added_repeated
#                                print("Case 4")
#                            cds_features[i][4] = cds_features[i][4] + "-tagged"
#                            if cds_features[i][5]:
#                                cds_features[i][5] = cds_features[i][5] + "-tagged"
#                        else:
#                            print("This case should not happen, end of gene should not come before start.")



def save_new_gff(new_gff_file, cds_features, sequences):
    with open(new_gff_file, 'w') as f:
        f.write("##gff-version 3\n")
        for feature in cds_features:
            feature_line = feature[0] + "\t" + "SJR" + "\t" + "gene" + "\t" + \
                    str(feature[1]) + "\t" + str(feature[2]) + "\t.\t" + \
                    feature[3] + "\t.\t" + "ID=" + feature[4] + ";Name=" + feature[4]
            if feature[5]:
                feature_line = feature_line + ";gene=" + feature[5]
            f.write(feature_line + "\n")
        f.write("###\n")
        f.write("##FASTA\n")
        for genome_contig, genome_seq in sequences.items():
            f.write(">" + genome_contig + "\n")
            for seq in genome_seq:
                f.write(seq + "\n")
#                    cds_features.append([contig, start, end, orientation, cds_name, gene_name])
            


if __name__ == "__main__":

    # Prepare input arguments
    # Create the argument parser
    parser = argparse.ArgumentParser(description="This script requires user-supplied variables.")

    # Add required arguments
    parser.add_argument("--orig_genome_gff", required=True, help="Path and name of original genome annotation gff file.")
    parser.add_argument("--modifications_gff",  required=True, help="Path and name of gff file with engineered genome modifications.")
    parser.add_argument("--final_genome_gff", required=True, help="Path and name of output gff file.")

    # Parse arguments
    args = parser.parse_args()

    # Access the arguments
    orig_gff = args.orig_genome_gff
    modifications_gff = args.modifications_gff
    final_gff = args.final_genome_gff

    # If no arguments are supplied, argparse will exit with an error message.
    print(f"Original genome annotation gff file: {orig_gff}")
    print(f"GFF file with engineered genome modifications: {modifications_gff}")
    print(f"Output gff file: {final_gff}")
    

    # Run main annotation
    cds_features, sequences = parse_gff(orig_gff)
    add_features, add_seqs  = parse_gff(modifications_gff)

    for add_contig, add_seq in add_seqs.items():
        print("------------------------------------")
        print("------------------------------------")
        print("Processing",add_contig)
        oneline_add_seq = "".join(add_seq)
        
        list_of_matches = [] #index of first homo_n characters in genome, index of last homo_n characters in genome, genome_contigs, orientation
        for genome_contig, genome_seq in sequences.items():
            oneline_genome_seq = "".join(genome_seq)
            if (oneline_add_seq[:homo_n] in oneline_genome_seq) and (oneline_add_seq[-homo_n:] in oneline_genome_seq):
                list_of_matches.append([oneline_genome_seq.find(oneline_add_seq[:homo_n]), \
                                        oneline_genome_seq.find(oneline_add_seq[-homo_n:]), \
                                        genome_contig, \
                                        "+"])
            elif (reverse_complement(oneline_add_seq[:homo_n]) in oneline_genome_seq) and (reverse_complement(oneline_add_seq[-homo_n:]) in oneline_genome_seq):
                list_of_matches.append([oneline_genome_seq.find(reverse_complement(oneline_add_seq[:homo_n])), \
                                        oneline_genome_seq.find(reverse_complement(oneline_add_seq[-homo_n:])), \
                                        genome_contig, \
                                        "-"])
                                        
        if len(list_of_matches) > 1:
        
            print("Warning, sequence was found more than once in the genome.")
            
        elif len(list_of_matches) < 1:
        
            print("Warning, sequence was not found in the genome.")
            
        else:
        
            match = list_of_matches[0]
            oneline_genome_seq = "".join(sequences[match[2]])
            if match[3] == "+":
                match_homopiece_1_left = match[0]
                match_homopiece_2_left = match[1]
            else:
                match_homopiece_1_left = match[1]
                match_homopiece_2_left = match[0]
                oneline_add_seq = reverse_complement(oneline_add_seq)
            
            match_found = 0
            
            # if it is an inserting plasmid
            if (match[3] == "+" and match[0] == match[1] + homo_n) or (match[3] == "-" and match[0] + homo_n == match[1]):
                print("'" + add_contig + "'", "appears to be an inserting plasmid.")
                match_found = 1

                num_repeats = int(add_contig.split(" ")[2].split("=")[1])
                oneline_add_seq_repeated = "".join([oneline_add_seq]*num_repeats)
                len_seq_added_for_changes          = len(oneline_add_seq_repeated)
                len_seq_added_for_additions        = len(oneline_add_seq)
                
                new_oneline_genome_seq = oneline_genome_seq[:match_homopiece_1_left] + oneline_add_seq_repeated + oneline_genome_seq[match_homopiece_1_left:]
                if len(new_oneline_genome_seq) - len(oneline_genome_seq) != len_seq_added_for_changes:
                    print("Something wrong with length count.")
                
            # if it is a replacing plasmid
            elif (match[3] == "+" and match[0] + homo_n <= match[1]) or (match[3] == "-" and match[1] + homo_n <= match[0]):
                print("'" + add_contig + "'", "appears to be a replacing plasmid.")
                match_found = 1

                len_seq_added_for_changes   = len(oneline_add_seq) - (match_homopiece_2_left - match_homopiece_1_left + 1 + homo_n - 1)
                num_repeats = 1
                len_seq_added_for_additions = len(oneline_add_seq)
                
                new_oneline_genome_seq = oneline_genome_seq[:match_homopiece_1_left] + oneline_add_seq + oneline_genome_seq[match_homopiece_2_left + homo_n:]
                if len(new_oneline_genome_seq) - len(oneline_genome_seq) != len_seq_added_for_changes:
                    print("Something wrong with length count.")
                    
            else:
                print("Could not integrate sequence in genome.")
                
            if match_found:
                last_agreement, first_agreement, first_agreement_in_orig_genome, first_agreement_in_new_genome = \
                    compare_new_old_genome_seq(new_oneline_genome_seq, oneline_genome_seq)
                
                change_existing_cds_features(cds_features, last_agreement, first_agreement_in_orig_genome, len_seq_added_for_changes)
                
                add_new_cds_features(add_features, add_contig, num_repeats, match, match_homopiece_1_left, len_seq_added_for_additions, cds_features)

                sequences[match[2]] = [new_oneline_genome_seq[i:i+text_width] for i in range(0, len(new_oneline_genome_seq), text_width)] # new_oneline_genome_seq
        
        save_new_gff(final_gff, cds_features, sequences)
        
