# 30/04/21

import re

#################################################
# INSERTION OF GERMLINE MUTATIONS IN PRE - MIRNAS
#################################################

# From the database ADmiRE, we have download the supplementary materials that provide all the PRE-GERMLINE mutations
# We want to insert those mutations in the corrisponding PRE mirnas.

file_TCGA_Pre_Germ = open("#TCGA_PanCanAtlas_Pre_Germ.txt", "r")
file_TCGA_Pre_Germ = file_TCGA_Pre_Germ.readlines()

fasta_file = open("fall_human_hairpin_sequences_not_mutated", "r")
fasta_file = fasta_file.read()
regex_headers_and_sequences = "(>.*\n)([AGUC\n]+)"  # splitting the fasta file in groups of headers and sequences
headers_and_sequences_split = re.findall(regex_headers_and_sequences, fasta_file)

annotation_file = open("hg19_mirna.gtf", "r") 
annotation_file = annotation_file.readlines()

#

for row in file_TCGA_Pre_Germ[1:]:
    row = row.split("\t")
    chr = row[0]
    position_snp = int(row[1])
    mirna_name_snp = row[2]
    starting_base = row[4]
    len_starting_base = len(starting_base)
    final_base = row[5]

    # make the complementary SNPs:
    complementary_final_base = ""
    for base in final_base:
        if (base == "A"): complementary_base = "T"
        elif (base == "T"): complementary_base = "A"
        elif (base == "C"): complementary_base ="G"
        elif (base == "G"): complementary_base = "C"
        elif (base == "-"): complementary_base = "-"
        complementary_final_base = complementary_final_base + complementary_base

    len_final_base = len(final_base)
    mutation_type = row[16]
    mutation = row[17]

    for line in annotation_file:
        line = line.split("\t")
        chr_annotation = line[0][3:]
        start_coordinate_annotation = int(line[3])
        end_coordinate_annotation = int(line[4])
        strand_annotation = line[6]
        name_annotation = line[8].split(";")[0][13:].replace("\"", "")

        if chr == chr_annotation and mirna_name_snp == name_annotation:

            for element in headers_and_sequences_split:
                header = element[0].replace("\n", "")
                header_split_name = header.split(" ")[0][5:]  # selecting only the name of the entire header
                sequence = element[1].replace("\n", "")
                reversed_sequence = sequence[::-1]

                if header_split_name == mirna_name_snp:  # matching the name in the fasta file with TCGA_Pre_Germ file

                    relative_position_snp = position_snp - start_coordinate_annotation + 1
                    # the relative position provided above is correct, but since in python we start counting from 0 we have to substruct 1 to the relative position
                    relative_position_snp = relative_position_snp - 1

                    if len_starting_base == 1 and final_base == "-":  # FOR SINGLE BASE DELETION (as G>-)

                        if strand_annotation == "+":
                            modified_sequence = sequence[:relative_position_snp] + sequence[relative_position_snp + 1:]
                            new_header = ">" + mirna_name_snp + "\t" + mutation_type + ", " + mutation
                            print(new_header)
                            print(modified_sequence)

                        else: # if the strand is negative, we have to use the reverde sequence, to add the complementary mutation and to reverse again the sequence
                            modified_sequence = reversed_sequence[:relative_position_snp] + reversed_sequence[relative_position_snp + 1:]
                            new_header = ">" + mirna_name_snp + "\t" + mutation_type + ", " + mutation
                            print(new_header)
                            print(modified_sequence[::-1])

                    elif len_starting_base > 1 and final_base == "-":  # FOR MULTIPLE DELETIONS (as ACGT>-)

                        if strand_annotation == "+":
                            modified_sequence = sequence[:relative_position_snp] + sequence[relative_position_snp + len_starting_base:]
                            new_header = ">" + mirna_name_snp + "\t" + mutation_type + ", " + mutation
                            print(new_header)
                            print(modified_sequence)

                        else: # if the strand is negative, we have to use the reverde sequence, to add the complementary mutation and to reverse again the sequence
                            modified_sequence = reversed_sequence[:relative_position_snp] + reversed_sequence[
                                                                                   relative_position_snp + len_starting_base:]
                            new_header = ">" + mirna_name_snp + "\t" + mutation_type + ", " + mutation
                            print(new_header)
                            print(modified_sequence[::-1])

                    elif starting_base == "-" and len_final_base == 1:  # FOR SINGLE BASE INSERTION (as ->G)

                        if strand_annotation == "+":
                            modified_sequence = sequence[:relative_position_snp] + final_base + sequence[relative_position_snp:]
                            new_header = ">" + mirna_name_snp + "\t" + mutation_type + ", " + mutation
                            print(new_header)
                            print(modified_sequence)

                        else: # if the strand is negative, we have to use the reverde sequence, to add the complementary mutation and to reverse again the sequence
                            modified_sequence = reversed_sequence[:relative_position_snp] + complementary_final_base + reversed_sequence[
                                                                                                relative_position_snp:]
                            new_header = ">" + mirna_name_snp + "\t" + mutation_type + ", " + mutation
                            print(new_header)
                            print(modified_sequence[::-1])

                    elif starting_base == "-" and len_final_base > 1:  # FOR MULTIPLE INSERTIONS (as ->ACGT)

                        if strand_annotation == "+":
                            modified_sequence = sequence[:relative_position_snp] + final_base +sequence[relative_position_snp:]
                            new_header = ">" + mirna_name_snp + "\t" + mutation_type + ", " + mutation
                            print(new_header)
                            print(modified_sequence)

                        else: # if the strand is negative, we have to use the reverde sequence, to add the complementary mutation and to reverse again the sequence
                            modified_sequence = reversed_sequence[:relative_position_snp] + complementary_final_base + reversed_sequence[
                                                                                                relative_position_snp:]
                            new_header = ">" + mirna_name_snp + "\t" + mutation_type + ", " + mutation
                            print(new_header)
                            print(modified_sequence[::-1])

                    elif len_starting_base == 1 and final_base != "-":  # FOR SINGLE BASE EXCHANGE (as C>T)

                        if strand_annotation == "+":
                            modified_sequence = sequence[:relative_position_snp] + final_base + sequence[
                                                                                              relative_position_snp + 1:]
                            new_header = ">" + mirna_name_snp + "\t" + mutation_type + ", " + mutation
                            print(new_header)
                            print(modified_sequence)

                        else: # if the strand is negative, we have to use the reverde sequence, to add the complementary mutation and to reverse again the sequence
                            modified_sequence = reversed_sequence[:relative_position_snp] + complementary_final_base + reversed_sequence[
                                                                                                relative_position_snp + 1:]
                            new_header = ">" + mirna_name_snp + "\t" + mutation_type + ", " + mutation
                            print(new_header)
                            print(modified_sequence[::-1])

                    elif len_starting_base > 1 and final_base != "-":  # FOR MULTIPLE BASE EXCHANGE (as TC>TTA)

                        if strand_annotation == "+":
                            modified_sequence = sequence[:relative_position_snp] + final_base + sequence[
                                                                                                relative_position_snp + len_starting_base:]
                            new_header = ">" + mirna_name_snp + "\t" + mutation_type + ", " + mutation
                            print(new_header)
                            print(modified_sequence)

                        else:  # if the strand is negative, we have to use the reverde sequence, to add the complementary mutation and to reverse again the sequence
                            modified_sequence = reversed_sequence[
                                                :relative_position_snp] + complementary_final_base + reversed_sequence[
                                                                                                     relative_position_snp + len_starting_base:]
                            new_header = ">" + mirna_name_snp + "\t" + mutation_type + ", " + mutation
                            print(new_header)
                            print(modified_sequence[::-1])



