# 21/05/21

# We want to select 1918 pseudo mirnas sequences from the ones proposed by the paper of triplet svm.
# They should not be the same of the ones present in iMcRNA paper, since they will be used as independent dataset to test the models

import re

pseudo_file_triplet_svm = open("8494_hairpins_over_fe_15_bp_18_from_cds.txt", "r")
pseudo_file_triplet_svm = pseudo_file_triplet_svm.read()

pseudo_imcRNA_tool = open("iMcRNA_1612_false_premiRNAs_modified", "r")
pseudo_imcRNA_tool = pseudo_imcRNA_tool.read()

regex_headers_and_sequences = "(>.*\r\n)([AGUC\r\n]+)"  # splitting the fasta file in groups of headers and sequences
pseudo_imcRNA_tool_split = re.findall(regex_headers_and_sequences, pseudo_imcRNA_tool)

pseudo_file_triplet_svm_split = re.findall(regex_headers_and_sequences, pseudo_file_triplet_svm)

headers_svm_file_list = []
headers_imcrna_file_list = []

for tuple_svm in pseudo_file_triplet_svm_split:
    header_svm_file = tuple_svm[0]
    header_svm_file = header_svm_file.split("\t")[0]
    headers_svm_file_list.append(header_svm_file)
    #sequence_svm_file = tuple_svm[1]

for tuple_imcrna in pseudo_imcRNA_tool_split:
    header_imcrna_file = tuple_imcrna[0].replace("\r\n", "")
    headers_imcrna_file_list.append(header_imcrna_file)

#print(headers_svm_file_list)
#print(headers_imcrna_file_list)

headers_not_in_common = list(set(headers_svm_file_list) - set(headers_imcrna_file_list))
#print(headers_not_in_common)

outfile = open("1918_pseudo_mirna_sequences_triplet_svm", "w")
count = 0
for tuple_svm in pseudo_file_triplet_svm_split:
    header_svm_file = tuple_svm[0]
    header_svm_file = header_svm_file.split("\t")[0]
    sequence_svm_file = tuple_svm[1].replace("\n", "")

    for element in headers_not_in_common:

        if header_svm_file == element:
            count = count + 1
            if count <= 1918:
                outfile.write(header_svm_file + "\n")
                outfile.write(sequence_svm_file + "\n")


# Save the file: 1918_pseudo_mirna_sequences_triplet_svm




