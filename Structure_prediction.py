# 27/01/21
# 1/03/21

import re

###################################################
# Reading and management of the RNAfold output file
###################################################

file_output_rnafold = open("rnafold_output_all_human_hairpin_sequences_not_mutated", "r")  # opening the RNAfold output file
RNAfold_output = file_output_rnafold.readlines()  # reading the file by rows

row_to_print = " "
#rows = header + " \n " + sequence + " \n " + dot_bracket_mfe + " \t " + mfe + " \t " + efe + " \t " + centroid

# use regular expression to recognize every different energy values from the others
regex_header = "(>.*)"  # for HEADER
regex_sequence = "([AUCG]+)"  # for SEQUENCE
regex_dot_bracket_notation = "([\.\(\)]+) [ \(.*\)]"  # for DOT BRACKET NOTATION
regex_MFE = "( \(.*\))"  # for MFE
regex_ensembleFE = "( \[.*\])"  # for ENSEMBLE FREE ENERGY
regex_centroid = "( \{.*\})"  # for CENTROID STRUCTURE

for row in RNAfold_output:
    HEADER = re.findall(regex_header, row)  #findall function will find all the elements of the corresponding regex
    SEQUENCE = re.findall(regex_sequence, row)
    DOT_BRACKET_NOTATION_MFE = re.findall(regex_dot_bracket_notation, row)
    MFE = re.findall(regex_MFE, row)
    ENSEMBLE_FE = re.findall(regex_ensembleFE, row)
    CENTROID_STR = re.findall(regex_centroid, row)

    if (HEADER):
        HEADER = HEADER[0]  # to consider only the string
        HEADER = HEADER.lower()
        row_to_print = HEADER

    if (SEQUENCE):
        SEQUENCE = SEQUENCE[0]
        row_to_print = row_to_print + " \n " + SEQUENCE

    if (DOT_BRACKET_NOTATION_MFE) and (MFE):
        DOT_BRACKET_NOTATION_MFE = DOT_BRACKET_NOTATION_MFE[0]
        row_to_print = row_to_print + " \n " + DOT_BRACKET_NOTATION_MFE
        MFE = MFE[0]
        row_to_print = row_to_print + " \t " + MFE

    if (ENSEMBLE_FE):
        ENSEMBLE_FE = ENSEMBLE_FE[0]
        row_to_print = row_to_print + " \t " + ENSEMBLE_FE

    if (CENTROID_STR):
        CENTROID_STR = CENTROID_STR[0]
        row_to_print = row_to_print + " \t " + CENTROID_STR
        print(row_to_print)


#save the table obtained
