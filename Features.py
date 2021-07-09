# miRNA FEATURES
# 6/02/21
# 16/02/21
# 9/03/21

import re  # module to exploit regular expressions

import textwrap  # module to format strings

import difflib  # module to compute and work with differences between string

# The aim of this script is to analyze some features of the output obtained running miRNA sequences in RNAfold.
# RNAfold is a tool for structure prediction and its typical output is the sequence, the so-called "dot-bracket notation" and the energy values.
# The dot-bracket notation is the most common and simple representation for RNA secondary structures.
# It is characterized by dots and brackets: in particular matching pairs of parenthesis '()' stand for base pairs and dots '..' stand for unpaired bases.
# It means that the base corresponding to the '(' will pair with the base corresponding to the ')', both looking the notation from the centre or from the ends.

# Working with this script we deal with a complication: some miRNA are characterized by "open" structures that cause mistakes in the evaluation of some features.
# In these situations, we do not consider these problematic structures and so we shorten the sequence and the dot-bracket notation as well.
# That's why you will find 2 versions of both sequence and dot-bracket notation: the untouched one (the entire sequence or notation including the problematic structures) and the truncated one (which is the one that we use for the analysis of the main features).

# Analyze the following features:

##########################
# DINUCLEOTIDE FREQUENCIES
##########################
# Starting from the sequences split every 2 bases, measure how many times we see A, G, C and U bases combined in couples (4^2 = 16 combinations).
# - AA, AG, AC, AU;
# - GA, GG, GC, GU;
# - CA, CG, CC, CU;
# - UA, UG, UC, UU.
# Then measure the frequency of every dinucletide, so the occurence of every dinucleotide should be divided by the total number of bases (the total length of the sequence), divided by 2.

###########################
# TRINUCLEOTIDE FREQUENCIES
###########################
# Starting from the sequences split every 3 bases, measure how many times we see A, G, C and U bases combined in triplets (4^3 = 64 combinations).
# - AAA, AAG, AAC, AAU, AGA, AGG, AGC, AGU, ACA, ACG, ACC, ACU, AUA, AUG, AUC, AUU;
# - GAA, GAG, GAC, GAU, GGA, GGG, GGC, GGU, GCA, GCG, GCC, GCU, GUA, GUG, GUC, GUU;
# - CAA, CAG, CAC, CAU, CGA, CGG, CGC, CGU, CCA, CCG, CCC, CCU, CUA, CUG, CUC, CUU;
# - UAA, UAG, UAC, UAU, UGA, UGG, UGC, UGU, UCA, UCG, UCC, UCU, UUA, UUG, UUC, UUU.
# Then measure the frequency of every trinucletide, so the occurence of every trinucleotide should be divided by the total number of bases (the total length of the sequence), divided by 3.

# N.B.: all the combinations of the 4 bases has been calculated in the following manner:
# import itertools
# bases = "AGCU"
# all_combinations = [''.join(i) for i in itertools.product(bases, bases, bases)]

##################################
# FREQUENCIES OF TRIPLET ELEMENTS
#################################
# In the predicted secondary structure (dot-bracket notation), there are only two statuses for each nucleotide, paired or unpaired, indicated by brackets '((' or '))' and dots ('.'), respectively.
# The left bracket '(' means that the paired nucleotide is located near the 5'-end and can be paired with another nucleotide at the 3'-end, which is indicated by a right bracket ')'.
# We decide to consider separately the paired condition at 5' and 3': thus for any 3 adjacent nucleotides there are 3^3 = 27 possible structure compositions:
# (((, ((), ((., ()(, ()), ()., (.(, (.), (.., )((, )(), )(., ))(, ))), ))., ).(, ).), ).., .((, .(), .(., .)(, .)), .)., ..(, ..), ...
# Also in this case we want to measure the frequency of every triplet element, so the occurence of every triplet element should be divided by the total length of the dot-bracket notation, that correponds to the total length of the sequence, divided by 3.

# N.B.: all the combinations has been calculated in the following manner:
#import itertools
#triplet_elements = "()."
#all_combinations = [''.join(i) for i in itertools.product(triplet_elements, triplet_elements, triplet_elements)]

##################################
# TOTAL NUMBER OF NOT PAIRED BASES
##################################
# In the dot-bracket notation, the total number of not paired bases should be equal to the total number of dots.

#############################################################
# TOTAL NUMBER OF LOOPS IN THE UNTOUCHED DOT-BRACKET NOTATION
#############################################################
# In general in the dot-bracket notation the loop is defined as a certain number of consecutive not paired bases (dots), contained in a certain number of matching pairs of parenthesis '()'
# There are different ways of defining a loop: as (.{3,}) ; ((.{3,})) ; (.{4,}) ; ((.{4,})) ; (.{5,}) ; ((.{5,}))
# We consider loops the ones with the more stringent definition: structures with at least 3 dots, contained in at least 1 matching pair of parenthesis '()'.
# Moreover we have to be sure of distinguish real loops with internal loop or bulges.
# Measure how many loops are present in the untouched dot-bracket notation.

##################################################
# TOTAL NUMBER OF LOOPS DEFINING THEM AS ((.{3,}))
##################################################
# Measure how many loops are present in the untouched dot-bracket notation, defining them as ((.{3,}))

################################################
# TOTAL NUMBER OF LOOPS DEFINING THEM AS (.{4,})
################################################
# Measure how many loops are present in the untouched dot-bracket notation, defining them as (.{4,})

##################################################
# TOTAL NUMBER OF LOOPS DEFINING THEM AS ((.{4,}))
##################################################
# Measure how many loops are present in the untouched dot-bracket notation, defining them as ((.{4,}))

################################################
# TOTAL NUMBER OF LOOPS DEFINING THEM AS (.{5,})
################################################
# Measure how many loops are present in the untouched dot-bracket notation, defining them as (.{5,})

##################################################
# TOTAL NUMBER OF LOOPS DEFINING THEM AS ((.{5,}))
##################################################
# Measure how many loops are present in the untouched dot-bracket notation, defining them as ((.{5,}))

#################################
# INDEXES OF THE LOOP OF INTEREST
#################################
# If the untouched dot-bracket notation has more than 1 loop, the loop of interest is the one whose structure is the longest respect to the others related to the other loops. So, we meausure the indexes of the loop of interest.
# Instead if the untouched dot-bracket notation has only 1 loop, the loop of interest correspond to the only one loop present.
# Measure the starting and ending indexes of the loop, taking into account only the consecutive not paired bases (dots) and not the matching pairs of parenthesis in which it is included.

####################################
# CENTRALITY OF THE LOOP OF INTEREST
####################################
# The definition of loop of interest is the same as before.
# Measure the position of the loop of interest (considering only the starting index), divided by the total length of the structure.

###################################
# DIMENSION OF THE LOOP OF INTEREST
###################################
# In the dot-bracket notation, the dimension of the loop of interest corresponds to the maximum number of consecutive not paired bases.
# If the number of loop is 1, we will measure the dimension of the only one loop present.

##############################
# TOTAL NUMBER OF PAIRED BASES
##############################
# It is equal to the total number of bases that composes the sequence (the total length of the sequence) minus the total number of not paired bases.

#############################
# FREQUENCY OF PAIRED COUPLES
#############################
# The possible pairings are: AU, UA, GU, UG , CG and GC.
# GU and UG are not canonical pairings that RNAfold already takes into account.
# Measure how many times we see AU, UA, GU, UG, CG and GC. Then divide the number of occurences by the total number of paired couples.
# The total number of paired couples can be obtained summing up all the occurences of the paired couples or using the total number of paired bases.

###############
# STEM SEQUENCE
###############
# The stem sequence corresponds to the sequence of consecutive paired bases.
# Of course, a stem is present only if a loop is present since they together form the so called stem-loop.

########
# BULGES
########
# Bulges are unpaired stretches of nucleotides.
# They are characterized by at least 1 not paired base inside the stem, so, they can be recognized as at least 1 dot inside 2 opened or closed brackets (in fact they can be present both at 5' or 3' of the sequence).
# Measure the total number of bulges present and their indexes.

################
# INTERNAL LOOPS
################
# Internal loops are found where the double stranded RNA separates due to no-base pairing between the nucleotides.
# Internal loops differ from stem loops as they occur in middle of a stretch of double stranded RNA.
# Using regular expression to found bulges, we found also internal loop. So, the challenge here is to recognize real bulges from internal loops.
# Analyze if the bases at the end of the not paired bases (dots) are paired: if the indexes of the corrisponding paired bases are consecutive means that it is a real bulge, otherwise it is an internal loop.
# Measure the total number of internal loops present and their indexes.

###############
# ENERGY VALUES
###############
# We report the following energy values: minimum free energy, ensemble free energy and centroid structure

#================================================================================================================================================================

RNAfold_output_file = open("table_rnafold_output_pseudo_miRNA_shuffled_sequences_3_dots_1_parentesis", "r")  # the input of the script is the summurizing table obtained from the RNAfold output
RNAfold_output = RNAfold_output_file.readlines()  # read the file by lines


# Regular expression are used to:

# 1- define which are the headers, the sequences, the dot-bracket notations and the energy values
regex_header = "(>.*)"  # for HEADER
regex_sequence = "([AUCG]+)"  # for SEQUENCE
regex_dot_bracket = "([\.\(\)]+) \t  \(.*\) \t  \[.*\] \t  \{.*\}"  # for DOT-BRACKET NOTATION
regex_energy_values = "[\.\(\)]+ \t ( \(.*\) \t  \[.*\] \t  \{.*\})"  # for respectively MINIMUM FREE ENERGY, ENSEMBLE FREE ENERGY and CENTROID STRUCTURE

# 2- define loops of interest
regex_loop_of_interest = "(\({1}\.{3,}\){1})+"  # at least 3 dots, inside 1 matching pair of parenthesis: '(...)'

# 3- other definitions for loops
regex_loop_3_dots_2_parentesis = "(\({2}\.{3,}\){2})+"  # at least 3 dots, inside 2 matching pairs of parenthesis: '((...))'
regex_loop_4_dots_1_parentesis = "(\({1}\.{4,}\){1})+"  # at least 4 dots, inside 1 matching pair of parenthesis: '(....)'
regex_loop_4_dots_2_parentesis = "(\({2}\.{4,}\){2})+"  # at least 4 dots, inside 2 matching pairs of parenthesis: '((....))'
regex_loop_5_dots_1_parentesis = "(\({1}\.{5,}\){1})+"  # at least 5 dots, inside 1 matching pair of parenthesis: '(.....)'
regex_loop_5_dots_2_parentesis = "(\({2}\.{5,}\){2})+"  # at least 5 dots, inside 2 matching pairs of parenthesis: '((....))'

# 4- define putative bulges
regex_putative_bulge_5_prime = "\(\.{1,}\("  # at least 1 dot, inside 1 not matching pair of parenthesis (both opened or closed parentesis): '(.(' or ').)'
regex_putative_bulge_3_prime = "\)\.{1,}\)"
regex_putative_bulge = regex_putative_bulge_5_prime + "|" + regex_putative_bulge_3_prime  # '|' is the inclusive or

#================================================================================================================================================================

warning_list = []  # list that will store all the miRNA headers with problems

for line in RNAfold_output:
    header = re.findall(regex_header, line)  # findall() function should find all the headers present in the file lines: in this case only one.
    sequence = re.findall(regex_sequence, line)
    dots_brackets = re.findall(regex_dot_bracket, line)
    energy_values = re.findall(regex_energy_values, line)

    # The output of regular expression is a list that should be empty if the element has not been found or full if the element has been found.
    # Since for every line the regular expression will found only one element, the string of the element is in position 0 of the list.

    if (header):  # it is equal to say 'if header != []:'
        header_str = header[0]

    if (sequence):
        # we distinguish two versions of the sequence string: the untouched one that will keep always the sequence as it is and another that during the scritp can be truncated, if some problematic structures or not paired bases at the ends are found.
        untouched_sequence = sequence[0]
        sequence_str = sequence[0]

    if (dots_brackets):
        # also here we distinguish two versions of the dot-bracket notation string
        untouched_dots_brackets = dots_brackets[0]
        dots_brackets_str = dots_brackets[0]

    if (energy_values):
        energy_values_str = energy_values[0]

        # First of all, find the number of loops in the untouched dot-bracket notation:
        find_loops = re.findall(regex_loop_of_interest, untouched_dots_brackets)

        #==================================  Working with sequences without loops

        if len(find_loops) == 0:  # if we have 0 loop

            # Even if we have not loops we continue to analyze the features.
            # If a feature cannot be calculated because needs the presence of the loop, we assign to it a default value of '0'.
            # We consider both the sequence and the dot-bracket notation as it is whitout cutting it. So, 'untouched_sequence' and 'untouched_dots_brackets' are interchangeable.

            total_length_seq = len(untouched_sequence)

            ### DINUCLEOTIDE FREQUENCIES

            sequence_split_in_two = textwrap.wrap(sequence_str, 2)  # split the sequence in couples

            total_count_AA = sequence_split_in_two.count("AA")  # find the number of AA occurences
            frequency_AA = total_count_AA / float(
                total_length_seq / 2)  # float() function returns a floating point number

            total_count_AG = sequence_split_in_two.count("AG")
            frequency_AG = total_count_AG / float(total_length_seq / 2)

            total_count_AC = sequence_split_in_two.count("AC")
            frequency_AC = total_count_AC / float(total_length_seq / 2)

            total_count_AU = sequence_split_in_two.count("AU")
            frequency_AU = total_count_AU / float(total_length_seq / 2)

            total_count_GA = sequence_split_in_two.count("GA")
            frequency_GA = total_count_GA / float(total_length_seq / 2)

            total_count_GG = sequence_split_in_two.count("GG")
            frequency_GG = total_count_GG / float(total_length_seq / 2)

            total_count_GC = sequence_split_in_two.count("GC")
            frequency_GC = total_count_GC / float(total_length_seq / 2)

            total_count_GU = sequence_split_in_two.count("GU")
            frequency_GU = total_count_GU / float(total_length_seq / 2)

            total_count_CA = sequence_split_in_two.count("CA")
            frequency_CA = total_count_CA / float(total_length_seq / 2)

            total_count_CG = sequence_split_in_two.count("CG")
            frequency_CG = total_count_CG / float(total_length_seq / 2)

            total_count_CC = sequence_split_in_two.count("CC")
            frequency_CC = total_count_CC / float(total_length_seq / 2)

            total_count_CU = sequence_split_in_two.count("CU")
            frequency_CU = total_count_CU / float(total_length_seq / 2)

            total_count_UA = sequence_split_in_two.count("UA")
            frequency_UA = total_count_UA / float(total_length_seq / 2)

            total_count_UG = sequence_split_in_two.count("UG")
            frequency_UG = total_count_UG / float(total_length_seq / 2)

            total_count_UC = sequence_split_in_two.count("UC")
            frequency_UC = total_count_UC / float(total_length_seq / 2)

            total_count_UU = sequence_split_in_two.count("UU")
            frequency_UU = total_count_UU / float(total_length_seq / 2)

            ### TRINUCLEOTIDE FREQUENCIES

            sequence_split_in_three = textwrap.wrap(sequence_str, 3)

            total_count_AAA = sequence_split_in_three.count("AAA")
            frequency_AAA = total_count_AAA / float(total_length_seq / 3)

            total_count_AAG = sequence_split_in_three.count("AAG")
            frequency_AAG = total_count_AAG / float(total_length_seq / 3)

            total_count_AAC = sequence_split_in_three.count("AAC")
            frequency_AAC = total_count_AAC / float(total_length_seq / 3)

            total_count_AAU = sequence_split_in_three.count("AAU")
            frequency_AAU = total_count_AAU / float(total_length_seq / 3)

            total_count_AGA = sequence_split_in_three.count("AGA")
            frequency_AGA = total_count_AGA / float(total_length_seq / 3)

            total_count_AGG = sequence_split_in_three.count("AGG")
            frequency_AGG = total_count_AGG / float(total_length_seq / 3)

            total_count_AGC = sequence_split_in_three.count("AGC")
            frequency_AGC = total_count_AGC / float(total_length_seq / 3)

            total_count_AGU = sequence_split_in_three.count("AGU")
            frequency_AGU = total_count_AGU / float(total_length_seq / 3)

            total_count_ACA = sequence_split_in_three.count("ACA")
            frequency_ACA = total_count_ACA / float(total_length_seq / 3)

            total_count_ACG = sequence_split_in_three.count("ACG")
            frequency_ACG = total_count_ACG / float(total_length_seq / 3)

            total_count_ACC = sequence_split_in_three.count("ACC")
            frequency_ACC = total_count_ACC / float(total_length_seq / 3)

            total_count_ACU = sequence_split_in_three.count("ACU")
            frequency_ACU = total_count_ACU / float(total_length_seq / 3)

            total_count_AUA = sequence_split_in_three.count("AUA")
            frequency_AUA = total_count_AUA / float(total_length_seq / 3)

            total_count_AUG = sequence_split_in_three.count("AUG")
            frequency_AUG = total_count_AUG / float(total_length_seq / 3)

            total_count_AUC = sequence_split_in_three.count("AUC")
            frequency_AUC = total_count_AUC / float(total_length_seq / 3)

            total_count_AUU = sequence_split_in_three.count("AUU")
            frequency_AUU = total_count_AUU / float(total_length_seq / 3)

            total_count_GAA = sequence_split_in_three.count("GAA")
            frequency_GAA = total_count_GAA / float(total_length_seq / 3)

            total_count_GAG = sequence_split_in_three.count("GAG")
            frequency_GAG = total_count_GAG / float(total_length_seq / 3)

            total_count_GAC = sequence_split_in_three.count("GAC")
            frequency_GAC = total_count_GAC / float(total_length_seq / 3)

            total_count_GAU = sequence_split_in_three.count("GAU")
            frequency_GAU = total_count_GAU / float(total_length_seq / 3)

            total_count_GGA = sequence_split_in_three.count("GGA")
            frequency_GGA = total_count_GGA / float(total_length_seq / 3)

            total_count_GGG = sequence_split_in_three.count("GGG")
            frequency_GGG = total_count_GGG / float(total_length_seq / 3)

            total_count_GGC = sequence_split_in_three.count("GGC")
            frequency_GGC = total_count_GGC / float(total_length_seq / 3)

            total_count_GGU = sequence_split_in_three.count("GGU")
            frequency_GGU = total_count_GGU / float(total_length_seq / 3)

            total_count_GCA = sequence_split_in_three.count("GCA")
            frequency_GCA = total_count_GCA / float(total_length_seq / 3)

            total_count_GCG = sequence_split_in_three.count("GCG")
            frequency_GCG = total_count_GCG / float(total_length_seq / 3)

            total_count_GCC = sequence_split_in_three.count("GCC")
            frequency_GCC = total_count_GCC / float(total_length_seq / 3)

            total_count_GCU = sequence_split_in_three.count("GCU")
            frequency_GCU = total_count_GCU / float(total_length_seq / 3)

            total_count_GUA = sequence_split_in_three.count("GUA")
            frequency_GUA = total_count_GUA / float(total_length_seq / 3)

            total_count_GUG = sequence_split_in_three.count("GUG")
            frequency_GUG = total_count_GUG / float(total_length_seq / 3)

            total_count_GUC = sequence_split_in_three.count("GUC")
            frequency_GUC = total_count_GUC / float(total_length_seq / 3)

            total_count_GUU = sequence_split_in_three.count("GUU")
            frequency_GUU = total_count_GUU / float(total_length_seq / 3)

            total_count_CAA = sequence_split_in_three.count("CAA")
            frequency_CAA = total_count_CAA / float(total_length_seq / 3)

            total_count_CAG = sequence_split_in_three.count("CAG")
            frequency_CAG = total_count_CAG / float(total_length_seq / 3)

            total_count_CAC = sequence_split_in_three.count("CAC")
            frequency_CAC = total_count_CAC / float(total_length_seq / 3)

            total_count_CAU = sequence_split_in_three.count("CAU")
            frequency_CAU = total_count_CAU / float(total_length_seq / 3)

            total_count_CGA = sequence_split_in_three.count("CGA")
            frequency_CGA = total_count_CGA / float(total_length_seq / 3)

            total_count_CGG = sequence_split_in_three.count("CGG")
            frequency_CGG = total_count_CGG / float(total_length_seq / 3)

            total_count_CGC = sequence_split_in_three.count("CGC")
            frequency_CGC = total_count_CGC / float(total_length_seq / 3)

            total_count_CGU = sequence_split_in_three.count("CGU")
            frequency_CGU = total_count_CGU / float(total_length_seq / 3)

            total_count_CCA = sequence_split_in_three.count("CCA")
            frequency_CCA = total_count_CCA / float(total_length_seq / 3)

            total_count_CCG = sequence_split_in_three.count("CCG")
            frequency_CCG = total_count_CCG / float(total_length_seq / 3)

            total_count_CCC = sequence_split_in_three.count("CCC")
            frequency_CCC = total_count_CCC / float(total_length_seq / 3)

            total_count_CCU = sequence_split_in_three.count("CCU")
            frequency_CCU = total_count_CCU / float(total_length_seq / 3)

            total_count_CUA = sequence_split_in_three.count("CUA")
            frequency_CUA = total_count_CUA / float(total_length_seq / 3)

            total_count_CUG = sequence_split_in_three.count("CUG")
            frequency_CUG = total_count_CUG / float(total_length_seq / 3)

            total_count_CUC = sequence_split_in_three.count("CUC")
            frequency_CUC = total_count_CUC / float(total_length_seq / 3)

            total_count_CUU = sequence_split_in_three.count("CUU")
            frequency_CUU = total_count_CUU / float(total_length_seq / 3)

            total_count_UAA = sequence_split_in_three.count("UAA")
            frequency_UAA = total_count_UAA / float(total_length_seq / 3)

            total_count_UAG = sequence_split_in_three.count("UAG")
            frequency_UAG = total_count_UAG / float(total_length_seq / 3)

            total_count_UAC = sequence_split_in_three.count("UAC")
            frequency_UAC = total_count_UAC / float(total_length_seq / 3)

            total_count_UAU = sequence_split_in_three.count("UAU")
            frequency_UAU = total_count_UAU / float(total_length_seq / 3)

            total_count_UGA = sequence_split_in_three.count("UGA")
            frequency_UGA = total_count_UGA / float(total_length_seq / 3)

            total_count_UGG = sequence_split_in_three.count("UGG")
            frequency_UGG = total_count_UGG / float(total_length_seq / 3)

            total_count_UGC = sequence_split_in_three.count("UGC")
            frequency_UGC = total_count_UGC / float(total_length_seq / 3)

            total_count_UGU = sequence_split_in_three.count("UGU")
            frequency_UGU = total_count_UGU / float(total_length_seq / 3)

            total_count_UCA = sequence_split_in_three.count("UCA")
            frequency_UCA = total_count_UCA / float(total_length_seq / 3)

            total_count_UCG = sequence_split_in_three.count("UCG")
            frequency_UCG = total_count_UCG / float(total_length_seq / 3)

            total_count_UCC = sequence_split_in_three.count("UCC")
            frequency_UCC = total_count_UCC / float(total_length_seq / 3)

            total_count_UCU = sequence_split_in_three.count("UCU")
            frequency_UCU = total_count_UCU / float(total_length_seq / 3)

            total_count_UUA = sequence_split_in_three.count("UUA")
            frequency_UUA = total_count_UUA / float(total_length_seq / 3)

            total_count_UUG = sequence_split_in_three.count("UUG")
            frequency_UUG = total_count_UUG / float(total_length_seq / 3)

            total_count_UUC = sequence_split_in_three.count("UUC")
            frequency_UUC = total_count_UUC / float(total_length_seq / 3)

            total_count_UUU = sequence_split_in_three.count("UUU")
            frequency_UUU = total_count_UUU / float(total_length_seq / 3)

            # From now on working with 'dots_brackets_str': if more than 1 loop is present, it corresponds with the truncated dot-bracket notation overwritten

            ### FREQUENCIES OF TRIPLET ELEMENTS

            dots_brackets_split_in_three = textwrap.wrap(dots_brackets_str, 3)

            find_bracket_bracket_bracket_5_prime = dots_brackets_split_in_three.count("(((")
            find_bracket_bracket_bracket_3_prime = dots_brackets_split_in_three.count(")))")
            total_count_bracket_bracket_bracket = find_bracket_bracket_bracket_5_prime + find_bracket_bracket_bracket_3_prime
            frequency_bracket_bracket_bracket = total_count_bracket_bracket_bracket / float(total_length_seq / 3)

            find_bracket_bracket_dot_5_prime = dots_brackets_split_in_three.count("((.")
            find_bracket_bracket_dot_3_prime = dots_brackets_split_in_three.count(")).")
            total_count_bracket_bracket_dot = find_bracket_bracket_dot_5_prime + find_bracket_bracket_dot_3_prime
            frequency_bracket_bracket_dot = total_count_bracket_bracket_dot / float(total_length_seq / 3)

            find_bracket_dot_bracket_5_prime = dots_brackets_split_in_three.count("(.(")
            find_bracket_dot_bracket_3_prime = dots_brackets_split_in_three.count(").)")
            total_count_bracket_dot_bracket = find_bracket_dot_bracket_5_prime + find_bracket_dot_bracket_3_prime
            frequency_bracket_dot_bracket = total_count_bracket_dot_bracket / float(total_length_seq / 3)

            find_bracket_dot_dot_5_prime = dots_brackets_split_in_three.count("(..")
            find_bracket_dot_dot_3_prime = dots_brackets_split_in_three.count(")..")
            total_count_bracket_dot_dot = find_bracket_dot_dot_5_prime + find_bracket_dot_dot_3_prime
            frequency_bracket_dot_dot = total_count_bracket_dot_dot / float(total_length_seq / 3)

            find_dot_bracket_bracket_5_prime = dots_brackets_split_in_three.count(".((")
            find_dot_bracket_bracket_3_prime = dots_brackets_split_in_three.count(".))")
            total_count_dot_bracket_bracket = find_dot_bracket_bracket_5_prime + find_dot_bracket_bracket_3_prime
            frequency_dot_bracket_bracket = total_count_dot_bracket_bracket / float(total_length_seq / 3)

            find_dot_bracket_dot_5_prime = dots_brackets_split_in_three.count(".(.")
            find_dot_bracket_dot_3_prime = dots_brackets_split_in_three.count(".).")
            total_count_dot_bracket_dot = find_dot_bracket_dot_5_prime + find_dot_bracket_dot_3_prime
            frequency_dot_bracket_dot = total_count_dot_bracket_dot / float(total_length_seq / 3)

            find_dot_dot_bracket_5_prime = dots_brackets_split_in_three.count("..(")
            find_dot_dot_bracket_3_prime = dots_brackets_split_in_three.count("..)")
            total_count_dot_dot_bracket = find_dot_dot_bracket_5_prime + find_dot_dot_bracket_3_prime
            frequency_dot_dot_bracket = total_count_dot_dot_bracket / float(total_length_seq / 3)

            total_count_dot_dot_dot = dots_brackets_split_in_three.count("...")
            frequency_dot_dot_dot = total_count_dot_dot_dot / float(total_length_seq / 3)

            #

            total_count_bracket_bracket_bracket_2 = dots_brackets_split_in_three.count("(()")
            frequency_bracket_bracket_bracket_2 = total_count_bracket_bracket_bracket_2 / float(total_length_seq / 3)

            total_count_bracket_bracket_bracket_3 = dots_brackets_split_in_three.count("()(")
            frequency_bracket_bracket_bracket_3 = total_count_bracket_bracket_bracket_3 / float(total_length_seq / 3)

            total_count_bracket_bracket_bracket_4 = dots_brackets_split_in_three.count("())")
            frequency_bracket_bracket_bracket_4 = total_count_bracket_bracket_bracket_4 / float(total_length_seq / 3)

            total_count_bracket_bracket_bracket_5 = dots_brackets_split_in_three.count(")((")
            frequency_bracket_bracket_bracket_5 = total_count_bracket_bracket_bracket_5 / float(total_length_seq / 3)

            total_count_bracket_bracket_bracket_6 = dots_brackets_split_in_three.count(")()")
            frequency_bracket_bracket_bracket_6 = total_count_bracket_bracket_bracket_6 / float(total_length_seq / 3)

            total_count_bracket_bracket_bracket_7 = dots_brackets_split_in_three.count("))(")
            frequency_bracket_bracket_bracket_7 = total_count_bracket_bracket_bracket_7 / float(total_length_seq / 3)

            total_count_bracket_bracket_dot_2 = dots_brackets_split_in_three.count("().")
            frequency_bracket_bracket_dot_2 = total_count_bracket_bracket_dot_2 / float(total_length_seq / 3)

            total_count_bracket_bracket_dot_3 = dots_brackets_split_in_three.count(")(.")
            frequency_bracket_bracket_dot_3 = total_count_bracket_bracket_dot_3 / float(total_length_seq / 3)

            total_count_bracket_dot_bracket_2 = dots_brackets_split_in_three.count("(.)")
            frequency_bracket_dot_bracket_2 = total_count_bracket_dot_bracket_2 / float(total_length_seq / 3)

            total_count_bracket_dot_bracket_3 = dots_brackets_split_in_three.count(").(")
            frequency_bracket_dot_bracket_3 = total_count_bracket_dot_bracket_3 / float(total_length_seq / 3)

            total_count_dot_bracket_bracket_2 = dots_brackets_split_in_three.count(".()")
            frequency_dot_bracket_bracket_2 = total_count_dot_bracket_bracket_2 / float(total_length_seq / 3)

            total_count_dot_bracket_bracket_3 = dots_brackets_split_in_three.count(".)(")
            frequency_dot_bracket_bracket_3 = total_count_dot_bracket_bracket_3 / float(total_length_seq / 3)

            ### TOTAL NUMBER OF NOT PAIRED BASES

            numb_not_paired_bases = dots_brackets_str.count(
                '.')  # not paired bases correspond to dots. Count() function counts how many dots there are

            ### TOTAL NUMBER OF LOOPS IN THE UNTOUCHED DOT-BRACKET NOTATION
            # We obtain the total number of loops from the list that contains all the loops found in the untouched dot-bracket notation ('find_loops' list)
            numb_loops = len(
                find_loops)  # len() function count how many elements there are in the list: it corresponds to how many loops there are

            ### TOTAL NUMBER OF LOOPS DEFINING THEM AS ((.{3,}))
            # Defining loops as '((...))', find the number of loops in the untouched dot-bracket notation:
            find_loop_3_dots_2_parentesis = re.findall(regex_loop_3_dots_2_parentesis, untouched_dots_brackets)
            numb_loops_3_dots_2_parentesis = len(find_loop_3_dots_2_parentesis)

            ### TOTAL NUMBER OF LOOPS DEFINING THEM AS (.{4,})
            # Defining loops as '(....)', find the number of loops in the untouched dot-bracket notation:
            find_loop_4_dots_1_parentesis = re.findall(regex_loop_4_dots_1_parentesis, untouched_dots_brackets)
            numb_loops_4_dots_1_parentesis = len(find_loop_4_dots_1_parentesis)

            ### TOTAL NUMBER OF LOOPS DEFINING THEM AS ((.{4,}))
            # Defining loops as '((....))', find the number of loops in the untouched dot-bracket notation:
            find_loop_4_dots_2_parentesis = re.findall(regex_loop_4_dots_2_parentesis, untouched_dots_brackets)
            numb_loops_4_dots_2_parentesis = len(find_loop_4_dots_2_parentesis)

            ### TOTAL NUMBER OF LOOPS DEFINING THEM AS (.{5,})
            # Defining loops as '(.....)', find the number of loops in the untouched dot-bracket notation:
            find_loop_5_dots_1_parentesis = re.findall(regex_loop_5_dots_1_parentesis, untouched_dots_brackets)
            numb_loops_5_dots_1_parentesis = len(find_loop_5_dots_1_parentesis)

            ### TOTAL NUMBER OF LOOPS DEFINING THEM AS ((.{5,}))
            # Defining loops as '((....))', find the number of loops in the untouched dot-bracket notation:
            find_loop_5_dots_2_parentesis = re.findall(regex_loop_5_dots_2_parentesis, untouched_dots_brackets)
            numb_loops_5_dots_2_parentesis = len(find_loop_5_dots_2_parentesis)

            ###  INDEXES OF THE LOOP OF INTEREST
            # Since we do not have a loop, we assign a default value of '0'.
            starting_index_loop_of_interest = 0  # default value
            length_loop_of_interest = 0  # default value
            ending_index_loop_of_interest = 0  # default value

            ###  CENTRALITY OF THE LOOP OF INTEREST

            loop_centrality = starting_index_loop_of_interest / float(total_length_seq)  # default value

            ### DIMENSION OF THE LOOP OF INTEREST

            # length_loop_of_interest

            ### TOTAL NUMBER OF PAIRED BASES

            numb_paired_bases = total_length_seq - numb_not_paired_bases  # it corresponds to the total number of bases minus the number of bases not paired

            ### FREQUENCY OF PAIRED COUPLES
            # We have to start from the ends of the dot-bracket notation and proceed till the centre, analyzing if the bases are paired.

            initial_index = 0  # it is the index used to count the positions from the start of the dot-bracket notation
            final_index = total_length_seq - 1  # it is the index used to count the positions from the end of the dot-bracket notation

            couples = []  # list that will contain all the couples found

            indexes_paired_bases = []  # list that will contain all the indexes of all the couples found

            for element in dots_brackets_str:

                if dots_brackets_str[initial_index] == "(" and dots_brackets_str[final_index] == ")":  # there is a pairing

                    couple_of_bases = sequence_str[initial_index] + sequence_str[final_index]  # extract the first base that corresponds to '(' and the last base that corresponds to ')'
                    couples.append(couple_of_bases)

                    # save also the indexes of the paired bases: they should be used later to recognize bulges from internal loops!
                    indexes_to_be_added = []  # sub list that will contain the indexes of the 2 bases that compose the pairing
                    indexes_to_be_added.append(initial_index)  # index of the base at 5'
                    indexes_to_be_added.append(final_index)  # index of the base at 3'
                    indexes_paired_bases.append(indexes_to_be_added)

                    # to move closer to the centre, the initial index should be increased of 1, instead the final index should be decreased of 1
                    initial_index = initial_index + 1
                    final_index = final_index - 1

                elif dots_brackets_str[initial_index] == "." and dots_brackets_str[final_index] == ")":  # there isn't a pairing

                    # the initial index should be increased of 1 to move closer to the centre, instead the final index remains fixed
                    initial_index = initial_index + 1

                elif dots_brackets_str[initial_index] == "(" and dots_brackets_str[final_index] == ".":  # there isn't a pairing

                    # the initial index remains fixed, instead the final index should be decrease of 1 to move closer to the centre
                    final_index = final_index - 1

                elif dots_brackets_str[initial_index] == "." and dots_brackets_str[final_index] == ".":  # there isn't a pairing

                    # to move closer to the centre the initial index should be increased of 1, instead the final index should be decreased of 1
                    initial_index = initial_index + 1
                    final_index = final_index - 1

            # Now we have to divide all the base couples that have been put in 'couples' list:

            all_CG_couples = ""  # it is an empty string that will contain all the couples CG
            all_GC_couples = ""
            all_UA_couples = ""
            all_AU_couples = ""
            all_UG_couples = ""
            all_GU_couples = ""

            for couple in couples:  # put every couple of bases found in the corresponding string
                if couple == "CG":
                    all_CG_couples = all_CG_couples + couple

                elif couple == "GC":
                    all_GC_couples = all_GC_couples + couple

                elif couple == "UA":
                    all_UA_couples = all_UA_couples + couple

                elif couple == "AU":
                    all_AU_couples = all_AU_couples + couple

                elif couple == "UG":
                    all_UG_couples = all_UG_couples + couple

                elif couple == "GU":
                    all_GU_couples = all_GU_couples + couple

            count_CG = len(all_CG_couples) / 2  # count how many CG occurences there are in the 'all_CG_couples' string. We use len() function to measure the bases and then divide by 2 to define how many couples there are
            count_GC = len(all_GC_couples) / 2
            count_UA = len(all_UA_couples) / 2
            count_AU = len(all_AU_couples) / 2
            count_UG = len(all_UG_couples) / 2
            count_GU = len(all_GU_couples) / 2

            total_numb_pairings = count_CG + count_GC + count_UA + count_AU + count_UG + count_GU  # measure the total number of pairings found

            if total_numb_pairings <= 0:
                warning_list.append(header_str)
                continue  # skip it

            else:

                frequency_couple_CG = count_CG / float(total_numb_pairings)  # measure the number of occurence of CG over the total number of pairings
                frequency_couple_GC = count_GC / float(total_numb_pairings)
                frequency_couple_UA = count_UA / float(total_numb_pairings)
                frequency_couple_AU = count_AU / float(total_numb_pairings)
                frequency_couple_UG = count_UG / float(total_numb_pairings)
                frequency_couple_GU = count_GU / float(total_numb_pairings)

            ### STEM SEQUENCE
            # Since we do not have a loop, we do not have neither a loop.

            reversed_stem_sequence_5_prime = "NA"  # default value

            reversed_stem_sequence_3_prime = "NA"  # default value

            ### DEFINE PUTATIVE BULGES

            putative_bulges = re.finditer(regex_putative_bulge, dots_brackets_str)  # finditer() function returns a tuple with an iterator yielding match-object instances over all non-overlapping matches for the RE pattern in string. The string is scanned left-to-right, and matches are returned in the order found.

            if (putative_bulges):  # if at least 1 putative bulge has been found

                numb_bulges = 0  # it is the iterator that should count how many bulges there are
                numb_internal_loop = 0  # it is the iterator that should count how many internal loop there are

                # Analyze every putative bulge found in order to recognize real bulges from internal loops
                for putative_bulge in putative_bulges:
                    # matched_putative_bulge = putative_bulge.group()  # group() function returns the string matched by the RE
                    matched_starting_index = putative_bulge.start()  # star() function returns the starting position of the match. It corresponds with the starting index of the putative bulge
                    matched_ending_index = putative_bulge.end() - 1  # end() function returns the ending position of the match. We have to decrease the number of 1 in order to obtain the real ending index

                    # Search if the bases around the putative bulge are paired in the corresponding position of the other filament
                    for sublist in indexes_paired_bases:  # 'indexes_paired_bases' is a list that contain all the indexes of the bases that compose all the pairings

                        # for putative bulges at 5':
                        if sublist[0] == matched_starting_index:
                            paired_base_starting_index = sublist[1]
                        if sublist[0] == matched_ending_index:
                            paired_base_ending_index = sublist[1]

                        # for putative bulges at 3':
                        if sublist[1] == matched_starting_index:
                            paired_base_starting_index = sublist[0]
                        if sublist[1] == matched_ending_index:
                            paired_base_ending_index = sublist[0]

                    ### DEFINE REAL BULGES
                    # If the indexes of the paired bases are consecutive (so, the absolute value of the subtraction of them is 1), it means that the putative bulge is a real bulge. Otherwise it is an internal loop.

                    if abs(paired_base_starting_index - paired_base_ending_index) == 1:  # abs() function is used to return the absolute value

                        bulge = putative_bulge

                        numb_bulges = numb_bulges + 1  # to count the number of bulges: if putative_bulge == bulge, increase the number of bulges

                        ### BULGES INDEXES

                        # bulge_starting_index = bulge.start() + 1  # to know the starting index of the first unpaired base (dot), we have to increase the number of 1
                        # bulge_ending_index = bulge.end() - 2  # to know the ending index of the last unpaired base (dot), we have to decrease the number of 2
                        # bulge_indexes = str(bulge_starting_index) + ":" + str(bulge_ending_index)

                    ### DEFINE INTERNAL LOOP

                    else:
                        internal_loop = putative_bulge

                        numb_internal_loop = numb_internal_loop + 1  # to count the number of internal loop: if putative_bulge == internal_loop, increase the number of internal loop

                        ### INTERNAL LOOP INDEXES

                        # internal_loop_starting_index = internal_loop.start() + 1
                        # internal_loop_ending_index = internal_loop.end() - 2
                        # internal_loop_indexes = str(internal_loop_starting_index) + ":" + str(internal_loop_ending_index)

                ### TOTAL NUMBER OF BULGES

                # numb_bulges

                ### TOTAL NUMBER OF INTERNAL LOOP

                # numb_internal_loop/2  # the unpaired bases of the internal loop are count once for every filament. So, we have to divide the number of internal loop by 2 in order to obtain the real number of internal loops present.

                # Working with 'energy_values_str'

                ### ENERGY VALUES OF THE UNTOUCHED SEQUENCE

                # energy_values_str

        # ================================== Working with sequences with at least 1 loop

        complete_list_loops = []  # these are lists that will contain the structures and some other important information of the loop(s) to be used for the analysis of the features
        loop_length_information = []

        count = 0  # it allows to count how many iteration the following for-loop does

        if len(find_loops) > 0:  # if we have at least 1 loop!

            # Before starting, we have to remove any possible "open" loop in the dot-bracket notation or dots (not paired bases) at the ends of dot-bracket notation.

            for loop_found in find_loops[::-1]:  # inverting the 'find_loops' list with the command [::-1]: in this way we can start the analysis of the loops from the 3' to 5' of the sequence (and not vice versa)

                # see the position of the loop
                len_loop = loop_found.count(".")  # count the length of the loop (how many dots are present in the loop): it is necessary in order to measure the final index of the loop
                starting_index_loop = dots_brackets_str.index(loop_found) + 1  # it is the index of the first dot that composes the loop. Add + 1 because we do not have to take into account '(' that enclose the loop
                final_index_loop = starting_index_loop + (len_loop - 1)  # it is the index of the last dot that composes the loop

                # divide the notation in 2 parts according to the position of the loop found
                dots_brackets_str_3_prime = dots_brackets_str[final_index_loop+1:]
                dots_brackets_str_5_prime = dots_brackets_str[:starting_index_loop]

                # going from the 3' and then 5' after the loop and search the presence of respectively a '(' and a ')': it means that there there is a possible open loop or another unwanted structure that we want to delete
                final_index_dot_bracket_notation_3_prime = len(dots_brackets_str_3_prime) - 1

                for component_3_prime in dots_brackets_str_3_prime:
                    if component_3_prime == "(":
                        termination_condition_3_prime_index = dots_brackets_str_3_prime.index(component_3_prime)  # it is the position where we will start deleting the structure
                        dots_brackets_str_3_prime = dots_brackets_str_3_prime[:termination_condition_3_prime_index]  # overwrite the notation without the problematic structure
                        final_index_dot_bracket_notation_3_prime = len(dots_brackets_str_3_prime) - 1
                        break

                notation_solved_3_prime = dots_brackets_str_3_prime
                while dots_brackets_str_3_prime[final_index_dot_bracket_notation_3_prime] == ".":  # if present, remove the dots at the end of the notation
                    dots_brackets_str_3_prime = dots_brackets_str_3_prime[:final_index_dot_bracket_notation_3_prime]  # overwrite the notation
                    final_index_dot_bracket_notation_3_prime = final_index_dot_bracket_notation_3_prime - 1  # define the new final index of the overwritten notation

                for component_5_prime in dots_brackets_str_5_prime[::-1]:  # inverting the part of the notation in order to analyze it from 3' to 5' (and not vice versa)
                    if component_5_prime == ")":
                        termination_condition_5_prime_index = dots_brackets_str_5_prime[::-1].index(component_5_prime)
                        dots_brackets_str_5_prime = dots_brackets_str_5_prime[::-1][:termination_condition_5_prime_index]
                        dots_brackets_str_5_prime = dots_brackets_str_5_prime[::-1]   # re-inverting the part of the notation
                        break

                notation_solved_5_prime = dots_brackets_str_5_prime
                while dots_brackets_str_5_prime[0] == ".":  # remove dots in this case at the beginning of the notation (because the notation has been already re-inverted)
                    dots_brackets_str_5_prime = dots_brackets_str_5_prime[1:]

                dots_brackets_str = dots_brackets_str_5_prime + "." * len_loop + dots_brackets_str_3_prime  # it is the "NEW" DOTS BRACKETS NOTATION (whithout any problematic structures)
                starting_index_new_dots_bracket_notation = untouched_dots_brackets.index(dots_brackets_str)
                ending_index_new_dots_bracket_notation_ending = starting_index_new_dots_bracket_notation + len(dots_brackets_str)
                sequence_str = untouched_sequence[starting_index_new_dots_bracket_notation:ending_index_new_dots_bracket_notation_ending]  # it is the "NEW" SEQUENCE
                total_length_seq = len(sequence_str)  # total number of bases present in the sequence

                loop_to_be_chosen = []  # to create sub lists of the list 'complete_list_loops': all the information of every loop (if present more than 1 loop) are contained in different sublists
                loop_to_be_chosen.append(loop_found)
                loop_to_be_chosen.append(dots_brackets_str)
                loop_to_be_chosen.append(starting_index_new_dots_bracket_notation)
                loop_to_be_chosen.append(ending_index_new_dots_bracket_notation_ending)
                loop_to_be_chosen.append(sequence_str)
                loop_to_be_chosen.append(total_length_seq)

                complete_list_loops.append(loop_to_be_chosen)  # the list contains all the structures found with their information (relative dots brackets notation, starting and ending index, sequence and length)

                loop_length_information.append(total_length_seq)  # here we store only the length of all the selected structures

                # now we want to overwrite our dot-bracket notation without the structure already solved
                notation_solved = notation_solved_5_prime + "." * len_loop + notation_solved_3_prime
                differences = [notation for notation in difflib.ndiff(untouched_dots_brackets, notation_solved) if notation[0] != ' ']  # find the part of the notation not already solved, substracting from the complete dot-bracket notation the part already solved
                dots_brackets_str = " "
                dots_brackets_str = dots_brackets_str.join(differences).replace("-", "").replace(" ", "")  # difflib.ndiff() function worked with list, so now we have to re-convert the list into string. Overwrite the new dot-bracket notation
                count = count + 1

                if count == len(find_loops):  # it means that we have analyzed all the loops and their notations.
                    break

        # We have to consider only the longest structure as input for the analysis of the following features

            max_length = max(loop_length_information)  # find which is the longest structure

        # N.B.: of couse all these pieces of code are more meaniful in the case of presence of more than one loop.. but in any case work also in the case of only one loop

            for every_structure in complete_list_loops:
                # search and selecting in the length information (that is the 5th element of every sublist) only the structure that have the corresponding 'max_length'
                if every_structure[5] == max_length:
                    loop_of_interest = every_structure[0]
                    dots_brackets_str = every_structure[1]
                    starting_index_new_dots_bracket_notation = every_structure[2]
                    ending_index_new_dots_bracket_notation_ending = every_structure[3]
                    sequence_str = every_structure[4]
                    total_length_seq = every_structure[5]
                    # these are respectively the loop, the dot bracket notation, the sequence string and length of the longest structure that will be used for the following analysis of the features


        ###  INDEXES OF THE LOOP OF INTEREST
        # We report the indexes of the only one loop considered. They are the indexes of the loop in the dot-bracket notation considered ('dots_brackets_str') and not in the original one ('untouched_dots_brackets')
            starting_index_loop_of_interest = dots_brackets_str.index(loop_of_interest) + 1  # it is the index of the first dot that composes the loop. Add + 1 because we do not have to take into account '(' that enclose the loop
            length_loop_of_interest = loop_of_interest.count(".")  # count the length of the loop (how many dots are present in the loop): it is necessary in order to measure the final index of the loop
            ending_index_loop_of_interest = starting_index_loop_of_interest + (length_loop_of_interest - 1)  # it is the index of the last dot that composes the loop


        # We have to start from the loop and proceed till the ends, analyzing if the bases are paired.

            initial_index = starting_index_loop_of_interest - 1  # it is the starting index of the first '(' of the loop , used to count the positions from the 5' of the loop
            final_index = ending_index_loop_of_interest + 1  # it is the ending index of the last ')' of the loop, used to count the positions from the 3' end of the loop

            couples = []  # list that will contain all the couples found

            indexes_paired_bases = []  # list that will contain all the indexes of all the couples found

            # for element in dots_brackets_str:
            while initial_index > 0 and final_index < len(dots_brackets_str):

                if dots_brackets_str[initial_index] == "(" and dots_brackets_str[final_index] == ")":  # there is a pairing

                    couple_of_bases = sequence_str[initial_index] + sequence_str[final_index]  # extract the first base that corresponds to '(' and the last base that corresponds to ')'
                    couples.append(couple_of_bases)

                    # save also the indexes of the paired bases: they should be used later to recognize bulges from internal loops!
                    indexes_to_be_added = []  # sub list that will contain the indexes of the 2 bases that compose the pairing
                    indexes_to_be_added.append(initial_index)  # index of the base at 5'
                    indexes_to_be_added.append(final_index)  # index of the base at 3'
                    indexes_paired_bases.append(indexes_to_be_added)

                    # to move farther from the centre, the initial index should be decreased of 1, instead the final index should be increased of 1
                    initial_index = initial_index - 1
                    final_index = final_index + 1

                elif dots_brackets_str[initial_index] == "." and dots_brackets_str[final_index] == ")":  # there isn't a pairing

                    # the initial index should be decreased of 1 to move farther from the centre, instead the final index remains fixed
                    initial_index = initial_index - 1

                elif dots_brackets_str[initial_index] == "(" and dots_brackets_str[final_index] == ".":  # there isn't a pairing

                    # the initial index remains fixed, instead the final index should be increased of 1 to move farther from the centre
                    final_index = final_index + 1

                elif dots_brackets_str[initial_index] == "." and dots_brackets_str[final_index] == ".":  # there isn't a pairing

                    # farther from the centre the initial index should be decreased of 1, instead the final index should be increased of 1
                    initial_index = initial_index - 1
                    final_index = final_index + 1

            # Now make the structure simmetric: using the initial_index and final_index, overwrite the dot-bracket notation, the sequence and re-measure the lenght
            dots_brackets_str = dots_brackets_str[initial_index:final_index]
            sequence_str = sequence_str[initial_index:final_index]
            total_length_seq = len(sequence_str)


            # From now working with 'sequence_str': if more than 1 loop is present, it corresponds with the truncated sequence overwritten

            ### DINUCLEOTIDE FREQUENCIES

            sequence_split_in_two = textwrap.wrap(sequence_str, 2)  # split the sequence in couples

            total_count_AA = sequence_split_in_two.count("AA")  # find the number of AA occurences
            frequency_AA = total_count_AA / float(total_length_seq / 2)  # float() function returns a floating point number

            total_count_AG = sequence_split_in_two.count("AG")
            frequency_AG = total_count_AG/float(total_length_seq / 2)

            total_count_AC = sequence_split_in_two.count("AC")
            frequency_AC = total_count_AC/float(total_length_seq / 2)

            total_count_AU = sequence_split_in_two.count("AU")
            frequency_AU = total_count_AU/float(total_length_seq / 2)

            total_count_GA = sequence_split_in_two.count("GA")
            frequency_GA = total_count_GA/float(total_length_seq / 2)

            total_count_GG = sequence_split_in_two.count("GG")
            frequency_GG = total_count_GG/float(total_length_seq / 2)

            total_count_GC = sequence_split_in_two.count("GC")
            frequency_GC = total_count_GC/float(total_length_seq / 2)

            total_count_GU = sequence_split_in_two.count("GU")
            frequency_GU = total_count_GU/float(total_length_seq / 2)

            total_count_CA = sequence_split_in_two.count("CA")
            frequency_CA = total_count_CA/float(total_length_seq / 2)

            total_count_CG = sequence_split_in_two.count("CG")
            frequency_CG = total_count_CG/float(total_length_seq / 2)

            total_count_CC = sequence_split_in_two.count("CC")
            frequency_CC = total_count_CC/float(total_length_seq / 2)

            total_count_CU = sequence_split_in_two.count("CU")
            frequency_CU = total_count_CU/float(total_length_seq / 2)

            total_count_UA = sequence_split_in_two.count("UA")
            frequency_UA = total_count_UA/float(total_length_seq / 2)

            total_count_UG = sequence_split_in_two.count("UG")
            frequency_UG = total_count_UG/float(total_length_seq / 2)

            total_count_UC = sequence_split_in_two.count("UC")
            frequency_UC = total_count_UC/float(total_length_seq / 2)

            total_count_UU = sequence_split_in_two.count("UU")
            frequency_UU = total_count_UU/float(total_length_seq / 2)

            ### TRINUCLEOTIDE FREQUENCIES

            sequence_split_in_three = textwrap.wrap(sequence_str, 3)

            total_count_AAA = sequence_split_in_three.count("AAA")
            frequency_AAA = total_count_AAA / float(total_length_seq / 3)

            total_count_AAG = sequence_split_in_three.count("AAG")
            frequency_AAG = total_count_AAG / float(total_length_seq / 3)

            total_count_AAC = sequence_split_in_three.count("AAC")
            frequency_AAC = total_count_AAC / float(total_length_seq / 3)

            total_count_AAU = sequence_split_in_three.count("AAU")
            frequency_AAU = total_count_AAU / float(total_length_seq / 3)

            total_count_AGA = sequence_split_in_three.count("AGA")
            frequency_AGA = total_count_AGA / float(total_length_seq / 3)

            total_count_AGG = sequence_split_in_three.count("AGG")
            frequency_AGG = total_count_AGG / float(total_length_seq / 3)

            total_count_AGC = sequence_split_in_three.count("AGC")
            frequency_AGC = total_count_AGC / float(total_length_seq / 3)

            total_count_AGU = sequence_split_in_three.count("AGU")
            frequency_AGU = total_count_AGU / float(total_length_seq / 3)

            total_count_ACA = sequence_split_in_three.count("ACA")
            frequency_ACA = total_count_ACA / float(total_length_seq / 3)

            total_count_ACG = sequence_split_in_three.count("ACG")
            frequency_ACG = total_count_ACG / float(total_length_seq / 3)

            total_count_ACC = sequence_split_in_three.count("ACC")
            frequency_ACC = total_count_ACC / float(total_length_seq / 3)

            total_count_ACU = sequence_split_in_three.count("ACU")
            frequency_ACU = total_count_ACU / float(total_length_seq / 3)

            total_count_AUA = sequence_split_in_three.count("AUA")
            frequency_AUA = total_count_AUA / float(total_length_seq / 3)

            total_count_AUG = sequence_split_in_three.count("AUG")
            frequency_AUG = total_count_AUG / float(total_length_seq / 3)

            total_count_AUC = sequence_split_in_three.count("AUC")
            frequency_AUC = total_count_AUC / float(total_length_seq / 3)

            total_count_AUU = sequence_split_in_three.count("AUU")
            frequency_AUU = total_count_AUU / float(total_length_seq / 3)

            total_count_GAA = sequence_split_in_three.count("GAA")
            frequency_GAA = total_count_GAA / float(total_length_seq / 3)

            total_count_GAG = sequence_split_in_three.count("GAG")
            frequency_GAG = total_count_GAG / float(total_length_seq / 3)

            total_count_GAC = sequence_split_in_three.count("GAC")
            frequency_GAC = total_count_GAC / float(total_length_seq / 3)

            total_count_GAU = sequence_split_in_three.count("GAU")
            frequency_GAU = total_count_GAU / float(total_length_seq / 3)

            total_count_GGA = sequence_split_in_three.count("GGA")
            frequency_GGA = total_count_GGA / float(total_length_seq / 3)

            total_count_GGG = sequence_split_in_three.count("GGG")
            frequency_GGG = total_count_GGG / float(total_length_seq / 3)

            total_count_GGC = sequence_split_in_three.count("GGC")
            frequency_GGC = total_count_GGC / float(total_length_seq / 3)

            total_count_GGU = sequence_split_in_three.count("GGU")
            frequency_GGU = total_count_GGU / float(total_length_seq / 3)

            total_count_GCA = sequence_split_in_three.count("GCA")
            frequency_GCA = total_count_GCA / float(total_length_seq / 3)

            total_count_GCG = sequence_split_in_three.count("GCG")
            frequency_GCG = total_count_GCG / float(total_length_seq / 3)

            total_count_GCC = sequence_split_in_three.count("GCC")
            frequency_GCC = total_count_GCC / float(total_length_seq / 3)

            total_count_GCU = sequence_split_in_three.count("GCU")
            frequency_GCU = total_count_GCU / float(total_length_seq / 3)

            total_count_GUA = sequence_split_in_three.count("GUA")
            frequency_GUA = total_count_GUA / float(total_length_seq / 3)

            total_count_GUG = sequence_split_in_three.count("GUG")
            frequency_GUG = total_count_GUG / float(total_length_seq / 3)

            total_count_GUC = sequence_split_in_three.count("GUC")
            frequency_GUC = total_count_GUC / float(total_length_seq / 3)

            total_count_GUU = sequence_split_in_three.count("GUU")
            frequency_GUU = total_count_GUU / float(total_length_seq / 3)

            total_count_CAA = sequence_split_in_three.count("CAA")
            frequency_CAA = total_count_CAA / float(total_length_seq / 3)

            total_count_CAG = sequence_split_in_three.count("CAG")
            frequency_CAG = total_count_CAG / float(total_length_seq / 3)

            total_count_CAC = sequence_split_in_three.count("CAC")
            frequency_CAC = total_count_CAC / float(total_length_seq / 3)

            total_count_CAU = sequence_split_in_three.count("CAU")
            frequency_CAU = total_count_CAU / float(total_length_seq / 3)

            total_count_CGA = sequence_split_in_three.count("CGA")
            frequency_CGA = total_count_CGA / float(total_length_seq / 3)

            total_count_CGG = sequence_split_in_three.count("CGG")
            frequency_CGG = total_count_CGG / float(total_length_seq / 3)

            total_count_CGC = sequence_split_in_three.count("CGC")
            frequency_CGC = total_count_CGC / float(total_length_seq / 3)

            total_count_CGU = sequence_split_in_three.count("CGU")
            frequency_CGU = total_count_CGU / float(total_length_seq / 3)

            total_count_CCA = sequence_split_in_three.count("CCA")
            frequency_CCA = total_count_CCA / float(total_length_seq / 3)

            total_count_CCG = sequence_split_in_three.count("CCG")
            frequency_CCG = total_count_CCG / float(total_length_seq / 3)

            total_count_CCC = sequence_split_in_three.count("CCC")
            frequency_CCC = total_count_CCC / float(total_length_seq / 3)

            total_count_CCU = sequence_split_in_three.count("CCU")
            frequency_CCU = total_count_CCU / float(total_length_seq / 3)

            total_count_CUA = sequence_split_in_three.count("CUA")
            frequency_CUA = total_count_CUA / float(total_length_seq / 3)

            total_count_CUG = sequence_split_in_three.count("CUG")
            frequency_CUG = total_count_CUG / float(total_length_seq / 3)

            total_count_CUC = sequence_split_in_three.count("CUC")
            frequency_CUC = total_count_CUC / float(total_length_seq / 3)

            total_count_CUU = sequence_split_in_three.count("CUU")
            frequency_CUU = total_count_CUU / float(total_length_seq / 3)

            total_count_UAA = sequence_split_in_three.count("UAA")
            frequency_UAA = total_count_UAA / float(total_length_seq / 3)

            total_count_UAG = sequence_split_in_three.count("UAG")
            frequency_UAG = total_count_UAG / float(total_length_seq / 3)

            total_count_UAC = sequence_split_in_three.count("UAC")
            frequency_UAC = total_count_UAC / float(total_length_seq / 3)

            total_count_UAU = sequence_split_in_three.count("UAU")
            frequency_UAU = total_count_UAU / float(total_length_seq / 3)

            total_count_UGA = sequence_split_in_three.count("UGA")
            frequency_UGA = total_count_UGA / float(total_length_seq / 3)

            total_count_UGG = sequence_split_in_three.count("UGG")
            frequency_UGG = total_count_UGG / float(total_length_seq / 3)

            total_count_UGC = sequence_split_in_three.count("UGC")
            frequency_UGC = total_count_UGC / float(total_length_seq / 3)

            total_count_UGU = sequence_split_in_three.count("UGU")
            frequency_UGU = total_count_UGU / float(total_length_seq / 3)

            total_count_UCA = sequence_split_in_three.count("UCA")
            frequency_UCA = total_count_UCA / float(total_length_seq / 3)

            total_count_UCG = sequence_split_in_three.count("UCG")
            frequency_UCG = total_count_UCG / float(total_length_seq / 3)

            total_count_UCC = sequence_split_in_three.count("UCC")
            frequency_UCC = total_count_UCC / float(total_length_seq / 3)

            total_count_UCU = sequence_split_in_three.count("UCU")
            frequency_UCU = total_count_UCU / float(total_length_seq / 3)

            total_count_UUA = sequence_split_in_three.count("UUA")
            frequency_UUA = total_count_UUA / float(total_length_seq / 3)

            total_count_UUG = sequence_split_in_three.count("UUG")
            frequency_UUG = total_count_UUG / float(total_length_seq / 3)

            total_count_UUC = sequence_split_in_three.count("UUC")
            frequency_UUC = total_count_UUC / float(total_length_seq / 3)

            total_count_UUU = sequence_split_in_three.count("UUU")
            frequency_UUU = total_count_UUU / float(total_length_seq / 3)

            # From now on working with 'dots_brackets_str': if more than 1 loop is present, it corresponds with the truncated dot-bracket notation overwritten

            ### FREQUENCIES OF TRIPLET ELEMENTS

            dots_brackets_split_in_three = textwrap.wrap(dots_brackets_str, 3)

            find_bracket_bracket_bracket_5_prime = dots_brackets_split_in_three.count("(((")
            find_bracket_bracket_bracket_3_prime = dots_brackets_split_in_three.count(")))")
            total_count_bracket_bracket_bracket = find_bracket_bracket_bracket_5_prime + find_bracket_bracket_bracket_3_prime
            frequency_bracket_bracket_bracket = total_count_bracket_bracket_bracket / float(total_length_seq / 3)

            find_bracket_bracket_dot_5_prime = dots_brackets_split_in_three.count("((.")
            find_bracket_bracket_dot_3_prime = dots_brackets_split_in_three.count(")).")
            total_count_bracket_bracket_dot = find_bracket_bracket_dot_5_prime + find_bracket_bracket_dot_3_prime
            frequency_bracket_bracket_dot = total_count_bracket_bracket_dot / float(total_length_seq / 3)

            find_bracket_dot_bracket_5_prime = dots_brackets_split_in_three.count("(.(")
            find_bracket_dot_bracket_3_prime = dots_brackets_split_in_three.count(").)")
            total_count_bracket_dot_bracket = find_bracket_dot_bracket_5_prime + find_bracket_dot_bracket_3_prime
            frequency_bracket_dot_bracket = total_count_bracket_dot_bracket / float(total_length_seq / 3)

            find_bracket_dot_dot_5_prime = dots_brackets_split_in_three.count("(..")
            find_bracket_dot_dot_3_prime = dots_brackets_split_in_three.count(")..")
            total_count_bracket_dot_dot = find_bracket_dot_dot_5_prime + find_bracket_dot_dot_3_prime
            frequency_bracket_dot_dot = total_count_bracket_dot_dot / float(total_length_seq / 3)

            find_dot_bracket_bracket_5_prime = dots_brackets_split_in_three.count(".((")
            find_dot_bracket_bracket_3_prime = dots_brackets_split_in_three.count(".))")
            total_count_dot_bracket_bracket = find_dot_bracket_bracket_5_prime + find_dot_bracket_bracket_3_prime
            frequency_dot_bracket_bracket = total_count_dot_bracket_bracket / float(total_length_seq / 3)

            find_dot_bracket_dot_5_prime = dots_brackets_split_in_three.count(".(.")
            find_dot_bracket_dot_3_prime = dots_brackets_split_in_three.count(".).")
            total_count_dot_bracket_dot = find_dot_bracket_dot_5_prime + find_dot_bracket_dot_3_prime
            frequency_dot_bracket_dot = total_count_dot_bracket_dot / float(total_length_seq / 3)

            find_dot_dot_bracket_5_prime = dots_brackets_split_in_three.count("..(")
            find_dot_dot_bracket_3_prime = dots_brackets_split_in_three.count("..)")
            total_count_dot_dot_bracket = find_dot_dot_bracket_5_prime + find_dot_dot_bracket_3_prime
            frequency_dot_dot_bracket = total_count_dot_dot_bracket / float(total_length_seq / 3)

            total_count_dot_dot_dot = dots_brackets_split_in_three.count("...")
            frequency_dot_dot_dot = total_count_dot_dot_dot / float(total_length_seq / 3)

            #

            total_count_bracket_bracket_bracket_2 = dots_brackets_split_in_three.count("(()")
            frequency_bracket_bracket_bracket_2 = total_count_bracket_bracket_bracket_2 / float(total_length_seq / 3)

            total_count_bracket_bracket_bracket_3 = dots_brackets_split_in_three.count("()(")
            frequency_bracket_bracket_bracket_3 = total_count_bracket_bracket_bracket_3 / float(total_length_seq / 3)

            total_count_bracket_bracket_bracket_4 = dots_brackets_split_in_three.count("())")
            frequency_bracket_bracket_bracket_4 = total_count_bracket_bracket_bracket_4 / float(total_length_seq / 3)

            total_count_bracket_bracket_bracket_5 = dots_brackets_split_in_three.count(")((")
            frequency_bracket_bracket_bracket_5 = total_count_bracket_bracket_bracket_5 / float(total_length_seq / 3)

            total_count_bracket_bracket_bracket_6 = dots_brackets_split_in_three.count(")()")
            frequency_bracket_bracket_bracket_6 = total_count_bracket_bracket_bracket_6 / float(total_length_seq / 3)

            total_count_bracket_bracket_bracket_7 = dots_brackets_split_in_three.count("))(")
            frequency_bracket_bracket_bracket_7 = total_count_bracket_bracket_bracket_7 / float(total_length_seq / 3)

            total_count_bracket_bracket_dot_2 = dots_brackets_split_in_three.count("().")
            frequency_bracket_bracket_dot_2 = total_count_bracket_bracket_dot_2 / float(total_length_seq / 3)

            total_count_bracket_bracket_dot_3 = dots_brackets_split_in_three.count(")(.")
            frequency_bracket_bracket_dot_3 = total_count_bracket_bracket_dot_3 / float(total_length_seq / 3)

            total_count_bracket_dot_bracket_2 = dots_brackets_split_in_three.count("(.)")
            frequency_bracket_dot_bracket_2 = total_count_bracket_dot_bracket_2 / float(total_length_seq / 3)

            total_count_bracket_dot_bracket_3 = dots_brackets_split_in_three.count(").(")
            frequency_bracket_dot_bracket_3 = total_count_bracket_dot_bracket_3 / float(total_length_seq / 3)

            total_count_dot_bracket_bracket_2 = dots_brackets_split_in_three.count(".()")
            frequency_dot_bracket_bracket_2 = total_count_dot_bracket_bracket_2 / float(total_length_seq / 3)

            total_count_dot_bracket_bracket_3 = dots_brackets_split_in_three.count(".)(")
            frequency_dot_bracket_bracket_3 = total_count_dot_bracket_bracket_3 / float(total_length_seq / 3)

            ### TOTAL NUMBER OF NOT PAIRED BASES

            numb_not_paired_bases = dots_brackets_str.count('.')  # not paired bases correspond to dots. Count() function counts how many dots there are

            ### TOTAL NUMBER OF LOOPS IN THE UNTOUCHED DOT-BRACKET NOTATION
            # We obtain the total number of loops from the list that contains all the loops found in the untouched dot-bracket notation ('find_loops' list)
            numb_loops = len(find_loops)  # len() function count how many elements there are in the list: it corresponds to how many loops there are

            ### TOTAL NUMBER OF LOOPS DEFINING THEM AS ((.{3,}))
            # Defining loops as '((...))', find the number of loops in the untouched dot-bracket notation:
            find_loop_3_dots_2_parentesis = re.findall(regex_loop_3_dots_2_parentesis, untouched_dots_brackets)
            numb_loops_3_dots_2_parentesis = len(find_loop_3_dots_2_parentesis)

            ### TOTAL NUMBER OF LOOPS DEFINING THEM AS (.{4,})
            # Defining loops as '(....)', find the number of loops in the untouched dot-bracket notation:
            find_loop_4_dots_1_parentesis = re.findall(regex_loop_4_dots_1_parentesis, untouched_dots_brackets)
            numb_loops_4_dots_1_parentesis = len(find_loop_4_dots_1_parentesis)

            ### TOTAL NUMBER OF LOOPS DEFINING THEM AS ((.{4,}))
            # Defining loops as '((....))', find the number of loops in the untouched dot-bracket notation:
            find_loop_4_dots_2_parentesis = re.findall(regex_loop_4_dots_2_parentesis, untouched_dots_brackets)
            numb_loops_4_dots_2_parentesis = len(find_loop_4_dots_2_parentesis)

            ### TOTAL NUMBER OF LOOPS DEFINING THEM AS (.{5,})
            # Defining loops as '(.....)', find the number of loops in the untouched dot-bracket notation:
            find_loop_5_dots_1_parentesis = re.findall(regex_loop_5_dots_1_parentesis, untouched_dots_brackets)
            numb_loops_5_dots_1_parentesis = len(find_loop_5_dots_1_parentesis)

            ### TOTAL NUMBER OF LOOPS DEFINING THEM AS ((.{5,}))
            # Defining loops as '((....))', find the number of loops in the untouched dot-bracket notation:
            find_loop_5_dots_2_parentesis = re.findall(regex_loop_5_dots_2_parentesis, untouched_dots_brackets)
            numb_loops_5_dots_2_parentesis = len(find_loop_5_dots_2_parentesis)

            ###  INDEXES OF THE LOOP OF INTEREST
            # We report the indexes of the only one loop considered. They are the indexes of the loop in the dot-bracket notation considered ('dots_brackets_str') and not in the original one ('untouched_dots_brackets')
            starting_index_loop_of_interest = dots_brackets_str.index(loop_of_interest) + 1  # it is the index of the first dot that composes the loop. Add + 1 because we do not have to take into account '(' that enclose the loop
            length_loop_of_interest = loop_of_interest.count(".")  # count the length of the loop (how many dots are present in the loop): it is necessary in order to measure the final index of the loop
            ending_index_loop_of_interest = starting_index_loop_of_interest + (length_loop_of_interest - 1)  # it is the index of the last dot that composes the loop

            ###  CENTRALITY OF THE LOOP OF INTEREST

            loop_centrality = starting_index_loop_of_interest/float(total_length_seq)

            ### DIMENSION OF THE LOOP OF INTEREST

            # length_loop_of_interest

            ### TOTAL NUMBER OF PAIRED BASES

            numb_paired_bases = total_length_seq - numb_not_paired_bases  # it corresponds to the total number of bases minus the number of bases not paired

            ### FREQUENCY OF PAIRED COUPLES
            # We have already calculated the occurences of the paired couples: they are all together strored in 'couples' list.
            # Now we have to divide them according to the single couple in 6 different sublists.

            all_CG_couples = ""  # it is an empty string that will contain all the couples CG
            all_GC_couples = ""
            all_UA_couples = ""
            all_AU_couples = ""
            all_UG_couples = ""
            all_GU_couples = ""

            for couple in couples:  # put every couple of bases found in the corresponding string
                if couple == "CG":
                    all_CG_couples = all_CG_couples + couple

                elif couple == "GC":
                    all_GC_couples = all_GC_couples + couple

                elif couple == "UA":
                    all_UA_couples = all_UA_couples + couple

                elif couple == "AU":
                    all_AU_couples = all_AU_couples + couple

                elif couple == "UG":
                    all_UG_couples = all_UG_couples + couple

                elif couple == "GU":
                    all_GU_couples = all_GU_couples + couple

            count_CG = len(all_CG_couples)/2  # count how many CG occurences there are in the 'all_CG_couples' string. We use len() function to measure the bases and then divide by 2 to define how many couples there are
            count_GC = len(all_GC_couples)/2
            count_UA = len(all_UA_couples)/2
            count_AU = len(all_AU_couples)/2
            count_UG = len(all_UG_couples)/2
            count_GU = len(all_GU_couples)/2

            total_numb_pairings = count_CG + count_GC + count_UA + count_AU + count_UG + count_GU  # measure the total number of pairings found

            if total_numb_pairings <= 0:
                warning_list.append(header_str)
                continue  # skip it

            else:

                frequency_couple_CG = count_CG / float(total_numb_pairings)  # measure the number of occurence of CG over the total number of pairings
                frequency_couple_GC = count_GC / float(total_numb_pairings)
                frequency_couple_UA = count_UA / float(total_numb_pairings)
                frequency_couple_AU = count_AU / float(total_numb_pairings)
                frequency_couple_UG = count_UG / float(total_numb_pairings)
                frequency_couple_GU = count_GU / float(total_numb_pairings)

            ### STEM SEQUENCE
            # We have to start to analyze the stem sequence from the ends of the loop and move farther from it.

            stem_sequence_5_prime = ""  # string that will contain the stem sequence at 5'
            stem_sequence_3_prime = ""

            five_prime_index = starting_index_loop_of_interest - 1  # it corresponds to 5' index of the loop. To move farther from the loop, it should be decreased of 1.
            three_prime_index = ending_index_loop_of_interest + 1  # it corresponds to 3' index of the loop. To move farther from the loop, it should be increased of 1.


            while five_prime_index != - 1 and three_prime_index != total_length_seq:  # this condition allows to span the entire sequence and stop at the end of the sequence

                if dots_brackets_str[five_prime_index] == "(" and dots_brackets_str[three_prime_index] == ")":  # there is a pairing

                    # Extraction of the base both at 5' and 3' and put respectively in 'stem_sequence_5_prime' string and 'stem_sequence_3_prime' string.
                    stem_sequence_5_prime = stem_sequence_5_prime + sequence_str[five_prime_index]
                    stem_sequence_3_prime = stem_sequence_3_prime + sequence_str[three_prime_index]

                    # to move farther from the loop, the 5' index should be decreased of 1, instead the 3' index should be increased of 1
                    five_prime_index = five_prime_index - 1
                    three_prime_index = three_prime_index + 1

                elif dots_brackets_str[five_prime_index] == "." and dots_brackets_str[three_prime_index] == ")":  # there isn't a pairing

                    # the 5' index should be decreased of 1 to move farther from the loop, instead the 3' index remains fixed
                    five_prime_index = five_prime_index - 1

                elif dots_brackets_str[five_prime_index] == "(" and dots_brackets_str[three_prime_index] == ".":  # there isn't a pairing

                    # the 5' index remains fixed, instead the 3' index should be increased of 1 to move farther from the loop
                    three_prime_index = three_prime_index + 1

                elif dots_brackets_str[five_prime_index] == "." and dots_brackets_str[three_prime_index] == ".":  # there isn't a pairing

                    five_prime_index = five_prime_index - 1
                    three_prime_index = three_prime_index + 1

                elif dots_brackets_str[five_prime_index] == ")" or dots_brackets_str[three_prime_index] == "(":  # the stem is finished

                    five_prime_index = five_prime_index - 1
                    three_prime_index = three_prime_index + 1

            # starting from the ends of the loop, the stem sequences already extracted (both the one at 5' and 3') should be reversed
            reversed_stem_sequence_5_prime = ''.join(reversed(stem_sequence_5_prime))

            reversed_stem_sequence_3_prime = ''.join(reversed(stem_sequence_3_prime))

            ### DEFINE PUTATIVE BULGES

            putative_bulges = re.finditer(regex_putative_bulge, dots_brackets_str)  # finditer() function returns a tuple with an iterator yielding match-object instances over all non-overlapping matches for the RE pattern in string. The string is scanned left-to-right, and matches are returned in the order found.

            if (putative_bulges):  # if at least 1 putative bulge has been found

                numb_bulges = 0  # it is the iterator that should count how many bulges there are
                numb_internal_loop = 0  # it is the iterator that should count how many internal loop there are

                # Analyze every putative bulge found in order to recognize real bulges from internal loops
                for putative_bulge in putative_bulges:
                    #matched_putative_bulge = putative_bulge.group()  # group() function returns the string matched by the RE
                    matched_starting_index = putative_bulge.start() # star() function returns the starting position of the match. It corresponds with the starting index of the putative bulge
                    matched_ending_index = putative_bulge.end() - 1 # end() function returns the ending position of the match. We have to decrease the number of 1 in order to obtain the real ending index

                    # Search if the bases around the putative bulge are paired in the corresponding position of the other filament
                    for sublist in indexes_paired_bases:  # 'indexes_paired_bases' is a list that contain all the indexes of the bases that compose all the pairings

                        # for putative bulges at 5':
                        if sublist[0] == matched_starting_index:
                            paired_base_starting_index = sublist[1]
                        if sublist[0] == matched_ending_index:
                            paired_base_ending_index = sublist[1]

                        # for putative bulges at 3':
                        if sublist[1] == matched_starting_index:
                            paired_base_starting_index = sublist[0]
                        if sublist[1] == matched_ending_index:
                            paired_base_ending_index = sublist[0]

                    ### DEFINE REAL BULGES
                    # If the indexes of the paired bases are consecutive (so, the absolute value of the subtraction of them is 1), it means that the putative bulge is a real bulge. Otherwise it is an internal loop.

                    if abs(paired_base_starting_index - paired_base_ending_index) == 1:  # abs() function is used to return the absolute value

                        bulge = putative_bulge

                        numb_bulges = numb_bulges + 1  # to count the number of bulges: if putative_bulge == bulge, increase the number of bulges

                        ### BULGES INDEXES

                        # bulge_starting_index = bulge.start() + 1  # to know the starting index of the first unpaired base (dot), we have to increase the number of 1
                        # bulge_ending_index = bulge.end() - 2  # to know the ending index of the last unpaired base (dot), we have to decrease the number of 2
                        # bulge_indexes = str(bulge_starting_index) + ":" + str(bulge_ending_index)

                    ### DEFINE INTERNAL LOOP

                    else:
                        internal_loop = putative_bulge

                        numb_internal_loop = numb_internal_loop + 1  # to count the number of internal loop: if putative_bulge == internal_loop, increase the number of internal loop

                        ### INTERNAL LOOP INDEXES

                        # internal_loop_starting_index = internal_loop.start() + 1
                        # internal_loop_ending_index = internal_loop.end() - 2
                        # internal_loop_indexes = str(internal_loop_starting_index) + ":" + str(internal_loop_ending_index)

                ### TOTAL NUMBER OF BULGES

                # numb_bulges

                ### TOTAL NUMBER OF INTERNAL LOOP

                # numb_internal_loop/2  # the unpaired bases of the internal loop are count once for every filament. So, we have to divide the number of internal loop by 2 in order to obtain the real number of internal loops present.

                # Working with 'energy_values_str'

                ### ENERGY VALUES OF THE UNTOUCHED SEQUENCE

                # energy_values_str

#========================
# PRINT ALL THE FEATURES:

            print("\n")
            print(header_str)
            print("DINUCLEOTIDE FREQUENCIES:")
            print("Frequency of AA = " + str(frequency_AA))  # str() function coverts a integer into a string (it is not possible to concatenate strings with integers)
            print("Frequency of AG = " + str(frequency_AG))
            print("Frequency of AC = " + str(frequency_AC))
            print("Frequency of AU = " + str(frequency_AU))
            print("Frequency of GA = " + str(frequency_GA))
            print("Frequency of GG = " + str(frequency_GG))
            print("Frequency of GC = " + str(frequency_GC))
            print("Frequency of GU = " + str(frequency_GU))
            print("Frequency of CA = " + str(frequency_CA))
            print("Frequency of CG = " + str(frequency_CG))
            print("Frequency of CC = " + str(frequency_CC))
            print("Frequency of CU = " + str(frequency_CU))
            print("Frequency of UA = " + str(frequency_UA))
            print("Frequency of UG = " + str(frequency_UG))
            print("Frequency of UC = " + str(frequency_UC))
            print("Frequency of UU = " + str(frequency_UU))
            print("TRINUCLEOTIDE FREQUENCIES:")
            print("Frequency of AAA = " + str(frequency_AAA))
            print("Frequency of AAG = " + str(frequency_AAG))
            print("Frequency of AAC = " + str(frequency_AAC))
            print("Frequency of AAU = " + str(frequency_AAU))
            print("Frequency of AGA = " + str(frequency_AGA))
            print("Frequency of AGG = " + str(frequency_AGG))
            print("Frequency of AGC = " + str(frequency_AGC))
            print("Frequency of AGU = " + str(frequency_AGU))
            print("Frequency of ACA = " + str(frequency_ACA))
            print("Frequency of ACG = " + str(frequency_ACG))
            print("Frequency of ACC = " + str(frequency_ACC))
            print("Frequency of ACU = " + str(frequency_ACU))
            print("Frequency of AUA = " + str(frequency_AUA))
            print("Frequency of AUG = " + str(frequency_AUG))
            print("Frequency of AUC = " + str(frequency_AUC))
            print("Frequency of AUU = " + str(frequency_AUU))
            print("Frequency of GAA = " + str(frequency_GAA))
            print("Frequency of GAG = " + str(frequency_GAG))
            print("Frequency of GAC = " + str(frequency_GAC))
            print("Frequency of GAU = " + str(frequency_GAU))
            print("Frequency of GGA = " + str(frequency_GGA))
            print("Frequency of GGG = " + str(frequency_GGG))
            print("Frequency of GGC = " + str(frequency_GGC))
            print("Frequency of GGU = " + str(frequency_GGU))
            print("Frequency of GCA = " + str(frequency_GCA))
            print("Frequency of GCG = " + str(frequency_GCG))
            print("Frequency of GCC = " + str(frequency_GCC))
            print("Frequency of GCU = " + str(frequency_GCU))
            print("Frequency of GUA = " + str(frequency_GUA))
            print("Frequency of GUG = " + str(frequency_GUG))
            print("Frequency of GUC = " + str(frequency_GUC))
            print("Frequency of GUU = " + str(frequency_GUU))
            print("Frequency of CAA = " + str(frequency_CAA))
            print("Frequency of CAG = " + str(frequency_CAG))
            print("Frequency of CAC = " + str(frequency_CAC))
            print("Frequency of CAU = " + str(frequency_CAU))
            print("Frequency of CGA = " + str(frequency_CGA))
            print("Frequency of CGG = " + str(frequency_CGG))
            print("Frequency of CGC = " + str(frequency_CGC))
            print("Frequency of CGU = " + str(frequency_CGU))
            print("Frequency of CCA = " + str(frequency_CCA))
            print("Frequency of CCG = " + str(frequency_CCG))
            print("Frequency of CCC = " + str(frequency_CCC))
            print("Frequency of CCU = " + str(frequency_CCU))
            print("Frequency of CUA = " + str(frequency_CUA))
            print("Frequency of CUG = " + str(frequency_CUG))
            print("Frequency of CUC = " + str(frequency_CUC))
            print("Frequency of CUU = " + str(frequency_CUU))
            print("Frequency of UAA = " + str(frequency_UAA))
            print("Frequency of UAG = " + str(frequency_UAG))
            print("Frequency of UAC = " + str(frequency_UAC))
            print("Frequency of UAU = " + str(frequency_UAU))
            print("Frequency of UGA = " + str(frequency_UGA))
            print("Frequency of UGG = " + str(frequency_UGG))
            print("Frequency of UGC = " + str(frequency_UGC))
            print("Frequency of UGU = " + str(frequency_UGU))
            print("Frequency of UCA = " + str(frequency_UCA))
            print("Frequency of UCG = " + str(frequency_UCG))
            print("Frequency of UCC = " + str(frequency_UCC))
            print("Frequency of UCU = " + str(frequency_UCU))
            print("Frequency of UUA = " + str(frequency_UUA))
            print("Frequency of UUG = " + str(frequency_UUG))
            print("Frequency of UUC = " + str(frequency_UUC))
            print("Frequency of UUU = " + str(frequency_UUU))
            print("FREQUENCIES OF TRIPLET ELEMENTS:")
            print("Frequency of '(((' and ')))' = " + str(frequency_bracket_bracket_bracket))
            print("Frequency of '((.' and ')).' = " + str(frequency_bracket_bracket_dot))
            print("Frequency of '(.(' and ').)' = " + str(frequency_bracket_dot_bracket))
            print("Frequency of '(..' and ')..' = " + str(frequency_bracket_dot_dot))
            print("Frequency of '.((' and '.))' = " + str(frequency_dot_bracket_bracket))
            print("Frequency of '.(.' and '.).' = " + str(frequency_dot_bracket_dot))
            print("Frequency of '..(' and '..)' = " + str(frequency_dot_dot_bracket))
            print("Frequency of '...' = " + str(frequency_dot_dot_dot))
            print("Frequency of '(()' = " + str(frequency_bracket_bracket_bracket_2))
            print("Frequency of '()(' = " + str(frequency_bracket_bracket_bracket_3))
            print("Frequency of '())' = " + str(frequency_bracket_bracket_bracket_4))
            print("Frequency of ')((' = " + str(frequency_bracket_bracket_bracket_5))
            print("Frequency of ')()' = " + str(frequency_bracket_bracket_bracket_6))
            print("Frequency of '))(' = " + str(frequency_bracket_bracket_bracket_7))
            print("Frequency of '().' = " + str(frequency_bracket_bracket_dot_2))
            print("Frequency of ')(.' = " + str(frequency_bracket_bracket_dot_3))
            print("Frequency of '(.)' = " + str(frequency_bracket_dot_bracket_2))
            print("Frequency of ').(' = " + str(frequency_bracket_dot_bracket_3))
            print("Frequency of '.()' = " + str(frequency_dot_bracket_bracket_2))
            print("Frequency of '.)(' = " + str(frequency_dot_bracket_bracket_3))
            print("Total number of not paired bases = " + str(numb_not_paired_bases))
            print("Total number of loops in the untouched dot-bracket notation = " + str(numb_loops))
            print("Total number of loops, defining them as ((.{3,})) = " + str(numb_loops_3_dots_2_parentesis))
            print("Total number of loops, defining them as (.{4,}) = " + str(numb_loops_4_dots_1_parentesis))
            print("Total number of loops, defining them as ((.{4,})) = " + str(numb_loops_4_dots_2_parentesis))
            print("Total number of loops, defining them as (.{5,}) = " + str(numb_loops_5_dots_1_parentesis))
            print("Total number of loops, defining them as ((.{5,})) = " + str(numb_loops_5_dots_2_parentesis))
            print("Indexes of the loop of interest = " + str(starting_index_loop_of_interest) + ":" + str(ending_index_loop_of_interest))
            print("Centrality of the loop of interest = " + str(loop_centrality))
            print("Dimension of the loop of interest = " + str(length_loop_of_interest))
            print("Total number of paired bases = " + str(numb_paired_bases))
            print("FREQUENCY OF PAIRED COUPLES:")
            print("Frequency of CG  = " + str(frequency_couple_CG))
            print("Frequency of GC  = " + str(frequency_couple_GC))
            print("Frequency of UA  = " + str(frequency_couple_UA))
            print("Frequency of AU  = " + str(frequency_couple_AU))
            print("Frequency of UG  = " + str(frequency_couple_UG))
            print("Frequency of GU  = " + str(frequency_couple_GU))
            print("Stem sequence at 5' = " + reversed_stem_sequence_5_prime)
            print("Stem sequence at 3' = " + reversed_stem_sequence_3_prime)
            # if numb_bulges != 0:
            #     print("Bulge indexes = " + bulge_indexes)
            # if numb_internal_loop/2 != 0:
            #     print("Internal loop indexes = " + internal_loop_indexes)
            print("Total number of bulges = " + str(numb_bulges))
            print("Total number of internal loops = " + str(numb_internal_loop / 2))
            print("Energy values of the untouched sequence = " + energy_values_str)

# Save all the features printed in a file


# To know which sequences have not been considered:
#print("\n")
#print("Warnings:")
print(warning_list)

