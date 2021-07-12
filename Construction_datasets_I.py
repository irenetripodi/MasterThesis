
# 24/02/21
# 5/03/21

import random  # module used to make random items

import re


#######################################
# CONSTRUCTION OF REAL miRNAS SEQUENCES
#######################################
fasta_mirbase = open("hairpin.fa", "r")
fasta_mirbase = fasta_mirbase.readlines()


hsa_fasta_mirbase = []
for line in fasta_mirbase:
    if line.startswith(">hsa"):
        hsa_fasta_mirbase.append(line)
#print(hsa_fasta_mirbase)


#########################################
# CONSTRUCTION OF PSEUDO miRNAS SEQUENCES
#########################################

# We generate two types of pseudo-miRNA:
#######################################################################################
### 1- STARTING FROM THE SEQUENCE OF ALL THE MIRNA NOT MUTATED AND SHUFFLE THE SEQUENCE
#######################################################################################
regex_header = "(>.*)"  # for header
regex_sequence = "([AUCG]+)"  # for sequence

# First of all selecting the headers of the sequences used and present in 'features_all_human_hairpin_sequences_not_mutated_TRUNCATED_SEQUENCES'
file_features = open("features_all_human_hairpin_sequences_not_mutated_TRUNCATED_SEQUENCES")
file_features = file_features.readlines()

sequences_used = []

for line in file_features:
    header_features = re.findall(regex_header, line)
    if (header_features):
        header_used = line.split(" ")[0]
        sequences_used.append(header_used)

number_sequences_used = len(sequences_used)

# Secondly, search in 'all_human_hairpin_sequences_not_mutated' the headers of the sequences used, present in 'sequences_used'

miRNA_sequences_not_mutated_file = open("all_human_hairpin_sequences_not_mutated_TRUNCATED_SEQUENCES", "r")
miRNA_sequences_not_mutated = miRNA_sequences_not_mutated_file.readlines()

for row in miRNA_sequences_not_mutated:
    header = re.findall(regex_header, row)
    sequence = re.findall(regex_sequence, row)

    if (header):
        header_str = header[0]
        header_name = header_str.split(" ")[0]

    if (sequence):
        sequence_str = sequence[0]
        sequence_list = list(sequence_str)  # firstly we have to convert the string to a list of characters, shuffle it, then join the result again.
        random.shuffle(sequence_list)  # shuffle() function takes a sequence and returns the sequence in a random order
        sequence_shuffled = ''.join(sequence_list)

        for elem in sequences_used:
            if header_name == elem:
                print(header_str)
                print(sequence_shuffled)
                break

        print(header_str)
        print(sequence_shuffled)

# Save the file as "pseudo_mirna_shuffled_sequences"

##############################################################################################################
### 2- STARTING FROM THE TRANSCRIPT SEQUENCE AND, LOOKING AT THE ANNOTATION, EXTRACT THE CODING EXON SEQUENCES
##############################################################################################################

regex_header = "(>.*)"  # for header
regex_sequence = "([AUCG]+)"  # for sequence
regex_sequence_ATCG = "([ATCG]+)"  # for sequence

##### FIRSTLY WE CREATE A SMALLER FILE OF THE TRANSCRIPTS DISCARDING ALL THE NON CODING RNAs:

# transcript_file = open("GRCh38_latest_rna.fna copia", "r")
# transcript = transcript_file.read()
#
# # recognition of sequences from headers
# regex_trascript = "(>.*\n)([AGTC\n]*)"  #  splitting the trascript file in groups of headers and sequences
# transcript_split = re.findall(regex_trascript, transcript)
# #print(transcript_split)  # creating a structure characterized by: row 0 (column 0 --> first header; column 1 --> first sequence); row 1 (column 0 --> second header; column 1 --> second sequence);
#
# #reorder transcript_split
# mRNA_list = []
# mRNA_str = ""
# #microRNA_list = []
# #non_coding_RNA_list = []
# #long_non_coding_RNA_list = []
#
# iterator = 0
#
# for line in transcript_split:
#       header_ = str(line[0]).replace("\n", "")  # deleting \n in the headers
#       header_split = header_.split(" ")
#       sequence_ = str(line[1]).replace("\n", "")  # deleting \n in the sequences
#       if "mRNA" in header_split:
#           iterator = iterator + 1
#           item = []
#           item.append(header_)
#           item.append(sequence_)
#           mRNA_list.append(item)
#           mRNA_str = mRNA_str + header_ + "\n" + sequence_ + "\n"
#
#       # if "microRNA" in header_split:
#       #     item = []
#       #     item.append(header_)
#       #     item.append(sequence_)
#       #     microRNA_list.append(item)
#       #
#       # if "non-coding RNA" in header_split:
#       #     item = []
#       #     item.append(header_)
#       #     item.append(sequence_)
#       #     non_coding_RNA_list.append(item)
#       #
#       # if "long non-coding RNA" in header_split:
#       #     item = []
#       #     item.append(header_)
#       #     item.append(sequence_)
#       #     long_non_coding_RNA_list.append(item)
#
#       if iterator == 1561:
#           print(mRNA_str)
#           break
#
# # Save the file with all the transcripts coding for proteins --> the file is called transcipt_coding_for_proteins

##### SECONDLY, WE FIND IN THE ANNOTATION FILE ALL THE CODING EXONS

annotation_file = open("GRCh38_latest_genomic.gff copia", "r")  # annotation file
annotation = annotation_file.readlines()

total_annotation_info = []
mRNA_found = False

for every_row in annotation:
    if not every_row.startswith("#"):
        every_row = every_row.split("\t")
        type = every_row[2]

        if type == "mRNA":
            mRNA_found = True
            annotation_info = []
            list = []
            list.append(type)
            start_position = every_row[3]
            list.append(start_position)
            end_position = every_row[4]
            list.append(end_position)
            info_column_split = every_row[8].split(";")

            for element in info_column_split:
                if element.startswith("ID"):
                    ID = element
                    list.append(ID)

                if element.startswith("Parent"):
                    parent = element
                    list.append(parent)

                if element.startswith("transcript_id"):
                    transcript_ID = element
                    list.append(transcript_ID)
            annotation_info.append(list)


        if type == "exon" and mRNA_found == True:
            list = []
            list.append(type)
            start_position = every_row[3]
            list.append(start_position)
            end_position = every_row[4]
            list.append(end_position)
            info_column_split = every_row[8].split(";")

            for element in info_column_split:
                if element.startswith("ID"):
                    ID = element
                    list.append(ID)

                if element.startswith("Parent"):
                    parent = element
                    list.append(parent)

                if element.startswith("transcript_id"):
                    transcript_ID = element
                    list.append(transcript_ID)
            annotation_info.append(list)

        if type == "CDS" and mRNA_found == True:
            list = []
            list.append(type)
            start_position = every_row[3]
            list.append(start_position)
            end_position = every_row[4]
            list.append(end_position)
            info_column_split = every_row[8].split(";")

            for element in info_column_split:
                if element.startswith("ID"):
                    ID = element
                    list.append(ID)

                if element.startswith("Parent"):
                    parent = element
                    list.append(parent)

                if element.startswith("transcript_id"):
                    transcript_ID = element
                    list.append(transcript_ID)
            annotation_info.append(list)
            mRNA_found = False
            total_annotation_info.append(annotation_info)

#print(total_annotation_info)


##### DEFINING EXONS INDEXES AND THEIR CDS

new_annotation_list = []
transcript_id_list = []

for sublist in total_annotation_info:
    transcript_id_list = []
    transcript_id_list.append(sublist[0][5].split("=")[1])
    exons_list = []
    cds_list = []
    for sub_list in sublist:

        if sub_list[0] == "exon":
            exon_item = []

            coord_start = sub_list[1]
            coord_end = sub_list[2]
            id_sublist = sub_list[3].split("=")[1]
            exon_item.append(id_sublist)
            exon_item.append(coord_start)
            exon_item.append(coord_end)
            exons_list.append(exon_item)

        if sub_list[0] == "CDS":
            cds_item = []

            coord_start = sub_list[1]
            coord_end = sub_list[2]
            id_sublist = sub_list[3].split("=")[1]
            cds_item.append(id_sublist)
            cds_item.append(coord_start)
            cds_item.append(coord_end)
            cds_list.append(cds_item)

    #transcript_id_list.append(transcript_name)
    transcript_id_list.append(exons_list)
    transcript_id_list.append(cds_list)
    new_annotation_list.append(transcript_id_list)

#print(new_annotation_list)


##### CALCULATE THE RELATIVE POSITIONS OF THE CODING EXONS

cds_exons_to_extract = []

for item in new_annotation_list:
    transcript_id = item[0]
    cont = 0
    for exon in item[1]:
        coord_exon_start = int(exon[1])
        coord_exon_end = int(exon[2])
        size = coord_exon_end - coord_exon_start + 1
        index_exon_start = cont
        index_exon_end = index_exon_start + size - 1
        cont = index_exon_end + 1
        exon.append(index_exon_start)
        exon.append(index_exon_end)

    for cds in item[2]:
        cds_indexes_list = []
        coord_cds_start = int(cds[1])
        coord_cds_end = int(cds[2])

        for exon in item[1]:
            cds_indexes = []
            coord_exon_start = int(exon[1])
            coord_exon_end = int(exon[2])

            if coord_cds_start > coord_exon_end: # my cds does not start in mu exon
                continue

            elif coord_cds_start >= coord_exon_start and coord_cds_end <= coord_exon_end:  # my cds is enclosed in my exon
                index_cds_start = index_exon_start + (coord_cds_start - coord_exon_start)
                index_cds_end = index_exon_end - (coord_exon_end - coord_cds_end)
                cds_indexes.append(transcript_id)  # append the transcript id
                cds_indexes.append(exon[0])  # append the name of the exon whose cds coordinates are expressed below
                cds_indexes.append(index_cds_start)
                cds_indexes.append(index_cds_end)
                cds_exons_to_extract.append(cds_indexes)

            else:
                if coord_cds_start >= coord_exon_start and coord_cds_start <= coord_exon_end and coord_cds_end > coord_exon_end: # my cds starts in my exon but does not end in my exon
                    index_cds_start = index_exon_start + (coord_cds_start - coord_exon_start)
                    index_cds_end = index_exon_end

                if coord_cds_start < coord_exon_start and coord_cds_end > coord_exon_end:  # my cds starts before my exon and ends after my exon (my exon is totally enclosed in my cds)
                    index_cds_start = index_exon_start
                    index_cds_end = index_exon_end

                if coord_cds_start < coord_exon_start and coord_cds_end >= coord_exon_start and coord_cds_end <= coord_exon_end:  # my cds does not start in my exon but ends in my exon
                    index_cds_start = index_exon_start
                    index_cds_end = index_cds_start + (coord_cds_end - coord_exon_start)

                cds_indexes.append(transcript_id)
                cds_indexes.append(exon[0])  # append the name of the exon whose cds coordinates are expressed below
                cds_indexes.append(index_cds_start)
                cds_indexes.append(index_cds_end)
                cds_exons_to_extract.append(cds_indexes)

#print(cds_exons_to_extract)


##### EXTRACT THE CODING EXONS IN THE TRANSCRIPTS

filtered_transcript_file = open("transcipt_coding_for_proteins", "r")  # transcript file
filtered_transcript = filtered_transcript_file.read()
# recognition of sequences from headers
regex_filtered_trascript = "(>.*\n)([AGTC\n]*)"  #  splitting the filtered trascript file in groups of headers and sequences
filtered_transcript_split = re.findall(regex_filtered_trascript, filtered_transcript)

for row in filtered_transcript_split:
     header_transcript = str(row[0]).replace("\n", "")  # deleting \n in the headers
     sequence_transcript = str(row[1]).replace("\n", "")  # deleting \n in the sequences

     for component in cds_exons_to_extract:
         transcript_id = component[0].replace("\n", "")
         exon_name = component[1]
         index_coding_exon_start = component[2]
         index_coding_exon_end = component[3]

         if header_transcript[1:].startswith(transcript_id):
            coding_exon_header = ">" + " " + exon_name
            coding_exon_sequence = sequence_transcript[index_coding_exon_start: index_coding_exon_end+1]
            #print(coding_exon_header)
            #print(coding_exon_sequence)

# Save the sequences in a file --> coding_exons_sequences


##### CREATE A LIST THAT WILL CONTAIN ALL THE HEADERS OF THE MIRNAs USED AND THEIR LENGTHS

# First of all selecting the headers of the sequences used and present in 'features_all_human_hairpin_sequences_not_mutated_TRUNCATED_SEQUENCES'
file_features = open("features_all_human_hairpin_sequences_not_mutated_3_dots_1_parentesis", "r")
file_features = file_features.readlines()

sequences_used = []

for line in file_features:
    header_features = re.findall(regex_header, line)
    sequence_features = re.findall(regex_sequence, line)
    if (header_features):
        header_used = line.split(" ")[0]
        sequences_used.append(header_used)

number_sequences_used = len(sequences_used)

# Then, selecting also the length of the sequences used present in the list 'sequences_used'

list_miRNAs_used = []

table_rnafold_file = open("table_rnafold_output_all_human_hairpin_sequences_not_mutated", "r")
table_rnafold_file = table_rnafold_file.readlines()

for line_rnafold in table_rnafold_file:
    header_rnafold = re.findall(regex_header, line_rnafold)
    sequence_rnafold = re.findall(regex_sequence, line_rnafold)

    if (header_rnafold):
        header_rnafold_used = header_rnafold[0].split(" ")[0]

    if (sequence_rnafold):
        len_seq = len(sequence_rnafold[0])

        for elem in sequences_used:
            if elem == header_rnafold_used:
                list_miRNAs_used.append(elem)
                list_miRNAs_used.append(len_seq)
                break

#print(list_miRNAs_used)  # it is the lisf of reference that we will use: it gives information of the headers used and the lenght of the miRNAs used

#####  EXRACT PART OF THE SEQUENCE FROM CODING EXONS SEQUENCES
# We have to select from the file generated 'coding_exons_sequences' the same number of sequences that we use for real miRNAs (1917).
# Every coding exon sequence (pseudo-miRNA) should be coupled randomly with one sequence of real miRNA. The so-called 'list_miRNAs_used' contains the headers and in the 2 positin the lenght of the sequence.
# We have to extract the sequence of coding exons with the same length of the real miRNA coupled with the pseudo-miRNA. The starting and ending position for the extraction are random, we have only to be sure that the lenght is the same .

pseudo_miRNAs_coding_exons_file = open("coding_exons_sequences", "r")
pseudo_miRNAs_coding_exons = pseudo_miRNAs_coding_exons_file.readlines()

iter = 0
pseudo_mirna_already_extracted =  []

for header_used in list_miRNAs_used:
    header_mirna_used = list_miRNAs_used[iter]
    lenght_mirna_used = list_miRNAs_used[iter + 1]

    for line_pseudo_mirna in pseudo_miRNAs_coding_exons:
        header_pseudo_mirna = re.findall(regex_header, line_pseudo_mirna)
        sequence_pseudo_mirna = re.findall(regex_sequence_ATCG, line_pseudo_mirna)

        if (header_pseudo_mirna):
            header_pseudo_mirna_str = header_pseudo_mirna[0]

        if (sequence_pseudo_mirna):
            sequence_pseudo_mirna_str = sequence_pseudo_mirna[0]
            len_sequence_pseudo_mirna = len(sequence_pseudo_mirna_str)

            if len_sequence_pseudo_mirna > lenght_mirna_used and header_pseudo_mirna_str not in pseudo_mirna_already_extracted:
                # generate randomly an initial and final position that are between the real positions of the real miRNA
                NEW_INITIAL_POS = random.sample(range(0, (len_sequence_pseudo_mirna-lenght_mirna_used)), 1)  # random() function output is a list
                NEW_INITIAL_POS = int(map(int, NEW_INITIAL_POS)[0])  # Convert list in string and the string into integer
                NEW_FINAL_POS = NEW_INITIAL_POS + lenght_mirna_used
                check_NEW_LENGHT = NEW_FINAL_POS - NEW_INITIAL_POS  # the real mirna lenght and the pseudo mirna lenght must be the same

                if check_NEW_LENGHT == lenght_mirna_used:
                    extracted_pseudo_miRNA_sequence = sequence_pseudo_mirna_str[NEW_INITIAL_POS:NEW_FINAL_POS]
                    new_header = header_pseudo_mirna_str + ", " + str(NEW_INITIAL_POS) + ":" + str(NEW_FINAL_POS) + ", " + "pseudo-miRNA " + header_mirna_used[1:]
                    pseudo_mirna_already_extracted.append(header_pseudo_mirna_str)  # the coding exon header (whose sequence is already used to extract the pseudo mirna seq) cannot be used again
                    break


    iter = iter + 2
    print(new_header)
    print(extracted_pseudo_miRNA_sequence)


    # if iter > 3834:
    #     break




# Save the sequences in a file --> pseudo_mirna_coding_exons_sequences





