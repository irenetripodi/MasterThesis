# 24/05/21

import re

#### Adding the columns in TCGA_PanCanAtlas_Pre_Germ file, in order to know the classification of the sequences, before and after the mutation added

file_TCGA_Pre_Germ = open("TCGA_PanCanAtlas_Pre_Germ.txt", "r")  
file_TCGA_Pre_Germ = file_TCGA_Pre_Germ.readlines()

file_classification_BEFORE_mutation_RF_model = open("classification_before_germ_mutation_model_RF_tripletsvm.txt", "r")
classification_BEFORE_mutation_RF_model = file_classification_BEFORE_mutation_RF_model.readlines()

file_classification_AFTER_mutation_RF_model = open("classification_after_germ_mutation_model_RF_tripletsvm.txt", "r")
classification_AFTER_mutation_RF_model = file_classification_AFTER_mutation_RF_model.readlines()

regex_row_before = "[0-9]+[ ]+([A-Z0-9-]+)[ ]+([A-Za-z \n-]+)"
regex_row_after = "[0-9]+[ ]+([A-Z0-9->]+)[ ]+GERMLINE,[ ]+([A-Z0-9_-]+)[ ]+([A-Za-z \n-]+)"

outputfile = open("TCGA_PanCanAtlas_Pre_Germ_CLASSIFICATION_MODEL_RF_tripletsvm.txt", 'w')

for row in file_TCGA_Pre_Germ[1:]:
    row = row.replace("\n", "").replace("\r", "")
    row_split = row.split("\t")
    mirna_name_snp = row_split[2].upper()
    mutation = row_split[17]

    for line_before in classification_BEFORE_mutation_RF_model[1:]:  #  BEFORE THE MUTATION
        line_split_before = re.findall(regex_row_before, line_before)
        name_mirna_before = line_split_before[0][0]
        classification_before = line_split_before[0][1].replace("\n", "")

        if mirna_name_snp == name_mirna_before:
            row = row + "\t" + classification_before


            for line_after in classification_AFTER_mutation_RF_model[1:]:  # AFTER THE MUTATION
                line_split_after = re.findall(regex_row_after, line_after)
                name_mirna_after = line_split_after[0][0][1:]
                mutation_after = line_split_after[0][1]
                classification_after = line_split_after[0][2].replace("\n", "")

                if mirna_name_snp == name_mirna_after and mutation == mutation_after:
                    row = row + "\t" + classification_after
                    #print(row)
                    outputfile.write(row + "\n")
                    break
            break


outputfile.close()