# 31/05/21

from tabulate import tabulate

import scipy.stats as stats

##############################
# COMPARISON OF MIRNAS MUTATED
##############################

### Loading the files: Classified as miRNA before the mutation and Classified as pseudo-miRNA after the mutation
# 1
PRE_SOM_mirna_before_and_pseudo_after = open("PRE_SOM_miRNA_before_and_pseudo_after", "r")
PRE_SOM_mirna_before_and_pseudo_after = PRE_SOM_mirna_before_and_pseudo_after.readlines()

# 2
PRE_GERM_mirna_before_and_pseudo_after = open("PRE_GERM_miRNA_before_and_pseudo_after", "r")
PRE_GERM_mirna_before_and_pseudo_after = PRE_GERM_mirna_before_and_pseudo_after.readlines()

# 3
GNOMAD_WES_mirna_before_and_pseudo_after = open("GNOMADWES_miRNA_before_and_pseudo_after", "r")
GNOMAD_WES_mirna_before_and_pseudo_after = GNOMAD_WES_mirna_before_and_pseudo_after.readlines()


### Loading the files with all the mutations, not divided according to the diruptive or not mutations
# 1
TCGA_PRE_SOM = open("TCGA_PanCanAtlas_Pre_Som_CLASSIFICATION_MODEL_RF_tripletsvm.txt", "r")
TCGA_PRE_SOM = TCGA_PRE_SOM.readlines()

# 2
TCGA_PRE_GERM = open("TCGA_PanCanAtlas_Pre_Germ_CLASSIFICATION_MODEL_RF_tripletsvm.txt", "r")
TCGA_PRE_GERM = TCGA_PRE_GERM.readlines()

# 3
GNOMAD_WES = open("gnomAD_WES_CLASSIFICATION_MODEL_RF_tripletsvm.txt", "r")
GNOMAD_WES = GNOMAD_WES.readlines()

#=======================================================================================================================

###################
# SOMATIC MUTATIONS
###################
damaging_somatic_mutations = [] # mirna before the mutation, pseudo after the mutation
all_somatic_mutations = []

for line in PRE_SOM_mirna_before_and_pseudo_after:
    line = line.split("\t")
    mirna = line[2]
    damaging_somatic_mutations.append(mirna)
dictionary_damaging_somatic_mutations = {i:damaging_somatic_mutations.count(i) for i in damaging_somatic_mutations}
#print("Damaging somatic mutations: " + str(dictionary_damaging_somatic_mutations))

for line in TCGA_PRE_SOM:
    line = line.split("\t")
    mirna = line[2]
    all_somatic_mutations.append(mirna)
dictionary_all_somatic_mutations = {i:all_somatic_mutations.count(i) for i in all_somatic_mutations}
#print("All somatic mutations: " + str(dictionary_all_somatic_mutations))

####################
# GERMLINE MUTATIONS
####################
damaging_germline_mutations = [] # mirna before the mutation, pseudo after the mutation
all_germline_mutations = []

for line in PRE_GERM_mirna_before_and_pseudo_after:
    line = line.split("\t")
    mirna = line[2]
    damaging_germline_mutations.append(mirna)
dictionary_damaging_germline_mutations = {i:damaging_germline_mutations.count(i) for i in damaging_germline_mutations}
#print("Damaging germline mutations: " + str(dictionary_damaging_germline_mutations))


for line in TCGA_PRE_GERM:
    line = line.split("\t")
    mirna = line[2]
    all_germline_mutations.append(mirna)
dictionary_all_germline_mutations = {i:all_germline_mutations.count(i) for i in all_germline_mutations}
#print("All germline mutations: " + str(dictionary_all_germline_mutations))

######################
# GNOMAD WES MUTATIONS
######################
damaging_gnomad_wes_mutations = [] # mirna before the mutation, pseudo after the mutation
all_gnomad_wes_mutations = []

for line in GNOMAD_WES_mirna_before_and_pseudo_after:
    line = line.split("\t")
    mirna = line[19]
    damaging_gnomad_wes_mutations.append(mirna)
dictionary_damaging_gnomad_wes_mutations = {i:damaging_gnomad_wes_mutations.count(i) for i in damaging_gnomad_wes_mutations}
#print("Damaging gnomad-wes mutations: " + str(dictionary_damaging_gnomad_wes_mutations))


for line in GNOMAD_WES:
    line = line.split("\t")
    mirna = line[19]
    all_gnomad_wes_mutations.append(mirna)
dictionary_all_gnomad_wes_mutations = {i:all_gnomad_wes_mutations.count(i) for i in all_gnomad_wes_mutations}
#print("All gnomad-wes mutations: " + str(dictionary_all_gnomad_wes_mutations))

#print(len(dictionary_all_gnomad_wes_mutations))
#=======================================================================================================================
table = []
table_header = ['MIRNA NAME', 'TYPE', 'TOT MUTATIONS', 'DAMAGING MUTATIONS']
table.append(table_header)

for X_all_gnomad_wes, Y_all_gnomad_wes in dictionary_all_gnomad_wes_mutations.items():
    mirna_name_all_gnomad_wes = X_all_gnomad_wes
    numb_all_gnomad_wes = Y_all_gnomad_wes
    if dictionary_damaging_gnomad_wes_mutations.has_key(mirna_name_all_gnomad_wes):
        numb_damaging_gnomad_wes = dictionary_damaging_gnomad_wes_mutations[mirna_name_all_gnomad_wes]
    else:
        numb_damaging_gnomad_wes = 0
    row_table_gnomad_wes = [mirna_name_all_gnomad_wes, 'GnomAD-WES', numb_all_gnomad_wes-numb_damaging_gnomad_wes, numb_damaging_gnomad_wes]
    table.append(row_table_gnomad_wes)

    #
    if dictionary_all_germline_mutations.has_key(mirna_name_all_gnomad_wes):
        numb_all_germline = dictionary_all_germline_mutations[mirna_name_all_gnomad_wes]
    else:
        numb_all_germline = 0
    if dictionary_damaging_germline_mutations.has_key(mirna_name_all_gnomad_wes):
        numb_damaging_germline = dictionary_damaging_germline_mutations[mirna_name_all_gnomad_wes]
    else:
        numb_damaging_germline = 0
    row_table_germline = [mirna_name_all_gnomad_wes, 'Germ', numb_all_germline-numb_damaging_germline, numb_damaging_germline]
    table.append(row_table_germline)
    oddsratio, pvalue = stats.fisher_exact([[numb_all_gnomad_wes-numb_damaging_gnomad_wes, numb_damaging_gnomad_wes], [numb_all_germline-numb_damaging_germline, numb_damaging_germline]])
    print(mirna_name_all_gnomad_wes + " " + str(pvalue) + "\n")

print(tabulate(table))

#=======================================================================================================================
table = []
table_header = ['MIRNA NAME', 'TYPE', 'TOT MUTATIONS', 'DAMAGING MUTATIONS']
table.append(table_header)

for X_all_gnomad_wes, Y_all_gnomad_wes in dictionary_all_gnomad_wes_mutations.items():
    mirna_name_all_gnomad_wes = X_all_gnomad_wes
    numb_all_gnomad_wes = Y_all_gnomad_wes
    if dictionary_damaging_gnomad_wes_mutations.has_key(mirna_name_all_gnomad_wes):
        numb_damaging_gnomad_wes = dictionary_damaging_gnomad_wes_mutations[mirna_name_all_gnomad_wes]
    else:
        numb_damaging_gnomad_wes = 0
    row_table_gnomad_wes = [mirna_name_all_gnomad_wes, 'GnomAD-WES', numb_all_gnomad_wes-numb_damaging_gnomad_wes, numb_damaging_gnomad_wes]
    table.append(row_table_gnomad_wes)

    #

    if dictionary_all_somatic_mutations.has_key(mirna_name_all_gnomad_wes):
        numb_all_somatic = dictionary_all_somatic_mutations[mirna_name_all_gnomad_wes]
    else:
        numb_all_somatic = 0
    if dictionary_damaging_somatic_mutations.has_key(mirna_name_all_gnomad_wes):
        numb_damaging_somatic = dictionary_damaging_somatic_mutations[mirna_name_all_gnomad_wes]
    else:
        numb_damaging_somatic = 0
    row_table_somatic = [mirna_name_all_gnomad_wes, 'Som', numb_all_somatic-numb_damaging_somatic, numb_damaging_somatic]
    table.append(row_table_somatic)
    oddsratio, pvalue = stats.fisher_exact([[numb_all_gnomad_wes - numb_damaging_gnomad_wes, numb_damaging_gnomad_wes],[numb_all_somatic-numb_damaging_somatic, numb_damaging_somatic]])
    print(mirna_name_all_gnomad_wes + " " + str(pvalue) + "\n")

print(tabulate(table))

