# ML
# 23/03/21



feature_table_mirna_tris = read.table("table_features_all_human_hairpin_sequences_not_mutated_tris.txt", sep = "\t", header = T, check.names = F)
#typeof(feature_table_mirna_tris)
feature_table_pseudo_mirna_tris = read.table("table_features_pseudo_miRNA_shuffled_tris.txt", sep = "\t", header = T, check.names = F)
# N.B.: the variables have the name '_tris' because the table_features are obtained running the RNAfold structure with the '_tris' version of the features script 

feature_table_tris = rbind(feature_table_mirna_tris, feature_table_pseudo_mirna_tris)


###############
# NORMALIZATION
###############
# Normalization is needed to make data more comparable

column_numb_not_paired_bases = feature_table_tris[, 102:102]  
#typeof(column_numb_not_paired_bases)  # it is a vector of integers
#View(column_numb_not_paired_bases)
mean_numb_not_paired_bases = mean(column_numb_not_paired_bases)
new_column_numb_not_paired_bases = c()
for (row in column_numb_not_paired_bases) {
  averaged_row = row/mean_numb_not_paired_bases
  new_column_numb_not_paired_bases = append(new_column_numb_not_paired_bases, averaged_row)
}

column_numb_loop = feature_table_tris[, 103:103]
#View(column_numb_loop)
mean_numb_loop = mean(column_numb_loop)
new_column_numb_loop = c()
for (row in column_numb_loop) {
  averaged_row = row/mean_numb_loop
  new_column_numb_loop = append(new_column_numb_loop, averaged_row)
}

column_numb_loop_3_dots_2_parentesis = feature_table_tris[, 104:104]
#View(column_numb_loop_3_dots_2_parentesis)
mean_column_numb_loop_3_dots_2_parentesis = mean(column_numb_loop_3_dots_2_parentesis)
new_column_numb_loop_3_dots_2_parentesis = c()
for (row in column_numb_loop_3_dots_2_parentesis) {
  averaged_row = row/mean_column_numb_loop_3_dots_2_parentesis
  new_column_numb_loop_3_dots_2_parentesis = append(new_column_numb_loop_3_dots_2_parentesis, averaged_row)
}

column_numb_loop_4_dots_1_parentesis = feature_table_tris[, 105:105]
#View(column_numb_loop_4_dots_1_parentesis)
mean_column_numb_loop_4_dots_1_parentesis = mean(column_numb_loop_4_dots_1_parentesis)
new_column_numb_loop_4_dots_1_parentesis = c()
for (row in column_numb_loop_4_dots_1_parentesis) {
  averaged_row = row/mean_column_numb_loop_4_dots_1_parentesis
  new_column_numb_loop_4_dots_1_parentesis = append(new_column_numb_loop_4_dots_1_parentesis, averaged_row)
}

column_numb_loop_4_dots_2_parentesis = feature_table_tris[, 106:106]
#View(column_numb_loop_4_dots_2_parentesis)
mean_column_numb_loop_4_dots_2_parentesis = mean(column_numb_loop_4_dots_2_parentesis)
new_column_numb_loop_4_dots_2_parentesis = c()
for (row in column_numb_loop_4_dots_2_parentesis) {
  averaged_row = row/mean_column_numb_loop_4_dots_2_parentesis
  new_column_numb_loop_4_dots_2_parentesis = append(new_column_numb_loop_4_dots_2_parentesis, averaged_row)
}

column_numb_loop_5_dots_1_parentesis = feature_table_tris[, 107:107]
#View(column_numb_loop_5_dots_1_parentesis)
mean_column_numb_loop_5_dots_1_parentesis = mean(column_numb_loop_5_dots_1_parentesis)
new_column_numb_loop_5_dots_1_parentesis = c()
for (row in column_numb_loop_5_dots_1_parentesis) {
  averaged_row = row/mean_column_numb_loop_5_dots_1_parentesis
  new_column_numb_loop_5_dots_1_parentesis = append(new_column_numb_loop_5_dots_1_parentesis, averaged_row)
}

column_numb_loop_5_dots_2_parentesis = feature_table_tris[, 108:108]
#View(column_numb_loop_5_dots_2_parentesis)
mean_column_numb_loop_5_dots_2_parentesis = mean(column_numb_loop_5_dots_2_parentesis)
new_column_numb_loop_5_dots_2_parentesis = c()
for (row in column_numb_loop_5_dots_2_parentesis) {
  averaged_row = row/mean_column_numb_loop_5_dots_2_parentesis
  new_column_numb_loop_5_dots_2_parentesis = append(new_column_numb_loop_4_dots_2_parentesis, averaged_row)
}

column_centrality = feature_table_tris[, 109:109]
#View(column_centrality)
mean_centrality = mean(column_centrality)
new_column_centrality = c()
for (row in column_centrality) {
  averaged_row = row/mean_centrality
  new_column_centrality = append(new_column_centrality, averaged_row)
}                                  

column_dim_loop = feature_table_tris[, 110:110]
#View(column_dim_loop)
mean_dim_loop = mean(column_dim_loop)
new_column_dim_loop = c()
for (row in column_dim_loop) {
  averaged_row = row/mean_dim_loop
  new_column_dim_loop = append(new_column_dim_loop, averaged_row)
}

column_numb_paired_bases = feature_table_tris[, 111:111]
#View(column_numb_paired_bases)
mean_numb_paired_bases = mean(column_numb_paired_bases)
new_column_numb_paired_bases = c()
for (row in column_numb_paired_bases) {
  averaged_row = row/mean_numb_paired_bases
  new_column_numb_paired_bases = append(new_column_numb_paired_bases, averaged_row)
}

column_numb_bulges = feature_table_tris[, 118:118]
#View(column_numb_bulges)
mean_numb_bulges = mean(column_numb_bulges)
new_column_numb_bulges = c()
for (row in column_numb_bulges) {
  averaged_row = row/mean_numb_bulges 
  new_column_numb_bulges = append(new_column_numb_bulges, averaged_row)
}

column_numb_internal_loop = feature_table_tris[, 119:119]
#View(column_numb_internal_loop)
mean_numb_internal_loop = mean(column_numb_internal_loop)
new_column_numb_internal_loop = c()
for (row in column_numb_internal_loop) {
  averaged_row = row/mean_numb_internal_loop
  new_column_numb_internal_loop = append(new_column_numb_internal_loop, averaged_row)
}

column_minimum_free_energy = feature_table_tris[, 120:120]
#View(column_minimum_free_energy)
mean_minimum_free_energy = mean(column_minimum_free_energy)
new_column_minimum_free_energy = c()
for (row in column_minimum_free_energy) {
  averaged_row = row/mean_minimum_free_energy
  new_column_minimum_free_energy = append(new_column_minimum_free_energy, averaged_row)
}

column_ensemble_free_energy = feature_table_tris[, 121:121]
#View(column_ensemble_free_energy)
mean_ensemble_free_energy = mean(column_ensemble_free_energy)
new_column_ensemble_free_energy = c()
for (row in column_ensemble_free_energy) {
  averaged_row = row/mean_ensemble_free_energy
  new_column_ensemble_free_energy = append(new_column_ensemble_free_energy, averaged_row)
}

column_centroid_structure = feature_table_tris[, 122:122]
#View(column_centroid_structure)
mean_centroid_structure = mean(column_centroid_structure)
new_column_centroid_structure = c()
for (row in column_centroid_structure) {
  averaged_row = row/mean_centroid_structure
  new_column_centroid_structure = append(new_column_centroid_structure, averaged_row)
}

# replacing the values with the normalized ones
normalized_feature_table = feature_table_tris[, 1:123] 

normalized_feature_table[, 102:102] <- new_column_numb_not_paired_bases
normalized_feature_table[, 103:103] <- new_column_numb_loop
normalized_feature_table[, 104:104] <- new_column_numb_loop_3_dots_2_parentesis
normalized_feature_table[, 105:105] <- new_column_numb_loop_4_dots_1_parentesis
normalized_feature_table[, 106:106] <- new_column_numb_loop_4_dots_2_parentesis
normalized_feature_table[, 107:107] <- new_column_numb_loop_5_dots_1_parentesis
normalized_feature_table[, 108:108] <- new_column_numb_loop_5_dots_2_parentesis
normalized_feature_table[, 109:109] <- new_column_centrality
normalized_feature_table[, 110:110] <- new_column_dim_loop
normalized_feature_table[, 111:111] <- new_column_numb_paired_bases
normalized_feature_table[, 118:118] <- new_column_numb_bulges
normalized_feature_table[, 119:119] <- new_column_numb_internal_loop
normalized_feature_table[, 120:120] <- new_column_minimum_free_energy
normalized_feature_table[, 121:121] <- new_column_ensemble_free_energy
normalized_feature_table[, 122:122] <- new_column_centroid_structure

normalized_feature_table = feature_table_tris[, 2:123]  # do not consider the miRNA name column 
normalized_feature_table_with_names = feature_table_tris[, 1:123]  # consider also the miRNA name column

##########
# BOXPLOTS
##########

boxplot(normalized_feature_table[,1]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of AA", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA"), ylab = "Distribution of the frequency of AA", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,2]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of AG", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of AG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,3]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of AC", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA"), ylab = "Distribution of the frequency of AC", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,4]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of AU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of AU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,5]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of GA", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of GA", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,6]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of GG", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of GG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,7]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of GC", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of GC", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,8]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of GU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of GU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,9]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of CA", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of CA", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,10]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of CG", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of CG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,11]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of CC", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of CC", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,12]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of CU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of CU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,13]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of UA", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of UA", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,14]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of UG", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of UG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,15]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of UC", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of UC", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,16]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of UU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of UU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,17]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of AAA", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of AAA", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,18]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of AAG", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled", 
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of AAG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,19]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of AAC", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of AAC", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,20]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of AAU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of AAU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,21]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of AGA", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of AGA", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,22]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of AGG", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of AGG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,23]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of AGC", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of AGC", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,24]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of AGU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of AGU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,25]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of ACA", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of ACA", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,26]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of ACG",sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled", 
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of ACG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,27]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of ACC", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of ACC", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,28]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of ACU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of ACU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,29]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of AUA", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of AUA", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,30]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of AUG", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of AUG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,31]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of AUC", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of AUC", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,32]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of AUU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of AUU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,33]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of GAA",sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled", 
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of GAA", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,34]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of GAG", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of GAG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,35]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of GAC", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of GAC", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,36]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of GAU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of GAU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,37]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of GGA", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of GGA", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,38]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of GGG", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of GGG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,39]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of GGC", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of GGC", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,40]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of GGU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of GGU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,41]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of GCA", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of GCA", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,42]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of GCG", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of GCG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,43]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of GCC", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of GCC", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,44]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of GCU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of GCU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,45]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of GUA", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of GUA", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,46]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of GUG", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of GUG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,47]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of GUC", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of GUC", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,48]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of GUU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA"), ylab = "Distribution of the frequency of GUU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,49]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of CAA", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of CAA", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,50]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of CAG", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of CAG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,51]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of CAC",sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled", 
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of CAC", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,52]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of CAU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of CAU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,53]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of CGA", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of CGA", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,54]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of CGG", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of CGG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,55]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of CGC", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of CGC", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,56]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of CGU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of CGU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,57]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of CCA", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of CCA", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,58]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of CCG", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of CCG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,59]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of CCC", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of CCC", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,60]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of CCU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of CCU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,61]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of CUA", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of CUA", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,62]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of CUG", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of CUG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,63]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of CUC", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of CUC", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,64]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of CUU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of CUU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,65]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of UAA", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of UAA", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,66]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of UAG", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of UAG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,67]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of UAC", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of UAC", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,68]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of UAU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of UAU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,69]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of UGA", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of UGA", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,70]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of UGG", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of UGG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,71]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of UGC", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of UGC", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,72]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of UGU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of UGU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,73]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of UCA", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of UCA", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,74]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of UCG", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of UCG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,75]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of UCC", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of UCC", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,76]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of UCU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of UCU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,77]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of UUA", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of UUA", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,78]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of UUG", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of UUG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,79]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of UUC", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of UUC", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,80]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of UUU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of UUU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,81]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of ((( and )))", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA"), ylab = "Distribution of the frequency of ((( and )))", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,82]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of ((. and )).", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of ((. and )).", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,83]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of (.( and ).)", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of (.( and ).)", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,84]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of (.. and )..", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of (.. and )..", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,85]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of .(( and .))", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of .(( and .))", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,86]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of .(. and .).", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of .(. and .).", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,87]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of ..( and ..)", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of ..( and ..)", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,88]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of ...", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of ...", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,89]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of (()", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of (()", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,90]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of ()(", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of ()(", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,91]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of ())", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of ())", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,92]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of )((", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of )((", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,93]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of )()", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of )()", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,94]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of ))(", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of ))(", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,95]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of ().", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of ().", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,96]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of )(.", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of )(.", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,97]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of (.)", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of (.)", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,98]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of ).(", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of ).(", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,99]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of .()", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of .()", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,100]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of .)(", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of .)(", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,101]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Total number of not paired bases", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the total number of not paired bases", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,102]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Total number of loops ", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA"), ylab = "Distribution of the total number of loops", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,103]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Total number of loops, defining them as ((.{3,}))", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the total number of loops, defining them as ((.{3,}))", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,104]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of Total number of loops, defining them as (.{4,})", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the total number of loops, defining them as (.{4,})", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,105]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Total number of loops, defining them as ((.{4,}))", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the total number of loops, defining them as ((.{4,}))", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,106]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Total number of loops, defining them as (.{5,})", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the total number of loops, defining them as (.{5,})", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,107]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Total number of loops, defining them as ((.{5,}))", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the total number of loops, defining them as ((.{5,}))", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,108]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Centrality of the loop of interest", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the centrality of the loop of interest", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,109]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Dimension of the loop of interest", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the dimension of the loop of interest", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,110]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Total number of paired bases", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the total number of paired bases", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,111]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of paired CG", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of paired CG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,112]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of paired GC", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of paired GC", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,113]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of paired UA", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of paired UA", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,114]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of paired AU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of paired AU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,115]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of paired UG", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of paired UG", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,116]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Frequency of paired GU", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the frequency of paired GU", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,117]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Total number of bulges", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the total number of bulges", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,118]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Total number of internal loops", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the total number of internal loops", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,119]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Minimum free energy", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the minimum free energy", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,120]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Ensamble free energy", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the ensamble free energy", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)

boxplot(normalized_feature_table[,121]~ normalized_feature_table$Class, col = c("green", "red"), main = "BOXPLOT: Centroid structure energy", sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled",
        names = c("pre-miRNAs","pseudo pre-miRNA", "our reshuffled pseudo"), ylab = "Distribution of the centroid structure energy", xlab = "",  
        cex.axis = 0.7, cex.main = 1, outline = FALSE, las = 1, varwidth = TRUE)


# for every feature the WILCOX TEST is measured

iterator = 1
colnames(normalized_feature_table)
for (column in normalized_feature_table[,1:121])
{
  name_feature = colnames(normalized_feature_table[iterator])
 distribution_mirna = normalized_feature_table[1:1918,iterator]
 distribution_pseudo_mirnas = normalized_feature_table[1919:3836,iterator]
 result_wilcox.test = wilcox.test(distribution_mirna,distribution_pseudo_mirnas)
 #print(name_feature)
 #print(result_wilcox.test$p.value)
 print(c(name_feature, result_wilcox.test$p.value))
 iterator = iterator + 1
}


