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

######
# PCA
#####
numeric_class_mirna_pseudo = c(rep(c(1), 1918), rep(c(2), 1918))  # 1=real mirna; 2=pseudo
normalized_feature_table[,122] <- numeric_class_mirna_pseudo
PCA_features = prcomp(normalized_feature_table[,-122]) # do not considering the numeric class
attributes(PCA_features)
M = PCA_normalized_feature_table$x
R = PCA_normalized_feature_table$rotation

#PCA_biplot = biplot(PCA_features)

plot_first_and_second_feature = plot(M[,1],M[,2],col=rep(c("green","red"),table(numeric_class_mirna_pseudo)),main="PCA FEATURES",sub = "Human pre-miRNAs from miRBase and pseudo pre-miRNAs reshuffled", xlab="PC1",ylab="PC2", cex.main=2, cex.lab=1.5, pch=20)
# Now use arrows() function to draw arrows. To do that, we require the starting point (0,0) of the plot and the final point, corresponding to the tip of the arrow.
# Starting and ending points need to be provided in the form of a vector for both the X and Y axes, separetely. 
# Our arrows correspond to the rotation matrix (i.e eigenvectors) of out PCA: the rotation matrix contains only the direction of the components.
# Moreover, to get the "final" coordinate we need to multiply the rotation (direction) by the maximun eigenvalue on that direction. 
#max_X = R[,1]*max(M[,1])
#max_Y = R[,2]*max(M[,2])
#arrows(0,0,max_X,max_Y,col="purple") 
legend(-230,60,legend=c("pre-miRNAs","pseudo pre-miRNA"),col=c("green","red"),pch=20,cex=0.8, bty = "n")
#The text function is similar to arrows: we just need to specify the final coordinates and the text (labels parameter) that you want to print (in our case the names of the features)
#text(max_X,max_Y,labels=rownames(R),pos=1,offset=0.5) 


