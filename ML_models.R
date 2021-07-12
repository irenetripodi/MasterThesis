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

prova = normalized_feature_table = feature_table_tris[, 80:90]
View(prova)


#saveRDS(normalized_feature_table, "normalized_feature_table.rds")  # save it
#saveRDS(normalized_feature_table_with_names, "normalized_feature_table_with_names.rds")  

################
# DATA SPLITTING
################
# The following step consists in the partion of the dataset, in order to divide it in 2 groups: 
# - training group, that is the one used to train the tool 
# - testing group, that is the one used then to validate your method and to see if the tool works well also with data not already seen

#install.packages('caret')
library(caret)

# DATA SPLITTING
training_index = createDataPartition(normalized_feature_table$Class, p = 0.8, list = F, times = 1)
# p is the percentage of data that goes to training
# times is the number of partitions to create
# list = FALSE avoids returning the data as a list

training_group <- normalized_feature_table[ training_index,]
#head(training_group)
saveRDS(training_group, "training_group_RF_reshuffled.rds")  # save it
testing_group  <- normalized_feature_table[-training_index,]
saveRDS(testing_group, "testing_group_RF_reshuffled.rds")  # save it

#training_group_with_names <- normalized_feature_table_with_names[ training_index,]
testing_group_with_names <- normalized_feature_table_with_names[-training_index,]
##########
# TRAINING
##########
# In this step allows to train the tool with the specific part of the dataset selected. 
# The train function can be used to:
# - evaluate, using resampling, the effect of model tuning parameters on performance
# - choose the “optimal” model across these parameters
# - estimate model performance from a training set

?trainControl
fit_control <- trainControl(method = "repeatedcv", number = 10, repeats = 10)
fit_control
# method tells the resampling method
# number	tells the number of folds or number of resampling iterations
# repeats: for repeated k-fold cross-validation only: the number of complete sets of folds to compute


?train
training_tris = train(Class ~ ., data = training_group, method = "rf", trControl = fit_control, tuneGrid = expand.grid(mtry=c(5,10:30,50,70,100)))
# trControl	indicates the list of values that define how this function acts; tuneGrid provides a data frame with possible tuning values. 
training_tris  # the Accuracy obtained is quite good 

# Kappa values represents the degree of accuracy and reliability in a statistical classification
# it tells how many data (the %) for example coming from "Conserved in both chimpazee and mouse" have been classified correctly in "Conserved in both chimpazee and mouse" class.  
confusionMatrix(training_tris)

attributes(training_tris)
training_tris$finalModel  # per essere più precisi si potrebbe salvare solo finalModel (i risultati sono più compatti)

# To see how much the feature are important for the determination of the prediction of the class use:  
training_tris$finalModel$importance 
# We can sort the feature according values (with decreasing order) in order to see which are the features more important:
sorted_features_tris <- training_tris$finalModel$importance[order(training_tris$finalModel$importance, decreasing = T),]
sorted_features_tris
# Appear more important the mismacths against structural elements to distinguish which miRNA are conserved and where and which not!


#########
# TESTING
#########
# The last step is to test the tool with data not already seen (testing group discarded before)
?predict
testing_tris <- predict(training_tris, testing_group)
testing_tris  # the result is a list of the classes predicted for the testing group

confusionMatrix(testing_tris, testing_group$Class) 



################
# SAVE THE MODEL
################

#save(training_tris, file= "training_tris.Rdata")

# Now we can save our best model in a file so that we can load it later and make other predictions.
model <- training_tris  # model with all the training data
saveRDS(model, "model.rds")  # save it
saveRDS(testing_tris, "model_testing.rds")  
training_tris$finalModel  # per essere più precisi si potrebbe salvare solo finalModel (i risultati sono più compatti)
