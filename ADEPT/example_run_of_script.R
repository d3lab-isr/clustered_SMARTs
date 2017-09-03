## First make sure you are in the same working directory as the data set
data_set = read.csv("sample_data_ADEPT")

# Next, make sure you have read in the function ADEPT_regression from the 
# adept_estimation.R file

## regressing using all covariates
regression_output1 = ADEPT_regression(data_set, "Y", "A1_ind", "response", "A2_ind", covariates=c("covariate_1", "covariate_2", "covariate_3"), "cluster.id")  ## note these names match the names in the data set

## regressing using no covariates
regression_output2 = ADEPT_regression(data_set, "Y", "A1_ind", "response", "A2_ind", covariates=NULL, "cluster.id")  ## note these names match the names in the data set
