## First make sure you are in the same working directory as the data set
data_set = read.csv("sample_data_prototypical")

## regressing using all covariates
regression_output1 = prototypical_regression(data_set, "response", "A1", "R", "A2", covariates=c("cov1", "cov2"), "clust")  ## note these names match the names in the data set