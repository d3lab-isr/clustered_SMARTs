## This script generate data from a prototypical cluster randomized SMART
## with no covariates

library(MASS)

##Set means and variance values
p_1 = .3 ##Probability of response to treatment 1
p_2 = .2 ##Probability of response to treatment -1
cell_mu = c(10,8,2,4,-2,3)  ##ordered vector of means for each cell of SMART
cell_var = c(100,100,100,100,100,100) ##ordered vector of variances for each cell of SMART
cell_cor = c(.1, .2, .15, .12, .11, .07) ##ordered vector of correlations for each cell of SMART 
cell_cov = cell_var*cell_cor 

## Return table of treatment means and treatment variances
q_1 = 1-p_1
q_2 = 1- p_2
treat_mu = c(cell_mu[1]*p_1 + cell_mu[2]*q_1, cell_mu[1]*p_1 + cell_mu[3]*q_1,
cell_mu[4]*p_2 + cell_mu[5]*q_2, cell_mu[4]*p_2 + cell_mu[6]*q_2) #Treatment means

treat_var = c(cell_var[1]*p_1 + cell_var[2]*q_1, cell_var[1]*p_1 + cell_var[3]*q_1, 	   cell_var[4]*p_2 + cell_var[5]*q_2, cell_var[4]*p_2 + cell_var[6]*q_2)  ##Treatment Variances

mean_diff_term = c((cell_mu[1]-cell_mu[2])^2*p_1*q_1, (cell_mu[1]-cell_mu[3])^2*p_1*q_1,
(cell_mu[4]-cell_mu[5])^2*p_2*q_2, (cell_mu[4]-cell_mu[6])^2*p_2*q_2)

treat_var = treat_var + mean_diff_term

treat_cov = c(cell_cov[1]*p_1 + cell_cov[2]*q_1, cell_cov[1]*p_1 + cell_cov[3]*q_1, 	   cell_cov[4]*p_2 + cell_cov[5]*q_2, cell_cov[4]*p_2 + cell_cov[6]*q_2)  ##Treatment Covariances

treat_cov = treat_cov + mean_diff_term
treat_cor = treat_cov/treat_var

treat_summary = data.frame(treat_mu, treat_var, treat_cor)
row.names(treat_summary) = c("A1=1,A2=1", "A1=1,A2=-1", "A1=-1,A2=1", "A1=-1,A2=-1" )

m = 5

## Create the data set
N = 200 ##This is the number of clusters
treat1_data = matrix(nrow = 0, ncol = 3)
treat2_data = matrix(nrow = 0, ncol = 3)
treat3_data = matrix(nrow = 0, ncol = 3)
treat4_data = matrix(nrow = 0, ncol = 3)
m_vec = rep(m,N) #rpois(N, m) + 1 
unique_m = unique(m_vec)
max_m = max(unique_m)
p_A1 = .5 ##Probability of receiving treatment A1
p_A2 = .5 ##Probability of receiving treatment A2 for non-respnoders
full_data = matrix(nrow = 0, ncol = 5)
cell_cor_mat = list()
cell_chol_cor_mat = list()
cell_cor_mat[[max_m]] = list()
cell_chol_cor_mat[[max_m]] = list()
for(j in 1:6){mat = matrix(cell_cor[j], max_m, max_m)
	diag(mat) = 1
	cell_cor_mat[[max_m]][[j]] = mat
	cell_chol_cor_mat[[max_m]][[j]] = t(chol(mat))}
if(length(unique_m) > 1){
for(i in unique_m){
if(i==max_m){}else{
cell_cor_mat[[i]] = list()
cell_chol_cor_mat[[i]] = list()
	for(j in 1:6){
		cell_cor_mat[[i]][[j]] = (cell_cor_mat[[max_m]][[j]])[1:i,1:i]
		cell_chol_cor_mat[[i]][[j]] = (cell_chol_cor_mat[[max_m]][[j]])[1:i,1:i]}}}}

for(i in 1:N){
	m = m_vec[i]
	A1 = 2*rbinom(1, 1, p_A1) + 1
	if(A1 == 1){
		R = rbinom(1, 1, p_1)
		
		## Generate data for A_1 = 1 and R = 1
		if(R == 1){
			A2 = 1
			Y = rnorm(m, 0, sd = sqrt(cell_var[1]))
			Y = cell_chol_cor_mat[[m]][[1]]%*%Y
			Y = Y + cell_mu[1]
			temp_mat = cbind(rep(i,m), rep(1,m), rep(1,m), rep(1,m), Y)
			full_data = rbind(full_data, temp_mat)
			treat1_data = rbind(treat1_data, temp_mat[,c(1, 3, 5)])
			treat2_data = rbind(treat2_data, temp_mat[,c(1, 3, 5)])
			
		} else{
			A2 = 2*rbinom(1, 1, p_A2) -1
			
			## Generate data for A_1 = 1 and R = 0 and A_2 = 1
			if(A2 == 1){
				Y = rnorm(m, 0, sd = sqrt(cell_var[2]))
				Y = cell_chol_cor_mat[[m]][[2]]%*%Y
				Y = Y + cell_mu[2]
				temp_mat = cbind(rep(i,m), rep(1,m), rep(0,m), rep(1,m), Y)
				full_data = rbind(full_data, temp_mat)
				treat1_data = rbind(treat1_data, temp_mat[,c(1, 3, 5)])
			
			## Generate data for A_1 = 1 and R = 0 and A_2 = -1	
			} else{
				Y = rnorm(m, 0, sd = sqrt(cell_var[3]))
				Y = cell_chol_cor_mat[[m]][[3]]%*%Y
				Y = Y + cell_mu[3]
				temp_mat = cbind(rep(i,m), rep(1,m), rep(0,m), rep(-1,m), Y)
				full_data = rbind(full_data, temp_mat)
				treat2_data = rbind(treat2_data, temp_mat[,c(1, 3, 5)])		
			}
		}
		
	} else{
		R = rbinom(1, 1, p_2)
		
		## Generate data for A_1 = -1 and R = 1
		if(R == 1){
			A2 = 1
			Y = rnorm(m, 0, sd = sqrt(cell_var[4]))
			Y = cell_chol_cor_mat[[m]][[4]]%*%Y
			Y = Y + cell_mu[4]
			temp_mat = cbind(rep(i,m), rep(-1,m), rep(1,m), rep(1,m), Y)
			full_data = rbind(full_data, temp_mat)
			treat3_data = rbind(treat3_data, temp_mat[,c(1, 3, 5)])
			treat4_data = rbind(treat4_data, temp_mat[,c(1, 3, 5)])
			
		} else{
			A2 = 2*rbinom(1, 1, p_A2) -1
			
			## Generate data for A_1 = -1 and R = 0 and A_2 = 1
			if(A2 == 1){
				Y = rnorm(m, 0, sd = sqrt(cell_var[4]))
				Y = cell_chol_cor_mat[[m]][[5]]%*%Y
				Y = Y + cell_mu[5]
				temp_mat = cbind(rep(i,m), rep(-1,m), rep(0,m), rep(1,m), Y)
				full_data = rbind(full_data, temp_mat)
				treat3_data = rbind(treat3_data, temp_mat[,c(1, 3, 5)])
							
			## Generate data for A_1 = -1 and R = 0 and A_2 = -1	
			} else{
				Y = rnorm(m, 0, sd = sqrt(cell_var[6]))
				Y = cell_chol_cor_mat[[m]][[6]]%*%Y
				Y = Y + cell_mu[6]
				temp_mat = cbind(rep(i,m), rep(-1,m), rep(0,m), rep(-1,m), Y)
				full_data = rbind(full_data, temp_mat)
				treat4_data = rbind(treat4_data, temp_mat[,c(1, 3, 5)])
			}
		}
	}
}
full_data_frame = as.data.frame(full_data)
names(full_data_frame) = c("clus_id", "A1", "R", "A2", "Y")

## full_data_frame contains the generated data set
