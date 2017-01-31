### Comments for Tim ###

## check this with each prototypical sim,
## make it possible to only do independent working correlation?
## make it possible to do different weightings, make sure to change in estimation and r and sigma estiamtion
## make sure when doing t-statistic, I subtract 2 times cov
## check for when I have more than one covariate and no covariates
## run code on data set with no preset values like d which came from generating the data
## get rid of all silly comments
## could functionize even further if i wanted, like computing v hat getting mean, getting se
## need to commment in functions and in standard error estimation, comments should be relatively the same

###					###

## Currently this does the prototypical regression highlighted in the paper, with fixed weights of 2 and 4.


prototypical_regression = function(dat, outcome, A1, R, A2, covariates = NULL, cluster) {
	## dat is a data frame where all data is stored
	
	## outcome is a character (i.e. string) indicating the name of the outcome in dat
	
	## A1 is a character (i.e. string) indicating the name of the A1 in dat
	## note the A1 column should be coded +/-1 in dat (as in paper)
	
	## R is a character (i.e. string) indicating the name of the responders in dat
	## note the A1 column should be coded 1-responder/0-nonresponder
	## in dat (as in paper)
	
	## A2 is a character (i.e. string) indicating the name of the A2 in dat
	## note the A2 column should be coded +/-1 in dat (as in paper)
	
	## covariates is a character vector indicating the names of the covariates
	## in dat.  Note that our function fits a linear model in terms of the covariates
	## (i.e. the marginal mean model in the paper), extensions to non-linear models
	## can be done by transforming the covariates prior to analysis. Leave blank
	## if there are no covariates
	
	## cluster is a character indicating the name of the cluster identifiers
	## in dat. Cluster identifiers should be coded numerically
	
	## Begin code
	
	## Get preliminary information
	clus_id = dat[,cluster] ## get the id's for each cluster
	d = length(covariates)  ## get the dimension of covariates
	tot_dim = d+5 ## get the total number of dimensions
	N = length(unique(clus_id)) ## get the total number of clusters 
	
	
	## For case where I have covariates. In this case specifying V matters and hence
	## we do V estimation iteration highlighted in the paper
	
	
	## Organize data into a data frame which is order by cluster id, A1, R, A2, covariates
	## outcome
	full_data_frame = dat[,c(cluster, A1, R, A2, covariates, outcome)]
	cov_names = NULL
	if(d>0){
	for(i in 1:d){cov_names[i] = paste("X", i, sep="")}}
	names(full_data_frame) = c("clus_id", "A1", "R", "A2", cov_names, "Y")
	
	## This function takes our data set, and a specified exchangeable V (which is specified
	## through r-correlation, and sigma-variace) and solves the estimating equation.
	## It returns our estimates of the covariates and new estimates for sigma and R 
	## using the method highlighted in the paper
	## It also returns the B matrix, which will be used for standard error estimation
	mean_and_V_estimate = function(full_data_frame, r, sigma){
		## full_data_frame is our nicely ordered data frame
		
		## r is the pre-specified correlations for each regime
		## sigma is the the pre-specified variances for each regime
		## r and sigma together completely specify V (V = sigma*exch(r))
		## r and sigma are length 3 vectors with each cell specifying 
		## the sigma and r for DTR (1,1), (1,-1), (-1,1), and (-1,-1) respectively
	
	
		tot_dim = dim(full_data_frame)[2] ## get total dimension
		d = tot_dim - 5 ## get covariates dimension
		clus_id = full_data_frame$clus_id ## get the cluster ids
		
		## Now that V is specified, we have everything we need to solve the estimating
		## equation in the paper.  We do this using linear algebra. Specifically, we 
		## rearrange the estimating equation as a solution to B^(-1)A
	
		## A is the vector which represent sum(W*I*D^T*V*Y) in paper, the vector that 
		## does not involve any covariates
		## B is the matrix which represents sum(W*I*D^T*V*D) in paper, this
		## is the matrix multiplied by the covariates (beta, eta) in the equation
		## To solve the estimating equation, you merely need to find B^(-1)*A
		A = rep(0, d+4)
		B = matrix(0, d+4, d+4)
		
	
		
		
		## Compute A and B by iteratively summing over the contribution of each cluster
		for(i in unique(clus_id)){
			
			## obtain data set for cluster i
			cur_data = full_data_frame[clus_id==i,] ## data set
			cur_resp = cur_data$Y ## outcomes
			cur_m = length(cur_resp) ## cluster size
			ind_vec = cur_data[1,2:4]  ## A1, R, A2, constant across all ind in cluster
			if(d>0) {cov_vec = as.matrix(cur_data[,5:(tot_dim-1)]) ##Obtain our covariates
			} else{cov_vec = NULL}
			
			## Now we go through all possibilities of A1, R, and A2 and use this information
			## to figure out which regime our cluster is consistent with
			
			## A1 = 1 and R = 1, these clusters are consistent with 1,1 and 1,-1
			if(ind_vec[1]==1 && ind_vec[2] == 1){
				
				## We calculate the cluster's contribution makes to A and B due 
				## due to consistency with DTR 1,1
				
				## calculate V^-1 in estimating equation
				cur_cor = r[1]
				cur_sigma = sigma[1]
				outer_term = 1/(cur_sigma*(1+(cur_m-2)*cur_cor + (1-cur_m)*cur_cor^2))
				diag_part = ((cur_m-2)*cur_cor + 1)*outer_term
				off_diag = -cur_cor*outer_term 
				v_inv = matrix(off_diag, cur_m, cur_m) 
				diag(v_inv) = diag_part
				
				## calculate D in estimating equation
				cur_d = cbind(rep(1,cur_m), rep(1,cur_m), rep(1,cur_m), rep(1,cur_m), cov_vec)
				
				## add in contribution to A and B
				A = A + 2*t(cur_d)%*%v_inv%*%cur_resp
				B = B + 2*t(cur_d)%*%v_inv%*%cur_d
				
				
				## We calculate the cluster's contribution makes to A and B due 
				## due to consistency with DTR 1,-1
				
				## calculate V^-1 in estimating equation
				cur_cor = r[2]
				cur_sigma = sigma[2]
				outer_term = 1/(cur_sigma*(1+(cur_m-2)*cur_cor + (1-cur_m)*cur_cor^2))
				diag_part = ((cur_m-2)*cur_cor + 1)*outer_term
				off_diag = -cur_cor*outer_term 
				v_inv = matrix(off_diag, cur_m, cur_m)
				diag(v_inv) = diag_part
				
				## calculate D in estimating equation
				cur_d = cbind(rep(1,cur_m), rep(1,cur_m), rep(-1,cur_m), rep(-1,cur_m), cov_vec)
				
				## add in contribution to A and B
				A = A + 2*t(cur_d)%*%v_inv%*%cur_resp
				B = B + 2*t(cur_d)%*%v_inv%*%cur_d
			}
			
			
			
			## A1 = 1 and R = 0
			if(ind_vec[1]==1 && ind_vec[2] == 0){
				
				## A2 = 1, these clusters are consistent with 1,1 only
				if(ind_vec[3] == 1){
					
					## We calculate the cluster's contribution makes to A and B due 
					## due to consistency with DTR 1,1
					
					## calculate V^-1 in estimating equation
					cur_cor = r[1]
					cur_sigma = sigma[1]
					outer_term = 1/(cur_sigma*(1+(cur_m-2)*cur_cor + (1-cur_m)*cur_cor^2))
					diag_part = ((cur_m-2)*cur_cor + 1)*outer_term
					off_diag = -cur_cor*outer_term 
					v_inv = matrix(off_diag, cur_m, cur_m)
					diag(v_inv) = diag_part
					
					## calculate D in estimating equation
					cur_d = cbind(rep(1,cur_m), rep(1,cur_m), rep(1,cur_m), rep(1,cur_m), cov_vec)
				
					## add in contribution to A and B
					A = A + 4*t(cur_d)%*%v_inv%*%cur_resp
					B = B + 4*t(cur_d)%*%v_inv%*%cur_d
					
					
				## A2 = -1, these clusters are consistent with 1,-1 only
				} else{
					
					## We calculate the cluster's contribution makes to A and B due 
					## due to consistency with DTR 1,-1
					
					## calculate V^-1 in estimating equation
					cur_cor = r[2]
					cur_sigma = sigma[2]
					outer_term = 1/(cur_sigma*(1+(cur_m-2)*cur_cor + (1-cur_m)*cur_cor^2))
					diag_part = ((cur_m-2)*cur_cor + 1)*outer_term
					off_diag = -cur_cor*outer_term 
					v_inv = matrix(off_diag, cur_m, cur_m)
					diag(v_inv) = diag_part
					
					## calculate D in estimating equation
					cur_d = cbind(rep(1,cur_m), rep(1,cur_m), rep(-1,cur_m), rep(-1,cur_m), cov_vec)
				
					## add in contribution to A and B
					A = A + 4*t(cur_d)%*%v_inv%*%cur_resp
					B = B + 4*t(cur_d)%*%v_inv%*%cur_d
					}
			}
			
			
			## A1 = 1 and R = 1, these clusters are consistent with -1,1 and -1,-1
			if(ind_vec[1]==-1 && ind_vec[2] == 1){
				
				## We calculate the cluster's contribution makes to A and B due 
				## due to consistency with DTR -1,1
				
				## calculate V^-1 in estimating equation
				cur_cor = r[3]
				cur_sigma = sigma[3]
				outer_term = 1/(cur_sigma*(1+(cur_m-2)*cur_cor + (1-cur_m)*cur_cor^2))
				diag_part = ((cur_m-2)*cur_cor + 1)*outer_term
				off_diag = -cur_cor*outer_term 
				v_inv = matrix(off_diag, cur_m, cur_m) 
				diag(v_inv) = diag_part
				
				## calculate D in estimating equation
				cur_d = cbind(rep(1,cur_m), rep(-1,cur_m), rep(1,cur_m), rep(-1,cur_m), cov_vec)
				
				## add in contribution to A and B
				A = A + 2*t(cur_d)%*%v_inv%*%cur_resp
				B = B + 2*t(cur_d)%*%v_inv%*%cur_d
				
				
				## We calculate the cluster's contribution makes to A and B due 
				## due to consistency with DTR -1,-1
				
				## calculate V^-1 in estimating equation
				cur_cor = r[4]
				cur_sigma = sigma[4]
				outer_term = 1/(cur_sigma*(1+(cur_m-2)*cur_cor + (1-cur_m)*cur_cor^2))
				diag_part = ((cur_m-2)*cur_cor + 1)*outer_term
				off_diag = -cur_cor*outer_term 
				v_inv = matrix(off_diag, cur_m, cur_m)
				diag(v_inv) = diag_part
				
				## calculate D in estimating equation
				cur_d = cbind(rep(1,cur_m), rep(-1,cur_m), rep(-1,cur_m), rep(1,cur_m), cov_vec)
				
				## add in contribution to A and B
				A = A + 2*t(cur_d)%*%v_inv%*%cur_resp
				B = B + 2*t(cur_d)%*%v_inv%*%cur_d
			}
			
			
			
			## A1 = -1 and R = 0
			if(ind_vec[1]== -1 && ind_vec[2] == 0){
				
				## A2 = 1, these clusters are consistent with -1,1 only
				if(ind_vec[3] == 1){
					
					## We calculate the cluster's contribution makes to A and B due 
					## due to consistency with DTR -1,1
					
					## calculate V^-1 in estimating equation
					cur_cor = r[3]
					cur_sigma = sigma[3]
					outer_term = 1/(cur_sigma*(1+(cur_m-2)*cur_cor + (1-cur_m)*cur_cor^2))
					diag_part = ((cur_m-2)*cur_cor + 1)*outer_term
					off_diag = -cur_cor*outer_term 
					v_inv = matrix(off_diag, cur_m, cur_m)
					diag(v_inv) = diag_part
					
					## calculate D in estimating equation
					cur_d = cbind(rep(1,cur_m), rep(-1,cur_m), rep(1,cur_m), rep(-1,cur_m), cov_vec)
				
					## add in contribution to A and B
					A = A + 4*t(cur_d)%*%v_inv%*%cur_resp
					B = B + 4*t(cur_d)%*%v_inv%*%cur_d
					
					
				## A2 = -1, these clusters are consistent with -1,-1 only
				} else{
					
					## We calculate the cluster's contribution makes to A and B due 
					## due to consistency with DTR -1,-1
					
					## calculate V^-1 in estimating equation
					cur_cor = r[4]
					cur_sigma = sigma[4]
					outer_term = 1/(cur_sigma*(1+(cur_m-2)*cur_cor + (1-cur_m)*cur_cor^2))
					diag_part = ((cur_m-2)*cur_cor + 1)*outer_term
					off_diag = -cur_cor*outer_term 
					v_inv = matrix(off_diag, cur_m, cur_m)
					diag(v_inv) = diag_part
					
					## calculate D in estimating equation
					cur_d = cbind(rep(1,cur_m), rep(-1,cur_m), rep(-1,cur_m), rep(1,cur_m), cov_vec)
				
					## add in contribution to A and B
					A = A + 4*t(cur_d)%*%v_inv%*%cur_resp
					B = B + 4*t(cur_d)%*%v_inv%*%cur_d
					}
			}
		}
		
		## Now that A and B are fully updated, we calculate our coefficient estimates
		## mean_val is a 4+d length vector containing estimates for
		##  beta_0, beta_1, beta_2, beta_3, eta_1,...eta_d, respectively.
		mean_val = solve(B, A)
		
		
		## Now we obtain updated for r and sigma based on our mean estimate
		
		## storage for appropriate sums needed to calculate r and sigma.
		## the var_sum and rho_sum are the numerators for our estimatates
		## and the denom_var and denom_rho are the denominators.  The 1,2,3,4
		## correspond, again to DTR (1,1), (1,-1), (-1,1), (-1,-1)
		var_sum1 = 0
		rho_sum1 = 0
		denom_var1 = 0
		denom_rho1 = 0
		
		var_sum2 = 0
		rho_sum2 = 0
		denom_var2 = 0
		denom_rho2 = 0
		
		var_sum3 = 0
		rho_sum3 = 0
		denom_var3 = 0
		denom_rho3 = 0
		
		var_sum4 = 0
		rho_sum4 = 0
		denom_var4 = 0
		denom_rho4 = 0
		
		## Iterate through the total number of clusters
		for(i in unique(clus_id)){
			
			## obtain data set for cluster i
			cur_data = full_data_frame[clus_id==i,]  ## data set
			cur_resp = cur_data$Y ## outcomes
			cur_m = length(cur_resp) ##cluster size
			ind_vec = cur_data[1,2:4] ## A1, R, A2, constant across all ind in cluster
			if(d>0) {cov_vec = as.matrix(cur_data[,5:(tot_dim-1)]) ##Obtain our covariates
			} else{cov_vec = NULL}
			
			## Now we go through all possibilities of A1, R, and A2 and use this information
			## to figure out which regime our cluster is consistent with
			
			## A1 = 1 and R = 1, these clusters are consistent with 1,1 and 1,-1
			if(ind_vec[1]==1 && ind_vec[2] == 1){
				
					## Residuals from DTR 1,1 and the residuals contribution to sigma and r
					## estimates
					cur_vec = cbind(rep(1,cur_m), rep(1,cur_m), rep(1,cur_m), rep(1,cur_m), cov_vec)
					cur_resid = cur_resp - cur_vec%*%mean_val
					inner_p = t(cur_resid)%*%cur_resid
		 			outer_p = cur_resid%*%t(cur_resid)
		 			var_sum1 = var_sum1 + inner_p
		 			rho_sum1 = rho_sum1 + sum(outer_p) - inner_p
		 			denom_var1 = denom_var1 + cur_m
		 			denom_rho1 = denom_rho1 + cur_m*(cur_m-1)
					
					## Residuals from DTR 1,-1 and the residuals contribution to sigma and r
					## estimates
					cur_vec = cbind(rep(1,cur_m), rep(1,cur_m), rep(-1,cur_m), rep(-1,cur_m), cov_vec)
					cur_resid = cur_resp - cur_vec%*%mean_val
					inner_p = t(cur_resid)%*%cur_resid
		 			outer_p = cur_resid%*%t(cur_resid)
		 			var_sum2 = var_sum2 + inner_p
		 			rho_sum2 = rho_sum2 + sum(outer_p) - inner_p
		 			denom_var2 = denom_var2 + cur_m
		 			denom_rho2 = denom_rho2 + cur_m*(cur_m-1)		
			}
				
				
			## A1 = 1 and R = 0
			if(ind_vec[1]==1 && ind_vec[2] == 0){
				
				## A2 = 1, these clusters are consistent with 1,1 only
				if(ind_vec[3] == 1){
					
					## Residuals from DTR 1,1 and the residuals contribution to sigma and r
					## estimates
					cur_vec = cbind(rep(1,cur_m), rep(1,cur_m), rep(1,cur_m), rep(1,cur_m), cov_vec)
					cur_resid = cur_resp - cur_vec%*%mean_val
					inner_p = t(cur_resid)%*%cur_resid
		 			outer_p = cur_resid%*%t(cur_resid)
		 			var_sum1 = var_sum1 + 2*inner_p
		 			rho_sum1 = rho_sum1 + 2*sum(outer_p) - 2*inner_p
		 			denom_var1 = denom_var1 + 2*cur_m
		 			denom_rho1 = denom_rho1 + 2*cur_m*(cur_m-1)
					
				## A2 = -1, these clusters are consistent with 1,-1 only
				} else{
					
					## Residuals from DTR 1,-1 and the residuals contribution to sigma and r
					## estimates
					cur_vec = cbind(rep(1,cur_m), rep(1,cur_m), rep(-1,cur_m), rep(-1,cur_m), cov_vec)
					cur_resid = cur_resp - cur_vec%*%mean_val
					inner_p = t(cur_resid)%*%cur_resid
		 			outer_p = cur_resid%*%t(cur_resid)
		 			var_sum2 = var_sum2 + 2*inner_p
		 			rho_sum2 = rho_sum2 + 2*sum(outer_p) - 2*inner_p
		 			denom_var2 = denom_var2 + 2*cur_m
		 			denom_rho2 = denom_rho2 + 2*cur_m*(cur_m-1)	
					}
			}
			
			
			## A1 = -1 and R = 1, these clusters are consistent with -1,1 and -1,-1
			if(ind_vec[1]== -1 && ind_vec[2] == 1){
				
					## Residuals from DTR -1,1 and the residuals contribution to sigma and r
					## estimates
					cur_vec = cbind(rep(1,cur_m), rep(-1,cur_m), rep(1,cur_m), rep(-1,cur_m), cov_vec)
					cur_resid = cur_resp - cur_vec%*%mean_val
					inner_p = t(cur_resid)%*%cur_resid
		 			outer_p = cur_resid%*%t(cur_resid)
		 			var_sum3 = var_sum3 + inner_p
		 			rho_sum3 = rho_sum3 + sum(outer_p) - inner_p
		 			denom_var3 = denom_var3 + cur_m
		 			denom_rho3 = denom_rho3 + cur_m*(cur_m-1)
					
					## Residuals from DTR -1,-1 and the residuals contribution to sigma and r
					## estimates
					cur_vec = cbind(rep(1,cur_m), rep(-1,cur_m), rep(-1,cur_m), rep(1,cur_m), cov_vec)
					cur_resid = cur_resp - cur_vec%*%mean_val
					inner_p = t(cur_resid)%*%cur_resid
		 			outer_p = cur_resid%*%t(cur_resid)
		 			var_sum4 = var_sum4 + inner_p
		 			rho_sum4 = rho_sum4 + sum(outer_p) - inner_p
		 			denom_var4 = denom_var4 + cur_m
		 			denom_rho4 = denom_rho4 + cur_m*(cur_m-1)		
			}
				
				
			## A1 = -1 and R = 0
			if(ind_vec[1]==-1 && ind_vec[2] == 0){
				
				## A2 = 1, these clusters are consistent with -1,1 only
				if(ind_vec[3] == 1){
					
					## Residuals from DTR -1,1 and the residuals contribution to sigma and r
					## estimates
					cur_vec = cbind(rep(1,cur_m), rep(-1,cur_m), rep(1,cur_m), rep(-1,cur_m), cov_vec)
					cur_resid = cur_resp - cur_vec%*%mean_val
					inner_p = t(cur_resid)%*%cur_resid
		 			outer_p = cur_resid%*%t(cur_resid)
		 			var_sum3 = var_sum3 + 2*inner_p
		 			rho_sum3 = rho_sum3 + 2*sum(outer_p) - 2*inner_p
		 			denom_var3 = denom_var3 + 2*cur_m
		 			denom_rho3 = denom_rho3 + 2*cur_m*(cur_m-1)
					
				## A2 = -1, these clusters are consistent with -1,-1 only
				} else{
					
					## Residuals from DTR -1,-1 and the residuals contribution to sigma and r
					## estimates
					cur_vec = cbind(rep(1,cur_m), rep(-1,cur_m), rep(-1,cur_m), rep(1,cur_m), cov_vec)
					cur_resid = cur_resp - cur_vec%*%mean_val
					inner_p = t(cur_resid)%*%cur_resid
		 			outer_p = cur_resid%*%t(cur_resid)
		 			var_sum4 = var_sum4 + 2*inner_p
		 			rho_sum4 = rho_sum4 + 2*sum(outer_p) - 2*inner_p
		 			denom_var4 = denom_var4 + 2*cur_m
		 			denom_rho4 = denom_rho4 + 2*cur_m*(cur_m-1)	
					}
			}
		}
		
		## We now calculate the our estimates of r and sigma
		var_est1 = var_sum1/(denom_var1-d-1)
		r_est1 = rho_sum1/denom_rho1/var_est1 
		
		var_est2 = var_sum2/(denom_var2-d-1)
		r_est2 = rho_sum2/denom_rho2/var_est2
		
		var_est3 = var_sum3/(denom_var3-d-1)
		r_est3 = rho_sum3/denom_rho3/var_est3
		
		var_est4 = var_sum4/(denom_var4-d-1)
		r_est4 = rho_sum4/denom_rho4/var_est4
		
		return(list(mean_val, c(var_est1, var_est2, var_est3, var_est4), c(r_est1, r_est2, r_est3, r_est4), B))
	}
		
	## We now iterate our estimator, this iteration is the exact same as the paper
	
	## The first time is with V specified as the identity matrix
	first_round = mean_and_V_estimate(full_data_frame = full_data_frame, r=c(0,0,0,0), sigma=c(1,1,1,1))
	
	## The second time uses the estimated values of r and sigma from the first round
	second_round = mean_and_V_estimate(full_data_frame = full_data_frame, r=first_round[[3]], sigma=first_round[[2]])
	
	## The third time used the estimated values of r and sigma from the second round
	third_round = mean_and_V_estimate(full_data_frame = full_data_frame, r=second_round[[3]], sigma=second_round[[2]])
	
	mean_val = third_round[[1]]   ## These are our final covariate estimates
	
	## We now obtain variance-covariance estimates.  First we need to keep values from our estimation
	## Note that we need to use the V used in estimation (i.e. second round V) 
	## when computing the standard errors. We do not care about the newest estimation of V. 
	## (i.e. returned from third round estimation)
	
	B = third_round[[4]]  ## This is the B matrix which will be used for standard errors
	r=second_round[[3]]  ## This is the r used for estimation
	sigma=second_round[[2]] ## This is the sigma used for estimation
	
	bread_inv = B/N  ## This is J in the paper, i.e. the outer matrix in covariance estimation
	
	##Find the meat i.e. the inner matrix by summing over the contribution of each cluster
	meat = matrix(0,d+4,d+4)
	for(i in unique(clus_id)){
		
		## obtain data set for cluster i
		cur_data = full_data_frame[clus_id==i,]  ## data set
		cur_resp = cur_data$Y ## outcomes
		cur_m = length(cur_resp) ## cluster size
		ind_vec = cur_data[1,2:4] ## A1, R, A2, constant across all ind in cluster
		if(d>0) {cov_vec = as.matrix(cur_data[,5:(tot_dim-1)]) ##Obtain our covariates
		} else{cov_vec = NULL}
			
		
		## Now we go through all possibilities of A1, R, and A2 and use this information
		## to figure out which regime our cluster is consistent with
			
		## A1 = 1 and R = 1, these clusters are consistent with 1,1 and 1,-1
		if(ind_vec[1]==1 && ind_vec[2] == 1){
			
				## We calculate the cluster's contribution makes to hat_U (defined in paper)
				## due to consistency with DTR 1,1
				
				## calculate V^-1 in estimating equation
				cur_cor = r[1]
				cur_sigma = sigma[1]
				outer_term = 1/(cur_sigma*(1+(cur_m-2)*cur_cor + (1-cur_m)*cur_cor^2))
				diag_part = ((cur_m-2)*cur_cor + 1)*outer_term
				off_diag = -cur_cor*outer_term 
				v_inv = matrix(off_diag, cur_m, cur_m)
				diag(v_inv) = diag_part			
				
				## calculate D in estimating equation	
				cur_d = cbind(rep(1,cur_m), rep(1,cur_m), rep(1,cur_m), rep(1,cur_m), cov_vec)
				
				## calculate residuals
				cur_resid = cur_resp - cur_d%*%mean_val
				
				## calculate hat_U contribution
				cur_u = 2*t(cur_d)%*%v_inv%*%cur_resid
				
				## We calculate the cluster's contribution makes to hat_U (defined in paper)
				## due to consistency with DTR 1,-1
				
				## calculate V^-1 in estimating equation
				cur_cor = r[2]
				cur_sigma = sigma[2]
				outer_term = 1/(cur_sigma*(1+(cur_m-2)*cur_cor + (1-cur_m)*cur_cor^2))
				diag_part = ((cur_m-2)*cur_cor + 1)*outer_term
				off_diag = -cur_cor*outer_term 
				v_inv = matrix(off_diag, cur_m, cur_m)
				diag(v_inv) = diag_part			
				
				## calculate D in estimating equation	
				cur_d = cbind(rep(1,cur_m), rep(1,cur_m), rep(-1,cur_m), rep(-1,cur_m), cov_vec)
				
				## calculate residuals
				cur_resid = cur_resp - cur_d%*%mean_val
				
				## calculate hat_U contribution
				cur_u = cur_u + 2*t(cur_d)%*%v_inv%*%cur_resid
				
				## add in this hat_U's contribution to the meat
				meat = meat + cur_u%*%t(cur_u)
		}
			
		## A1 = 1 and R = 0
		if(ind_vec[1]==1 && ind_vec[2] == 0){
			
			## A2 = 1, these clusters are consistent with 1,1 only
			if(ind_vec[3] == 1){
				
				## We calculate the cluster's contribution makes to hat_U (defined in paper)
				## due to consistency with DTR 1,1
				
				## calculate V^-1 in estimating equation
				cur_cor = r[1]
				cur_sigma = sigma[1]
				outer_term = 1/(cur_sigma*(1+(cur_m-2)*cur_cor + (1-cur_m)*cur_cor^2))
				diag_part = ((cur_m-2)*cur_cor + 1)*outer_term
				off_diag = -cur_cor*outer_term 
				v_inv = matrix(off_diag, cur_m, cur_m)
				diag(v_inv) = diag_part				
				
				## calculate D in estimating equation
				cur_d = cbind(rep(1,cur_m), rep(1,cur_m), rep(1,cur_m), rep(1,cur_m), cov_vec)
				
				## calculate residuals
				cur_resid = cur_resp - cur_d%*%mean_val
				
				## calculate hat_U contribution
				cur_u = 4*t(cur_d)%*%v_inv%*%cur_resid
				
				## add in this hat_U's contribution to the meat
				meat = meat + cur_u%*%t(cur_u)
				
			## A2 = -1, these clusters are consistent with 1,-1 only
			} else{
				
				## We calculate the cluster's contribution makes to hat_U (defined in paper)
				## due to consistency with DTR 1,-1
				
				## calculate V^-1 in estimating equation
				cur_cor = r[2]
				cur_sigma = sigma[2]
				outer_term = 1/(cur_sigma*(1+(cur_m-2)*cur_cor + (1-cur_m)*cur_cor^2))
				diag_part = ((cur_m-2)*cur_cor + 1)*outer_term
				off_diag = -cur_cor*outer_term 
				v_inv = matrix(off_diag, cur_m, cur_m)
				diag(v_inv) = diag_part				
				
				## calculate D in estimating equation
				cur_d = cbind(rep(1,cur_m), rep(1,cur_m), rep(-1,cur_m), rep(-1,cur_m), cov_vec)
				
				## calculate residuals
				cur_resid = cur_resp - cur_d%*%mean_val
				
				## calculate hat_U contribution
				cur_u = 4*t(cur_d)%*%v_inv%*%cur_resid
				
				## add in this hat_U's contribution to the meat
				meat = meat + cur_u%*%t(cur_u)
				}
		}
		
		## A1 = -1 and R = 1, these clusters are consistent with -1,1 and -1,-1
		if(ind_vec[1]== -1 && ind_vec[2] == 1){
			
				## We calculate the cluster's contribution makes to hat_U (defined in paper)
				## due to consistency with DTR -1,1
				
				## calculate V^-1 in estimating equation
				cur_cor = r[3]
				cur_sigma = sigma[3]
				outer_term = 1/(cur_sigma*(1+(cur_m-2)*cur_cor + (1-cur_m)*cur_cor^2))
				diag_part = ((cur_m-2)*cur_cor + 1)*outer_term
				off_diag = -cur_cor*outer_term 
				v_inv = matrix(off_diag, cur_m, cur_m)
				diag(v_inv) = diag_part			
				
				## calculate D in estimating equation	
				cur_d = cbind(rep(1,cur_m), rep(-1,cur_m), rep(1,cur_m), rep(-1,cur_m), cov_vec)
				
				## calculate residuals
				cur_resid = cur_resp - cur_d%*%mean_val
				
				## calculate hat_U contribution
				cur_u = 2*t(cur_d)%*%v_inv%*%cur_resid
				
				## We calculate the cluster's contribution makes to hat_U (defined in paper)
				## due to consistency with DTR -1,-1
				
				## calculate V^-1 in estimating equation
				cur_cor = r[4]
				cur_sigma = sigma[4]
				outer_term = 1/(cur_sigma*(1+(cur_m-2)*cur_cor + (1-cur_m)*cur_cor^2))
				diag_part = ((cur_m-2)*cur_cor + 1)*outer_term
				off_diag = -cur_cor*outer_term 
				v_inv = matrix(off_diag, cur_m, cur_m)
				diag(v_inv) = diag_part			
				
				## calculate D in estimating equation	
				cur_d = cbind(rep(1,cur_m), rep(-1,cur_m), rep(-1,cur_m), rep(1,cur_m), cov_vec)
				
				## calculate residuals
				cur_resid = cur_resp - cur_d%*%mean_val
				
				## calculate hat_U contribution
				cur_u = cur_u + 2*t(cur_d)%*%v_inv%*%cur_resid
				
				## add in this hat_U's contribution to the meat
				meat = meat + cur_u%*%t(cur_u)
		}
			
		## A1 = -1 and R = 0
		if(ind_vec[1]== -1 && ind_vec[2] == 0){
			
			## A2 = 1, these clusters are consistent with -1,1 only
			if(ind_vec[3] == 1){
				
				## We calculate the cluster's contribution makes to hat_U (defined in paper)
				## due to consistency with DTR -1,1
				
				## calculate V^-1 in estimating equation
				cur_cor = r[3]
				cur_sigma = sigma[3]
				outer_term = 1/(cur_sigma*(1+(cur_m-2)*cur_cor + (1-cur_m)*cur_cor^2))
				diag_part = ((cur_m-2)*cur_cor + 1)*outer_term
				off_diag = -cur_cor*outer_term 
				v_inv = matrix(off_diag, cur_m, cur_m)
				diag(v_inv) = diag_part				
				
				## calculate D in estimating equation
				cur_d = cbind(rep(1,cur_m), rep(-1,cur_m), rep(1,cur_m), rep(-1,cur_m), cov_vec)
				
				## calculate residuals
				cur_resid = cur_resp - cur_d%*%mean_val
				
				## calculate hat_U contribution
				cur_u = 4*t(cur_d)%*%v_inv%*%cur_resid
				
				## add in this hat_U's contribution to the meat
				meat = meat + cur_u%*%t(cur_u)
				
			## A2 = -1, these clusters are consistent with -1,-1 only
			} else{
				
				## We calculate the cluster's contribution makes to hat_U (defined in paper)
				## due to consistency with DTR -1,-1
				
				## calculate V^-1 in estimating equation
				cur_cor = r[4]
				cur_sigma = sigma[4]
				outer_term = 1/(cur_sigma*(1+(cur_m-2)*cur_cor + (1-cur_m)*cur_cor^2))
				diag_part = ((cur_m-2)*cur_cor + 1)*outer_term
				off_diag = -cur_cor*outer_term 
				v_inv = matrix(off_diag, cur_m, cur_m)
				diag(v_inv) = diag_part				
				
				## calculate D in estimating equation
				cur_d = cbind(rep(1,cur_m), rep(-1,cur_m), rep(-1,cur_m), rep(1,cur_m), cov_vec)
				
				## calculate residuals
				cur_resid = cur_resp - cur_d%*%mean_val
				
				## calculate hat_U contribution
				cur_u = 4*t(cur_d)%*%v_inv%*%cur_resid
				
				## add in this hat_U's contribution to the meat
				meat = meat + cur_u%*%t(cur_u)
				}
		}	
	}
	
	## Now with the the meat and bread, we calculate the variance covariace matrix
	meat = meat/N
	bread = solve(bread_inv)
	var_cov_mat = bread%*%meat%*%bread/N
		
			
	## Now we make the output beautiful.  We end up returning 
	## 4 things: the coefficient estimates-coefficients, their var_cov 
	## matrix-var_cov_coefficients, the mean under each regime estimates-DTR_mean_estimates, 
	## and the var_cov matrix of those estimates-var_cov_DTR_mean.  We obtain the 
	## mean under each regime using a transformation matrix
	
	eta_names = NULL
	if(d>0){for(i in 1:d){eta_names[i] = paste("eta", i, sep="")}}
	row.names(mean_val) = NULL
	mean_val_return = data.frame(mean_val)
	names(mean_val_return) = "estimates"
	row.names(mean_val_return) = c("beta_0", "beta_1", "beta_2", "beta_3", eta_names)
	
	row.names(var_cov_mat) = NULL
	var_cov_mat_return = data.frame(var_cov_mat)
	names(var_cov_mat_return) = c("beta_0", "beta_1", "beta_2", "beta_3", eta_names)
	row.names(var_cov_mat_return) = c("beta_0", "beta_1", "beta_2", "beta_3", eta_names)
	
	trans_mat = matrix(c(1,1,1,1,1,1,-1,-1,1,-1,1,-1,1,-1,-1,1), nrow=4, ncol=4)
		
	dtr_mean = data.frame(trans_mat%*%mean_val[1:4,])
	names(dtr_mean) = "DTR mean estimates"
	row.names(dtr_mean) = c("(1,1)", "(1,-1)", "(-1,1)", "(-1,-1)")
	
	dtr_mean_var_cov_mat = data.frame(trans_mat%*%(var_cov_mat[1:4,1:4])%*%(trans_mat))
	names(dtr_mean_var_cov_mat) = c("(1,1)", "(1,-1)", "(-1,1)", "(-1,-1)")
	row.names(dtr_mean_var_cov_mat) = c("(1,1)", "(1,-1)", "(-1,1)", "(-1,-1)")

	return_list = list(mean_val_return, var_cov_mat_return, dtr_mean, dtr_mean_var_cov_mat)
	names(return_list) = c("coefficients", "var_cov_coefficients", "DTR_mean_estimates", "var_cov_DTR_mean")
	return(return_list)  
}