## Here is a function to calculate either the sample size (N) or power (beta)
## for an prototypical style cluster randomized SMART
## Formula 6 and the second Formula 7 in the paper

prototypical_sample_size = function(N = NULL, beta = NULL, alpha, m, delta, ICC, p_1, p_neg1, Cor_Y_X = NULL){
  ## N is the total number of clusters
  ## beta is the type 2 error (1 - power)
  ## alpha is the type 1 error
  ## delta is the standardized effect size
  ## ICC is the inter cluster correlation (rho with no covariates and rho* with covariates)
  ## p_1 is P(R = 1 | A_1 = 1)
  ## p_neg1 is P(R = 1 | A_1 = -1)
  ## Cor_Y_X is correlation between Y and the cluster level covariate X. It should only
  ## be defined when there is a covariate. If it is left undefined the non-covariate formula
  ## will be used
  
  ## To use the formula, either N or beta should be left blank, but not both.
  ## When N is blank, the sample size necessary to obtain that power will be returned
  ## When beta is blank, the power will be returned for a given sample size.
  
  z_alpha = qnorm(1-alpha/2)  ## Get the standard normal quanitle
  cov_string = "with a cluster level covariate"
  if(is.null(Cor_Y_X)){  ## determines if we wil use the formula with no covariates
    Cor_Y_X = 0
    cov_string = "with no covariates"
  }
  
  ## the case where we solve for N
  if(is.null(N)){
    z_beta = qnorm(1-beta)
    
    first_term = 4*(z_beta+z_alpha)^2/(m*delta^2)
    second_term = (1+(m-1)*ICC)
    third_term = 1+(1-p_1 + 1-p_neg1)/2
    
    # Sample size
    sample_size = ceiling(first_term*second_term*third_term*(1-Cor_Y_X^2))
    
    cat("The sample size required for a prototypical design", cov_string, "is", toString(sample_size))
    
    
    ## the case where we solve for the power
  } else if(is.null(beta)){
    second_term = (1+(m-1)*ICC)
    third_term = 1+(1-p_1 + 1-p_neg1)/2
    
    # The power value by back solving the equation
    power_val = pnorm(sqrt(N/second_term/third_term/(1-Cor_Y_X^2)/4*delta^2*m) - z_alpha)
    cat("The power for a prototypical design", cov_string, "is", power_val)
    
  } else {
    cat("Either beta or N must be left undefined")
  }
}

## Example Run
prototypical_sample_size(beta = .2, alpha = .05, m = 10, delta = .282, ICC = .01, p_1 = .2, p_neg1 = .3, Cor_Y_X = .2)