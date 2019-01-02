# VCF_cost.R
# VCF cost functions for likelihood
# Author: Cait Harrigan

#Gaussian likelihood maximization 
gaussian_mll <- function(phis, mean_square_phis){
  # mean_phis read from counts file
  # mean_square_phis read from quadratic phis file
  # Score a segment using likelihood under normal
  n <- length(phis)
  LL <- sum( (n/2)*log(2*pi), (n/2)*log(mean_square_phis), 
             1/(2*mean_square_phis) * sum( (phis - mean(phis))^2 ) )
  return(LL)
}

# [END]