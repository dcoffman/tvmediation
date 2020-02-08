
smoothest <- function(t.seq, t.coeff, t.est, deltat) {
  # Compute local polynomial estimation using rule of thumb for bandwidth selection
  #
  # Args:
  #   t.seq     -->   a vector of time points at each obs
  #   t.coeff   -->   estimated coefficients
  #   t.est     -->   time points to make the estimation
  #   deltat    -->   half the time between two time points
  #
  # Returns:
  #  smoothest  --> list containing ********
  #
  ##
  
  # Equations 4 & 5
  bw_alpha1 <- thumbBw(t.seq[-1], t.coeff[1, ], deg = 1, kernel = gaussK)
  bw_beta2 <- thumbBw(t.seq[-1], t.coeff[3, ], deg = 1, kernel = gaussK)
  hat.alpha.1 = locPolSmootherC(t.seq[-1], t.coeff[1, ], t.est - deltat,
                                bw_alpha1, deg = 1, kernel = gaussK)$beta0
  hat.beta.2 = locPolSmootherC(t.seq[-1], t.coeff[3, ], t.est, bw_beta2,
                               deg = 1, kernel = gaussK)$beta0
  
  return(list(bw_alpha1 = bw_alpha1, bw_beta2 = bw_beta2,
              hat.alpha.1 = hat.alpha.1, hat.beta.2 = hat.beta.2,
              est.M = hat.alpha.1*hat.beta.2))
}
