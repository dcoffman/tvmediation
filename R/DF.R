
DF <- function(t.seq, t, deltat, temp.coeff, bwa, bwb) {
  # Calculates derivative for asymptotic standard deviations
  #
  # Args:
  #   t.seq       -->   a vector of time points at for each obs
  #   t           -->   t.est at time j
  #   deltat      -->   half the time between two time points
  #   temp.coeff  -->   t.coeff
  #   bwa         -->   smoothed alpha
  #   bwb         -->   smoothed beta
  #
  # Returns:
  #   DF          -->
  #
  ##

  dydx <- NULL
  hat.alpha.1 <- locPolSmootherC(t.seq[-1], temp.coeff[1, ], t - deltat,
                                 bwa, deg = 1, kernel = gaussK)$beta0
  hat.beta.2 <- locPolSmootherC(t.seq[-1], temp.coeff[3, ], t,
                                bwb, deg = 1, kernel = gaussK)$beta0
  common.factor <- c(hat.beta.2, 0, hat.alpha.1)
  #browser()
  for (m in 1:(length(t.seq) - 1)) {
   weight1 <- localLinearWt(bwa, t-deltat, t.seq, m)
   weight2 <- localLinearWt(bwb, t, t.seq, m)
   chuck.three <- common.factor * c(weight1, 0, weight2)
   dydx <- c(dydx, chuck.three)
  }
  return(dydx)

}
