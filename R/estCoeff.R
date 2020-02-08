
estCoeff <- function(newMO.j.est) {
  
  # Estimate coefficients.
  #
  # Args:
  #   newMO.j.est  -->   a list containing mean centered mediators ($M) and outcomes ($Y)
  #
  # Returns:
  #   coeff.est    -->   estimated coefficients
  #
  ##
  
  coeff.est <- MASS::ginv(t(newMO.j.est$M) %*% (newMO.j.est$M)) %*% t(newMO.j.est$M) %*% (newMO.j.est$Y)
  return(coeff.est)
}
