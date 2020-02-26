
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
  
  ### Current method: Ridge Regression
  sym_newMO <- t(newMO.j.est$M)%*%(newMO.j.est$M)
  diag(sym_newMO) <- diag(sym_newMO) + .001
  coeff.est <- solve(sym_newMO)%*%t(newMO.j.est$M)%*%(newMO.j.est$Y)
  
  ### Original way that rendered singular matrices which were non-invertible
  # coeff.est <- solve(t(newMO.j.est$M)%*%(newMO.j.est$M))%*%t(newMO.j.est$M)%*%(newMO.j.est$Y)
  
  ### Generalized Moore-Penrose Inverse method
  # coeff.est <- MASS::ginv(t(newMO.j.est$M) %*% (newMO.j.est$M)) %*% t(newMO.j.est$M) %*% (newMO.j.est$Y)
  
  return(coeff.est)
}
