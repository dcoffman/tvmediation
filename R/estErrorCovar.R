
estErrorCovar <- function(M.new.j, M.new.k,
                          Y.new.j, Y.new.k,
                          nomissing.index.j, nomissing.index.k) {
  # Estimate covariance of the error term at time t_j and t_k for every j & k when j != k
  #
  #
  # Args:
  #   mediator.new.j      -->   mean centered mediator at time j
  #   mediator.new.k      -->   mean centered mediator at time k of j
  #   outcome.new.j       -->   mean centered outcomes at time j
  #   outcome.new.k       -->   mean centered outcomes at time k of j
  #   nomissing.index.j   -->   index of missing values within mediators j
  #   nomissing.index.k   -->   index of missing values within mediators k
  #
  # Returns:
  #   est.error.covar     -->   matrix of covariance estimates
  #
  ##

  N.j = length(Y.new.j)
  N.k = length(Y.new.k)

  Pj <- (M.new.j) %*% solve(t(M.new.j) %*% (M.new.j)) %*% t(M.new.j)
  Pk <- (M.new.k) %*% solve(t(M.new.k) %*% (M.new.k)) %*% t(M.new.k)

  I.minus.Pj <- (diag(1, N.j) - Pj)
  I.minus.Pk <- (diag(1, N.k) - Pk)
  ej = I.minus.Pj %*% Y.new.j
  ek = I.minus.Pk %*% Y.new.k

  ind.j = which(nomissing.index.j == 1)
  ind.k = which(nomissing.index.k == 1)
  check.index <- na.omit(cbind(1:length(ind.j), match(ind.j, (ind.k))))
  Mat = matrix(0, nrow = length(ind.j), ncol = length(ind.k))
  Mat[check.index] = 1

  result <- sum(diag(ej %*% t(ek))) / sum(diag(I.minus.Pj %*% (Mat) %*% t(I.minus.Pk)))
  return(result)

}
