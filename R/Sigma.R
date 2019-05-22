
Sigma <- function(M.new.j, M.new.k,
                  nomissing.index.j, nomissing.index.k ) {
  # Estimate covariance matrix
  #
  # Args:
  #   mediator.new.j      -->   mean centered mediator at time j
  #   mediator.new.k      -->   mean centered mediator at time k of j
  #   nomissing.index.j   -->   index of missing values within mediators j
  #   nomissing.index.k   -->   index of missing values within mediators k
  #
  # Returns:
  #   Sigma               -->
  #
  ##

  ind.j <- which(nomissing.index.j == 1)
  ind.k <- which(nomissing.index.k == 1)
  check.index <- na.omit(cbind(1:length(ind.j), match(ind.j, (ind.k))))
  Mat <- matrix(0, nrow = length(ind.j), ncol = length(ind.k))
  Mat[check.index] <- 1
  temp.matrix <- solve(t(M.new.j) %*% (M.new.j)) %*%
    t(M.new.j) %*% Mat %*% (M.new.k) %*%
    solve(t(M.new.k) %*% (M.new.k))

  return(temp.matrix)

}
