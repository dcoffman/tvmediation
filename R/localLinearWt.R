
localLinearWt <- function(bw, t, t.seq, j) {
  # Return the weight of the local polynomial regression
  #
  # Args:
  #   bw            -->   smoothed alpha or beta
  #   t             -->   t.est at time j or t-deltat
  #   t.seq         -->   a vector of time points at each obs
  #   j             -->   counter
  #
  # Returns:
  #   localLinearWt -->   weights for the second step of the estimation procedure
  #
  ##
  
  # Second step of the estimation procedure, find the weight
  C <- matrix(1, nrow = 2, ncol = length(t.seq) - 1)
  C[2,] <- t - t.seq[-1]
  
  W <- matrix(0, nrow = length(t.seq) - 1 , ncol = length(t.seq) - 1)
  diag(W) <- kern(t, bw, t.seq[-1])
  
  weight <- solve((C) %*% W %*% t(C)) %*% C[ ,j] %*% diag(W)[j]
  return(weight[1])
  
}
