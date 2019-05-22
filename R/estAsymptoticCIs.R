
estAsymptoticCIs <- function(trt, t.seq, M, Y, t.est, t.coeff, est.smooth, deltat) {
  # Estimate asymptotic confidence intervals and standard deviation
  #
  # Args:
  #   trt         -->   a vector with treatment values
  #   t.seq       -->   a vector of time points at each obs
  #   M           -->   matrix of mediator values
  #   Y           -->   matrix of outcome outcomes
  #   t.est       -->   time points to make the estimation
  #   t.coeff     -->   estimated coefficients
  #   est.smooth  -->   smoothed coefficients
  #   deltat      -->   half the time between two time points
  #
  # Returns:
  #   hat.alpha.1 -->  estiamted treatment effect on mediator
  #   hat.beta.2  -->  estiamted mediation effect on outcome
  #   sd          -->  estimated standard deviation of the mediation effect
  #
  ##

  start.time <- Sys.time()
  print("Beginning asymptotic process.")
  temp.col <- NULL
  for (j in 2:length(t.seq)) {
    temp.row    <- NULL
    for (k in 2:length(t.seq)) {
      # create empty vector, store raw mediator.
      temp.M.j <- NULL
      temp.M.k <- NULL
      temp.M.j <- cbind(M[j - 1, ], M[j, ])
      temp.M.k <- cbind(M[k - 1, ], M[k, ])
      # Derive centered Mediators and Outcomes
      newMO.j.est <- newMediatorOutcome(trt, temp.M.j, Y[j - 1, ])
      newMO.k.est <- newMediatorOutcome(trt, temp.M.k, Y[k - 1, ])
      # compute matrix for multiplying with covariance matrix
      temp.sigma <- Sigma(newMO.j.est$M, newMO.k.est$M,
                          newMO.j.est$nomissingIndex, newMO.k.est$nomissingIndex)
      # Estimate big covariance matirx
      if (j < length(t.seq) || k < length(t.seq)) {
        cat(sprintf("estimating error covariance between time %.1f and %.1f", t.seq[j], t.seq[k]), "\r")
      } else {
        print(sprintf("estimating error covariance between time %.1f and %.1f", t.seq[j], t.seq[k]))
      }
      est.covar.error <- estErrorCovar(newMO.j.est$M, newMO.k.est$M,
                                       newMO.j.est$Y, newMO.k.est$Y,
                                       newMO.j.est$nomissingIndex, newMO.k.est$nomissingIndex)
      # multiply sigma and covar matrix, store results.
      results <- est.covar.error * temp.sigma
      temp.row <- cbind(temp.row, results)
    }
    temp.col <- rbind(temp.col, temp.row)
  }

  # DERIVE SD's
  sd = NULL
  for (j in 1:length(t.est)) {
    t <- t.est[j]
    dydx <- DF(t.seq, t, deltat, t.coeff, est.smooth$bw_alpha1, est.smooth$bw_beta2)
    #sd.temp <- sqrt(t(dydx) %*% temp.col %*% dydx)  # GENERATES NaNs B/C THE PRODUCT CAN BE LESS THAN ZERO
    sd.temp <- t(dydx) %*% temp.col %*% dydx
    sd <- cbind(sd, sd.temp)
  }

  # DERIVE CI's
  CI.upper <- est.smooth$true.M + (1.96 * sd)
  CI.lower <- est.smooth$true.M - (1.96 * sd)

  # GENERATE OUTPUT
  result <- list(sd = sd, CI.upper = CI.upper, CI.lower = CI.lower)

  end.time <- Sys.time()
  total.time <- end.time - start.time
  print(sprintf("Process complete. Elapsed time = %.3f secs", as.numeric(total.time, units = "secs")))
  return(result)

}
