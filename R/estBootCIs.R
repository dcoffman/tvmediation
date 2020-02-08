
estBootCIs <- function(trt, t.seq, M, Y, t.est, deltat, replicates) {
  # Estimate bootstrap confidence intervals
  #
  # Args:
  #   trt         -->   a vector with treatment values
  #   t.seq       -->   a vector of time points at each obs
  #   M           -->   matrix of mediator values
  #   Y           -->   matrix of outcome outcomes
  #   t.est       -->   time points to make the estimation
  #   deltat      -->   half the time between two time points
  #   replicates  -->   Number of replicates for bootstrapping confidence intervals.
  #
  # Returns:
  #   boot.se     -->  bootstrapped standard error for mediation effect
  #   CI.upper    -->  upper limit of percentile bootstrapped CI 
  #   CI.lower    -->  lower limit of percentile bootstrapped CI 
  #
  ##
  
  start.time <- Sys.time()
  print("Beginning bootstrap.")
  N <- length(trt)
  storage.boot <- matrix(0, nrow = replicates, ncol = length(t.est))
  
  for(k1 in 1:replicates) {
    if (k1 < replicates) {
      cat(sprintf("Boostrapping iteration %03d", k1), " \r")
    } else {
      print(sprintf("Boostrapping iteration %03d", k1))
    }
    set.seed(k1^2)
    index.sample <- sample(1:N, N, replace = TRUE)
    t.coeff.boot <- NULL
    for (j in 2:length(t.seq)) {
      # create empty vector, store raw mediator.
      temp.M.boot <- NULL
      temp.Y.boot <- NULL
      temp.M.boot  <- cbind(M[j - 1, index.sample], M[j, index.sample])
      temp.Y.boot <- Y[ , index.sample]
      temp.trt.boot <- trt[index.sample]
      # Derive centered Mediators and Outcomes
      newMO.j.est.boot <- newMediatorOutcome(temp.trt.boot, temp.M.boot, temp.Y.boot[j - 1, ])
      # Estimate coefficients, then store them.
      coeff.est.boot <- estCoeff(newMO.j.est.boot)
      t.coeff.boot <- cbind(t.coeff.boot, coeff.est.boot)  # store coeff estimates at t.seq
    }
    # Equations 4 & 5
    est.smooth.boot <- smoothest(t.seq, t.coeff.boot, t.est, deltat)
    storage.boot[k1, ] <- est.smooth.boot$est.M
  }
  
  # COMPUTE 2.5% AND 97.5% QUANTILES FOR EACH COLUMN
  CI.all <- apply(storage.boot, 2, quantile, probs = c(0.025, 0.975))
  
  # GENERATE OUTPUT
  results <- list(boot.se=apply(storage.boot, 2, sd), CI.upper = CI.all[2, ], CI.lower = CI.all[1, ])
  
  end.time <- Sys.time()
  total.time <- end.time - start.time
  print(sprintf("Process complete. Elapsed time = %.3f secs", as.numeric(total.time, units = "secs")))
  return(results)
  
}
