##########################################################################
#### Bootstrap for computing the CI for coefficients alpha1 and beta2 ####
##########################################################################

bootci_coeff_2trt <- function(trt, t.seq, M, Y, t.est, deltat, replicates) {
  
  start.time <- Sys.time()
  print("Beginning bootstrap for coefficient CIs.")
  N <- length(trt)
  storage.boot.1 <- matrix(0, nrow = replicates, ncol = length(t.est))
  storage.boot.2 <- matrix(0, nrow = replicates, ncol = length(t.est))
  
  for(c1 in 1:replicates) {
    if (c1 < replicates) {
      cat(sprintf("Boostrapping iteration %03d", c1), " \r")
    } else {
      print(sprintf("Boostrapping iteration %03d", c1))
    }
    set.seed(c1^2)
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
    storage.boot.1[c1, ] <- est.smooth.boot$hat.alpha.1
    storage.boot.2[c1, ] <- est.smooth.boot$hat.beta.2
  }
  
  # COMPUTE 2.5% AND 97.5% QUANTILES FOR EACH COLUMN
  CI.alpha1 <- apply(storage.boot.1, 2, quantile, probs = c(0.025, 0.975))
  CI.beta2 <- apply(storage.boot.2, 2, quantile, probs = c(0.025, 0.975))
  
  # GENERATE OUTPUT
  results <- list(CI.upper.alpha = CI.alpha1[2, ], CI.lower.alpha = CI.alpha1[1, ],
                  CI.upper.beta = CI.beta2[2, ], CI.lower.beta = CI.beta2[1, ])
  
  end.time <- Sys.time()
  total.time <- end.time - start.time
  print(sprintf("Process complete. Elapsed time = %.3f secs", as.numeric(total.time, units = "secs")))
  return(results)
}
