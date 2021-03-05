#' Bootstrapping samples to estimate mediation effects confidence intervals for continuous 
#' outcome and two treatment (exposure) groups.
#' 
#' Part of the set of internal functions for estimating bootstrapped confidence intervals 
#' for continuous outcome and two treatment groups when user argument \code{CI="boot"}.
#' 
#' @param trt         a vector indicating treatment group
#' @param t.seq       a vector of time points at each obs
#' @param M           matrix of mediator values
#' @param Y           matrix of outcome values
#' @param t.est       time points at which to make the estimation
#' @param deltat      a small constant which controls the time-lag of the effect of the 
#'                    mediator on the outcome, half the time between two time points
#' @param replicates  number of replicates for bootstrapping confidence intervals.
#' 
#' @return {boot.se}{bootstrapped standard error for the estimated mediation effect}
#' @return {CI.upper}{percentile bootstrapped CI upper limit for the estimated mediation effect}
#' @return {CI.lower}{percentile bootstrapped CI lower limit for the estimated mediation effect}
#' 
#' 

estBootCIs <- function(trt, t.seq, M, Y, t.est, deltat, replicates) {
 
  start.time <- Sys.time()
  print("Beginning bootstrap for mediation effect CIs.")
  N <- length(trt)
  storage.boot <- matrix(0, nrow = replicates, ncol = length(t.est))
  
  for(k1 in 1:replicates) {
    if (k1 < replicates) {
      cat(sprintf("Boostrapping iteration %03d", k1), " \r")
    } else {
      print(sprintf("Boostrapping iteration %03d", k1))
    }
    index.sample <- sample(1:N, N, replace = TRUE)
    t.coeff.boot <- NULL
    for (j in 2:length(t.seq)) {
      # create empty vector, store raw mediator.
      temp.M.boot <- NULL
      temp.Y.boot <- NULL
      temp.M.boot  <- cbind(M[j - 1, index.sample], M[j, index.sample])
      temp.Y.boot <- Y[, index.sample]
      temp.trt.boot <- trt[index.sample]
      # Derive centered Mediators and Outcomes
      newMO.j.est.boot <- newMediatorOutcome(temp.trt.boot, temp.M.boot, temp.Y.boot[j-1,])
      # Estimate coefficients, then store them.
      coeff.est.boot <- estCoeff(newMO.j.est.boot)
      t.coeff.boot <- cbind(t.coeff.boot, coeff.est.boot)  # store coeff estimates at t.seq
    }
    # Equations 4 & 5
    bw_alpha <- locpol::thumbBw(t.seq[-1], t.coeff.boot[1, ], deg = 1, kernel = locpol::gaussK)
    bw_beta  <- locpol::thumbBw(t.seq[-1], t.coeff.boot[3, ], deg = 1, kernel = locpol::gaussK)
    hat.alpha <- locpol::locPolSmootherC(t.seq[-1], t.coeff.boot[1,], t.est - deltat, bw_alpha,
                                        deg = 1, kernel = locpol::gaussK)$beta0
    hat.beta  <- locpol::locPolSmootherC(t.seq[-1], t.coeff.boot[3,], t.est, bw_beta, 
                                        deg = 1, kernel = locpol::gaussK)$beta0
    est.smooth.boot <- list(bw_alpha = bw_alpha, bw_beta = bw_beta,
                            hat.alpha = hat.alpha, hat.beta = hat.beta,
                            est.M = hat.alpha*hat.beta)
    storage.boot[k1,] <- est.smooth.boot$est.M
  }
  
  # COMPUTE 2.5% AND 97.5% QUANTILES FOR EACH COLUMN
  CI.all <- apply(storage.boot, 2, quantile, probs = c(0.025, 0.975))
  
  # GENERATE OUTPUT
  results <- list(boot.se=apply(storage.boot, 2, sd), CI.upper = CI.all[2, ], 
                  CI.lower = CI.all[1, ])
  
  end.time <- Sys.time()
  total.time <- end.time - start.time
  print(sprintf("Process complete. Elapsed time = %.3f secs", 
                as.numeric(total.time, units = "secs")))
  return(results)
  
}
