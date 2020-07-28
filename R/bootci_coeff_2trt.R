#' Bootstrapping samples to estimate confidence intervals for coefficients of continuous outcome and two treatment (exposure) groups
#' 
#' Part of the set of internal functions for estimating bootstrapped confidence intervals for the coefficients of the mediation model for continuous outcome and two treatment (exposure) groups.
#' 
#' @param trt         a vector with treatment values
#' @param t.seq       a vector of time points for each observation
#' @param M           matrix of mediator values in wide format
#' @param Y           matrix of outcome outcomes in wide format
#' @param t.est       time points to make the estimation                              Default = t.seq
#' @param deltat      a small constant which controls the time-lag of the effect of the mediator on the outcome.
#' @param replicates  number of replicates for bootstrapping confidence intervals.    Default = 1000
#' 
#' @return \item{CI.upper.alpha}{upper limit of confidence intervals for coefficient hat.alpha}
#' @return \item{CI.lower.alpha}{lower limit of confidence intervals for coefficient hat.alpha}
#' @return \item{CI.upper.gamma}{upper limit of confidence intervals for coefficient hat.gamma}
#' @return \item{CI.lower.gamma}{lower limit of confidence intervals for coefficient hat.gamma}
#' @return \item{CI.upper.beta}{upper limit of confidence intervals for coefficient hat.beta}
#' @return \item{CI.lower.beta}{lower limit of confidence intervals for coefficient hat.beta}
#' @return \item{CI.upper.tao}{upper limit of confidence intervals for coefficient hat.tao}
#' @return \item{CI.lower.tao}{lower limit of confidence intervals for coefficient hat.tao}
#' 
#' @export
#' 

##########################################################################
#### Bootstrap for computing the CI for coefficients alpha1 and beta2 ####
##########################################################################

bootci_coeff_2trt <- function(trt, t.seq, M, Y, t.est, deltat, replicates) {
  
  start.time <- Sys.time()
  print("Beginning bootstrap for coefficient CIs.")
  N <- length(trt)
  storage.boot.1 <- matrix(0, nrow = replicates, ncol = length(t.est))
  storage.boot.2 <- matrix(0, nrow = replicates, ncol = length(t.est))
  storage.boot.3 <- matrix(0, nrow = replicates, ncol = length(t.est))
  storage.boot.4 <- matrix(0, nrow = replicates, ncol = length(t.est))
  
  
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
      # Estimate coefficients alpha, beta and gamma
      coeff.est.boot <- estCoeff(newMO.j.est.boot)
      
      # Steps to derive the total effect coefficient tao
      temp.X.new <- scale(temp.trt.boot, center = TRUE, scale = FALSE)
      temp.Y.new <- scale(temp.Y.boot[j - 1, ], center = TRUE, scale = FALSE)
      nomissing.X.new <- complete.cases(temp.X.new)
      nomissing.Y.new <- complete.cases(temp.Y.new)
      nomissing.index.new <- nomissing.X.new * nomissing.Y.new
      
      temp.X.new <- temp.X.new[which(nomissing.index.new == 1),]
      temp.Y.new <- temp.Y.new[which(nomissing.index.new == 1)]
      temp.sym_newMO <- t(temp.X.new)%*%(temp.X.new)
      temp.coeff.tao <- solve(temp.sym_newMO)%*%t(temp.X.new)%*%(temp.Y.new)
      
      # Store the coefficients
      t.coeff.all <- rbind(coeff.est.boot, temp.coeff.tao)
      t.coeff.boot <- cbind(t.coeff.boot, t.coeff.all)  # store coeff estimates at t.seq
    }
    # Equations 4 & 5
    est.smooth.boot <- smoothest(t.seq, t.coeff.boot, t.est, deltat)
    storage.boot.1[c1, ] <- est.smooth.boot$hat.alpha
    storage.boot.2[c1, ] <- est.smooth.boot$hat.gamma
    storage.boot.3[c1, ] <- est.smooth.boot$hat.beta
    storage.boot.4[c1, ] <- est.smooth.boot$hat.tao
    
  }
  
  # COMPUTE 2.5% AND 97.5% QUANTILES FOR EACH COLUMN
  CI.alpha <- apply(storage.boot.1, 2, quantile, probs = c(0.025, 0.975))
  CI.gamma <- apply(storage.boot.2, 2, quantile, probs = c(0.025, 0.975))
  CI.beta <- apply(storage.boot.3, 2, quantile, probs = c(0.025, 0.975))
  CI.tao <- apply(storage.boot.4, 2, quantile, probs = c(0.025, 0.975))
  
  # GENERATE OUTPUT
  results <- list(CI.upper.alpha = CI.alpha[2, ], CI.lower.alpha = CI.alpha[1, ],
                  CI.upper.gamma = CI.gamma[2, ], CI.lower.gamma = CI.gamma[1, ],
                  CI.upper.beta = CI.beta[2, ], CI.lower.beta = CI.beta[1, ],
                  CI.upper.tao = CI.tao[2, ], CI.lower.tao = CI.tao[1, ])
  
  end.time <- Sys.time()
  total.time <- end.time - start.time
  print(sprintf("Process complete. Elapsed time = %.3f secs", as.numeric(total.time, units = "secs")))
  return(results)
}
