#' Bootstrap samples to estimate confidence intervals for coefficients for a continuous outcome and three treatment groups
#' 
#' Part of the set of internal functions for estimating bootstrapped confidence intervals for the coefficients of the mediation model for a continuous outcome and three treatment groups
#' 
#' @param T1                a vector indicating assignment to treatment 1
#' @param T2                a vector indicating assignment to treatment 2
#' @param t.seq             a vector of time points for each observation
#' @param mediator          matrix of mediator values in wide format
#' @param outcome           matrix of outcome values in wide format
#' @param t.est             time points at which to make the estimation. Default = t.seq.
#' @param original.coeff    a list of the estimated coefficients.
#' @param boot.sample       number of replicates for bootstrapping confidence intervals. Default = 1000.
#' 
#' @return \item{alw1}{CI lower limit for estimated Treatment 1 effect on mediator}
#' @return \item{aup1}{CI upper limit for estimated Treatment 1 effect on mediator}
#' @return \item{alw2}{CI lower limit for estimated Treatment 2 effect on mediator}
#' @return \item{aup2}{CI upper limit for estimated Treatment 2 effect on mediator}
#' @return \item{glw1}{CI lower limit for estimated Treatment 1 direct effect on outcome}
#' @return \item{gup1}{CI upper limit for estimated Treatment 1 direct effect on outcome}
#' @return \item{glw2}{CI lower limit for estimated Treatment 2 direct effect on outcome}
#' @return \item{gup2}{CI upper limit for estimated Treatment 2 direct effect on outcome}
#' @return \item{tlw1}{CI lower limit for estimated Treatment 1 total effect on outcome}
#' @return \item{tup1}{CI upper limit for estimated Treatment 1 total effect on outcome}
#' @return \item{tlw2}{CI lower limit for estimated Treatment 2 total effect on outcome}
#' @return \item{tup2}{CI upper limit for estimated Treatment 2 total effect on outcome}
#' @return \item{blw}{CI lower limit for estimated effect of mediator on outcome}
#' @return \item{bup}{CI upper limit for estimated effect of mediator on outcome}
#' 
#' @export
#' 

bootci_coeff_3trt <- function(T1, T2, t.seq, mediator, outcome, t.est, original.coeff, boot.sample = 1000)
{
  original.alpha1 <- (original.coeff$hat.alpha1)
  original.alpha2 <- (original.coeff$hat.alpha2)
  original.gamma1 <- (original.coeff$hat.gamma1)
  original.gamma2 <- (original.coeff$hat.gamma2)
  original.tau1 <- (original.coeff$hat.tau1)
  original.tau2 <- (original.coeff$hat.tau2)
  original.beta <- (original.coeff$hat.beta)
  
  N <- length(T1)
  x <- mediator
  y <- outcome
  
  storage.boot1 <- matrix(0, nrow=boot.sample, ncol=length(t.est))
  storage.boot2 <- matrix(0, nrow=boot.sample, ncol=length(t.est))
  storage.boot3 <- matrix(0, nrow=boot.sample, ncol=length(t.est))
  storage.boot4 <- matrix(0, nrow=boot.sample, ncol=length(t.est))
  storage.boot5 <- matrix(0, nrow=boot.sample, ncol=length(t.est))
  storage.boot6 <- matrix(0, nrow=boot.sample, ncol=length(t.est))
  storage.boot7 <- matrix(0, nrow=boot.sample, ncol=length(t.est))
  
  start.time <- Sys.time()
  print(paste("Beginning bootstrap for coefficient CIs."))
  
  for(c1 in 1:boot.sample)
  {
    set.seed(c1^2)
    if (c1 < boot.sample) {
      cat(sprintf("Boostrapping iteration %03d", c1), " \r")
    } else {
      print(sprintf("Boostrapping iteration %03d", c1))
    }
    
    index.sample <- sample(1:N, N, replace=TRUE)
    x.boot <- x[,index.sample]
    y.boot <- y[,index.sample]
    T1.boot <- T1[index.sample]
    T2.boot <- T2[index.sample]
    deltat.boot <- max(diff(t.seq))/2
    
    result.boot <- tvmcurve_3trt(T1.boot, T2.boot, t.seq, x.boot, y.boot, t.est)
    storage.boot1[c1,] <- result.boot$hat.alpha1
    storage.boot2[c1,] <- result.boot$hat.alpha2
    storage.boot3[c1,] <- result.boot$hat.gamma1
    storage.boot4[c1,] <- result.boot$hat.gamma2
    storage.boot5[c1,] <- result.boot$hat.tau1
    storage.boot6[c1,] <- result.boot$hat.tau2
    storage.boot7[c1,] <- result.boot$hat.beta
  }
  
  orig.se1.all <- apply(storage.boot1, 2, sd)
  orig.se2.all <- apply(storage.boot2, 2, sd)
  orig.se3.all <- apply(storage.boot3, 2, sd)
  orig.se4.all <- apply(storage.boot4, 2, sd)
  orig.se5.all <- apply(storage.boot5, 2, sd)
  orig.se6.all <- apply(storage.boot6, 2, sd)
  orig.se7.all <- apply(storage.boot7, 2, sd)
  
  t.length <- dim(storage.boot1)[2]
  alpha.lower.1 <- rep(0,t.length)
  alpha.upper.1 <- rep(0,t.length)
  alpha.lower.2 <- rep(0,t.length)
  alpha.upper.2 <- rep(0,t.length)
  gamma.lower.1 <- rep(0,t.length)
  gamma.upper.1 <- rep(0,t.length)
  gamma.lower.2 <- rep(0,t.length)
  gamma.upper.2 <- rep(0,t.length)
  tau.lower.1 <- rep(0,t.length)
  tau.upper.1 <- rep(0,t.length)
  tau.lower.2 <- rep(0,t.length)
  tau.upper.2 <- rep(0,t.length)
  beta.lower <- rep(0,t.length)
  beta.upper <- rep(0,t.length)
  
  for(c2 in 1:(dim(storage.boot1)[2])) # run through different time points 
  {
    temp1 <- storage.boot1[,c2]
    temp2 <- storage.boot2[,c2]
    temp3 <- storage.boot3[,c2]
    temp4 <- storage.boot4[,c2]
    temp5 <- storage.boot5[,c2]
    temp6 <- storage.boot6[,c2]
    temp7 <- storage.boot7[,c2]
    
    orig.est1 <- original.alpha1[c2]
    orig.est2 <- original.alpha2[c2]
    orig.est3 <- original.gamma1[c2]
    orig.est4 <- original.gamma2[c2]
    orig.est5 <- original.tau1[c2]
    orig.est6 <- original.tau2[c2]
    orig.est7 <- original.beta[c2]
    
    ###percentile bootstrap method 
    alpha.lower.1[c2] <- quantile(temp1, 0.975)
    alpha.upper.1[c2] <- quantile(temp1, 0.025)
    
    alpha.lower.2[c2] <- quantile(temp2, 0.975)
    alpha.upper.2[c2] <- quantile(temp2, 0.025)
    
    gamma.lower.1[c2] <- quantile(temp3, 0.975)
    gamma.upper.1[c2] <- quantile(temp3, 0.025)
    
    gamma.lower.2[c2] <- quantile(temp4, 0.975)
    gamma.upper.2[c2] <- quantile(temp4, 0.025)
    
    tau.lower.1[c2] <- quantile(temp5, 0.975)
    tau.upper.1[c2] <- quantile(temp5, 0.025)
    
    tau.lower.2[c2] <- quantile(temp6, 0.975)
    tau.upper.2[c2] <- quantile(temp6, 0.025)
    
    beta.lower[c2] <- quantile(temp7, 0.975)
    beta.upper[c2] <- quantile(temp7, 0.025)
    
  }
  
  coeff_CI <- list(alw1 = alpha.lower.1, aup1 = alpha.upper.1, alw2 = alpha.lower.2, aup2 = alpha.upper.2,
                   glw1 = gamma.lower.1, gup1 = gamma.upper.1, glw2 = gamma.lower.2, gup2 = gamma.upper.2,
                   tlw1 = tau.lower.1, tup1 = tau.upper.1, tlw2 = tau.lower.2, tup2 = tau.upper.2,
                   blw = beta.lower, bup = beta.upper)
  
  end.time <- Sys.time()
  total.time <- end.time - start.time
  print(sprintf("Process complete. Elapsed time = %.3f secs", 
                as.numeric(total.time, units = "secs")))
  
  return(coeff_CI)
  
}

