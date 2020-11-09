#' Bootstrap samples to estimate confidence intervals for the mediation effect for a binary outcome
#' 
#' Part of the set of internal functions for estimating bootstrapped confidence intervals for the mediation effect for a binary outcome when user argument \code{CI="boot"}.
#' 
#' @param treatment    a vector indicating treatment group
#' @param t.seq        a vector of unique time points for each observation
#' @param m            matrix of mediator values in wide format
#' @param outcome      matrix of outcome values in wide format
#' @param coeff_data   a merged dataset of indirect and direct effects and CIs estimated from \code{bootci_coeff_binary}
#' @param replicates   number of replicates for bootstrapping CIs. Default = 1000.
#' 
#' @return \item{timeseq}{time points of estimation}
#' @return \item{alpha_hat}{time-varying treatment effect on the mediator}
#' @return \item{CI.lower.a}{CI lower limit for estimated coefficient alpha_hat}
#' @return \item{CI.upper.a}{CI upper limit for estimated coefficient alpha_hat}
#' @return \item{gamma_hat}{time-varying treatment effect on the outcome (direct effect)}
#' @return \item{CI.lower.g}{CI lower limit for estimated coefficient gamma_hat}
#' @return \item{CI.upper.g}{CI upper limit for estimated coefficient gamma_hat}
#' @return \item{beta_hat}{time-varying effect of the mediator on the outcome}
#' @return \item{CI.lower.b}{CI lower limit for estimated coefficient beta_hat}
#' @return \item{CI.upper.b}{CI upper limit for estimated coefficient beta_hat}
#' @return \item{tau_hat}{time-varying treatment effect on outcome (total effect)}
#' @return \item{CI.lower.t}{CI lower limit for estimated coefficient tau_hat}
#' @return \item{CI.upper.t}{CI upper limit for estimated coefficient tau_hat}
#' @return \item{medEffect}{time varying mediation effect}
#' @return \item{CI.lower}{CI lower limit for medEffect}
#' @return \item{CI.upper}{CI upper limit for medEffect}
#' 
#' @export
#' 

bootci_tvmb <- function(treatment, t.seq, m, outcome, coeff_data, replicates = 1000){
  
  set.seed(27)
  reps <- replicates
  
  start.time <- Sys.time()
  print("Beginning bootstrap for mediation effect CIs.")
  
  n <- length(treatment)
  nm <- nrow(outcome)
  
  #matrix for indirect effects for each individual at each time
  IE <- matrix(NA, nrow = reps, ncol = nm-1)
  
  for(i in 1:reps){
    
    if (i < reps) {
      cat(sprintf("Boostrapping iteration %03d", i), " \r")
    } else {
      print(sprintf("Boostrapping iteration %03d", i))
    }
    
    #get indexes to use for bootstrap sample
    index1 <- sample(1:n, size=n, replace=TRUE)
    
    aAllTemp <- vector()
    bAllTemp <- vector()
    
    #fit m ~ x for first t.seq time points and extract slope
    for(k in 1:nm){
      fit1 <- lm((m[k,index1]) ~ treatment[index1], na.action=na.omit)
      aAllTemp <- append(aAllTemp, fit1$coefficients[[2]])
    }
    
    for(j in 2:nm){
      
      fit2 <- glm(outcome[j,index1] ~ treatment[index1] + m[j-1,index1],
                  family="binomial", na.action=na.omit)
      bHat <- fit2$coefficients[[3]]
      gHat <- fit2$coefficients[[2]]
      
      sd2 <- sqrt(gHat^2*var(treatment[index1], na.rm=TRUE) + 
                  bHat^2*var(m[(j-1),index1], na.rm=TRUE) + 
                  2*gHat*bHat*cov(treatment[index1],m[(j-1),index1], use="complete.obs") + 
                    (pi^2/3))
      
      #append standardized coefficient
      bAllTemp <- append(bAllTemp, bHat/sd2)
    }  
    
    #calculate mediation effect for an individual at each time point
    
    t.seq.b <- t.seq
    t.seq.b <- t.seq.b[-1]
    
    t.seq.b2 <- t.seq
    t.seq.b2 <- t.seq.b2[-1]
    
    coeff_a_temp <- cbind(t.seq, aAllTemp)
    coeff_b_temp <- cbind(t.seq.b2, bAllTemp)
    coeff_dat <- merge(coeff_a_temp, coeff_b_temp, by.x = "t.seq", by.y = "t.seq.b2",
                       all.x = TRUE)
    
    #b(t)*a(t-1) because `a` starts at t=1 while `b` starts at t=2
    for(l in 1:nrow(coeff_dat)){
      if(!is.na(coeff_dat$bAllTemp[l])){
        coeff_dat$medProd[l] <- coeff_dat$bAllTemp[l]*coeff_dat$aAllTemp[l-1]
      }
    }
    
    #calculate smoothed mediation effects
    medProdTemp <- coeff_dat$medProd
    medProdTemp <- medProdTemp[which(!is.na(medProdTemp))]
    
    smooth <- loess(medProdTemp ~ t.seq.b2[1:length(t.seq.b2)], span = 0.5, degree=1)
    
    pred <- predict(smooth, t.seq[2:nm])
    IE[i,] <- pred
  }
  
  #calculate and smooth 2.5 and 97.5 quantiles from bootstrapping
  quantiles <- matrix(NA, nrow=2, ncol=(nm-1))
  lower <- 0.025
  upper <- 1 - lower
  
  for(i in 1:(nm-1)){
    quantiles[1,i] <- quantile(IE[,i], c(lower), na.rm=TRUE)
    quantiles[2,i] <- quantile(IE[,i], c(upper), na.rm=TRUE)
  }
  
  smoothLow <- loess(quantiles[1,] ~ t.seq[2:nm], span = 0.1, degree=1)
  smoothUp <- loess(quantiles[2,] ~ t.seq[2:nm], span = 0.1, degree=1)
  
  CI_1 <- data.frame(cbind(t.seq.b2, smoothLow$fitted, smoothUp$fitted))
  names(CI_1) <- c("t.seq", "CI.lower", "CI.upper")
  
  IE_t <- t(IE)
  IE_t <- data.frame(cbind(t.seq.b, IE_t))
  
  final_CI <- merge(coeff_data, CI_1, all.x = TRUE)
  
  end.time <- Sys.time()
  total.time <- end.time - start.time
  print(sprintf("Process complete. Elapsed time = %.3f secs", 
                as.numeric(total.time, units = "secs")))
  
  final_list <- list(bootstrap_result = IE_t, all_results = final_CI)
  return(final_list)
}
