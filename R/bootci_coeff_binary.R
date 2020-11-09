#' Bootstrap samples to estimate confidence intervals for binary outcome coefficients
#' 
#' Internal function for estimating bootstrapped confidence intervals for the coefficients of the mediation model for a binary outcome.
#' 
#' @param treatment    a vector indicating treatment group
#' @param t.seq        a vector of unique time points for each observation
#' @param m            matrix of mediator values in wide format
#' @param outcome      matrix of outcome values in wide format
#' @param replicates   Number of replicates for bootstrapping confidence intervals. Default = 1000.
#' 
#' @return \item{t.seq}{time points of estimation}
#' @return \item{CI.lower.a}{CI lower limit for alpha_hat}
#' @return \item{CI.upper.a}{CI upper limit for alpha_hat}
#' @return \item{CI.lower.g}{CI lower limit for gamma_hat}
#' @return \item{CI.upper.g}{CI upper limit for gamma_hat}
#' @return \item{CI.lower.b}{CI lower limit for beta_hat}
#' @return \item{CI.upper.b}{CI upper limit for beta_hat}
#' @return \item{CI.lower.t}{CI lower limit for tau_hat}
#' @return \item{CI.upper.t}{CI upper limit for tau_hat}
#' 
#' @export
#' 

bootci_coeff_binary <- function(treatment, t.seq, m, outcome, replicates = 1000){
 
  set.seed(27)
  reps <- replicates
  
  start.time <- Sys.time()
  print("Beginning bootstrap for coefficient CIs.")
  
  n <- length(treatment)
  nm <- nrow(outcome)
  
  #matrices for total, direct and indirect effects for each individual at each time
  IE_a <- matrix(NA, nrow = reps, ncol = nm)
  IE_g <- matrix(NA, nrow = reps, ncol = nm)
  IE_b <- matrix(NA, nrow = reps, ncol = nm)
  IE_t <- matrix(NA, nrow = reps, ncol = nm)
  
  for(i in 1:reps){
    if (i < reps) {
      cat(sprintf("Boostrapping iteration %03d", i), " \r")
    } else {
      print(sprintf("Boostrapping iteration %03d", i))
    }
    
    #get indexes to use for bootstrap sample
    index1 <- sample(1:n, size=n, replace=TRUE)
    
    aAllTemp <- vector()
    gAllTemp <- vector()
    bAllTemp <- vector()
    tAllTemp <- vector()
    
    #fit mediator(m) ~ exposure(x) for first t.seq time points and extract slope
    for(k in 1:nm){
      fit1 <- lm((m[k,index1]) ~ treatment[index1], na.action=na.omit)
      aAllTemp <- append(aAllTemp, fit1$coefficients[[2]])
    }
    
    #fit outcome(y) ~ mediator(m) + exposure(x) for first t.seq time points and extract slope
    for(j in 2:nm){
      
      fit2 <- glm(outcome[j,index1] ~ treatment[index1] + m[j-1,index1],
                  family="binomial", na.action=na.omit)
      bHat <- fit2$coefficients[[3]]
      gHat <- fit2$coefficients[[2]]
      
      sd2 <- sqrt(gHat^2*var(treatment[index1], na.rm=TRUE) + 
                  bHat^2*var(m[(j-1),index1], na.rm=TRUE) + 
                  2*gHat*bHat*cov(treatment[index1],m[(j-1),index1], use="complete.obs") + 
                    (pi^2/3))
      
      # append standardized coefficient
      gAllTemp <- append(gAllTemp, gHat/sd2)
      bAllTemp <- append(bAllTemp, bHat/sd2)
    }
    
    # fit outcome(y) ~ exposure(x) for first t.seq time points and extract slope
    for(l in 2:nm){

      fit3 <- glm(outcome[l,index1] ~ treatment[index1], family="binomial", na.action=na.omit)
      tHat <- fit3$coefficients[[2]]

      sd1 <- sqrt(tHat^2*var(treatment[index1], na.rm=TRUE) + (pi^2/3))

      # append standardized coefficient to list of all coefficients
      tAllTemp <- append(tAllTemp, tHat/sd1)
    }
    
    t.seq.b <- t.seq
    t.seq.b <- t.seq.b[-1]
    
    t.seq.b2 <- t.seq
    t.seq.b2 <- t.seq.b2[-1]
    
    coeff_a_temp <- cbind(t.seq, aAllTemp)
    coeff_g_temp <- cbind(t.seq.b2, gAllTemp)
    coeff_b_temp <- cbind(t.seq.b2, bAllTemp)
    coeff_t_temp <- cbind(t.seq.b2, tAllTemp)
    
    coeff_dat1 <- merge(coeff_a_temp, coeff_b_temp, by.x = "t.seq", by.y = "t.seq.b2",
                        all.x = TRUE)
    coeff_dat2 <- merge(coeff_dat1, coeff_g_temp, by.x = "t.seq", by.y = "t.seq.b2",
                       all.x = TRUE)
    coeff_dat <- merge(coeff_dat2, coeff_t_temp, by.x = "t.seq", by.y = "t.seq.b2",
                       all.x = TRUE)
    
    smootha <- loess(aAllTemp ~ t.seq[1:length(t.seq)], span = 0.3, degree = 1)
    smoothg <- loess(gAllTemp ~ t.seq.b[1:length(t.seq.b2)], span = 0.2, degree = 1)
    smoothb <- loess(bAllTemp ~ t.seq.b[1:length(t.seq.b2)], span = 0.2, degree = 1)
    smootht <- loess(tAllTemp ~ t.seq.b[1:length(t.seq.b2)], span = 0.2, degree = 1)
    
    
    pred_a <- predict(smootha, t.seq[1:nm])
    IE_a[i,] <- pred_a
    
    pred_g <- predict(smoothg, t.seq[1:nm])
    IE_g[i,] <- pred_g

    pred_b <- predict(smoothb, t.seq[1:nm])
    IE_b[i,] <- pred_b
    
    pred_t <- predict(smootht, t.seq[1:nm])
    IE_t[i,] <- pred_t
  }
  
  #calculate and smooth 2.5 and 97.5 quantiles from bootstrapping
  quantiles_a <- matrix(NA, nrow=2, ncol=nm)
  quantiles_g <- matrix(NA, nrow=2, ncol=nm)
  quantiles_b <- matrix(NA, nrow=2, ncol=nm)
  quantiles_t <- matrix(NA, nrow=2, ncol=nm)
  lower <- 0.025
  upper <- 1 - lower
  
  for(i in 1:nm){
    quantiles_a[1,i] <- quantile(IE_a[,i], c(lower), na.rm=TRUE)
    quantiles_a[2,i] <- quantile(IE_a[,i], c(upper), na.rm=TRUE)
  }
  
  for(i in 1:nm){
    quantiles_g[1,i] <- quantile(IE_g[,i], c(lower), na.rm=TRUE)
    quantiles_g[2,i] <- quantile(IE_g[,i], c(upper), na.rm=TRUE)
  }
  
  for(i in 1:nm){
    quantiles_b[1,i] <- quantile(IE_b[,i], c(lower), na.rm=TRUE)
    quantiles_b[2,i] <- quantile(IE_b[,i], c(upper), na.rm=TRUE)
  }
  
  for(i in 1:nm){
    quantiles_t[1,i] <- quantile(IE_t[,i], c(lower), na.rm=TRUE)
    quantiles_t[2,i] <- quantile(IE_t[,i], c(upper), na.rm=TRUE)
  }
  
  smoothLow_a <- loess(quantiles_a[1,] ~ t.seq[1:nm], span = 0.1, degree=1)
  smoothUp_a <- loess(quantiles_a[2,] ~ t.seq[1:nm], span = 0.1, degree=1)
  
  smoothLow_g <- loess(quantiles_g[1,] ~ t.seq[1:nm], span = 0.1, degree=1)
  smoothUp_g <- loess(quantiles_g[2,] ~ t.seq[1:nm], span = 0.1, degree=1)
  
  smoothLow_b <- loess(quantiles_b[1,] ~ t.seq[1:nm], span = 0.1, degree=1)
  smoothUp_b <- loess(quantiles_b[2,] ~ t.seq[1:nm], span = 0.1, degree=1)
  
  smoothLow_t <- loess(quantiles_t[1,] ~ t.seq[1:nm], span = 0.1, degree=1)
  smoothUp_t <- loess(quantiles_t[2,] ~ t.seq[1:nm], span = 0.1, degree=1)
  
  #creating a dataframe with the time sequences, mediation effect and quantiles
  test_t1 <- data.frame(cbind(t.seq, smoothLow_a$fitted, smoothUp_a$fitted))
  test_t2 <- data.frame(cbind(t.seq.b2, smoothLow_g$fitted, smoothUp_g$fitted))
  test_t3 <- data.frame(cbind(t.seq.b2, smoothLow_b$fitted, smoothUp_b$fitted))
  test_t4 <- data.frame(cbind(t.seq.b2, smoothLow_t$fitted, smoothUp_t$fitted))
  
  names(test_t1) <- c("t.seq", "CI.lower.a", "CI.upper.a")
  names(test_t2) <- c("t.seq", "CI.lower.g", "CI.upper.g")
  names(test_t3) <- c("t.seq", "CI.lower.b", "CI.upper.b")
  names(test_t4) <- c("t.seq", "CI.lower.t", "CI.upper.t")
  
  coeff_all1 <- merge(test_t1, test_t2, all.x = TRUE)
  coeff_all2 <- merge(coeff_all1, test_t3, all.x = TRUE)
  coeff_all <- merge(coeff_all2, test_t4, all.x = TRUE)
  
  end.time <- Sys.time()
  total.time <- end.time - start.time
  print(sprintf("Process complete. Elapsed time = %.3f secs", 
                as.numeric(total.time, units = "secs")))
  
  return(coeff_all)
}
