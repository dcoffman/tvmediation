##### *********************************************************************** #####
##### Bootstrapping samples to estimate confidence intervals for coefficients #####
##### *********************************************************************** #####

bootci_coeff_binary <- function(treatment, t.seq, m, outcome, replicates = 500){
  #bootstrapping
  set.seed(27)
  reps = replicates
  
  start.time <- Sys.time()
  print("Beginning bootstrap for coefficient CIs.")
  
  n = length(treatment)
  nm = nrow(outcome)
  
  #matrix with indirect (mediation) effects for each individual at each time
  IE_a1 = matrix(NA, nrow = reps, ncol = nm)
  # IE_b1 = matrix(NA, nrow = reps, ncol = nm)
  IE_b2 = matrix(NA, nrow = reps, ncol = nm)
  
  #take 500 or replicates number of bootstrap samples
  for(i in 1:reps){
    #get indexes to use for bootstrap sample
    index1 = sample(1:n, size=n,replace=TRUE)
    
    a1AllTemp = vector()
    # b1AllTemp = vector()
    b2AllTemp = vector()
    
    #fit mediator(m)~exposure(x) for first t.seq time points and extract slope
    for(k in 1:nm){
      fit1 <- lm((m[k,index1]) ~ treatment[index1],na.action=na.omit)
      a1AllTemp = append(a1AllTemp,fit1$coefficients[[2]])
    }
    
    #fit outcome(y) ~ mediator(m) + exposure(x) for first t.seq time points and extract slope
    for(j in 2:nm){
      
      fit2 = glm(outcome[j,index1] ~ treatment[index1] + m[j-1,index1],family="binomial",na.action=na.omit)
      b2Hat = fit2$coefficients[[3]]
      b1Hat = fit2$coefficients[[2]]
      
      sd2 = sqrt(b1Hat^2*var(treatment[index1],na.rm=TRUE)+b2Hat^2*var(m[(j-1),index1],na.rm=TRUE)+2*b1Hat*b2Hat*cov(treatment[index1],m[(j-1),index1],use="complete.obs")+(pi^2/3))
      
      #append standardized coefficient
      # b1AllTemp = append(b1AllTemp,b1Hat/sd2)
      b2AllTemp = append(b2AllTemp,b2Hat/sd2)
    }  
    
    t.seq.b2 <- t.seq
    t.seq.b2 <- t.seq.b2[-1]
    
    coeff_a1_temp <- cbind(t.seq, a1AllTemp)
    # coeff_b1_temp <- cbind(t.seq.b2, b1AllTemp)
    coeff_b2_temp <- cbind(t.seq.b2, b2AllTemp)
    
    # coeff_dat1 <- merge(coeff_a1_temp, coeff_b1_temp, by.x = "t.seq", by.y = "t.seq.b2",
    #                     all.x = TRUE)
    # coeff_dat <- merge(coeff_dat1, coeff_b2_temp, by.x = "t.seq", by.y = "t.seq.b2",
    #                    all.x = TRUE)
    coeff_dat <- merge(coeff_a1_temp, coeff_b2_temp, by.x = "t.seq", by.y = "t.seq.b2",
                       all.x = TRUE)
    
    smootha1 = loess(a1AllTemp ~ t.seq[1:length(t.seq)], span = 0.3, degree=1)
    # smoothb1 = loess(b1AllTemp ~ t.seq[1:length(t.seq.b2)], span = 0.2, degree=1)
    smoothb2 = loess(b2AllTemp ~ t.seq.b[1:length(t.seq.b2)], span = 0.2,degree=1)
    
    
    pred_a1 = predict(smootha1,t.seq[1:nm])
    IE_a1[i,] = pred_a1
    
    # pred_b1 = predict(smoothb1,t.seq[1:nm])
    # IE_b1[i,] = pred_b1
    
    pred_b2 = predict(smoothb2,t.seq[1:nm])
    IE_b2[i,] = pred_b2
  }
  
  #calculate and smooth 2.5 and 97.5 quantiles from bootstrapping
  quantiles_a1 = matrix(NA, nrow=2,ncol=nm)
  # quantiles_b1 = matrix(NA, nrow=2,ncol=nm)
  quantiles_b2 = matrix(NA, nrow=2,ncol=nm)
  lower = 0.025
  upper = 1 - lower
  
  for(i in 1:nm){
    quantiles_a1[1,i] = quantile(IE_a1[,i], c(lower),na.rm=TRUE)
    quantiles_a1[2,i] = quantile(IE_a1[,i], c(upper),na.rm=TRUE)
  }
  
  # for(i in 1:nm){
  #   quantiles_b1[1,i] = quantile(IE_b1[,i], c(lower),na.rm=TRUE)
  #   quantiles_b1[2,i] = quantile(IE_b1[,i], c(upper),na.rm=TRUE)
  # }
  
  for(i in 1:nm){
    quantiles_b2[1,i] = quantile(IE_b2[,i], c(lower),na.rm=TRUE)
    quantiles_b2[2,i] = quantile(IE_b2[,i], c(upper),na.rm=TRUE)
  }
  
  smoothLow_a1 = loess(quantiles_a1[1,] ~ t.seq[1:nm], span = 0.1,degree=1)
  smoothUp_a1 = loess(quantiles_a1[2,] ~ t.seq[1:nm], span = 0.1,degree=1)
  
  # smoothLow_b1 = loess(quantiles_b1[1,] ~ t.seq[1:nm], span = 0.1,degree=1)
  # smoothUp_b1 = loess(quantiles_b1[2,] ~ t.seq[1:nm], span = 0.1,degree=1)
  
  smoothLow_b2 = loess(quantiles_b2[1,] ~ t.seq[1:nm], span = 0.1,degree=1)
  smoothUp_b2 = loess(quantiles_b2[2,] ~ t.seq[1:nm], span = 0.1,degree=1)
  
  #creating a dataframe with the time sequences, mediation effect and quantiles
  test_t1 <- data.frame(cbind(t.seq, smoothLow_a1$fitted, smoothUp_a1$fitted))
  # test_t2 <- data.frame(cbind(t.seq.b2, smoothLow_b1$fitted, smoothUp_b1$fitted))
  test_t3 <- data.frame(cbind(t.seq.b2, smoothLow_b2$fitted, smoothUp_b2$fitted))
  
  names(test_t1) <- c("t.seq", "CI.lower.a1", "CI.upper.a1")
  # names(test_t2) <- c("t.seq", "CI.lower.b1", "CI.upper.b1")
  names(test_t3) <- c("t.seq", "CI.lower.b2", "CI.upper.b2")
  
  # coeff_all <- merge(test_t1, test_t2, all.x = TRUE)
  coeff_all <- merge(test_t1, test_t3, all.x = TRUE)
  # coeff_all <- merge(coeff_all, test_t3, all.x = TRUE)
  
  end.time <- Sys.time()
  total.time <- end.time - start.time
  print(sprintf("Process complete. Elapsed time = %.3f secs", as.numeric(total.time, units = "secs")))
  
  return(coeff_all)
}
