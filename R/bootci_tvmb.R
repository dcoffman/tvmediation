##### ****************************************************** #####
##### Bootstrapping samples to estimate confidence intervals #####
##### ****************************************************** #####;

bootci_tvmb <- function(treatment, t.seq, m, outcome, coeff_data, replicates = 500){
  #bootstrapping
  set.seed(27)
  reps = replicates
  
  start.time <- Sys.time()
  print("Beginning bootstrap for coefficient CIs.")
  
  n = length(treatment)
  nm = nrow(outcome)
  
  #matrix with indirect (mediation) effects for each individual at each time
  IE = matrix(NA,nrow = reps,ncol = nm-1)
  
  #take 500 or replicates number of bootstrap samples
  for(i in 1:reps){
    #get indexes to use for bootstrap sample
    index1 = sample(1:n,size=n,replace=TRUE)
    
    a1AllTemp = vector()
    b2AllTemp = vector()
    
    #fit m~x for first t.seq time points and extract slope
    for(k in 1:nm){
      fit1 <- lm((m[k,index1]) ~ treatment[index1],na.action=na.omit)
      a1AllTemp = append(a1AllTemp,fit1$coefficients[[2]])
    }
    
    for(j in 2:nm){
      
      fit2 = glm(outcome[j,index1] ~ treatment[index1] + m[j-1,index1],family="binomial",na.action=na.omit)
      b2Hat = fit2$coefficients[[3]]
      b1Hat = fit2$coefficients[[2]]
      
      sd2 = sqrt(b1Hat^2*var(treatment[index1],na.rm=TRUE)+b2Hat^2*var(m[(j-1),index1],na.rm=TRUE)+2*b1Hat*b2Hat*cov(treatment[index1],m[(j-1),index1],use="complete.obs")+(pi^2/3))
      
      #append standardized coefficient
      b2AllTemp = append(b2AllTemp,b2Hat/sd2)
    }  
    
    #calculate mediation effect for an individual at time points 2-50
    #this is taking product of a1[t-1]*b2[t]
    
    t.seq.b <- t.seq
    t.seq.b <- t.seq.b[-1]
    
    t.seq.b2 <- t.seq
    t.seq.b2 <- t.seq.b2[-1]
    
    coeff_a_temp <- cbind(t.seq, a1AllTemp)
    coeff_b_temp <- cbind(t.seq.b2, b2AllTemp)
    coeff_dat <- merge(coeff_a_temp, coeff_b_temp, by.x = "t.seq", by.y = "t.seq.b2",
                       all.x = TRUE)
    
    #calculate mediation effects
    #really b2(t)*a1(t-1) because a1 starts at t=1 while b2 starts at t=2
    for(l in 1:nrow(coeff_dat)){
      if(!is.na(coeff_dat$b2AllTemp[l])){
        coeff_dat$medProd[l] = coeff_dat$b2AllTemp[l]*coeff_dat$a1AllTemp[l-1]
      }
    }
    
    #calculate smooth line for products
    medProdTemp <- coeff_dat$medProd
    medProdTemp <- medProdTemp[which(!is.na(medProd))]
    
    smooth = loess(medProdTemp ~ t.seq.b2[1:length(t.seq.b2)], span = 0.5,degree=1)
    
    pred = predict(smooth,t.seq[2:nm])
    IE[i,] = pred
  }
  
  #calculate and smooth 2.5 and 97.5 quantiles from bootstrapping
  quantiles = matrix(NA, nrow=2,ncol=(nm-1))
  lower = 0.025
  upper = 1 - lower
  
  for(i in 1:(nm-1)){
    quantiles[1,i] = quantile(IE[,i], c(lower),na.rm=TRUE)
    quantiles[2,i] = quantile(IE[,i], c(upper),na.rm=TRUE)
  }
  
  smoothLow = loess(quantiles[1,] ~ t.seq[2:nm], span = 0.1, degree=1)
  smoothUp = loess(quantiles[2,] ~ t.seq[2:nm], span = 0.1, degree=1)
  
  t.seq.b3 <- t.seq.b2[-1]
  CI_1 <- data.frame(cbind(t.seq.b3, smoothLow$fitted, smoothUp$fitted))
  
  #creating a dataframe with the time sequences, mediation effect and quantiles
  # test4_1 <- data.frame(cbind(t.seq.b2, medProd, smoothProd$fitted))
  # names(test4_1) <- c("t.seq", "medProd", "smoothProd")
  # test4 <- merge(test4_1, CI_1, by.x = "t.seq", by.y = "t.seq.b3", all.x = TRUE)
  
  names(CI_1) <- c("t.seq", "CI.low", "CI.upper")
  
  IE_t <- t(IE)
  IE_t <- data.frame(cbind(t.seq.b, IE_t))
  
  final_CI <- merge(coeff_data, CI_1, all.x = TRUE)
  # results <- final_CI
  # names(results) <- c("timeseq", "medEffect", "alpha1_hat", "beta2_hat", "CI.low", "CI.upper")
  
  end.time <- Sys.time()
  total.time <- end.time - start.time
  print(sprintf("Process complete. Elapsed time = %.3f secs", as.numeric(total.time, units = "secs")))
  
  final_list <- list(bootstrap_result = IE_t, all_results = final_CI)
  return(final_list)
}
