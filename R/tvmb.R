#' Start of the tvmediation function for binary outcome function
#' @export

tvmb <- function(treatment, t.seq, mediator, outcome, plot = FALSE, CI="boot", replicates = 500, verbose = "FALSE"){
  
  # Estimating time-varying mediation effect for binary outcome function
  #
  # Args:
  #   treatment   -->   a vector with treatment values
  #   t.seq       -->   a vector of unique time points for each observation
  #   mediator    -->   matrix of mediator values in wide format
  #   outcome     -->   matrix of outcome outcomes in wide format
  #   plot        -->   TRUE or FALSE for plotting mediation effect                     Default = "FALSE".
  #   CI          -->   "none" or "boot" method of deriving confidence intervals.       Default = "boot".
  #   replicates  -->   Number of replicates for bootstrapping confidence intervals.    Default = 500.
  #   verbose     -->   TRUE or FALSE for printing results to screen.                   Default = "FALSE"
  #
  # Returns:
  #   alpha1_hat       -->   estiamted treatment effect on mediator
  #   beta2_hat        -->   estiamted mediation effect on outcome
  #   medEffect        -->   time varying mediation effect
  #   CI.low           -->   Lower confidence intervals
  #   CI.upper         -->   Upper confidence intervals
  #
  # Optional Returns:
  #   plot1_a1         -->   plot for alpha1_hat across t.seq
  #   plot2_b2         -->   plot for beta2_hat across t.seq
  #   MedEff           -->   plot for medEffect across t.seq
  #   MedEff_CI        -->   plot for CIs of medEffect
  #   bootstrap        -->   plot for estimated medEffects from 
  #                          bootstrapped samples across t.seq
  ##

# Checking for any NA in the treatment vector and removing those indeces
# from treatment as well as outcome and mediator vectors

  index = vector()  
  index=which(!is.na(treatment))
  treatment = treatment[index]
  outcome = outcome[,index]
  m = mediator[,index]

  t.seq <- sort(unique(t.seq)) # check if the t.seq is unique or not

  # set n and nm
  n = length(treatment)
  nm = nrow(outcome)

  ##### *************************************************** #####
  ##### Estimating mediation effect product and difference  #####
  ##### *************************************************** #####

  # Checking for sparseness of the outcome matrix. If all observations
  # corresponding to each time point is missing the function will stop.

  index_ex = vector()

  for(i in 2:nm){
    check_pt <- c(is.na(outcome[i,]))
    num_TRUE <- sum(check_pt, na.rm = TRUE)
    if(num_TRUE == n){
      index_ex <- append(index_ex,i)
      print(paste("Observations missing for row",i,"of the outcome matrix."))
      next()
    }
  }

  if(length(index_ex) == 0){
 
  # Regression of exposure(treatment - y) on mediator(x)
  # fit a1 coefficient for each time point
  a1All = vector()
  
  for(i in 1:nm){
    fit1 = lm(m[i,] ~ treatment,na.action=na.omit)
    a1All = append(a1All,fit1$coefficients[[2]])
  }
  
  # create smoothing line
  smootha1 = loess(a1All ~ t.seq[1:nm], span = 0.3, degree=1)
  
  # creating a dataframe with the time sequences, a1 and smoothed coeffeicients
  test1 <- data.frame(cbind(t.seq, a1All, smootha1$fitted))
  names(test1)[3] <- "smootha1"
  
  
  # Regression of outcome(y) on mediator(x2) and exposure(treatment - x1)
  # fit b1 and b2 coefficients for each time point
  b2All = vector()
  b1All = vector()
  sd2Real = vector()
  
  for(i in 2:nm){
    fit2 = glm(outcome[i,] ~ treatment + m[(i-1),],family="binomial",na.action=na.omit)
    
    b2Hat = fit2$coefficients[[3]]
    b1Hat = fit2$coefficients[[2]]
    
    #calculate sd to standardize with
    sd2 = sqrt(b1Hat^2*var(treatment,na.rm=TRUE)+b2Hat^2*var(m[(i-1),],na.rm=TRUE)+2*b1Hat*b2Hat*cov(treatment,m[(i-1),],use="complete.obs")+(pi^2/3))
    
    #append standardized coefficient to list of all coefficients
    b1All = append(b1All,b1Hat/sd2)
    b2All = append(b2All,b2Hat/sd2) 
  }
  
  # smooth
    t.seq.b <- t.seq
    t.seq.b <- t.seq.b[-1]
    smoothb2 = loess(b2All ~ t.seq.b[1:length(t.seq.b)], span = 0.2,degree=1)

  
  # creating a dataframe with the time sequences, b2 and smoothened coeffeicients
    test2 <- data.frame(cbind(t.seq.b, b2All, smoothb2$fitted))
    names(test2)[3] <- "smoothb2"
  
  
  #regress y on x at each time point to use in difference method
    cAll = vector()

    for(i in 2:nm){
    
        fit3 = glm(outcome[i,]~treatment, family="binomial",na.action=na.omit)
        cHat = fit3$coefficients[[2]]
        
        sd1 = sqrt(cHat^2*var(treatment,na.rm=TRUE)+(pi^2/3))
        
        #append standardized coefficient to list of all coefficients
        cAll = append(cAll,cHat/sd1)
    }
    
    coeff_alpha <- cbind(t.seq, a1All)
    coeff_beta <- cbind(t.seq.b, b1All, b2All, cAll)
    coeff_data <- merge(coeff_alpha, coeff_beta, by.x = "t.seq", by.y = "t.seq.b",
                        all.x = TRUE)
  
  #calculate mediation effects
  #really b2(t)*a1(t-1) because a1 starts at t=1 while b2 starts at t=2
    for(i in 1:nrow(coeff_data)){
      if(!is.na(coeff_data$b1All[i])){
        coeff_data$medProd[i] = coeff_data$b2All[i]*coeff_data$a1All[i-1]
        coeff_data$medDif[i] = coeff_data$cAll[i] - coeff_data$b1All[i]
      }
    }
  
  #calculate smooth line for products
    medProd <- coeff_data$medProd
    medProd <- medProd[which(!is.na(medProd))]
    smoothProd = loess(medProd ~ t.seq.b[1:length(t.seq.b)], span = 0.3,degree=1)
  
  #calculate smooth line for differences
    medDif <- coeff_data$medDif
    medDif <- medDif[which(!is.na(medDif))]
    smoothDif = loess(medDif ~ t.seq.b[1:length(t.seq.b)], span = 0.3,degree=1)
  
  #creating a dataframe with the time sequences, mediation effects and smoothened coeffeicients
    test_a <- data.frame(cbind(t.seq.b, medProd, smoothProd$fitted))
    names(test_a)[2] <- "med_pt"
    names(test_a)[3] <- "smooth"
    test_a$type <- "Prod"
    
    test_b <- data.frame(cbind(t.seq.b, medDif, smoothDif$fitted))
    names(test_b)[2] <- "med_pt"
    names(test_b)[3] <- "smooth"
    test_b$type <- "Diff"
    
    test3 <- rbind(test_a, test_b)
    
  ##### ****************************************************** #####
  ##### Bootstrapping samples to estimate confidence intervals #####
  ##### ****************************************************** #####;
  
  if(CI == "boot"){
    
    #bootstrapping
    set.seed(27)
    reps = replicates
    
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
    
    smoothLow = loess(quantiles[1,] ~ t.seq[2:nm], span = 0.1,degree=1)
    smoothUp = loess(quantiles[2,] ~ t.seq[2:nm], span = 0.1,degree=1)
    
    quantiles1 <- t(quantiles)
    
    #creating a dataframe with the time sequences, mediation effect and quantiles
    test4 <- data.frame(cbind(t.seq.b2, medProd, smoothProd$fitted, quantiles1))
    names(test4) <- c("t.seq", "medProd", "smoothProd", "LowQnt", "UpQnt")
    
    IE_t <- t(IE)
    IE_t <- data.frame(cbind(t.seq.b, IE_t))
    
    final_dat <- merge(coeff_data, test4, all.x = TRUE)
    final_results <- final_dat %>%
                     select(-b1All, -cAll, -medDif, -smoothProd)
    names(final_results) <- c("timeseq", "medEffect", "alpha1_hat", "beta2_hat", "CI.low", "CI.upper")
  }else{
    final_dat <- coeff_data
    final_results <- final_dat %>%
      select(t.seq, a1All, b2All, medProd)
    names(final_results) <- c("timeseq", "alpha1_hat", "beta2_hat", "medEffect")
  }
  
  ##### ************************************************ #####
  ##### Creating plots of the mediation effect estimates #####
  ##### ************************************************ #####
  
  if(plot == "TRUE"){
    # First Plot: plotting calculated coefficients and smoothed coefficients using ggplot
    plot1_a1 <- ggplot(data = test1, aes(t.seq,a1All)) +
      geom_point() +
      labs(title = "Plotting the alpha coffecients",
           x = "Time Sequence",
           y = "Alpha1") +
      geom_line(data = test1, aes(t.seq,smootha1), color = "red", size = 0.75)
    
    # Second plot: plotting beta2 coeffeicients (calculated and smoothed) across
    # the time using ggplot
    plot2_b2 <- ggplot(data = test2, aes(t.seq.b, b2All)) +
      geom_point() +
      labs(title = "Plotting the beta2 coffecients",
           x = "Time Sequence",
           y = "Beta2") +
      geom_line(aes(t.seq.b, smoothb2), color = "red", size = 0.75)
    
    # Third plot: plotting the mediation effects across time using ggplot
    plot3 <- ggplot(data = test3, aes(t.seq.b, med_pt, color = as.factor(type))) +
      geom_point() +
      geom_line(aes(t.seq.b, smooth), size = 0.5) +
      labs(title = "Plotting the mediation effect",
           x = "Time Sequence",
           y = "Mediation Effect",
           color = "Effect Type") +
      scale_color_manual(labels = c("Prod", "Diff"),
                         values = c("indianred1", "blue"))
    if(CI == "boot"){
      # Fourth plot: plotting the mediation effect with 95% CIs
      plot4 <- ggplot(data = test4, aes(t.seq.b2, medProd)) +
        geom_point(shape = 1, size = 1.75, stroke = 1.75, color = "blue") +
        geom_line(aes(t.seq.b2, smoothProd), size = 0.8, color = "blue") +
        geom_line(aes(t.seq.b2, LowQnt, color = "red"), size = 1.25) +
        geom_line(aes(t.seq.b2, UpQnt, color = "red"), size = 1.25) +
        geom_line(aes(t.seq.b2, 0)) +
        labs(title = "Mediation Effect with 95% CIs (computed with bootstrap)",
             x = "Time Sequence",
             y = "Mediation Effect") + 
        theme(legend.position = "none")
      
      # Fifth plot: plotting the mediation effect from 500 bootstrap samples
      plot5 <- ggplot(data = IE_t, aes(t.seq.b, V2)) + 
        geom_line() + 
        labs(title = "Bootstrap result of Mediation Effect",
             x = "Time Sequence",
             y = "Mediation Effect")
      for (i in 3:ncol(IE_t)) {
        x <- data.frame(cbind(t.seq.b,IE_t[,i]))
        names(x)[2] <- "val"
        plot5 <- plot5 + geom_line(data = x, aes(t.seq.b, val))
      }
      
      plot_results <- list("plot1_a1" = plot1_a1,
                      "plot2_b2" = plot2_b2,
                      "MedEff" = plot3,
                      "MedEff_CI" = plot4,
                      "bootstrap" = plot5)
    }else{
      plot_results <- list("plot1_a1" = plot1_a1,
                           "plot2_b2" = plot2_b2,
                           "MedEff" = plot3)
  }
}
    
  #end of 'if' condition on plot
    
  if(verbose == "TRUE"){
    if(plot == "TRUE"){
      print(final_results)
      print(plot_results)
    }else{
      print(final_results)
    }
  }
  
  # Enclosing all the plots in a list object to return
    if(plot == "TRUE" & CI == "boot"){
      results <- list("Estimates" = final_results,
                      "plot1_a1" = plot1_a1,
                      "plot2_b2" = plot2_b2,
                      "MedEff" = plot3,
                      "MedEff_CI" = plot4,
                      "bootstrap" = plot5)
    }
    else if(plot == "TRUE" & CI != "boot"){
      results <- list("Estimates" = final_results,
                      "plot1_a1" = plot1_a1,
                      "plot2_b2" = plot2_b2,
                      "MedEff" = plot3)  
    }
    else{
      results <- list("Estimates" = final_results)
    }
    
    return(results)
   
}else{
      print("Function tvmb stopped running because of missing observations from the outcome matrix.")
}
  
} #end of function
