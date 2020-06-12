#' Time Varying Mediation Function: Binary Outcome and Two Treatment (Exposure) Groups
#' @export

tvmb <- function(treatment, t.seq, mediator, outcome, plot = FALSE, CI="boot", replicates = 1000, verbose = FALSE)
{
  
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
  #   timeseq          -->   time points of estimation
  #   alpha1_hat       -->   estiamted exposure effect on mediator (indirect effect)
  #   CI.lower.a1      -->   Lower confidence intervals for alpha1_hat
  #   CI.upper.a1      -->   Upper confidence intervals for alpha1_hat
  #   beta2_hat        -->   estiamted mediation effect on outcome (indirect effect)
  #   CI.lower.b2      -->   Lower confidence intervals for beta2_hat
  #   CI.upper.b2      -->   Upper confidence intervals for beta2_hat
  #   b1All            -->   estimated exposure effect on outcome (direct effect)
  #   cAll             -->   estimated effect of exposure on outcome (total effect)
  #   medDiff          -->   time varying mediation effect (difference term)
  #   medEffect        -->   time varying mediation effect (product term)
  #   CI.low           -->   Lower confidence intervals for medEffect
  #   CI.upper         -->   Upper confidence intervals for medEffect
  #
  # Optional Returns:
  #   plot1_a1         -->   plot for alpha1_hat with CIs across t.seq
  #   plot2_b2         -->   plot for beta2_hat with CIs across t.seq
  #   MedEff           -->   plot for mediation effects (difference and product) across t.seq
  #   MedEff_CI        -->   plot for CIs of medEffect
  #   bootstrap        -->   plot for estimated medEffects from bootstrapped samples across t.seq
  ##

  ## Testing the class type of the arguments passed in the function ##
  ctm <- class(mediator)
  cto <- class(outcome)
  ctt <- class(treatment)
  cts <- class(t.seq)
  flag <- 0
  
  if(ctm != "matrix"){
    print("Error: `mediator` is not of class type `matrix`.")
    flag <- flag + 1
  }
  if(cto != "matrix"){
    print("Error: `outcome` is not of class type `matrix`.")
    flag <- flag + 1
  }
  if(is.vector(treatment) != TRUE || ctt != "numeric"){
    print("Error: `treatment` is not of class type `numeric vector`.")
    flag <- flag + 1
  }
  if(is.vector(t.seq) != TRUE || cts != "numeric"){
    print("Error: `t.seq` is not of class type `numeric vector`.")
    flag <- flag + 1
  }
  if(flag == 0){
    # Checking for any NA in the treatment vector and removing those indeces
    # from treatment as well as outcome and mediator matrices
    if(CI == "boot" || CI == "none"){
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
        
        # Regression of exposure treatment(y) on mediator(x)
        # fit a1 coefficient for each time point
        a1All = vector()
        
        for(i in 1:nm){
          fit1 = lm(m[i,] ~ treatment,na.action=na.omit)
          a1All = append(a1All,fit1$coefficients[[2]])
        }
        
        # create smoothing line
        smootha1 = loess(a1All ~ t.seq[1:nm], span = 0.3, degree=1)
        
        # creating a dataframe with the time sequences, a1 and smoothed coefficients
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
        
        
        # creating a dataframe with the time sequences, b2 and smoothened coefficients
        test2 <- data.frame(cbind(t.seq.b, b2All, smoothb2$fitted))
        names(test2)[3] <- "smoothb2"
        
        
        #regress outcome(y) on exposure(x) at each time point to use in difference method
        cAll = vector()
        
        for(i in 2:nm){
          
          fit3 = glm(outcome[i,]~treatment, family="binomial",na.action=na.omit)
          cHat = fit3$coefficients[[2]]
          
          sd1 = sqrt(cHat^2*var(treatment,na.rm=TRUE)+(pi^2/3))
          
          #append standardized coefficient to list of all coefficients
          cAll = append(cAll,cHat/sd1)
        }
        
        
        ##### *********************************************************************** #####
        ##### Bootstrapping samples to estimate confidence intervals for coefficients #####
        ##### *********************************************************************** #####
        
        coeff_CI <- bootci_coeff_binary(treatment, t.seq, m, outcome, replicates)
        
        #*********************************************************************************#
        
        
        #### Formatting the results into a single dataframe ####
        
        coeff_alpha1 <- merge(test1, coeff_CI, by.x = "t.seq") %>%
          select(-CI.lower.b2, -CI.upper.b2)
        coeff_beta1 <- cbind(t.seq.b, b1All, cAll)
        coeff_beta2 <- merge(test2, coeff_CI, by.x = "t.seq.b", by.y = "t.seq", all.y = TRUE) %>%
          select(-CI.lower.a1, -CI.upper.a1)
        
        coeff_data1 <- merge(coeff_alpha1, coeff_beta2, by.x = "t.seq", by.y = "t.seq.b",
                             all.x = TRUE)
        
        coeff_data <- merge(coeff_data1, coeff_beta1, by.x = "t.seq", by.y = "t.seq.b",
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
        
        #creating a dataframe with the time sequences, mediation effects and smoothened coefficients
        test_a <- data.frame(cbind(t.seq.b, medProd, smoothProd$fitted))
        names(test_a)[2] <- "med_pt"
        names(test_a)[3] <- "smooth_Prod"
        test_a$type <- "Prod"
        
        test_b <- data.frame(cbind(t.seq.b, medDif, smoothDif$fitted))
        names(test_b)[2] <- "med_pt"
        names(test_b)[3] <- "smooth"
        test_b$type <- "Diff"
        
        # test3 <- rbind(test_a, test_b)
        
        coeff_data <- merge(coeff_data, test_a, by.x = "t.seq", by.y = "t.seq.b",
                            all.x = TRUE) %>%
          select(-med_pt, -type)
        names(coeff_data)[14] <- "smooth_medProd"
        
        coeff_data <- merge(coeff_data, test_b, by.x = "t.seq", by.y = "t.seq.b",
                            all.x = TRUE) %>%
          select(-med_pt, -type)
        names(coeff_data)[15] <- "smooth_medDif"
        
        
        ##### ****************************************************** #####
        ##### Bootstrapping samples to estimate confidence intervals #####
        ##### ****************************************************** #####;
        
        if(CI == "boot"){
          list_all <- bootci_tvmb(treatment, t.seq, m, outcome, coeff_data, replicates)
          IE_t <- list_all$bootstrap_result
          final_dat <- list_all$all_results
          final_dat1 <- final_dat %>%
            select(- c(a1All, b2All, medProd, medDif))
          final_results <- final_dat1[c(1:9, 11, 10, 12, 13)]
          names(final_results)[c(1, 2, 5, 10, 11)] <- c("timeseq", "alpha1_hat", "beta2_hat", "medDiff", "medEffect")
        }else{
          final_dat <- coeff_data
          final_dat1 <- final_dat %>%
            select(- c(a1All, b2All, medProd, medDif))
          final_results <- final_dat1[c(1:9, 11, 10)]
          names(final_results)[c(1, 2, 5, 10, 11)] <- c("timeseq", "alpha1_hat", "beta2_hat", "medDiff", "medEffect")
        }
        
        ##### ************************************************ #####
        ##### Creating plots of the mediation effect estimates #####
        ##### ************************************************ #####
        
        if(plot == TRUE){
          
          lt <- length(final_results$timeseq)
          l <- min(final_results$timeseq)
          u <- max(final_results$timeseq)
          
          if(u <= 1){
            if(lt <= 10){
              i <- 0.2
            }else if(lt>10 && lt<=20){
              i <- 0.15
            }else if(lt>20){
              i <- 0.25
            }
          }else if(u>1 && u<=30){
            i <- 2
          }else if(u>30 && u <=50){
            i <- 5
          }else if(u>50){
            i <- 10
          }
          
          # First Plot: plotting alpha1 coefficients (smoothed) using across
          # the time using ggplot
          plot1_a1 <- ggplot(data = final_results, aes(timeseq, alpha1_hat)) +
            geom_line(color = "red", size = 0.75) +
            geom_line(aes(timeseq, CI.lower.a1), color = "blue", size = 0.8, linetype = "dashed") +
            geom_line(aes(timeseq, CI.upper.a1), color = "blue", size = 0.8, linetype = "dashed") +
            labs(title = "Plotting the alpha coefficients",
                 x = "Time (in days)",
                 y = "Alpha1") +
            scale_x_continuous(breaks = seq(l, u, i))
          
          # Second plot: plotting beta2 coefficients (smoothed) across
          # the time using ggplot
          plot2_b2 <- ggplot(data = final_results, aes(timeseq, beta2_hat)) +
            geom_line(color = "red", size = 0.75) +
            geom_line(aes(timeseq, CI.lower.b2), color = "blue", size = 0.8, linetype = "dashed") +
            geom_line(aes(timeseq, CI.upper.b2), color = "blue", size = 0.8, linetype = "dashed") +
            labs(title = "Plotting the beta2 coefficients",
                 x = "Time (in days)",
                 y = "Beta2") +
            scale_x_continuous(breaks = seq(l, u, i))
          
          # Third plot: plotting the mediation effects across time using ggplot
          plot3_a <- ggplot(data = final_results, aes(timeseq, medDiff)) +
            geom_line(color = "black", size = 0.75) +
            labs(title = "Plotting the mediation (difference) effect",
                 x = "Time (in days)",
                 y = "Mediation Effect") +
            scale_x_continuous(breaks = seq(l, u, i))
          
          plot3_b <- ggplot(data = final_results, aes(timeseq, medEffect)) +
            geom_line(color = "red", size = 0.75) +
            labs(title = "Plotting the mediation (product) effect",
                 x = "Time (in days)",
                 y = "Mediation Effect") +
            scale_x_continuous(breaks = seq(l, u, i))
          
          plot3 <- ggarrange(plot3_a, plot3_b)
          
          if(CI == "boot"){
            # Fourth plot: plotting the mediation effect with 95% CIs
            plot4 <- ggplot(data = final_results, aes(timeseq, medEffect)) +
              geom_line(size = 1, color = "red") +
              geom_line(aes(timeseq, CI.low), color = "blue", size = 0.8, linetype = "dashed") +
              geom_line(aes(timeseq, CI.upper), color = "blue", size = 0.8, linetype = "dashed") +
              geom_line(aes(timeseq, 0)) +
              labs(title = "Mediation Effect with 95% CIs (computed with percentile bootstrap)",
                   x = "Time (in days)",
                   y = "Mediation Effect") + 
              theme(legend.position = "none") +
              scale_x_continuous(breaks = seq(l, u, i))
            
            # Fifth plot: plotting the mediation effect from 500 bootstrap samples
            plot5 <- ggplot(data = IE_t, aes(t.seq.b, V2)) + 
              geom_line() + 
              labs(title = "Bootstrap result of Mediation Effect",
                   x = "Time (in days)",
                   y = "Mediation Effect") +
              scale_x_continuous(breaks = seq(l, u, i))
            for (i in 2:ncol(IE_t)) {
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
        
        if(verbose == TRUE){
          if(plot == TRUE){
            print(final_results)
            print(plot_results)
          }else{
            print(final_results)
          }
        }
        
        # Enclosing all the plots in a list object to return
        if(plot == TRUE & CI == "boot"){
          results <- list("Estimates" = final_results,
                          "plot1_a1" = plot1_a1,
                          "plot2_b2" = plot2_b2,
                          "MedEff" = plot3,
                          "MedEff_CI" = plot4,
                          "bootstrap" = plot5)
        }
        else if(plot == TRUE & CI != "boot"){
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
    } else{
      print(paste("Error:Accepted values for CI are 'boot' and 'none';you have entered an unacceptable value for CI."))
    }#end of CI="boot" if condition
  }else{
    print(paste("tvmb() stopped execution due to unacceptable class type of function argument(s)."))
  }
  

} #end of function
