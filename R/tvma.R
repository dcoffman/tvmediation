#' Time-varying mediation function for continuous outcome and 2 treatment arms (exposure groups)
#' @export

tvma <- function(treatment, t.seq, mediator, outcome, t.est = t.seq, plot = FALSE, CI="boot", replicates = 1000, verbose = FALSE)
  {
      # Estimate time-varying mediation effect and bootstrap standard errors
      #
      ### ARGUMENTS:
      #   treatment   -->   a vector with treatment values
      #   t.seq       -->   a vector of time points for each observation
      #   mediator    -->   matrix of mediator values in wide format
      #   outcome     -->   matrix of outcome outcomes in wide format
      #   t.est       -->   time points to make the estimation                              Default = t.seq
      #   plot        -->   TRUE or FALSE for plotting mediation effect                     Default = "FALSE"
      #   CI          -->   "none" or "boot" method of deriving confidence intervals.       Default = "boot"
      #   replicates  -->   Number of replicates for bootstrapping confidence intervals.    Default = 1000
      #   verbose     -->   TRUE or FALSE for printing results to screen.                   Default = "FALSE"
      #
      # RETURNS:
      #
      ### Dataframe 'final_results' with the following elements
      #
      #   hat.alpha.1       -->   estimated treatment effect on mediator
      #   hat.beta.2        -->   estimated mediation effect on outcome
      #   est.M             -->   time varying mediation effect
      #   CI.upper.alpha    -->   Upper confidence intervals for coefficient alpha1
      #   CI.lower.alpha    -->   Lower confidence intervals for coefficient alpha1
      #   CI.upper.beta     -->   Upper confidence intervals for coefficient beta2
      #   CI.lower.beta     -->   Lower confidence intervals for coefficient beta2
      # 
      # Optional results based on user argument CI = "boot" passed
      #   boot.se.m         -->   estimated standard error of the mediation effect
      #   CI.upper          -->   Upper confidence intervals
      #   CI.lower          -->   Lower confidence intervals
      #  
      ### Optional plot results
      #
      #   Alpha_CI          -->   plot for hat.alpha.1 across t.seq with CI
      #   Beta_CI           -->   plot for hat.beta.2 across t.seq with CI
      #   MedEff            -->   plot for mediation.effect across t.seq
      #   MedEff_CI         -->   plot for CIs of mediation.effect across t.seq
      #
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
    if(CI == "boot" || CI == "none"){
      deltat <- max(diff(t.seq)) / 2  # half the time between two measures
      N <- length(treatment)
      t.coeff <- NULL
      
      for (j in 2:length(t.seq)){
        # create empty vector, store raw mediator.
        temp.mediator.j  <- NULL
        temp.mediator.j  <- cbind(mediator[j - 1, ], mediator[j, ])
        
        # Derive centered Mediators and Outcomes
        newMO.j.est <- newMediatorOutcome(treatment, temp.mediator.j, outcome[j - 1, ])
        
        # Estimate coefficients, then store them.
        coeff.est <- estCoeff(newMO.j.est)
        t.coeff <- cbind(t.coeff, coeff.est)  # store coeff estimates at t.seq
      }
      
      # EQUATIONS 4 & 5
      est.smooth <- smoothest(t.seq, t.coeff, t.est, deltat)
      
      ## calculating the CI for the coefficients alpha1 and beta2
      coeff_CI_2trt <- bootci_coeff_2trt(treatment, t.seq, mediator, outcome, t.est, deltat, replicates)
      
      test1 <- cbind(as.data.frame(est.smooth), as.data.frame(coeff_CI_2trt), t.est)
      
      # CALCULATE CONFIDENCE INTERVALS
      if(CI == "boot"){
        results_ci <- estBootCIs(treatment, t.seq, mediator, outcome, t.est, deltat, replicates)
        test2 <- cbind(as.data.frame(results_ci), t.est)
        
        final_dat <- merge(test1, test2, all.x = TRUE)
        final_results <- final_dat %>%
          select(-bw_alpha1, -bw_beta2)
        names(final_results)[9] <- c("boot.se.m")
      }
      else{
        final_results <- test1 %>%
          select(-bw_alpha1, -bw_beta2)
      }
      
      #### Plot construction ####
      if(plot == TRUE){
        
        l <- min(final_results$t.est)
        u <- max(final_results$t.est)
        
        if(u <= 1){
          i <- 0.2
        }else if(u>1 && u<=30){
          i <- 2
        }else if(u>30 && u <=50){
          i <- 5
        }else if(u>50){
          i <- 10
        }
        
        # First Plot: plotting alpha1 coefficients across time using ggplot
        plot1_a1 <- ggplot(data = final_results, aes(t.est,hat.alpha.1)) +
          geom_line(color = "red", size = 0.75) +
          geom_line(aes(t.est, CI.lower.alpha), size = 0.8, color = "blue", linetype = "dashed") +
          geom_line(aes(t.est, CI.upper.alpha), size = 0.8, color = "blue", linetype = "dashed") +
          labs(title = "Plotting the alpha1 coefficients",
               x = "Time (in days)",
               y = "Alpha1") +
          scale_x_continuous(breaks = seq(l, u, i))
        
        # Second plot: plotting beta2 coefficients across time using ggplot
        plot2_b2 <- ggplot(data = final_results, aes(t.est,hat.beta.2)) +
          geom_line(color = "red", size = 0.75) +
          geom_line(aes(t.est, CI.lower.beta), size = 0.8, color = "blue", linetype = "dashed") +
          geom_line(aes(t.est, CI.upper.beta), size = 0.8, color = "blue", linetype = "dashed") +
          labs(title = "Plotting the beta2 coefficients",
               x = "Time (in days)",
               y = "Beta2") +
          scale_x_continuous(breaks = seq(l, u, i))
        
        # Third plot: plotting the mediation effect of treatment arm
        plot3 <- ggplot(data = final_results, aes(t.est,est.M)) +
          geom_line(color = "red", size = 0.75) +
          labs(title = "Plotting the time-varying mediation effect",
               x = "Time (in days)",
               y = "Mediation Effect") +
          scale_x_continuous(breaks = seq(l, u, i))
        
        if(CI == "boot"){
          
          # Fourth plot: plotting the mediation effect of treatment with 95% CIs
          plot4 <- ggplot(data = final_results, aes(t.est, est.M)) +
            geom_line(size = 1, color = "red") +
            geom_line(aes(t.est, CI.lower), size = 0.8, color = "blue", linetype = "dashed") +
            geom_line(aes(t.est, CI.upper), size = 0.8, color = "blue", linetype = "dashed") +
            # geom_line(aes(t.est, 0)) +
            labs(title = "Mediation Effect with 95% CIs (computed with percentile bootstrap)",
                 x = "Time (in days)",
                 y = "Mediation Effect") + 
            theme(legend.position = "none") +
            scale_x_continuous(breaks = seq(l, u, i))   
          
          plot_results <- list("Alpha_CI" = plot1_a1,
                               "Beta_CI" = plot2_b2,
                               "MedEff" = plot3,
                               "MedEff_CI" = plot4)
        }else{
          plot_results <- list("Alpha_CI" = plot1_a1,
                               "Beta_CI" = plot2_b2,
                               "MedEff" = plot3)
        }
      }
      
      ## Print results to screen and return them VERBOSE condition ##
      if(verbose == TRUE){
        print("Time Varying Mediation Results:")
        if(plot == TRUE){
          print(final_results)
          print(plot_results)
        }else{
          print(final_results)
        }
      }
      
      ## Enclosing all the plots in a list object to return ##
      if(plot == TRUE & CI == "boot"){
        results <- list("Estimates" = final_results,
                        "Alpha_CI" = plot1_a1,
                        "Beta_CI" = plot2_b2,
                        "MedEff" = plot3,
                        "MedEff_CI" = plot4)
      }
      else if(plot == TRUE & CI != "boot"){
        results <- list("Estimates" = final_results,
                        "Alpha_CI" = plot1_a1,
                        "Beta_CI" = plot2_b2,
                        "MedEff" = plot3)  
      }
      else{
        results <- list("Estimates" = final_results)
      }
      
      return(results)
    }
    else{
      print(paste("Error:Accepted values for CI are 'boot' and 'none';you have entered an unacceptable value for CI."))
    }
  }else{
    print(paste("tvma() stopped execution due to unacceptable class type of function argument(s)."))
  }
}