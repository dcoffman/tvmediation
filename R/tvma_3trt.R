#' Time-varying mediation function for continuous outcome and 3 treatment arms (exposure groups)
#' @export

tvma_3trt <- function(NRT1, NRT2, t.seq, mediator, outcome, t.est = t.seq, plot = FALSE, CI="boot", replicates = 1000, verbose = FALSE)
  {
  
  # Estimate time-varying mediation effect and bootstrap standard errors
  # for continuous outcome and 3 treatment arms
  #
  # Args:
  #   NRT1        -->   a vector with treatment arm 1 values
  #   NRT2        -->   a vector with treatment arm 2 values
  #   t.seq       -->   a vector of time points for each observation
  #   mediator    -->   matrix of mediator values in wide format
  #   outcome     -->   matrix of outcome outcomes in wide format
  #   t.est       -->   time points to make the estimation                              Default = t.seq
  #   plot        -->   TRUE or FALSE for plotting mediation effect                     Default = "FALSE"
  #   CI          -->   "none" or "boot" method of deriving confidence intervals.       Default = "boot"
  #   replicates  -->   Number of replicates for bootstrapping confidence intervals.    Default = 1000
  #   verbose     -->   TRUE or FALSE for printing results to screen.                   Default = "FALSE"
  #
  # Returns:
  #
  ##### Dataframe 'final_results' with the following elements
  #
  #   hat.alpha1        -->   estimated NRT1 effect on mediator
  #   hat.alpha2        -->   estimated NRT2 effect on mediator
  #   hat.beta3         -->   estimated mediation effect on outcome
  #   hat.mediation1    -->   time varying mediation effect - NRT1 on outcome
  #   hat.mediation2    -->   time varying mediation effect - NRT2 on outcome
  ## Optional results based on user arguments passed
  #   CI.upper.NRT1     -->   Upper confidence intervals for NRT1
  #   CI.lower.NRT1     -->   Lower confidence intervals for NRT1
  #   CI.upper.NRT2     -->   Upper confidence intervals for NRT2
  #   CI.lower.NRT2     -->   Lower confidence intervals for NRT2
  #   SE_MedEff1        -->   estimated standard errors of the mediation effect for NRT1
  #   SE_MedEff2        -->   estimated standard errors of the mediation effect for NRT2
  #
  ##### Optional plot results
  #   plot1_a1         -->   plot for hat.alpha1 across t.seq
  #   plot2_a2         -->   plot for hat.alpha2 across t.seq
  #   plot3_b3         -->   plot for hat.beta3 across t.seq
  #   MedEff_NRT1      -->   plot for hat.mediation1 across t.seq
  #   MedEff_NRT2      -->   plot for hat.mediation2 across t.seq
  #   MedEff_CI_NRT1   -->   plot for CIs of medEffect for NRT1
  #   MedEff_CI_NRT2   -->   plot for CIs of medEffect for NRT2
  ##
  
  results_me <- tvmcurve_3trt(NRT1, NRT2, t.seq, mediator, outcome, t.est)
  alldata <- list(NRT1 = NRT1, NRT2 = NRT2, x = mediator, y = outcome, t.seq = t.seq)
  test1 <- cbind(as.data.frame(results_me), t.est)
  
  if(CI == "boot"){
    results_ci <- bootci_tvm_3trt(replicates, alldata, t.est)
    test2 <- cbind(as.data.frame(results_ci), t.est)
    
    final_dat <- merge(test1, test2, all.x = TRUE)
    final_results <- final_dat %>%
                     select(-orig.mediation1, -orig.mediation2)
    names(final_results)[c(1,7:12)] <- c("timeseq","CI.lower.NRT1", "CI.upper.NRT1","CI.lower.NRT2","CI.upper.NRT2","SE_MedEff1","SE_MedEff2")
  }else{
    final_results <- test1
    names(final_results)[6] <- c("timeseq")
    
  }
  
  if(plot == TRUE){
    
    # First Plot: plotting alpha1 coefficients across time using ggplot
    plot1_a1 <- ggplot(data = test1, aes(t.est,hat.alpha1)) +
                geom_line(color = "red", size = 0.75) +
                labs(title = "Plotting the alpha1 coffecients",
                     x = "Time Sequence",
                     y = "Alpha1")
    
    # Second plot: plotting alpha2 coeffeicients across time using ggplot
    plot2_a2 <- ggplot(data = test1, aes(t.est,hat.alpha2)) +
                geom_line(color = "red", size = 0.75) +
                labs(title = "Plotting the alpha2 coffecients",
                     x = "Time Sequence",
                     y = "Alpha2")
    
    # Third plot: plotting beta3 coeffeicients across time using ggplot
    plot3_b3 <- ggplot(data = test1, aes(t.est,hat.beta3)) +
                geom_line(color = "red", size = 0.75) +
                labs(title = "Plotting the beta3 coffecients",
                     x = "Time Sequence",
                     y = "Beta3")
    
    # Fourth plot: plotting the mediation effect of treatment arm1
    plot4 <- ggplot(data = test1, aes(t.est,hat.mediation1)) +
              geom_line(color = "red", size = 0.75) +
              labs(title = "Plotting the time-varying mediation effect",
                   x = "Time Sequence",
                   y = "Mediation Effect for NRT1")
      
    # Fifth plot: plotting the mediation effect of treatment arm2
    plot5 <- ggplot(data = test1, aes(t.est,hat.mediation2)) +
              geom_line(color = "red", size = 0.75) +
              labs(title = "Plotting the time-varying mediation effect",
                   x = "Time Sequence",
                   y = "Mediation Effect for NRT2")
    
    if(CI == "boot"){
      # Sixth plot: plotting the mediation effect of treatment arm1 with 95% CIs
      plot6 <- ggplot(data = test2, aes(t.est, orig.mediation1)) +
        geom_line(size = 1, color = "red") +
        geom_line(aes(t.est, plw1), size = 0.8, color = "blue", linetype = "dashed") +
        geom_line(aes(t.est, pup1), size = 0.8, color = "blue", linetype = "dashed") +
        geom_line(aes(t.est, 0)) +
        labs(title = "Mediation Effect with 95% CIs (computed with bootstrap)",
             x = "Time Sequence",
             y = "Mediation Effect for NRT1") + 
        theme(legend.position = "none")    
      
      # Seventh plot: plotting the mediation effect of treatment arm2 with 95% CIs
      plot7 <- ggplot(data = test2, aes(t.est, orig.mediation2)) +
        geom_line(size = 1, color = "red") +
        geom_line(aes(t.est, plw2), size = 0.8, color = "blue", linetype = "dashed") +
        geom_line(aes(t.est, pup2), size = 0.8, color = "blue", linetype = "dashed") +
        geom_line(aes(t.est, 0)) +
        labs(title = "Mediation Effect with 95% CIs (computed with bootstrap)",
             x = "Time Sequence",
             y = "Mediation Effect for NRT2") + 
        theme(legend.position = "none")
      
      plot_results <- list("plot1_a1" = plot1_a1,
                           "plot2_a2" = plot2_a2,
                           "plot3_b3" = plot3_b3,
                           "MedEff_NRT1" = plot4,
                           "MedEff_NRT2" = plot5,
                           "MedEff_CI_NRT1" = plot6,
                           "MedEff_CI_NRT2" = plot7)
    }else{
      plot_results <- list("plot1_a1" = plot1_a1,
                           "plot2_a2" = plot2_a2,
                           "plot3_b3" = plot3_b3,
                           "MedEff_NRT1" = plot4,
                           "MedEff_NRT2" = plot5)
    }
}
  
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
                    "plot2_a2" = plot2_a2,
                    "plot3_b3" = plot3_b3,
                    "MedEff_NRT1" = plot4,
                    "MedEff_NRT2" = plot5,
                    "MedEff_CI_NRT1" = plot6,
                    "MedEff_CI_NRT2" = plot7)
  }
  else if(plot == TRUE & CI != "boot"){
    results <- list("Estimates" = final_results,
                    "plot1_a1" = plot1_a1,
                    "plot2_a2" = plot2_a2,
                    "plot3_b3" = plot3_b3,
                    "MedEff_NRT1" = plot4,
                    "MedEff_NRT2" = plot5)  
  }
  else{
    results <- list("Estimates" = final_results)
  }
  
  return(results)
}