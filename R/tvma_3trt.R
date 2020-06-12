#' Time Varying Mediation Function: Continuous Outcome and Three Treatment (Exposure) Groups
#' @export

tvma_3trt <- function(NRT1, NRT2, t.seq, mediator, outcome, t.est = t.seq, plot = FALSE, CI="boot", replicates = 1000, verbose = FALSE)
  {
      # Estimate time-varying mediation effect and bootstrap standard errors
      # for continuous outcome and 3 treatment arms
      #
      ### ARGUMENTS:
      #
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
      ### RETURNS:
      #
      ### Dataframe 'final_results' with the following elements
      #
      #   hat.alpha1        -->   estimated NRT1 effect on mediator
      #   hat.alpha2        -->   estimated NRT2 effect on mediator
      #   hat.beta3         -->   estimated mediation effect on outcome
      #   hat.mediation1    -->   time varying mediation effect - NRT1 on outcome
      #   hat.mediation2    -->   time varying mediation effect - NRT2 on outcome
      #
      ### Optional results based on user argument CI = "boot" passed
      #
      #   CI.upper.NRT1     -->   Upper confidence intervals for NRT1
      #   CI.lower.NRT1     -->   Lower confidence intervals for NRT1
      #   CI.upper.NRT2     -->   Upper confidence intervals for NRT2
      #   CI.lower.NRT2     -->   Lower confidence intervals for NRT2
      #   SE_MedEff1        -->   estimated standard errors of the mediation effect for NRT1
      #   SE_MedEff2        -->   estimated standard errors of the mediation effect for NRT2
      #
      ### Optional plot results
      #
      #   plot1_a1         -->   plot for hat.alpha1 across t.seq with CI
      #   plot2_a2         -->   plot for hat.alpha2 across t.seq with CI
      #   plot3_b3         -->   plot for hat.beta3 across t.seq with CI
      #   MedEff_NRT1      -->   plot for hat.mediation1 across t.seq
      #   MedEff_NRT2      -->   plot for hat.mediation2 across t.seq
      #   MedEff_CI_NRT1   -->   plot for CIs of medEffect for NRT1
      #   MedEff_CI_NRT2   -->   plot for CIs of medEffect for NRT2
      #
      ##
      
      ## Testing the class type of the arguments passed in the function ##
      ctm <- class(mediator)
      cto <- class(outcome)
      ctt1 <- class(NRT1)
      ctt2 <- class(NRT2)
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
      if(is.vector(NRT1) != TRUE || ctt1 != "numeric"){
        print("Error: `NRT1` is not of class type `numeric vector`.")
        flag <- flag + 1
      }
      if(is.vector(NRT2) != TRUE || ctt2 != "numeric"){
        print("Error: `NRT2` is not of class type `numeric vector`.")
        flag <- flag + 1
      }
      if(is.vector(t.seq) != TRUE || cts != "numeric"){
        print("Error: `t.seq` is not of class type `numeric vector`.")
        flag <- flag + 1
      }
      if(flag == 0){
        if(CI == "boot" || CI == "none"){
          results_me <- tvmcurve_3trt(NRT1, NRT2, t.seq, mediator, outcome, t.est)
          alldata <- list(NRT1 = NRT1, NRT2 = NRT2, x = mediator, y = outcome, t.seq = t.seq)
          original.coeff <- list(hat.alpha1 = results_me$hat.alpha1, hat.alpha2 = results_me$hat.alpha2, hat.beta3 = results_me$hat.beta3)
          
          #### Computing CI for coefficients irrespective of user choice of bootstrapping ####
          bootcoeff_CI <- bootci_coeff_3trt(NRT1, NRT2, t.seq, mediator, outcome, t.est, original.coeff, replicates)
          
          test1 <- cbind(as.data.frame(results_me), as.data.frame(bootcoeff_CI), t.est)
          
          #### CI for mediation ####
          if(CI == "boot"){
            results_ci <- bootci_tvm_3trt(replicates, alldata, t.est)
            test2 <- cbind(as.data.frame(results_ci), t.est)
            
            final_dat <- merge(test1, test2, all.x = TRUE)
            final_results <- final_dat %>%
              select(-orig.mediation1, -orig.mediation2)
            names(final_results)[c(1,7:18)] <- c("timeseq",
                                                 "CI.lower.alpha1", "CI.upper.alpha1",
                                                 "CI.lower.alpha2", "CI.upper.alpha2",
                                                 "CI.lower.beta3", "CI.upper.beta3",
                                                 "CI.lower.NRT1", "CI.upper.NRT1",
                                                 "CI.lower.NRT2","CI.upper.NRT2",
                                                 "SE_MedEff1","SE_MedEff2")
          }
          else{
            final_results <- test1
            names(final_results)[6:12] <- c("CI.lower.alpha1", "CI.upper.alpha1",
                                            "CI.lower.alpha2", "CI.upper.alpha2",
                                            "CI.lower.beta3", "CI.upper.beta3",
                                            "timeseq")
          }
          
          #### Plot construction ####
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
            
            # First Plot: plotting alpha1 coefficients across time using ggplot
            plot1_a1 <- ggplot(data = test1, aes(t.est,hat.alpha1)) +
              geom_line(color = "red", size = 0.75) +
              geom_line(aes(t.est, alw1), size = 0.8, color = "blue", linetype = "dashed") +
              geom_line(aes(t.est, aup1), size = 0.8, color = "blue", linetype = "dashed") +
              labs(title = "Plotting the alpha1 coefficients",
                   x = "Time (in days)",
                   y = "Alpha1") +
              scale_x_continuous(breaks = seq(l, u, i))
            
            # Second plot: plotting alpha2 coefficients across time using ggplot
            plot2_a2 <- ggplot(data = test1, aes(t.est,hat.alpha2)) +
              geom_line(color = "red", size = 0.75) +
              geom_line(aes(t.est, alw2), size = 0.8, color = "blue", linetype = "dashed") +
              geom_line(aes(t.est, aup2), size = 0.8, color = "blue", linetype = "dashed") +
              labs(title = "Plotting the alpha2 coefficients",
                   x = "Time (in days)",
                   y = "Alpha2") +
              scale_x_continuous(breaks = seq(l, u, i))
            
            # Third plot: plotting beta3 coefficients across time using ggplot
            plot3_b3 <- ggplot(data = test1, aes(t.est,hat.beta3)) +
              geom_line(color = "red", size = 0.75) +
              geom_line(aes(t.est, blw3), size = 0.8, color = "blue", linetype = "dashed") +
              geom_line(aes(t.est, bup3), size = 0.8, color = "blue", linetype = "dashed") +
              labs(title = "Plotting the beta3 coefficients",
                   x = "Time (in days)",
                   y = "Beta3") +
              scale_x_continuous(breaks = seq(l, u, i))
            
            # Fourth plot: plotting the mediation effect of treatment arm1
            plot4 <- ggplot(data = test1, aes(t.est,hat.mediation1)) +
              geom_line(color = "red", size = 0.75) +
              labs(title = "Plotting the time-varying mediation effect (NRT1)",
                   x = "Time (in days)",
                   y = "Mediation Effect for NRT1") +
              scale_x_continuous(breaks = seq(l, u, i))
            
            # Fifth plot: plotting the mediation effect of treatment arm2
            plot5 <- ggplot(data = test1, aes(t.est,hat.mediation2)) +
              geom_line(color = "red", size = 0.75) +
              labs(title = "Plotting the time-varying mediation effect (NRT2)",
                   x = "Time (in days)",
                   y = "Mediation Effect for NRT2") +
              scale_x_continuous(breaks = seq(l, u, i))
            
            if(CI == "boot"){
              # Sixth plot: plotting the mediation effect of treatment arm1 with 95% CIs
              plot6 <- ggplot(data = test2, aes(t.est, orig.mediation1)) +
                geom_line(size = 1, color = "red") +
                geom_line(aes(t.est, plw1), size = 0.8, color = "blue", linetype = "dashed") +
                geom_line(aes(t.est, pup1), size = 0.8, color = "blue", linetype = "dashed") +
                geom_line(aes(t.est, 0)) +
                labs(title = "Mediation Effect (NRT1) with 95% CIs (computed with percentile bootstrap)",
                     x = "Time (in days)",
                     y = "Mediation Effect for NRT1") + 
                theme(legend.position = "none")  +
                scale_x_continuous(breaks = seq(l, u, i))  
              
              # Seventh plot: plotting the mediation effect of treatment arm2 with 95% CIs
              plot7 <- ggplot(data = test2, aes(t.est, orig.mediation2)) +
                geom_line(size = 1, color = "red") +
                geom_line(aes(t.est, plw2), size = 0.8, color = "blue", linetype = "dashed") +
                geom_line(aes(t.est, pup2), size = 0.8, color = "blue", linetype = "dashed") +
                geom_line(aes(t.est, 0)) +
                labs(title = "Mediation Effect (NRT2) with 95% CIs (computed with percentile bootstrap)",
                     x = "Time (in days)",
                     y = "Mediation Effect for NRT2") + 
                theme(legend.position = "none") +
                scale_x_continuous(breaks = seq(l, u, i))
              
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
          
          #### VERBOSE condition ####
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
          
        }else{
          print(paste("Error:Accepted values for CI are 'boot' and 'none';you have entered an unacceptable value for CI."))
        }
      }else{
        print(paste("tvma_3trt() stopped execution due to unacceptable class type of function argument(s)."))
      }
  }