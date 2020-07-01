#' Time Varying Mediation Analysis: Continuous Outcome and Two Treatment Arms (Exposure Groups)
#' 
#' Function to estimate the time-varying mediation effect and bootstrap standard errors, involving two treatment groups and continuous outcome.
#' 
#' @param treatment   a vector with treatment values
#' @param t.seq       a vector of time points for each observation
#' @param mediator    matrix of mediator values in wide format
#' @param outcome     matrix of outcome outcomes in wide format
#' @param t.est       a vector of time points to make the estimation                  Default = t.seq (OPTIONAL ARGUMENT)
#' @param plot        TRUE or FALSE for plotting mediation effect                     Default = "FALSE" (OPTIONAL ARGUMENT)
#' @param CI          "none" or "boot" method of deriving confidence intervals.       Default = "boot" (OPTIONAL ARGUMENT)
#' @param replicates  Number of replicates for bootstrapping confidence intervals.    Default = 1000 (OPTIONAL ARGUMENT)
#' @param verbose     TRUE or FALSE for printing results to screen.                   Default = "FALSE" (OPTIONAL ARGUMENT)
#' 
#' @return \item{hat.alpha.1}{estimated main treatment arm (exposure group) of interest effect on mediator (indirect effect component)}
#' @return \item{CI.lower.alpha1}{lower limit of confidence intervals for estimated coefficient hat.alpha.1}
#' @return \item{CI.upper.alpha1}{upper limit of confidence intervals for estimated coefficient hat.alpha.1}
#' @return \item{hat.beta.1}{estimated main treatment arm (exposure group) of interest effect on outcome (direct effect)}
#' @return \item{CI.lower.beta1}{lower limit of confidence intervals for estimated coefficient hat.beta.1}
#' @return \item{CI.upper.beta1}{upper limit of confidence intervals for estimated coefficient hat.beta.1}
#' @return \item{hat.beta.2}{estimated mediator effect on outcome (indirect effect component)}
#' @return \item{CI.lower.beta2}{lower limit of confidence intervals for estimated coefficient hat.beta.2}
#' @return \item{CI.upper.beta2}{upper limit of confidence intervals for estimated coefficient hat.beta.2}
#' @return \item{est.M}{time varying mediation effect - main treatment arm (exposure group) of interest on outcome}
#' @return \item{boot.se.m}{estimated standard error of est.M}
#' @return \item{CI.lower}{lower limit of confidence intervals of est.M}
#' @return \item{CI.upper}{upper limit of confidence intervals of est.M}
#' 
#' @section Plot Returns:
#' \enumerate{
#' \item{\code{Alpha1_CI }}{plot for hat.alpha.1 across t.est with CIs}
#' \item{\code{Beta1_CI }}{plot for hat.beta.1 across t.est with CIs}
#' \item{\code{Beta2_CI }}{plot for hat.beta.2 across t.est with CIs}
#' \item{\code{MedEff }}{plot for est.M across t.est}
#' \item{\code{MedEff_CI }}{plot for est.M with CIs across t.est}
#' }
#' 
#' @note
#' \enumerate{
#' \item{Currently supports 2 treatment options. Future releases may expand number of treatment options.}
#' \item{** IMPORTANT ** An alternate way of formatting the data and calling the function is documented in detail in the function tutorial for tvmb().}
#' }
#' 
#' @examples
#' data(smoker)
#' 
#' # REDUCE DATA SET TO ONLY 2 TREATMENT CONDITIONS (EXCLUDING COMBINATION NRT)
#' smoker.sub <- smoker[smoker$treatment != 4, ]
#' 
#' # GENERATE WIDE FORMATTED MEDIATORS
#' mediator <- LongToWide(smoker.sub$SubjectID,
#'                        smoker.sub$timeseq,
#'                        smoker.sub$NegMoodLst15min)
#' 
#' # GENERATE WIDE FORMATTED OUTCOMES
#' outcome <- LongToWide(smoker.sub$SubjectID,
#'                       smoker.sub$timeseq,
#'                       smoker.sub$cessFatig)
#' 
#' # GENERATE A BINARY TREATMENT VARIABLE
#' trt <- as.numeric(unique(smoker.sub[,c("SubjectID","varenicline")])[,2])-1
#' 
#' # GENERATE A VECTOR OF UNIQUE TIME POINTS
#' t.seq <- sort(unique(smoker.sub$timeseq))
#' 
#' # COMPUTE TIME VARYING MEDIATION ANALYSIS USING BOOTSTRAPPED CONFIDENCE INTERVALS
#' results <- tvma(trt, t.seq, mediator, outcome)
#' 
#' # COMPUTE TIME VARYING MEDIATION ANALYSIS FOR SPECIFIED POINTS IN TIME USING 500 REPLICATES
#' results <- tvma(trt, t.seq, mediator, outcome,
#'                 t.est = c(0.2, 0.4, 0.6, 0.8),
#'                 replicates = 500)
#' 
#' @references 
#' \enumerate{
#' \item{Fan, J. and Gijbels, I. (1996). Local polynomial modelling and its applications: monographs on statistics and applied probability 66 66. CRC Press.}
#' \item{Fan, J. and Zhang, W. (1999). Statistical estimation in varying coefficient models. The annals of Statistics 27 1491-1518.}
#' \item{Fan, J. and Zhang, W. (2000). Two-step estimation of functional linear models with applications to longitudinal data. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 62 303-322.}
#' }
#' 
#' @export
#' @importFrom stats complete.cases cov glm lm loess na.omit predict quantile sd var
#' @import dplyr
#' @import ggplot2
#' @import kedd
#' @import locpol

tvma <- function(treatment, t.seq, mediator, outcome, t.est = t.seq, plot = FALSE, CI="boot", replicates = 1000, verbose = FALSE)
  {
      # Estimate time-varying mediation effect and bootstrap standard errors
      #
      ### ARGUMENTS:
      #   treatment   -->   a vector with treatment values
      #   t.seq       -->   a vector of time points for each observation
      #   mediator    -->   matrix of mediator values in wide format
      #   outcome     -->   matrix of outcome outcomes in wide format
      #   t.est       -->   a vector of time points to make the estimation                  Default = t.seq
      #   plot        -->   TRUE or FALSE for plotting mediation effect                     Default = "FALSE"
      #   CI          -->   "none" or "boot" method of deriving confidence intervals.       Default = "boot"
      #   replicates  -->   Number of replicates for bootstrapping confidence intervals.    Default = 1000
      #   verbose     -->   TRUE or FALSE for printing results to screen.                   Default = "FALSE"
      #
      # RETURNS:
      #
      ### Dataframe 'final_results' with the following elements
      #
      #   hat.alpha.1       -->   estimated treatment effect on mediator (indirect effect component)
      #   CI.lower.alpha1   -->   lower limit of confidence intervals for coefficient hat.alpha.1
      #   CI.upper.alpha1   -->   upper limit of confidence intervals for coefficient hat.alpha.1
      #   hat.beta.1        -->   estimated treatment effect on outcome (direct effect)
      #   CI.lower.beta1    -->   lower limit of confidence intervals for coefficient hat.beta.1
      #   CI.upper.beta1    -->   upper limit of confidence intervals for coefficient hat.beta.1
      #   hat.beta.2        -->   estimated mediation effect on outcome (indirect effect component)
      #   CI.lower.beta2    -->   lower limit of confidence intervals for coefficient hat.beta.2
      #   CI.upper.beta2    -->   upper limit of confidence intervals for coefficient hat.beta.2
      #   est.M             -->   time varying mediation effect
  
      # 
      # Optional results based on user argument CI = "boot" passed
      #   boot.se.m         -->   estimated standard error of est.M
      #   CI.lower          -->   lower limit of confidence intervals of est.M
      #   CI.upper          -->   upper limit of confidence intervals of est.M
      #  
      ### Optional plot results
      #
      #   Alpha1_CI          -->   plot for hat.alpha.1 across t.est with CI
      #   Beta1_CI           -->   plot for hat.beta.1 across t.est with CI
      #   Beta2_CI           -->   plot for hat.beta.2 across t.est with CI
      #   MedEff             -->   plot for est.M across t.est
      #   MedEff_CI          -->   plot for est.M with CIs across t.est
      #
      ##
  
  ## Testing the class type of the arguments passed in the function ##
  ctm <- class(mediator)
  cto <- class(outcome)
  ctt <- class(treatment)
  cts <- class(t.seq)
  flag <- 0
  
  if(ctm[1] != "matrix"){
    print("Error: `mediator` is not of class type `matrix`.")
    flag <- flag + 1
  }
  if(cto[1] != "matrix"){
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
      
      index = vector()  
      index=which(!is.na(treatment))
      treatment = treatment[index]
      outcome = outcome[,index]
      mediator = mediator[,index]
      
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
          select(-bw_alpha1, -bw_beta1, -bw_beta2)
        final_results <- final_results[c("t.est","hat.alpha.1","CI.lower.alpha1","CI.upper.alpha1",
                                         "hat.beta.1", "CI.lower.beta1", "CI.upper.beta1",
                                         "hat.beta.2", "CI.lower.beta2", "CI.upper.beta2",
                                         "est.M", "boot.se", "CI.lower", "CI.upper")]
        
        names(final_results)[12] <- c("boot.se.m")
      }
      else{
        final_results <- test1 %>%
          select(-bw_alpha1, -bw_beta1, -bw_beta2)
        final_results <- final_results[c("t.est","hat.alpha.1","CI.lower.alpha1","CI.upper.alpha1",
                                         "hat.beta.1", "CI.lower.beta1", "CI.upper.beta1",
                                         "hat.beta.2", "CI.lower.beta2", "CI.upper.beta2",
                                         "est.M")]
      }
      
     
      #### Plot construction ####
      if(plot == TRUE){
        
        lt <- length(final_results$t.est)
        l <- min(final_results$t.est)
        u <- max(final_results$t.est)
        
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
        plot1_a1 <- ggplot(data = final_results, aes(t.est,hat.alpha.1)) +
          geom_line(color = "red", size = 0.75) +
          geom_line(aes(t.est, CI.lower.alpha1), size = 0.8, color = "blue", linetype = "dashed") +
          geom_line(aes(t.est, CI.upper.alpha1), size = 0.8, color = "blue", linetype = "dashed") +
          labs(title = "Plotting the alpha1 coefficients",
               x = "Time Sequence",
               y = "Alpha1") +
          scale_x_continuous(breaks = seq(l, u, i))
        
        # Second plot: plotting beta1 coefficients across time using ggplot
        plot2_b1 <- ggplot(data = final_results, aes(t.est,hat.beta.1)) +
          geom_line(color = "red", size = 0.75) +
          geom_line(aes(t.est, CI.lower.beta1), size = 0.8, color = "blue", linetype = "dashed") +
          geom_line(aes(t.est, CI.upper.beta1), size = 0.8, color = "blue", linetype = "dashed") +
          labs(title = "Plotting the beta1 coefficients",
               x = "Time Sequence",
               y = "Beta1") +
          scale_x_continuous(breaks = seq(l, u, i))
        
        # Third plot: plotting beta2 coefficients across time using ggplot
        plot3_b2 <- ggplot(data = final_results, aes(t.est,hat.beta.2)) +
          geom_line(color = "red", size = 0.75) +
          geom_line(aes(t.est, CI.lower.beta2), size = 0.8, color = "blue", linetype = "dashed") +
          geom_line(aes(t.est, CI.upper.beta2), size = 0.8, color = "blue", linetype = "dashed") +
          labs(title = "Plotting the beta2 coefficients",
               x = "Time Sequence",
               y = "Beta2") +
          scale_x_continuous(breaks = seq(l, u, i))
        
        # Fourth plot: plotting the mediation effect of treatment arm
        plot4 <- ggplot(data = final_results, aes(t.est,est.M)) +
          geom_line(color = "red", size = 0.75) +
          labs(title = "Plotting the time-varying mediation effect",
               x = "Time Sequence",
               y = "Mediation Effect") +
          scale_x_continuous(breaks = seq(l, u, i))
        
        if(CI == "boot"){
          
          # Fifth plot: plotting the mediation effect of treatment with 95% CIs
          plot5 <- ggplot(data = final_results, aes(t.est, est.M)) +
            geom_line(size = 1, color = "red") +
            geom_line(aes(t.est, CI.lower), size = 0.8, color = "blue", linetype = "dashed") +
            geom_line(aes(t.est, CI.upper), size = 0.8, color = "blue", linetype = "dashed") +
            # geom_line(aes(t.est, 0)) +
            labs(title = "Mediation Effect with 95% CIs (computed with percentile bootstrap)",
                 x = "Time Sequence",
                 y = "Mediation Effect") + 
            theme(legend.position = "none") +
            scale_x_continuous(breaks = seq(l, u, i))   
          
          plot_results <- list("Alpha1_CI" = plot1_a1,
                               "Beta1_CI" = plot2_b1,
                               "Beta2_CI" = plot3_b2,
                               "MedEff" = plot4,
                               "MedEff_CI" = plot5)
        }else{
          plot_results <- list("Alpha1_CI" = plot1_a1,
                               "Beta1_CI" = plot2_b1,
                               "Beta2_CI" = plot3_b2,
                               "MedEff" = plot4)
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
                        "Alpha1_CI" = plot1_a1,
                        "Beta1_CI" = plot2_b1,
                        "Beta2_CI" = plot3_b2,
                        "MedEff" = plot4,
                        "MedEff_CI" = plot5)
      }
      else if(plot == TRUE & CI != "boot"){
        results <- list("Estimates" = final_results,
                        "Alpha1_CI" = plot1_a1,
                        "Beta1_CI" = plot2_b1,
                        "Beta2_CI" = plot3_b2,
                        "MedEff" = plot4)  
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