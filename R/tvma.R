#' Time Varying Mediation Function: Continuous Outcome and Two Treatment Arms (Exposure Groups)
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
#' @return \item{hat.alpha}{estimated main treatment arm (exposure group) of interest effect on mediator (indirect effect component)}
#' @return \item{CI.lower.alpha}{lower limit of confidence intervals for estimated coefficient hat.alpha}
#' @return \item{CI.upper.alpha}{upper limit of confidence intervals for estimated coefficient hat.alpha}
#' @return \item{hat.gamma}{estimated main treatment arm (exposure group) of interest effect on outcome (direct effect)}
#' @return \item{CI.lower.gamma}{lower limit of confidence intervals for estimated coefficient hat.gamma}
#' @return \item{CI.upper.gamma}{upper limit of confidence intervals for estimated coefficient hat.gamma}
#' @return \item{hat.beta}{estimated mediator effect on outcome (indirect effect component)}
#' @return \item{CI.lower.beta}{lower limit of confidence intervals for estimated coefficient hat.beta}
#' @return \item{CI.upper.beta}{upper limit of confidence intervals for estimated coefficient hat.beta}
#' @return \item{hat.tao}{estimated main treatment arm (exposure group) of interest effect on outcome, excluding adjustment for mediator (total effect)}
#' @return \item{CI.lower.tao}{lower limit of confidence intervals for estimated coefficient hat.tao}
#' @return \item{CI.upper.tao}{upper limit of confidence intervals for estimated coefficient hat.tao}
#' @return \item{est.M}{time varying mediation effect - main treatment arm (exposure group) of interest on outcome}
#' @return \item{boot.se.m}{estimated standard error of est.M}
#' @return \item{CI.lower}{lower limit of confidence intervals of est.M}
#' @return \item{CI.upper}{upper limit of confidence intervals of est.M}
#' 
#' @section Plot Returns:
#' \enumerate{
#' \item{\code{Alpha_CI }}{plot for hat.alpha across t.est with CIs}
#' \item{\code{Gamma_CI }}{plot for hat.gamma across t.est with CIs}
#' \item{\code{Beta_CI }}{plot for hat.beta across t.est with CIs}
#' \item{\code{Tao_CI }}{plot for hat.tao across t.est with CIs}
#' \item{\code{MedEff }}{plot for est.M across t.est}
#' \item{\code{MedEff_CI }}{plot for est.M with CIs across t.est}
#' }
#' 
#' @note
#' \enumerate{
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
#' # results <- tvma(trt, t.seq, mediator, outcome)
#' 
#' # COMPUTE TIME VARYING MEDIATION ANALYSIS FOR SPECIFIED POINTS IN TIME USING 500 REPLICATES
#' # results <- tvma(trt, t.seq, mediator, outcome,
#' #                 t.est = c(0.2, 0.4, 0.6, 0.8),
#' #                 replicates = 500)
#' 
#' @references 
#' \enumerate{
#' \item{Fan, J. and Gijbels, I. (1996). Local polynomial modelling and its applications: monographs on statistics and applied probability 66 66. CRC Press.}
#' \item{Fan, J. and Zhang, W. (1999). Statistical estimation in varying coefficient models. The annals of Statistics 27 1491-1518.}
#' \item{Fan, J. and Zhang, W. (2000). Two-step estimation of functional linear models with applications to longitudinal data. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 62 303-322.}
#' \item{Cai, X., Piper, M. E., Li, R., & Coffman, D. L. (under review). Estimation and inference for the mediation effect in a time-varying mediation model. Journal of the Royal Statistics Society, Series C (Applied Statistics).}
#' \item{Baker, T. B., Piper, M. E., Stein, J. H., Smith, S. S., Bolt, D. M., Fraser, D. L., & Fiore, M. C. (2016). Effects of nicotine patch vs varenicline vs combination nicotine replacement therapy on smoking cessation at 26 weeks: A randomized clinical trial. JAMA, 315(4), 371-379. doi:10.1001/jama.2015.19284.}
#' \item{Efron, B.; Tibshirani, R. Bootstrap Methods for Standard Errors, Confidence Intervals, and Other Measures of Statistical Accuracy. Statist. Sci. 1 (1986), no. 1, 54--75. doi:10.1214/ss/1177013815.}
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
      #   hat.alpha        -->   estimated treatment effect on mediator (indirect effect component)
      #   CI.lower.alpha   -->   lower limit of confidence intervals for coefficient hat.alpha
      #   CI.upper.alpha   -->   upper limit of confidence intervals for coefficient hat.alpha
      #   hat.gamma        -->   estimated treatment effect on outcome (direct effect)
      #   CI.lower.gamma   -->   lower limit of confidence intervals for coefficient hat.gamma
      #   CI.upper.gamma   -->   upper limit of confidence intervals for coefficient hat.gamma
      #   hat.beta         -->   estimated mediation effect on outcome (indirect effect component)
      #   CI.lower.beta    -->   lower limit of confidence intervals for coefficient hat.beta
      #   CI.upper.beta    -->   upper limit of confidence intervals for coefficient hat.beta
      #   hat.tao          -->   estimated treatment effect on outcome, excluding adjustment for the mediator (total effect)
      #   CI.lower.tao     -->   lower limit of confidence intervals for coefficient hat.tao
      #   CI.upper.tao     -->   upper limit of confidence intervals for coefficient hat.tao
      #   est.M            -->   time varying mediation effect
  
      # 
      # Optional results based on user argument CI = "boot" passed
      #   boot.se.m         -->   estimated standard error of est.M
      #   CI.lower          -->   lower limit of confidence intervals of est.M
      #   CI.upper          -->   upper limit of confidence intervals of est.M
      #  
      ### Optional plot results
      #
      #   Alpha_CI          -->   plot for hat.alpha across t.est with CI
      #   Gamma_CI          -->   plot for hat.gamma across t.est with CI
      #   Beta_CI           -->   plot for hat.beta across t.est with CI
      #   Tao_CI            -->   plot for hat.tao across t.est with CI
      #   MedEff            -->   plot for est.M across t.est
      #   MedEff_CI         -->   plot for est.M with CIs across t.est
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
        
        # Estimate coefficients alpha, beta and gamma
        coeff.est <- estCoeff(newMO.j.est)
        
        # Steps to derive the total effect coefficient tao
        X.new <- scale(treatment, center = TRUE, scale = FALSE)
        Y.new <- scale(outcome[j - 1, ], center = TRUE, scale = FALSE)
        nomissing.X <- complete.cases(X.new)
        nomissing.Y <- complete.cases(Y.new)
        nomissing.index <- nomissing.X * nomissing.Y
        
        X.new <- X.new[which(nomissing.index == 1),]
        Y.new <- Y.new[which(nomissing.index == 1)]
        sym_newMO <- t(X.new)%*%(X.new)
        coeff.tao <- solve(sym_newMO)%*%t(X.new)%*%(Y.new)
        
        # Store the coefficients
        coeff.all <- rbind(coeff.est, coeff.tao)
        t.coeff <- cbind(t.coeff, coeff.all)  # store coeff estimates at t.seq
      }
      
      # EQUATIONS 4 & 5
      est.smooth <- smoothest(t.seq, t.coeff, t.est, deltat)
      
      ## calculating the CI for the coefficients alpha and beta
      coeff_CI_2trt <- bootci_coeff_2trt(treatment, t.seq, mediator, outcome, t.est, deltat, replicates)
      
      test1 <- cbind(as.data.frame(est.smooth), as.data.frame(coeff_CI_2trt), t.est)
      
      # CALCULATE CONFIDENCE INTERVALS
      if(CI == "boot"){
        results_ci <- estBootCIs(treatment, t.seq, mediator, outcome, t.est, deltat, replicates)
        test2 <- cbind(as.data.frame(results_ci), t.est)
        
        final_dat <- merge(test1, test2, all.x = TRUE)
        final_results <- final_dat %>%
          select(-bw_alpha, -bw_gamma, -bw_beta, -bw_tao)
        final_results <- final_results[c("t.est","hat.alpha","CI.lower.alpha","CI.upper.alpha",
                                         "hat.gamma", "CI.lower.gamma", "CI.upper.gamma",
                                         "hat.beta", "CI.lower.beta", "CI.upper.beta",
                                         "hat.tao", "CI.lower.tao", "CI.upper.tao",
                                         "est.M", "boot.se", "CI.lower", "CI.upper")]
        
        names(final_results)[15] <- c("boot.se.m")
      }
      else{
        final_results <- test1 %>%
          select(-bw_alpha1, -bw_beta1, -bw_beta2, -bw_tao)
        final_results <- final_results[c("t.est","hat.alpha","CI.lower.alpha","CI.upper.alpha",
                                         "hat.gamma", "CI.lower.gamma", "CI.upper.gamma",
                                         "hat.beta", "CI.lower.beta", "CI.upper.beta",
                                         "hat.tao", "CI.lower.tao", "CI.upper.tao",
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
        
        # First Plot: plotting alpha coefficients across time using ggplot
        plot1_a <- ggplot(data = final_results, aes(t.est,hat.alpha)) +
          geom_line(color = "red", size = 0.75) +
          geom_line(aes(t.est, CI.lower.alpha), size = 0.8, color = "blue", linetype = "dashed") +
          geom_line(aes(t.est, CI.upper.alpha), size = 0.8, color = "blue", linetype = "dashed") +
          labs(title = "Plotting the alpha coefficients",
               x = "Time Sequence",
               y = "Alpha") +
          scale_x_continuous(breaks = seq(l, u, i))
        
        # Second plot: plotting gamma coefficients across time using ggplot
        plot2_g <- ggplot(data = final_results, aes(t.est,hat.gamma)) +
          geom_line(color = "red", size = 0.75) +
          geom_line(aes(t.est, CI.lower.gamma), size = 0.8, color = "blue", linetype = "dashed") +
          geom_line(aes(t.est, CI.upper.gamma), size = 0.8, color = "blue", linetype = "dashed") +
          labs(title = "Plotting the gamma coefficients",
               x = "Time Sequence",
               y = "Gamma") +
          scale_x_continuous(breaks = seq(l, u, i))
        
        # Third plot: plotting beta coefficients across time using ggplot
        plot3_b <- ggplot(data = final_results, aes(t.est,hat.beta)) +
          geom_line(color = "red", size = 0.75) +
          geom_line(aes(t.est, CI.lower.beta), size = 0.8, color = "blue", linetype = "dashed") +
          geom_line(aes(t.est, CI.upper.beta), size = 0.8, color = "blue", linetype = "dashed") +
          labs(title = "Plotting the beta coefficients",
               x = "Time Sequence",
               y = "Beta") +
          scale_x_continuous(breaks = seq(l, u, i))
        
        # Fourth plot: plotting tao coefficients across time using ggplot
        plot4_t <- ggplot(data = final_results, aes(t.est,hat.tao)) +
          geom_line(color = "red", size = 0.75) +
          geom_line(aes(t.est, CI.lower.tao), size = 0.8, color = "blue", linetype = "dashed") +
          geom_line(aes(t.est, CI.upper.tao), size = 0.8, color = "blue", linetype = "dashed") +
          labs(title = "Plotting the tao coefficients",
               x = "Time Sequence",
               y = "Tao") +
          scale_x_continuous(breaks = seq(l, u, i))
        
        # Fifth plot: plotting the mediation effect of treatment arm
        plot5 <- ggplot(data = final_results, aes(t.est,est.M)) +
          geom_line(color = "red", size = 0.75) +
          labs(title = "Plotting the time-varying mediation effect",
               x = "Time Sequence",
               y = "Mediation Effect") +
          scale_x_continuous(breaks = seq(l, u, i))
        
        if(CI == "boot"){
          
          # Sixth plot: plotting the mediation effect of treatment with 95% CIs
          plot6 <- ggplot(data = final_results, aes(t.est, est.M)) +
            geom_line(size = 1, color = "red") +
            geom_line(aes(t.est, CI.lower), size = 0.8, color = "blue", linetype = "dashed") +
            geom_line(aes(t.est, CI.upper), size = 0.8, color = "blue", linetype = "dashed") +
            # geom_line(aes(t.est, 0)) +
            labs(title = "Mediation Effect with 95% CIs (computed with percentile bootstrap)",
                 x = "Time Sequence",
                 y = "Mediation Effect") + 
            theme(legend.position = "none") +
            scale_x_continuous(breaks = seq(l, u, i))   
          
          plot_results <- list("Alpha_CI" = plot1_a,
                               "Gamma_CI" = plot2_g,
                               "Beta_CI" = plot3_b,
                               "Tao_CI" = plot4_t,
                               "MedEff" = plot5,
                               "MedEff_CI" = plot6)
        }else{
          plot_results <- list("Alpha_CI" = plot1_a,
                               "Gamma_CI" = plot2_g,
                               "Beta_CI" = plot3_b,
                               "Tao_CI" = plot4_t,
                               "MedEff" = plot5)
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
                        "Alpha_CI" = plot1_a,
                        "Gamma_CI" = plot2_g,
                        "Beta_CI" = plot3_b,
                        "Tao_CI" = plot4_t,
                        "MedEff" = plot5,
                        "MedEff_CI" = plot6)
      }
      else if(plot == TRUE & CI != "boot"){
        results <- list("Estimates" = final_results,
                        "Alpha_CI" = plot1_a,
                        "Gamma_CI" = plot2_g,
                        "Beta_CI" = plot3_b,
                        "Tao_CI" = plot4_t,
                        "MedEff" = plot5)  
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