#' Time Varying Mediation Function: Continuous Outcome and Two Treatment Groups
#' 
#' Function to estimate the time-varying mediation effect and bootstrap standard 
#' errors for two treatment groups and a continuous outcome.
#' 
#' @param treatment   a vector indicating treatment group
#' @param t.seq       a vector of time points for each observation
#' @param mediator    matrix of mediator values in wide format
#' @param outcome     matrix of outcome values in wide format
#' @param t.est       a vector of time points at which to estimate. 
#'                    Default = t.seq (OPTIONAL ARGUMENT)
#' @param plot        TRUE or FALSE for producing plots. Default = "FALSE" 
#'                    (OPTIONAL ARGUMENT)
#' @param CI          "none" or "boot" method of deriving confidence intervals. 
#'                    Default = "boot" (OPTIONAL ARGUMENT)
#' @param replicates  number of replicates for bootstrapping confidence intervals. 
#'                    Default = 1000 (OPTIONAL ARGUMENT)
#' @param verbose     TRUE or FALSE for printing results to screen. 
#'                    Default = "FALSE" (OPTIONAL ARGUMENT)
#' 
#' @return \item{hat.alpha}{estimated time-varying treatment effect on mediator}
#' @return \item{CI.lower.alpha}{CI lower limit for estimated coefficient hat.alpha}
#' @return \item{CI.upper.alpha}{CI upper limit for estimated coefficient hat.alpha}
#' @return \item{hat.gamma}{estimated time-varying treatment effect on outcome (direct effect)}
#' @return \item{CI.lower.gamma}{CI lower limit for estimated coefficient hat.gamma}
#' @return \item{CI.upper.gamma}{CI upper limit for estimated coefficient hat.gamma}
#' @return \item{hat.beta}{estimated time-varying effect of the mediator on outcome}
#' @return \item{CI.lower.beta}{CI lower limit for estimated coefficient hat.beta}
#' @return \item{CI.upper.beta}{CI upper limit for estimated coefficient hat.beta}
#' @return \item{hat.tau}{estimated time-varying treatment effect on outcome (total effect)}
#' @return \item{CI.lower.tau}{CI lower limit for estimated coefficient hat.tau}
#' @return \item{CI.upper.tau}{CI upper limit for estimated coefficient hat.tau}
#' @return \item{est.M}{time varying mediation effect}
#' @return \item{boot.se.m}{estimated standard error for est.M}
#' @return \item{CI.lower}{CI lower limit for est.M}
#' @return \item{CI.upper}{CI upper limit for est.M}
#' 
#' @section Plot Returns:
#' \enumerate{
#' \item{\code{Alpha_CI }}{plot for hat.alpha with CIs over t.est}
#' \item{\code{Gamma_CI }}{plot for hat.gamma with CIs over t.est}
#' \item{\code{Beta_CI }}{plot for hat.beta with CIs over t.est}
#' \item{\code{Tau_CI }}{plot for hat.tau with CIs over t.est}
#' \item{\code{MedEff }}{plot for est.M over t.est}
#' \item{\code{MedEff_CI }}{plot for est.M with CIs over t.est}
#' }
#' 
#' @note
#' \enumerate{
#' \item{** IMPORTANT ** An alternate way of formatting the data and calling the 
#'     function is documented in detail in the tutorial for the tvmb() function.}
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
#' # COMPUTE TIME VARYING MEDIATION ANALYSIS FOR SPECIFIED POINTS IN TIME USING 250 REPLICATES
#' # results <- tvma(trt, t.seq, mediator, outcome,
#' #                 t.est = c(0.2, 0.4, 0.6, 0.8),
#' #                 replicates = 250)
#' 
#' @references 
#' \enumerate{
#' \item{Fan, J. and Gijbels, I. (1996). Local polynomial modelling and its 
#'      applications: Monographs on statistics and applied probability 66. 
#'      CRC Press.}
#' \item{Fan, J. and Zhang, W. (1999). Statistical estimation in varying 
#'       coefficient models. The Annals of Statistics, 27, 1491-1518.}
#' \item{Fan, J. and Zhang, W. (2000). Two-step estimation of functional linear 
#'       models with applications to longitudinal data. Journal of the Royal 
#'       Statistical Society: Series B (Statistical Methodology), 62, 303-322.}
#' \item{Cai, X., Piper, M. E., Li, R., & Coffman, D. L. (2020). 
#'       Estimation and inference for the mediation effect in a time-varying 
#'       mediation model. <https://arxiv.org/abs/2008.11797>}
#' \item{Baker, T. B., Piper, M. E., Stein, J. H., Smith, S. S., Bolt, D. M., 
#'       Fraser, D. L., & Fiore, M. C. (2016). Effects of nicotine patch vs 
#'       varenicline vs combination nicotine replacement therapy on smoking 
#'       cessation at 26 weeks: A randomized clinical trial. JAMA, 315(4), 
#'       371-379.}
#' \item{Efron, B. and Tibshirani, R. (1986). Bootstrap methods for standard 
#'       errors, confidence intervals, and other measures of statistical accuracy. 
#'       Statistical Science, 1, 54-75.}
#' }
#' 
#' @export
#' @importFrom stats complete.cases cov glm lm loess na.omit predict quantile sd var
#' @import dplyr
#' @import ggplot2
#' @import kedd
#' @import locpol
#' @import utf8
 
tvma <- function(treatment, t.seq, mediator, outcome, t.est = t.seq, plot = FALSE, CI="boot", replicates = 1000, verbose = FALSE)
  {
  ## Testing the class of the arguments passed in the function
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
      
      index <- vector()  
      index <- which(!is.na(treatment))
      treatment <- treatment[index]
      outcome <- outcome[,index]
      mediator <- mediator[,index]
      
      deltat <- max(diff(t.seq))/2  # half the time between two measures
      N <- length(treatment)
      t.coeff <- NULL
      
      for (j in 2:length(t.seq)){
        # create empty vector, store mediator
        temp.mediator.j  <- NULL
        temp.mediator.j  <- cbind(mediator[j-1, ], mediator[j, ])
        
        # Derive centered Mediators and Outcomes
        newMO.j.est <- newMediatorOutcome(treatment, temp.mediator.j, outcome[j-1,])
        
        # Estimate coefficients alpha, beta and gamma
        coeff.est <- estCoeff(newMO.j.est)
        
        # Estimate total effect coefficient
        X.new <- scale(treatment, center = TRUE, scale = FALSE)
        Y.new <- scale(outcome[j-1,], center = TRUE, scale = FALSE)
        nomissing.X <- complete.cases(X.new)
        nomissing.Y <- complete.cases(Y.new)
        nomissing.index <- nomissing.X*nomissing.Y
        
        X.new <- X.new[which(nomissing.index == 1),]
        Y.new <- Y.new[which(nomissing.index == 1)]
        sym_newMO <- t(X.new)%*%(X.new)
        coeff.tau <- solve(sym_newMO)%*%t(X.new)%*%(Y.new)
        
        # Store the coefficients
        coeff.all <- rbind(coeff.est, coeff.tau)
        t.coeff <- cbind(t.coeff, coeff.all)  # store coeff estimates at t.seq
      }
      
      est.smooth <- smoothest(t.seq, t.coeff, t.est, deltat)
      
      ## CIs for the coefficients alpha and beta
      coeff_CI_2trt <- bootci_coeff_2trt(treatment, t.seq, mediator, outcome, t.est, deltat, replicates)
      
      test1 <- cbind(as.data.frame(est.smooth), as.data.frame(coeff_CI_2trt), t.est)
      
      # CONFIDENCE INTERVALS for mediated effect
      if(CI == "boot"){
        results_ci <- estBootCIs(treatment, t.seq, mediator, outcome, t.est, deltat, replicates)
        test2 <- cbind(as.data.frame(results_ci), t.est)
        
        final_dat <- merge(test1, test2, all.x = TRUE)
        final_results <- final_dat %>%
          select(-bw_alpha, -bw_gamma, -bw_beta, -bw_tau)
        final_results <- final_results[c("t.est","hat.alpha","CI.lower.alpha","CI.upper.alpha",
                                         "hat.gamma", "CI.lower.gamma", "CI.upper.gamma",
                                         "hat.beta", "CI.lower.beta", "CI.upper.beta",
                                         "hat.tau", "CI.lower.tau", "CI.upper.tau",
                                         "est.M", "boot.se", "CI.lower", "CI.upper")]
        
        names(final_results)[15] <- c("boot.se.m")
      }
      else{
        final_results <- test1 %>%
          select(-bw_alpha1, -bw_beta1, -bw_beta2, -bw_tau)
        final_results <- final_results[c("t.est","hat.alpha","CI.lower.alpha","CI.upper.alpha",
                                         "hat.gamma", "CI.lower.gamma", "CI.upper.gamma",
                                         "hat.beta", "CI.lower.beta", "CI.upper.beta",
                                         "hat.tau", "CI.lower.tau", "CI.upper.tau",
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
        
        # First Plot: alpha coefficients over time
        plot1_a <- ggplot(data = final_results, aes(t.est, hat.alpha)) +
          geom_line(color = "red", size = 0.75) +
          geom_line(aes(t.est, CI.lower.alpha), size = 0.8, color = "blue", linetype = "dashed") +
          geom_line(aes(t.est, CI.upper.alpha), size = 0.8, color = "blue", linetype = "dashed") +
          labs(title = "Plot of the time-varying treatment effect on the mediator",
               x = "Time Sequence",
               y = "Alpha coefficient") +
          scale_x_continuous(breaks = seq(l, u, i))
        
        # Second plot: gamma coefficient over time
        plot2_g <- ggplot(data = final_results, aes(t.est, hat.gamma)) +
          geom_line(color = "red", size = 0.75) +
          geom_line(aes(t.est, CI.lower.gamma), size = 0.8, color = "blue", linetype = "dashed") +
          geom_line(aes(t.est, CI.upper.gamma), size = 0.8, color = "blue", linetype = "dashed") +
          labs(title = "Plot of the time-varying direct effect",
               x = "Time Sequence",
               y = "Gamma coefficient") +
          scale_x_continuous(breaks = seq(l, u, i))
        
        # Third plot: beta coefficient over time
        plot3_b <- ggplot(data = final_results, aes(t.est, hat.beta)) +
          geom_line(color = "red", size = 0.75) +
          geom_line(aes(t.est, CI.lower.beta), size = 0.8, color = "blue", linetype = "dashed") +
          geom_line(aes(t.est, CI.upper.beta), size = 0.8, color = "blue", linetype = "dashed") +
          labs(title = "Plot of the time-varying effect of the mediator on the outcome",
               x = "Time Sequence",
               y = "Beta coefficient") +
          scale_x_continuous(breaks = seq(l, u, i))
        
        # Fourth plot: tau coefficient over time
        plot4_t <- ggplot(data = final_results, aes(t.est, hat.tau)) +
          geom_line(color = "red", size = 0.75) +
          geom_line(aes(t.est, CI.lower.tau), size = 0.8, color = "blue", linetype = "dashed") +
          geom_line(aes(t.est, CI.upper.tau), size = 0.8, color = "blue", linetype = "dashed") +
          labs(title = "Plot of the time-varying total effect",
               x = "Time Sequence",
               y = "Tau cofficient") +
          scale_x_continuous(breaks = seq(l, u, i))
        
        # Fifth plot: mediation effect over time
        plot5 <- ggplot(data = final_results, aes(t.est, est.M)) +
          geom_line(color = "red", size = 0.75) +
          labs(title = "Plot of the time-varying mediation effect",
               x = "Time Sequence",
               y = "Mediation effect") +
          scale_x_continuous(breaks = seq(l, u, i))
        
        if(CI == "boot"){
          
          # Sixth plot: mediation effect with 95% CIs over time
          plot6 <- ggplot(data = final_results, aes(t.est, est.M)) +
            geom_line(size = 1, color = "red") +
            geom_line(aes(t.est, CI.lower), size = 0.8, color = "blue", linetype = "dashed") +
            geom_line(aes(t.est, CI.upper), size = 0.8, color = "blue", linetype = "dashed") +
            # geom_line(aes(t.est, 0)) +
            labs(title = "Time-varying mediation effect with 95% percentile bootstrap CIs",
                 x = "Time Sequence",
                 y = "Mediation effect") + 
            theme(legend.position = "none") +
            scale_x_continuous(breaks = seq(l, u, i))   
          
          plot_results <- list("Alpha_CI" = plot1_a,
                               "Gamma_CI" = plot2_g,
                               "Beta_CI" = plot3_b,
                               "Tau_CI" = plot4_t,
                               "MedEff" = plot5,
                               "MedEff_CI" = plot6)
        }else{
          plot_results <- list("Alpha_CI" = plot1_a,
                               "Gamma_CI" = plot2_g,
                               "Beta_CI" = plot3_b,
                               "Tau_CI" = plot4_t,
                               "MedEff" = plot5)
        }
      }
      
      ## Print results to screen
      if(verbose == TRUE){
        print("Time Varying Mediation Results:")
        if(plot == TRUE){
          print(final_results)
          print(plot_results)
        }else{
          print(final_results)
        }
      }
      
      ## Enclosing all the plots in a list object to return
      if(plot == TRUE & CI == "boot"){
        results <- list("Estimates" = final_results,
                        "Alpha_CI" = plot1_a,
                        "Gamma_CI" = plot2_g,
                        "Beta_CI" = plot3_b,
                        "Tau_CI" = plot4_t,
                        "MedEff" = plot5,
                        "MedEff_CI" = plot6)
      }
      else if(plot == TRUE & CI != "boot"){
        results <- list("Estimates" = final_results,
                        "Alpha_CI" = plot1_a,
                        "Gamma_CI" = plot2_g,
                        "Beta_CI" = plot3_b,
                        "Tau_CI" = plot4_t,
                        "MedEff" = plot5)  
      }
      else{
        results <- list("Estimates" = final_results)
      }
      
      return(results)
    }
    else{
      print(paste("Error: Accepted values for CI are 'boot' and 'none';
                  you have entered an unacceptable value for CI."))
    }
  }else{
    print(paste("tvma() stopped execution due to unacceptable class type of function argument(s)."))
  }
}