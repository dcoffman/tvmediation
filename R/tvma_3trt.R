#' Time Varying Mediation Function: Continuous Outcome and Three Treatment Groups
#' 
#' Function to estimate the time-varying mediation effect and bootstrap standard errors for three treatment groups and a continuous outcome.
#' 
#' @param T1          a vector indicating assignment to treatment 1
#' @param T2          a vector indicating assignment to treatment 2
#' @param t.seq       a vector of time points for each observation
#' @param mediator    matrix of mediator values in wide format
#' @param outcome     matrix of outcome values in wide format
#' @param t.est       a vector of time points at which to make the estimation. Default = t.seq. (OPTIONAL ARGUMENT)
#' @param plot        TRUE or FALSE for plotting mediation effect. Default = "FALSE". (OPTIONAL ARGUMENT)
#' @param CI          "none" or "boot" method of deriving confidence intervals. Default = "boot". (OPTIONAL ARGUMENT)
#' @param replicates  number of replicates for bootstrapping confidence intervals. Default = 1000. (OPTIONAL ARGUMENT)
#' @param grpname     name of the treatment arms (exposure groups) to be displayed in the results. Default = "T". (OPTIONAL ARGUMENT) 
#' @param verbose     TRUE or FALSE for printing results to screen. Default = "FALSE". (OPTIONAL ARGUMENT)
#' 
#' @return \item{hat.alpha1}{estimated Treatment 1 effect on mediator}
#' @return \item{CI.lower.alpha1}{CI lower limit for estimated coefficient hat.alpha1}
#' @return \item{CI.upper.alpha1}{CI upper limit for estimated coefficient hat.alpha1}
#' @return \item{hat.alpha2}{estimated Treatment 2 effect on mediator}
#' @return \item{CI.lower.alpha2}{CI lower limit for estimated coefficient hat.alpha2}
#' @return \item{CI.upper.alpha2}{CI upper limit for estimated coefficient hat.alpha2}
#' @return \item{hat.gamma1}{estimated Treatment 1 direct effect on outcome}
#' @return \item{CI.lower.gamma1}{CI lower limit for estimated coefficient hat.gamma1}
#' @return \item{CI.upper.gamma1}{CI upper limit for estimated coefficient hat.gamma1}
#' @return \item{hat.gamma2}{estimated Treatment 2 direct effect on outcome}
#' @return \item{CI.lower.gamma2}{CI lower limit for estimated coefficient hat.gamma2}
#' @return \item{CI.upper.gamma2}{CI upper limit for estimated coefficient hat.gamma2}
#' @return \item{hat.tau1}{estimated Treatment 1 total effect on outcome}
#' @return \item{CI.lower.tau1}{CI lower limit for estimated coefficient hat.tau1}
#' @return \item{CI.upper.tau1}{CI upper limit for estimated coefficient hat.tau1}
#' @return \item{hat.tau2}{estimated Treatment 2 total effect on outcome}
#' @return \item{CI.lower.tau2}{CI lower limit for estimated coefficient hat.tau2}
#' @return \item{CI.upper.tau2}{CI upper limit for estimated coefficient hat.tau2}
#' @return \item{hat.beta}{estimated mediator effect on outcome}
#' @return \item{CI.lower.beta}{CI lower limit for estimated coefficient hat.beta}
#' @return \item{CI.upper.beta}{CI upper limit for estimated coefficient hat.beta}
#' @return \item{hat.mediation1}{time varying mediation effect for Treatment 1 on outcome}
#' @return \item{SE_MedEff1}{estimated standard errors of hat.mediation1}
#' @return \item{CI.upper.T1}{CI upper limit for hat.mediation1}
#' @return \item{CI.lower.T1}{CI lower limit for hat.mediation1}
#' @return \item{hat.mediation2}{time varying mediation effect for Treatment 2 on outcome}
#' @return \item{SE_MedEff2}{estimated standard errors of hat.mediation2}
#' @return \item{CI.upper.T2}{CI upper limit for hat.mediation2}
#' @return \item{CI.lower.T2}{CI lower limit for hat.mediation2}
#' 
#' @section Plot Returns:
#' \enumerate{
#' \item{\code{plot1_a1 }}{plot for hat.alpha1 with CIs over t.est}
#' \item{\code{plot2_a2 }}{plot for hat.alpha2 with CIs over t.est}
#' \item{\code{plot3_g1 }}{plot for hat.gamma1 with CIs over t.est}
#' \item{\code{plot4_g2 }}{plot for hat.gamma2 with CIs over t.est}
#' \item{\code{plot5_t1 }}{plot for hat.tau1 with CIs over t.est}
#' \item{\code{plot6_t2 }}{plot for hat.tau2 with CIs over t.est}
#' \item{\code{plot7_b }}{plot for hat.beta with CIs over t.est}
#' \item{\code{MedEff_T1 }}{plot for hat.mediation1 over t.est}
#' \item{\code{MedEff_T2 }}{plot for hat.mediation2 over t.est}
#' \item{\code{MedEff_CI_T1 }}{plot for hat.mediation1 with CIs over t.est}
#' \item{\code{MedEff_CI_T2 }}{plot for hat.mediation2 with CIs over t.est}
#' }
#' 
#' @examples
#' data(smoker)
#' 
#' # GENERATE WIDE FORMATTED MEDIATORS
#' mediator <- LongToWide(smoker$SubjectID,
#'                        smoker$timeseq, 
#'                        smoker$NegMoodLst15min)
#' 
#' # GENERATE WIDE FORMATTED OUTCOMES
#' outcome <- LongToWide(smoker$SubjectID,
#'                       smoker$timeseq,
#'                       smoker$cessFatig)
#' 
#' # GENERATE TWO BINARY TREATMENT VARIABLES
#' NRT1 <- as.numeric(unique(smoker[,c("SubjectID","varenicline")])[,2])-1
#' NRT2 <- as.numeric(unique(smoker[,c("SubjectID","comboNRT")])[,2])-1
#' 
#' # GENERATE A VECTOR OF UNIQUE TIME POINTS
#' t.seq <- sort(unique(smoker$timeseq))
#' 
#' # COMPUTE TIME VARYING MEDIATION ANALYSIS USING BOOTSTRAPPED CONFIDENCE INTERVALS
#' # results <- tvma_3trt(NRT1, NRT2, t.seq, mediator, outcome)
#' 
#' # COMPUTE TIME VARYING MEDIATION ANALYSIS FOR SPECIFIED POINTS IN TIME USING 500 REPLICATES
#' # results <- tvma_3trt(NRT1, NRT2, t.seq, mediator, outcome,
#' #                      t.est = c(0.2, 0.4, 0.6, 0.8),
#' #                      replicates = 500)
#' 
#' @references 
#' \enumerate{
#' \item{Fan, J. and Gijbels, I. (1996). Local polynomial modelling and its applications: monographs on statistics and applied probability 66. CRC Press.}
#' \item{Fan, J. and Zhang, W. (1999). Statistical estimation in varying coefficient models. The Annals of Statistics, 27, 1491-1518.}
#' \item{Fan, J. and Zhang, W. (2000). Two-step estimation of functional linear models with applications to longitudinal data. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 62, 303-322.}
#' \item{Cai, X., Piper, M. E., Li, R., & Coffman, D. L. (under review). Estimation and inference for the mediation effect in a time-varying mediation model. Journal of the Royal Statistics Society, Series C (Applied Statistics).}
#' \item{Baker, T. B., Piper, M. E., Stein, J. H., Smith, S. S., Bolt, D. M., Fraser, D. L., & Fiore, M. C. (2016). Effects of nicotine patch vs varenicline vs combination nicotine replacement therapy on smoking cessation at 26 weeks: A randomized clinical trial. JAMA, 315(4), 371-379. doi:10.1001/jama.2015.19284.}
#' \item{Efron, B. and Tibshirani, R. (1986). Bootstrap methods for standard errors, confidence intervals, and other measures of statistical accuracy. Statistical Science, 1, 54-75. doi:10.1214/ss/1177013815.}
#' }
#' 
#' @export
#' 
#' @importFrom stats complete.cases cov glm lm loess na.omit predict quantile sd var
#' @import dplyr
#' @import ggplot2
#' @import kedd
#' @import locpol

tvma_3trt <- function(T1, T2, t.seq, mediator, outcome, t.est = t.seq, plot = FALSE, CI="boot", replicates = 1000, grpname = "T", verbose = FALSE)
  {
      ### Optional results based on user argument CI = "boot"
      #
      #   CI.upper.NRT1     -->   Upper confidence intervals for NRT1
      #   CI.lower.NRT1     -->   Lower confidence intervals for NRT1
      #   CI.upper.NRT2     -->   Upper confidence intervals for NRT2
      #   CI.lower.NRT2     -->   Lower confidence intervals for NRT2
      #   SE_MedEff1        -->   estimated standard errors of the mediation effect for NRT1
      #   SE_MedEff2        -->   estimated standard errors of the mediation effect for NRT2
      
      ## Testing the class of the arguments passed in the function
      ctm <- class(mediator)
      cto <- class(outcome)
      ctt1 <- class(T1)
      ctt2 <- class(T2)
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
      if(is.vector(T1) != TRUE || ctt1 != "numeric"){
        print("Error: `T1` is not of class type `numeric vector`.")
        flag <- flag + 1
      }
      if(is.vector(T2) != TRUE || ctt2 != "numeric"){
        print("Error: `T2` is not of class type `numeric vector`.")
        flag <- flag + 1
      }
      if(is.vector(t.seq) != TRUE || cts != "numeric"){
        print("Error: `t.seq` is not of class type `numeric vector`.")
        flag <- flag + 1
      }
      if(flag == 0){
        if(CI == "boot" || CI == "none"){
          results_me <- tvmcurve_3trt(T1, T2, t.seq, mediator, outcome, t.est)
          alldata <- list(T1 = T1, T2 = T2, x = mediator, y = outcome, t.seq = t.seq)
          original.coeff <- list(hat.alpha1 = results_me$hat.alpha1, hat.alpha2 = results_me$hat.alpha2,
                                 hat.gamma1 = results_me$hat.gamma1, hat.gamma2 = results_me$hat.gamma2, 
                                 hat.tau1 = results_me$hat.tau1, hat.tau2 = results_me$hat.tau2,
                                 hat.beta = results_me$hat.beta)
          
          #### Computing CI for coefficients irrespective of user choice of bootstrapping ####
          bootcoeff_CI <- bootci_coeff_3trt(T1, T2, t.seq, mediator, outcome, t.est, original.coeff, replicates)
          
          test1 <- cbind(as.data.frame(results_me), as.data.frame(bootcoeff_CI), t.est)
          
          #### CI for mediation ####
          if(CI == "boot"){
            results_ci <- bootci_tvm_3trt(replicates, alldata, t.est)
            test2 <- cbind(as.data.frame(results_ci), t.est)
            
            final_dat <- merge(test1, test2, all.x = TRUE)
            final_results <- final_dat %>%
              select(-orig.mediation1, -orig.mediation2)
            
            names(final_results)[c(1,11:24,29:30)] <- c("timeseq",
                                                 "CI.lower.alpha1", "CI.upper.alpha1",
                                                 "CI.lower.alpha2", "CI.upper.alpha2",
                                                 "CI.lower.gamma1", "CI.upper.gamma1",
                                                 "CI.lower.gamma2", "CI.upper.gamma2",
                                                 "CI.lower.tau1", "CI.upper.tau1",
                                                 "CI.lower.tau2", "CI.upper.tau2",
                                                 "CI.lower.beta", "CI.upper.beta",
                                                 "SE_MedEff1","SE_MedEff2")

            final_results <- final_results[c("timeseq",
                                             "hat.alpha1", "CI.lower.alpha1", "CI.upper.alpha1",
                                             "hat.alpha2", "CI.lower.alpha2", "CI.upper.alpha2",
                                             "hat.gamma1", "CI.lower.gamma1", "CI.upper.gamma1",
                                             "hat.gamma2", "CI.lower.gamma2", "CI.upper.gamma2",
                                             "hat.tau1", "CI.lower.tau1", "CI.upper.tau1",
                                             "hat.tau2", "CI.lower.tau2", "CI.upper.tau2",
                                             "hat.beta", "CI.lower.beta", "CI.upper.beta",
                                             "hat.mediation1", "SE_MedEff1", "plw1", "pup1",
                                             "hat.mediation2", "SE_MedEff2", "plw2", "pup2")]
            
            names(final_results)[c(25)]<- paste("CI.lower.",grpname,"1", sep = "")
            names(final_results)[c(26)]<- paste("CI.upper.",grpname,"1", sep = "")
            names(final_results)[c(29)]<- paste("CI.lower.",grpname,"2", sep = "")
            names(final_results)[c(30)]<- paste("CI.upper.",grpname,"2", sep = "")
          }
          else{
            final_results <- test1
            names(final_results)[10:24] <- c("CI.lower.alpha1", "CI.upper.alpha1",
                                             "CI.lower.alpha2", "CI.upper.alpha2",
                                             "CI.lower.gamma1", "CI.upper.gamma1",
                                             "CI.lower.gamma2", "CI.upper.gamma2",
                                             "CI.lower.tau1", "CI.upper.tau1",
                                             "CI.lower.tau2", "CI.upper.tau2",
                                             "CI.lower.beta", "CI.upper.beta",
                                             "timeseq")
            final_results <- final_results[c("timeseq",
                                             "hat.alpha1", "CI.lower.alpha1", "CI.upper.alpha1",
                                             "hat.alpha2", "CI.lower.alpha2", "CI.upper.alpha2",
                                             "hat.gamma1", "CI.lower.gamma1", "CI.upper.gamma1",
                                             "hat.gamma2", "CI.lower.gamma2", "CI.upper.gamma2",
                                             "hat.tau1", "CI.lower.tau1", "CI.upper.tau1",
                                             "hat.tau2", "CI.lower.tau2", "CI.upper.tau2",
                                             "hat.beta", "CI.lower.beta", "CI.upper.beta",
                                             "hat.mediation1", "hat.mediation2")]
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
            
            # First Plot: alpha1 coefficients smoothed across time
            plot1_a1 <- ggplot(data = final_results, aes(timeseq, hat.alpha1)) +
                                geom_line(color = "red", size = 0.75) +
                                geom_line(aes(timeseq, CI.lower.alpha1), size = 0.8, color = "blue", linetype = "dashed") +
                                geom_line(aes(timeseq, CI.upper.alpha1), size = 0.8, color = "blue", linetype = "dashed") +
                                labs(title = paste("Plot of the time-varying effect on the mediator (",grpname,"1)", sep = ""),
                                     x = "Time Sequence",
                                     y = "Alpha1 coefficient") +
                                scale_x_continuous(breaks = seq(l, u, i))
            
            # Second plot: alpha2 coefficients smoothed across time
            plot2_a2 <- ggplot(data = final_results, aes(timeseq, hat.alpha2)) +
                                geom_line(color = "red", size = 0.75) +
                                geom_line(aes(timeseq, CI.lower.alpha2), size = 0.8, color = "blue", linetype = "dashed") +
                                geom_line(aes(timeseq, CI.upper.alpha2), size = 0.8, color = "blue", linetype = "dashed") +
                                labs(title = paste("Plot of the time-varying effect on the mediator (",grpname,"2)", sep = ""),
                                     x = "Time Sequence",
                                     y = "Alpha2 coefficient") +
                                scale_x_continuous(breaks = seq(l, u, i))
            
            # Third plot: gamma1 coefficients smoothed across time
            plot3_g1 <- ggplot(data = final_results, aes(timeseq, hat.gamma1)) +
                                geom_line(color = "red", size = 0.75) +
                                geom_line(aes(timeseq, CI.lower.gamma1), size = 0.8, color = "blue", linetype = "dashed") +
                                geom_line(aes(timeseq, CI.upper.gamma1), size = 0.8, color = "blue", linetype = "dashed") +
                                labs(title = paste("Plot of the time-varying direct effect (",grpname,"1)", sep = ""),
                                     x = "Time Sequence",
                                     y = "Gamma1 coefficient") +
                                scale_x_continuous(breaks = seq(l, u, i))
            
            # Fourth plot: gamma2 coefficients smoothed across time
            plot4_g2 <- ggplot(data = final_results, aes(timeseq, hat.gamma2)) +
                                geom_line(color = "red", size = 0.75) +
                                geom_line(aes(timeseq, CI.lower.gamma2), size = 0.8, color = "blue", linetype = "dashed") +
                                geom_line(aes(timeseq, CI.upper.gamma2), size = 0.8, color = "blue", linetype = "dashed") +
                                labs(title = paste("Plot of the time-varying direct effect (",grpname,"2)", sep = ""),
                                     x = "Time Sequence",
                                     y = "Gamma2 coefficient") +
                                scale_x_continuous(breaks = seq(l, u, i))
            
            # Fifth plot: tau1 coefficients smoothed across time
            plot5_t1 <- ggplot(data = final_results, aes(timeseq, hat.tau1)) +
              geom_line(color = "red", size = 0.75) +
              geom_line(aes(timeseq, CI.lower.tau1), size = 0.8, color = "blue", linetype = "dashed") +
              geom_line(aes(timeseq, CI.upper.tau1), size = 0.8, color = "blue", linetype = "dashed") +
              labs(title = paste("Plot of the time-varying total effect (",grpname,"1)", sep = ""),
                   x = "Time Sequence",
                   y = "Tau1 coefficient") +
              scale_x_continuous(breaks = seq(l, u, i))
            
            # Sixth plot: tau2 coefficients smoothed across time
            plot6_t2 <- ggplot(data = final_results, aes(timeseq, hat.tau2)) +
              geom_line(color = "red", size = 0.75) +
              geom_line(aes(timeseq, CI.lower.tau2), size = 0.8, color = "blue", linetype = "dashed") +
              geom_line(aes(timeseq, CI.upper.tau2), size = 0.8, color = "blue", linetype = "dashed") +
              labs(title = paste("Plot of the time-varying total effect (",grpname,"2)", sep = ""),
                   x = "Time Sequence",
                   y = "Tau2 coefficient") +
              scale_x_continuous(breaks = seq(l, u, i))
            
            # Seventh plot: beta coefficients smoothed across time
            plot7_b <- ggplot(data = final_results, aes(timeseq, hat.beta)) +
                                geom_line(color = "red", size = 0.75) +
                                geom_line(aes(timeseq, CI.lower.beta), size = 0.8, color = "blue", linetype = "dashed") +
                                geom_line(aes(timeseq, CI.upper.beta), size = 0.8, color = "blue", linetype = "dashed") +
                                labs(title = "Plot of the time-varying effect of the mediator on the outcome",
                                     x = "Time Sequence",
                                     y = "Beta coefficient") +
                                scale_x_continuous(breaks = seq(l, u, i))
            
            # Eighth plot: T1 mediation effect smoothed across time
            plot8 <- ggplot(data = final_results, aes(timeseq, hat.mediation1)) +
                              geom_line(color = "red", size = 0.75) +
                              labs(title = paste("Time-varying mediation effect (",grpname,"1)", sep = ""),
                                   x = "Time Sequence",
                                   y = paste("Mediation Effect for",grpname,"1", sep = "")) +
                              scale_x_continuous(breaks = seq(l, u, i))
            
            # Ninth plot: T2 mediation effect smoothed across time
            plot9 <- ggplot(data = final_results, aes(timeseq, hat.mediation2)) +
                              geom_line(color = "red", size = 0.75) +
                              labs(title = paste("Time-varying mediation effect (",grpname,"2)", sep = ""),
                                   x = "Time Sequence",
                                   y = paste("Mediation Effect for",grpname,"2", sep = "")) +
                              scale_x_continuous(breaks = seq(l, u, i))
            
            if(CI == "boot"){
              
              CI.lower.T1 <- paste("CI.lower.",grpname,"1", sep = "")
              CI.upper.T1 <- paste("CI.upper.",grpname,"1", sep = "")
              CI.lower.T2 <- paste("CI.lower.",grpname,"2", sep = "")
              CI.upper.T2 <- paste("CI.upper.",grpname,"2", sep = "")
              
              # Tenth plot: T1 mediation effect with 95% CIs
              plot10 <- ggplot(data = final_results, aes(timeseq, hat.mediation1)) +
                                geom_line(size = 1, color = "red") +
                                geom_line(aes_string("timeseq", CI.lower.T1), size = 0.8, color = "blue", linetype = "dashed") +
                                geom_line(aes_string("timeseq", CI.upper.T1), size = 0.8, color = "blue", linetype = "dashed") +
                                geom_line(aes(timeseq, 0)) +
                                labs(title = paste("Mediation Effect (",grpname,"1) with 95% CIs (computed with percentile bootstrap)", sep = ""),
                                     x = "Time Sequence",
                                     y = paste("Mediation Effect for ", grpname, "1", sep = "")) + 
                                theme(legend.position = "none")  +
                                scale_x_continuous(breaks = seq(l, u, i))  
              
              # Eleventh plot: T2 mediation effect with 95% CIs
              plot11 <- ggplot(data = final_results, aes(timeseq, hat.mediation2)) +
                                geom_line(size = 1, color = "red") +
                                geom_line(aes_string("timeseq", CI.lower.T2), size = 0.8, color = "blue", linetype = "dashed") +
                                geom_line(aes_string("timeseq", CI.upper.T2), size = 0.8, color = "blue", linetype = "dashed") +
                                geom_line(aes(timeseq, 0)) +
                                labs(title = paste("Mediation Effect (",grpname,"2) with 95% CIs (computed with percentile bootstrap)", sep = ""),
                                     x = "Time Sequence",
                                     y = paste("Mediation Effect for ", grpname, "2", sep = "")) + 
                                theme(legend.position = "none")  +
                                scale_x_continuous(breaks = seq(l, u, i)) 
              
              plot_results <- list("plot1_a1" = plot1_a1,
                                   "plot2_a2" = plot2_a2,
                                   "plot3_g1" = plot3_g1,
                                   "plot4_g2" = plot4_g2,
                                   "plot5_t1" = plot5_t1,
                                   "plot6_t2" = plot6_t2,
                                   "plot7_b" = plot7_b,
                                   "MedEff_T1" = plot8,
                                   "MedEff_T2" = plot9,
                                   "MedEff_CI_T1" = plot10,
                                   "MedEff_CI_T2" = plot11)
            }else{
              plot_results <- list("plot1_a1" = plot1_a1,
                                   "plot2_a2" = plot2_a2,
                                   "plot3_g1" = plot3_g1,
                                   "plot4_g2" = plot4_g2,
                                   "plot5_t1" = plot5_t1,
                                   "plot6_t2" = plot6_t2,
                                   "plot7_b" = plot7_b,
                                   "MedEff_T1" = plot8,
                                   "MedEff_T2" = plot9)
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
                            "plot3_g1" = plot3_g1,
                            "plot4_g2" = plot4_g2,
                            "plot5_t1" = plot5_t1,
                            "plot6_t2" = plot6_t2,
                            "plot7_b" = plot7_b,
                            "MedEff_T1" = plot8,
                            "MedEff_T2" = plot9,
                            "MedEff_CI_T1" = plot10,
                            "MedEff_CI_T2" = plot11)
          }
          else if(plot == TRUE & CI != "boot"){
            results <- list("Estimates" = final_results,
                            "plot1_a1" = plot1_a1,
                            "plot2_a2" = plot2_a2,
                            "plot3_g1" = plot3_g1,
                            "plot4_g2" = plot4_g2,
                            "plot5_t1" = plot5_t1,
                            "plot6_t2" = plot6_t2,
                            "plot7_b" = plot7_b,
                            "MedEff_T1" = plot8,
                            "MedEff_T2" = plot9)  
          }
          else{
            results <- list("Estimates" = final_results)
          }
          
          return(results)
          
        }else{
          print(paste("Error: Accepted values for CI are 'boot' and 'none';
                      you have entered an unacceptable value for CI."))
        }
      }else{
        print(paste("tvma_3trt() stopped execution due to unacceptable class type of function argument(s)."))
      }
  }