#' Time Varying Mediation Function: Continuous Outcome and Three Treatment Arms (Exposure Groups)
#' 
#' Function to estimate the time-varying mediation effect and bootstrap standard errors, involving three treatment/exposure groups and continuous outcome.
#' 
#' @param NRT1        a vector with treatment arm 1 (exposure group 1) values
#' @param NRT2        a vector with treatment arm 2 (exposure group 2) values
#' @param t.seq       a vector of time points for each observation
#' @param mediator    matrix of mediator values in wide format
#' @param outcome     matrix of outcome outcomes in wide format
#' @param t.est       a vector of time points to make the estimation                                  Default = t.seq. (OPTIONAL ARGUMENT)
#' @param plot        TRUE or FALSE for plotting mediation effect                                     Default = "FALSE". (OPTIONAL ARGUMENT)
#' @param CI          "none" or "boot" method of deriving confidence intervals.                       Default = "boot". (OPTIONAL ARGUMENT)
#' @param replicates  number of replicates for bootstrapping confidence intervals.                    Default = 1000. (OPTIONAL ARGUMENT)
#' @param grpname     name of the treatment arms (exposure groups) to be displayed in the results.    Default = "NRT". (OPTIONAL ARGUMENT) 
#' @param verbose     TRUE or FALSE for printing results to screen.                                   Default = "FALSE". (OPTIONAL ARGUMENT)
#' 
#' @return \item{hat.alpha1}{estimated effect of treatment arm 1 (exposure group 1) on mediator (indirect effect component)}
#' @return \item{CI.lower.alpha1}{lower limit of confidence intervals for estimated coefficient hat.alpha1}
#' @return \item{CI.upper.alpha1}{upper limit of confidence intervals for estimated coefficient hat.alpha1}
#' @return \item{hat.alpha2}{estimated effect of treatment arm 2 (exposure group 2) on mediator (indirect effect component)}
#' @return \item{CI.lower.alpha2}{lower limit of confidence intervals for estimated coefficient hat.alpha2}
#' @return \item{CI.upper.alpha2}{upper limit of confidence intervals for estimated coefficient hat.alpha2}
#' @return \item{hat.gamma1}{estimated effect of treatment arm 1 (exposure group 1) on outcome (direct effect)}
#' @return \item{CI.lower.gamma1}{lower limit of confidence intervals for estimated coefficient hat.gamma1}
#' @return \item{CI.upper.gamma1}{upper limit of confidence intervals for estimated coefficient hat.gamma1}
#' @return \item{hat.gamma2}{estimated effect of treatment arm 2 (exposure group 2) on outcome (direct effect)}
#' @return \item{CI.lower.gamma2}{lower limit of confidence intervals for estimated coefficient hat.gamma2}
#' @return \item{CI.upper.gamma2}{upper limit of confidence intervals for estimated coefficient hat.gamma2}
#' @return \item{hat.beta}{estimated mediator effect on outcome (indirect effect component)}
#' @return \item{CI.lower.beta}{lower limit of confidence intervals for estimated coefficient hat.beta}
#' @return \item{CI.upper.beta}{upper limit of confidence intervals for estimated coefficient hat.beta}
#' @return \item{hat.mediation1}{time varying mediation effect - treatment arm 1 (exposure group 1) on outcome}
#' @return \item{SE_MedEff1}{estimated standard errors of hat.mediation1}
#' @return \item{CI.upper.NRT1}{upper limit of confidence intervals for hat.mediation1}
#' @return \item{CI.lower.NRT1}{lower limit of confidence intervals for hat.mediation1}
#' @return \item{hat.mediation2}{time varying mediation effect - treatment arm 2 (exposure group 2) on outcome}
#' @return \item{SE_MedEff2}{estimated standard errors of hat.mediation2}
#' @return \item{CI.upper.NRT2}{upper limit of confidence intervals for hat.mediation2}
#' @return \item{CI.lower.NRT2}{lower limit of confidence intervals for hat.mediation2}
#' 
#' @section Plot Returns:
#' \enumerate{
#' \item{\code{plot1_a1 }}{plot for hat.alpha1 across t.est with CIs}
#' \item{\code{plot2_a2 }}{plot for hat.alpha2 across t.est with CIs}
#' \item{\code{plot3_g1 }}{plot for hat.gamma1 across t.est with CIs}
#' \item{\code{plot4_g2 }}{plot for hat.gamma2 across t.est with CIs}
#' \item{\code{plot5_b }}{plot for hat.beta across t.est with CIs}
#' \item{\code{MedEff_NRT1 }}{plot for hat.mediation1 across t.est}
#' \item{\code{MedEff_NRT2 }}{plot for hat.mediation2 across t.est}
#' \item{\code{MedEff_CI_NRT1 }}{plot for hat.mediation1 with CIs across t.est}
#' \item{\code{MedEff_CI_NRT2 }}{plot for hat.mediation2 with CIs across t.est}
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
#' results <- tvma_3trt(NRT1, NRT2, t.seq, mediator, outcome)
#' 
#' # COMPUTE TIME VARYING MEDIATION ANALYSIS FOR SPECIFIED POINTS IN TIME USING 500 REPLICATES
#' results <- tvma_3trt(NRT1, NRT2, t.seq, mediator, outcome,
#'                      t.est = c(0.2, 0.4, 0.6, 0.8),
#'                      replicates = 500)
#' 
#' @references 
#' \enumerate{
#' \item{Fan, J. and Gijbels, I. (1996). Local polynomial modelling and its applications: monographs on statistics and applied probability 66 66. CRC Press.}
#' \item{Fan, J. and Zhang, W. (1999). Statistical estimation in varying coefficient models. The annals of Statistics 27 1491-1518.}
#' \item{Fan, J. and Zhang, W. (2000). Two-step estimation of functional linear models with applications to longitudinal data. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 62 303-322.}
#' }
#' 
#' @export
#' 
#' @importFrom stats complete.cases cov glm lm loess na.omit predict quantile sd var
#' @import dplyr
#' @import ggplot2
#' @import kedd
#' @import locpol

tvma_3trt <- function(NRT1, NRT2, t.seq, mediator, outcome, t.est = t.seq, plot = FALSE, CI="boot", replicates = 1000, grpname = "NRT", verbose = FALSE)
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
      #   t.est       -->   time points to make the estimation                                              Default = t.seq
      #   plot        -->   TRUE or FALSE for plotting mediation effect                                     Default = "FALSE"
      #   CI          -->   "none" or "boot" method of deriving confidence intervals.                       Default = "boot"
      #   replicates  -->   Number of replicates for bootstrapping confidence intervals.                    Default = 1000
      #   grpname     -->   Name of the treatment arms (exposure groups) to be displayed in the results.    Default = "NRT".
      #   verbose     -->   TRUE or FALSE for printing results to screen.                                   Default = "FALSE"
      #
      ### RETURNS:
      #
      ### Dataframe 'final_results' with the following elements
      #
      #   hat.alpha1        -->   estimated NRT1 (treatment) effect on mediator
      #   CI.lower.alpha1   -->   lower limit of 95% confidence intervals for hat.alpha1
      #   CI.upper.alpha1   -->   upper limit of 95% confidence intervals for hat.alpha1
      #   hat.alpha2        -->   estimated NRT2 (treatment) effect on mediator
      #   CI.lower.alpha2   -->   lower limit of 95% confidence intervals for hat.alpha2
      #   CI.upper.alpha2   -->   upper limit of 95% confidence intervals for hat.alpha2
      #   hat.gamma1        -->   estimated NRT1 (treatment) effect on outcome
      #   CI.lower.gamma1   -->   lower limit of 95% confidence intervals for hat.gamma1
      #   CI.upper.gamma1   -->   upper limit of 95% confidence intervals for hat.gamma1
      #   hat.gamma2        -->   estimated NRT2 (treatment) effect on outcome
      #   CI.lower.gamma2   -->   lower limit of 95% confidence intervals for hat.gamma2
      #   CI.upper.gamma2   -->   upper limit of 95% confidence intervals for hat.gamma2
      #   hat.beta          -->   estimated mediator effect on outcome
      #   CI.lower.beta     -->   lower limit of 95% confidence intervals for hat.beta
      #   CI.upper.beta     -->   upper limit of 95% confidence intervals for hat.beta
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
      # `plot1_a1` plot for `hat.alpha1` with 95% CIs across `timeseq`
      # `plot2_a2` plot for `hat.alpha2` with 95% CIs across `timeseq`
      # `plot3_g1` plot for `hat.beta1` with 95% CIs across `timeseq`
      # `plot4_g2` plot for `hat.beta2` with 95% CIs across `timeseq`
      # `plot5_b` plot for `hat.beta3` with 95% CIs across `timeseq`
      # `MedEff_NRT1` plot for `hat.mediation1` across `timeseq`
      # `MedEff_NRT2` plot for `hat.mediation2` across `timeseq`
      # `MedEff_CI_NRT1` plot for `hat.mediation1` with 95% CIs across `timeseq`
      # `MedEff_CI_NRT2` plot for `hat.mediation2` with 95% CIs across `timeseq`
      #
      ##
      
      ## Testing the class type of the arguments passed in the function ##
      ctm <- class(mediator)
      cto <- class(outcome)
      ctt1 <- class(NRT1)
      ctt2 <- class(NRT2)
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
          original.coeff <- list(hat.alpha1 = results_me$hat.alpha1, hat.alpha2 = results_me$hat.alpha2,
                                 hat.gamma1 = results_me$hat.gamma1, hat.gamma2 = results_me$hat.gamma2, hat.beta = results_me$hat.beta)
          
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
            
            names(final_results)[c(1,9:18,23:24)] <- c("timeseq",
                                                 "CI.lower.alpha1", "CI.upper.alpha1",
                                                 "CI.lower.alpha2", "CI.upper.alpha2",
                                                 "CI.lower.gamma1", "CI.upper.gamma1",
                                                 "CI.lower.gamma2", "CI.upper.gamma2",
                                                 "CI.lower.beta", "CI.upper.beta",
                                                 "SE_MedEff1","SE_MedEff2")

            final_results <- final_results[c("timeseq",
                                             "hat.alpha1", "CI.lower.alpha1", "CI.upper.alpha1",
                                             "hat.alpha2", "CI.lower.alpha2", "CI.upper.alpha2",
                                             "hat.gamma1", "CI.lower.gamma1", "CI.upper.gamma1",
                                             "hat.gamma2", "CI.lower.gamma2", "CI.upper.gamma2",
                                             "hat.beta", "CI.lower.beta", "CI.upper.beta",
                                             "hat.mediation1", "SE_MedEff1", "plw1", "pup1",
                                             "hat.mediation2", "SE_MedEff2", "plw2", "pup2")]
            
            names(final_results)[c(19)]<- paste("CI.lower.",grpname,"1", sep = "")
            names(final_results)[c(20)]<- paste("CI.upper.",grpname,"1", sep = "")
            names(final_results)[c(23)]<- paste("CI.lower.",grpname,"2", sep = "")
            names(final_results)[c(24)]<- paste("CI.upper.",grpname,"2", sep = "")
          }
          else{
            final_results <- test1
            names(final_results)[8:18] <- c( "CI.lower.alpha1", "CI.upper.alpha1",
                                             "CI.lower.alpha2", "CI.upper.alpha2",
                                             "CI.lower.gamma1", "CI.upper.gamma1",
                                             "CI.lower.gamma2", "CI.upper.gamma2",
                                             "CI.lower.beta", "CI.upper.beta",
                                             "timeseq")
            final_results <- final_results[c("timeseq",
                                             "hat.alpha1", "CI.lower.alpha1", "CI.upper.alpha1",
                                             "hat.alpha2", "CI.lower.alpha2", "CI.upper.alpha2",
                                             "hat.gamma1", "CI.lower.gamma1", "CI.upper.gamma1",
                                             "hat.gamma2", "CI.lower.gamma2", "CI.upper.gamma2",
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
            
            # First Plot: plotting alpha1 coefficients across time using ggplot
            plot1_a1 <- ggplot(data = final_results, aes(timeseq, hat.alpha1)) +
                                geom_line(color = "red", size = 0.75) +
                                geom_line(aes(timeseq, CI.lower.alpha1), size = 0.8, color = "blue", linetype = "dashed") +
                                geom_line(aes(timeseq, CI.upper.alpha1), size = 0.8, color = "blue", linetype = "dashed") +
                                labs(title = "Plotting the alpha1 coefficients",
                                     x = "Time Sequence",
                                     y = "Alpha1") +
                                scale_x_continuous(breaks = seq(l, u, i))
            
            # Second plot: plotting alpha2 coefficients across time using ggplot
            plot2_a2 <- ggplot(data = final_results, aes(timeseq, hat.alpha2)) +
                                geom_line(color = "red", size = 0.75) +
                                geom_line(aes(timeseq, CI.lower.alpha2), size = 0.8, color = "blue", linetype = "dashed") +
                                geom_line(aes(timeseq, CI.upper.alpha2), size = 0.8, color = "blue", linetype = "dashed") +
                                labs(title = "Plotting the alpha2 coefficients",
                                     x = "Time Sequence",
                                     y = "Alpha2") +
                                scale_x_continuous(breaks = seq(l, u, i))
            
            # Third plot: plotting gamma1 coefficients across time using ggplot
            plot3_g1 <- ggplot(data = final_results, aes(timeseq, hat.gamma1)) +
                                geom_line(color = "red", size = 0.75) +
                                geom_line(aes(timeseq, CI.lower.gamma1), size = 0.8, color = "blue", linetype = "dashed") +
                                geom_line(aes(timeseq, CI.upper.gamma1), size = 0.8, color = "blue", linetype = "dashed") +
                                labs(title = "Plotting the gamma1 coefficients",
                                     x = "Time Sequence",
                                     y = "Gamma1") +
                                scale_x_continuous(breaks = seq(l, u, i))
            
            # Fourth plot: plotting gamma2 coefficients across time using ggplot
            plot4_g2 <- ggplot(data = final_results, aes(timeseq, hat.gamma2)) +
                                geom_line(color = "red", size = 0.75) +
                                geom_line(aes(timeseq, CI.lower.gamma2), size = 0.8, color = "blue", linetype = "dashed") +
                                geom_line(aes(timeseq, CI.upper.gamma2), size = 0.8, color = "blue", linetype = "dashed") +
                                labs(title = "Plotting the gamma2 coefficients",
                                     x = "Time Sequence",
                                     y = "Gamma2") +
                                scale_x_continuous(breaks = seq(l, u, i))
            
            # Fifth plot: plotting beta coefficients across time using ggplot
            plot5_b <- ggplot(data = final_results, aes(timeseq, hat.beta)) +
                                geom_line(color = "red", size = 0.75) +
                                geom_line(aes(timeseq, CI.lower.beta), size = 0.8, color = "blue", linetype = "dashed") +
                                geom_line(aes(timeseq, CI.upper.beta), size = 0.8, color = "blue", linetype = "dashed") +
                                labs(title = "Plotting the beta coefficients",
                                     x = "Time Sequence",
                                     y = "Beta") +
                                scale_x_continuous(breaks = seq(l, u, i))
            
            # Sixth plot: plotting the mediation effect of treatment arm1
            plot6 <- ggplot(data = final_results, aes(timeseq, hat.mediation1)) +
                              geom_line(color = "red", size = 0.75) +
                              labs(title = paste("Plotting the time-varying mediation effect (",grpname,"1)", sep = ""),
                                   x = "Time Sequence",
                                   y = paste("Mediation Effect for",grpname,"1", sep = "")) +
                              scale_x_continuous(breaks = seq(l, u, i))
            
            # Seventh plot: plotting the mediation effect of treatment arm2
            plot7 <- ggplot(data = final_results, aes(timeseq, hat.mediation2)) +
                              geom_line(color = "red", size = 0.75) +
                              labs(title = paste("Plotting the time-varying mediation effect (",grpname,"2)", sep = ""),
                                   x = "Time Sequence",
                                   y = paste("Mediation Effect for",grpname,"2", sep = "")) +
                              scale_x_continuous(breaks = seq(l, u, i))
            
            if(CI == "boot"){
              
              CI.lower.NRT1 <- paste("CI.lower.",grpname,"1", sep = "")
              CI.upper.NRT1 <- paste("CI.upper.",grpname,"1", sep = "")
              CI.lower.NRT2 <- paste("CI.lower.",grpname,"2", sep = "")
              CI.upper.NRT2 <- paste("CI.upper.",grpname,"2", sep = "")
              
              # Eighth plot: plotting the mediation effect of treatment arm1 with 95% CIs
              plot8 <- ggplot(data = final_results, aes(timeseq, hat.mediation1)) +
                                geom_line(size = 1, color = "red") +
                                geom_line(aes_string("timeseq", CI.lower.NRT1), size = 0.8, color = "blue", linetype = "dashed") +
                                geom_line(aes_string("timeseq", CI.upper.NRT1), size = 0.8, color = "blue", linetype = "dashed") +
                                geom_line(aes(timeseq, 0)) +
                                labs(title = paste("Mediation Effect (",grpname,"1) with 95% CIs (computed with percentile bootstrap)", sep = ""),
                                     x = "Time Sequence",
                                     y = paste("Mediation Effect for ", grpname, "1", sep = "")) + 
                                theme(legend.position = "none")  +
                                scale_x_continuous(breaks = seq(l, u, i))  
              
              # Ninth plot: plotting the mediation effect of treatment arm2 with 95% CIs
              plot9 <- ggplot(data = final_results, aes(timeseq, hat.mediation2)) +
                                geom_line(size = 1, color = "red") +
                                geom_line(aes_string("timeseq", CI.lower.NRT2), size = 0.8, color = "blue", linetype = "dashed") +
                                geom_line(aes_string("timeseq", CI.upper.NRT2), size = 0.8, color = "blue", linetype = "dashed") +
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
                                   "plot5_b" = plot5_b,
                                   "MedEff_NRT1" = plot6,
                                   "MedEff_NRT2" = plot7,
                                   "MedEff_CI_NRT1" = plot8,
                                   "MedEff_CI_NRT2" = plot9)
            }else{
              plot_results <- list("plot1_a1" = plot1_a1,
                                   "plot2_a2" = plot2_a2,
                                   "plot3_g1" = plot3_g1,
                                   "plot4_g2" = plot4_g2,
                                   "plot5_b" = plot5_b,
                                   "MedEff_NRT1" = plot6,
                                   "MedEff_NRT2" = plot7)
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
                            "plot5_b" = plot5_b,
                            "MedEff_NRT1" = plot6,
                            "MedEff_NRT2" = plot7,
                            "MedEff_CI_NRT1" = plot8,
                            "MedEff_CI_NRT2" = plot9)
          }
          else if(plot == TRUE & CI != "boot"){
            results <- list("Estimates" = final_results,
                            "plot1_a1" = plot1_a1,
                            "plot2_a2" = plot2_a2,
                            "plot3_g1" = plot3_g1,
                            "plot4_g2" = plot4_g2,
                            "plot5_b" = plot5_b,
                            "MedEff_NRT1" = plot6,
                            "MedEff_NRT2" = plot7)  
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