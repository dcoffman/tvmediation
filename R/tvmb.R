#' Time Varying Mediation Function: Binary Outcome and Two Treatment Arms (Exposure Groups)
#' 
#' Function to estimate the time-varying mediation effect and bootstrap standard errors, involving two treatment groups and binary outcome.
#' 
#' @param treatment   a vector with treatment values
#' @param t.seq       a vector of unique time points for each observation
#' @param mediator    matrix of mediator values in wide format
#' @param outcome     matrix of outcome outcomes in wide format
#' @param plot        TRUE or FALSE for plotting mediation effect                     Default = "FALSE". (OPTIONAL ARGUMENT)
#' @param CI          "none" or "boot" method of deriving confidence intervals.       Default = "boot". (OPTIONAL ARGUMENT)
#' @param replicates  Number of replicates for bootstrapping confidence intervals.    Default = 1000. (OPTIONAL ARGUMENT)
#' @param verbose     TRUE or FALSE for printing results to screen.                   Default = "FALSE". (OPTIONAL ARGUMENT)
#' 
#' @return \item{timeseq}{time points of estimation}
#' @return \item{alpha_hat}{estimated main treatment arm (exposure group) of interest effect on mediator (indirect effect component)}
#' @return \item{CI.lower.a}{lower limit of confidence intervals for estimated coefficient alpha_hat}
#' @return \item{CI.upper.a}{upper limit of confidence intervals for estimated coefficient alpha_hat}
#' @return \item{gamma_hat}{estimated main treatment arm (exposure group) of interest effect on outcome (direct effect)}
#' @return \item{CI.lower.g}{lower limit of confidence intervals for estimated coefficient gamma_hat}
#' @return \item{CI.upper.g}{upper limit of confidence intervals for estimated coefficient gamma_hat}
#' @return \item{beta_hat}{estimated mediator effect on outcome (indirect effect component)}
#' @return \item{CI.lower.b}{lower limit of confidence intervals for estimated coefficient beta_hat}
#' @return \item{CI.upper.b}{upper limit of confidence intervals for estimated coefficient beta_hat}
#' @return \item{tao_hat}{estimated main treatment arm (exposure group) of interest effect on outcome, excluding adjustment for mediator (total effect)}
#' @return \item{CI.lower.t}{lower limit of confidence intervals for estimated coefficient tao_hat}
#' @return \item{CI.upper.t}{upper limit of confidence intervals for estimated coefficient tao_hat}
#' @return \item{medEffect}{time varying mediation effect - main treatment arm (exposure group) of interest on outcome (estimated by product method)}
#' @return \item{CI.lower}{lower limit of confidence intervals for medEffect}
#' @return \item{CI.upper}{upper limit of confidence intervals for medEffect}
#' 
#' @section Plot Returns:
#' \enumerate{
#' \item{\code{plot1_a }}{plot for alpha_hat across t.seq with CIs}
#' \item{\code{plot2_g }}{plot for gamma_hat across t.seq with CIs}
#' \item{\code{plot3_b }}{plot for beta_hat across t.seq with CIs}
#' \item{\code{plot4_t }}{plot for tao_hat across t.seq with CIs}
#' \item{\code{MedEff }}{plot for medEffect across t.seq}
#' \item{\code{MedEff_CI }}{plot for medEffect with CIs across t.seq}
#' \item{\code{bootstrap }}{plot for estimated medEffect(s) from bootstrapped samples across t.seq}
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
#'                       smoker.sub$smoke_status)
#' 
#' # GENERATE A BINARY TREATMENT VARIABLE
#' trt <- as.numeric(unique(smoker.sub[ , c("SubjectID","varenicline")])[ , 2])-1
#' 
#' # GENERATE A VECTOR OF UNIQUE TIME POINTS
#' t.seq <- sort(unique(smoker.sub$timeseq))
#' 
#' # COMPUTE TIME VARYING MEDIATION ANALYSIS USING BOOTSTRAPPED CONFIDENCE INTERVALS
#' results <- tvmb(trt, t.seq, mediator, outcome)
#' 
#' @references 
#' \enumerate{
#' \item{Fan, J. and Gijbels, I. (1996). Local polynomial modelling and its applications: monographs on statistics and applied probability 66 66. CRC Press.}
#' \item{Fan, J. and Zhang, W. (1999). Statistical estimation in varying coefficient models. The annals of Statistics 27 1491-1518.}
#' \item{Fan, J. and Zhang, W. (2000). Two-step estimation of functional linear models with applications to longitudinal data. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 62 303-322.}
#' \item{Baker, T. B., Piper, M. E., Stein, J. H., Smith, S. S., Bolt, D. M., Fraser, D. L., & Fiore, M. C. (2016). Effects of nicotine patch vs varenicline vs combination nicotine replacement therapy on smoking cessation at 26 weeks: A randomized clinical trial. JAMA, 315(4), 371-379. doi:10.1001/jama.2015.19284.}
#' \item{Efron, B.; Tibshirani, R. Bootstrap Methods for Standard Errors, Confidence Intervals, and Other Measures of Statistical Accuracy. Statist. Sci. 1 (1986), no. 1, 54--75. doi:10.1214/ss/1177013815.}
#' }
#' 
#' @export
#' 
#' @importFrom stats complete.cases cov glm lm loess na.omit predict quantile sd var
#' @import dplyr
#' @import ggplot2
#' @import kedd
#' @import locpol

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
  #   alpha_hat        -->   estimated exposure effect on mediator (indirect effect component)
  #   CI.lower.a       -->   lower limit of confidence intervals for alpha1_hat
  #   CI.upper.a       -->   upper limit of confidence intervals for alpha1_hat
  #   gamma_hat        -->   estimated exposure effect on outcome (direct effect)
  #   CI.lower.g       -->   lower limit of confidence intervals for gamma_hat
  #   CI.upper.g       -->   upper limit of confidence intervals for gamma_hat
  #   beta_hat         -->   estimated mediator effect on outcome (indirect effect component)
  #   CI.lower.b       -->   lower limit of confidence intervals for beta_hat
  #   CI.upper.b       -->   upper limit of confidence intervals for beta_hat
  #   tao_hat          -->   estimated effect of exposure on outcome, excluding adjustment for mediator (total effect)
  #   CI.lower.t       -->   lower limit of confidence intervals for tao_hat
  #   CI.upper.t       -->   upper limit of confidence intervals for tao_hat
  #   medDiff          -->   time varying mediation effect (difference term) ## excluded
  #   medEffect        -->   time varying mediation effect (product term)
  #   CI.lower         -->   lower limit of confidence intervals for medEffect
  #   CI.upper         -->   upper limit of confidence intervals for medEffect
  #
  # Optional Returns:
  #   plot1_a         -->   plot for alpha_hat with CIs across t.seq
  #   plot2_g         -->   plot for gamma_hat with CIs across t.seq
  #   plot3_b         -->   plot for beta_hat with CIs across t.seq
  #   plot4_t         -->   plot for tao_hat with CIs across t.seq
  #   MedEff          -->   plot for mediation effects (difference and product) across t.seq
  #   MedEff_CI       -->   plot for CIs of medEffect
  #   bootstrap       -->   plot for estimated medEffects from bootstrapped samples across t.seq
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
        # fit alpha coefficient for each time point
        aAll = vector()
        
        for(i in 1:nm){
          fit1 = lm(m[i,] ~ treatment,na.action=na.omit)
          aAll = append(aAll,fit1$coefficients[[2]])
        }
        
        # create smoothing line
        smootha = loess(aAll ~ t.seq[1:nm], span = 0.3, degree=1)
        
        # creating a dataframe with the time sequences, a1 and smoothed coefficients
        test1 <- data.frame(cbind(t.seq, aAll, smootha$fitted))
        names(test1)[3] <- "smootha"
        
        
        # Regression of outcome(y) on mediator(x2) and exposure(treatment - x1)
        # fit gamma and beta coefficients for each time point
        bAll = vector()
        gAll = vector()
        sd2Real = vector()
        
        for(i in 2:nm){
          fit2 = glm(outcome[i,] ~ treatment + m[(i-1),],family="binomial",na.action=na.omit)
          
          bHat = fit2$coefficients[[3]]
          gHat = fit2$coefficients[[2]]
          
          #calculate sd to standardize with
          sd2 = sqrt(gHat^2*var(treatment,na.rm=TRUE)+bHat^2*var(m[(i-1),],na.rm=TRUE)+2*gHat*bHat*cov(treatment,m[(i-1),],use="complete.obs")+(pi^2/3))
          
          #append standardized coefficient to list of all coefficients
          gAll = append(gAll,gHat/sd2)
          bAll = append(bAll,bHat/sd2) 
        }
        
        # smooth
        t.seq.b <- t.seq
        t.seq.b <- t.seq.b[-1]
        smoothg = loess(gAll ~ t.seq.b[1:length(t.seq.b)], span = 0.2, degree = 1)
        smoothb = loess(bAll ~ t.seq.b[1:length(t.seq.b)], span = 0.2, degree = 1)
        
        # regress outcome(y) on exposure(x) at each time point to use in difference method
        # fit tao coefficient for each time point
        tAll = vector()

        for(i in 2:nm){

          fit3 = glm(outcome[i,]~treatment, family="binomial",na.action=na.omit)
          tHat = fit3$coefficients[[2]]

          sd1 = sqrt(tHat^2*var(treatment,na.rm=TRUE)+(pi^2/3))

          # append standardized coefficient to list of all coefficients
          tAll = append(tAll,tHat/sd1)
        }

        # smooth
        smootht = loess(tAll ~ t.seq.b[1:length(t.seq.b)], span = 0.2, degree = 1)
        
        # creating a dataframe with the time sequences, and smoothed coefficients
        test2 <- data.frame(cbind(t.seq.b, gAll, smoothg$fitted, bAll, smoothb$fitted, tAll, smootht$fitted))
        names(test2)[3] <- "smoothg"
        names(test2)[5] <- "smoothb"
        names(test2)[7] <- "smootht"
        
        ##### *********************************************************************** #####
        ##### Bootstrapping samples to estimate confidence intervals for coefficients #####
        ##### *********************************************************************** #####
        
        coeff_CI <- bootci_coeff_binary(treatment, t.seq, m, outcome, replicates)
        
        #*********************************************************************************#
        
        
        #### Formatting the results into a single dataframe ####
        
        coeff_data1 <- merge(test1, test2, by.x = "t.seq", by.y = "t.seq.b", all.x = TRUE)
        coeff_data2 <- merge(coeff_data1, coeff_CI, by.x = "t.seq")
        coeff_data <- coeff_data2[c("t.seq","aAll","smootha", "CI.lower.a", "CI.upper.a",
                                            "gAll", "smoothg", "CI.lower.g", "CI.upper.g",
                                            "bAll", "smoothb", "CI.lower.b", "CI.upper.b",
                                            "tAll", "smootht", "CI.lower.t", "CI.upper.t")]
        
        
        #calculate mediation effects
        #really b(t)*a(t-1) because `a` starts at t=1 while `b` starts at t=2
        for(i in 1:nrow(coeff_data)){
          if(!is.na(coeff_data$gAll[i])){
            coeff_data$medProd[i] = coeff_data$bAll[i]*coeff_data$aAll[i-1]
            # coeff_data$medDif[i] = coeff_data$cAll[i] - coeff_data$b1All[i]
            ##exclusion of mediation effects estimate through difference method
          }
        }
        
        #calculate smooth line for mediation effect estimate via product method
        medProd <- coeff_data$medProd
        medProd <- medProd[which(!is.na(medProd))]
        smoothProd = loess(medProd ~ t.seq.b[1:length(t.seq.b)], span = 0.3,degree=1)
        
        # #calculate smooth line for differences
        # medDif <- coeff_data$medDif
        # medDif <- medDif[which(!is.na(medDif))]
        # smoothDif = loess(medDif ~ t.seq.b[1:length(t.seq.b)], span = 0.3,degree=1)
        
        #creating a dataframe with the time sequences, mediation effects and smoothened coefficients
        test_a <- data.frame(cbind(t.seq.b, medProd, smoothProd$fitted))
        names(test_a)[2] <- "med_pt"
        names(test_a)[3] <- "smooth_Prod"
        test_a$type <- "Prod"
        
        ##exclusion of mediation effects estimate through difference method
        # test_b <- data.frame(cbind(t.seq.b, medDif, smoothDif$fitted))
        # names(test_b)[2] <- "med_pt"
        # names(test_b)[3] <- "smooth"
        # test_b$type <- "Diff"
        
        # test3 <- rbind(test_a, test_b)
        
        coeff_data <- merge(coeff_data, test_a, by.x = "t.seq", by.y = "t.seq.b",
                            all.x = TRUE) %>%
                            select(-med_pt, -type)
        names(coeff_data)[19] <- "smooth_medProd"
        
        ##exclusion of mediation effects estimate through difference method
        # coeff_data <- merge(coeff_data, test_b, by.x = "t.seq", by.y = "t.seq.b",
        #                     all.x = TRUE) %>%
        #   select(-med_pt, -type)
        # names(coeff_data)[21] <- "smooth_medDif"
        
        
        ##### ****************************************************** #####
        ##### Bootstrapping samples to estimate confidence intervals #####
        ##### ****************************************************** #####;
        
        if(CI == "boot"){
          list_all <- bootci_tvmb(treatment, t.seq, m, outcome, coeff_data, replicates)
          IE_t <- list_all$bootstrap_result
          final_dat <- list_all$all_results
          final_dat1 <- final_dat %>%
            select(- c(aAll, gAll, bAll, tAll, medProd))
          final_results <- final_dat1
          names(final_results)[c(1, 2, 5, 8, 11, 14)] <- c("timeseq", "alpha_hat", "gamma_hat", "beta_hat", "tao_hat", "medEffect")
        }else{
          final_dat <- coeff_data
          final_dat1 <- final_dat %>%
            select(- c(aAll, gAll, bAll, tAll, medProd))
          final_results <- final_dat1
          names(final_results)[c(1, 2, 5, 8, 11, 14)] <- c("timeseq", "alpha_hat", "gamma_hat", "beta_hat", "tao_hat", "medEffect")
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
          
          # First Plot: plotting alpha1 coefficients (smoothed) using across the time using ggplot
          plot1_a <- ggplot(data = final_results, aes(timeseq, alpha_hat)) +
            geom_line(color = "red", size = 0.75) +
            geom_line(aes(timeseq, CI.lower.a), color = "blue", size = 0.8, linetype = "dashed") +
            geom_line(aes(timeseq, CI.upper.a), color = "blue", size = 0.8, linetype = "dashed") +
            labs(title = "Plotting the alpha coefficients",
                 x = "Time Sequence",
                 y = "Alpha") +
            scale_x_continuous(breaks = seq(l, u, i))
          
          # Second plot: plotting beta1 coefficients (smoothed) across the time using ggplot
          plot2_g <- ggplot(data = final_results, aes(timeseq, gamma_hat)) +
            geom_line(color = "red", size = 0.75) +
            geom_line(aes(timeseq, CI.lower.g), color = "blue", size = 0.8, linetype = "dashed") +
            geom_line(aes(timeseq, CI.upper.g), color = "blue", size = 0.8, linetype = "dashed") +
            labs(title = "Plotting the gamma coefficients",
                 x = "Time Sequence",
                 y = "Gamma") +
            scale_x_continuous(breaks = seq(l, u, i))
          
          # Third plot: plotting beta2 coefficients (smoothed) across the time using ggplot
          plot3_b <- ggplot(data = final_results, aes(timeseq, beta_hat)) +
            geom_line(color = "red", size = 0.75) +
            geom_line(aes(timeseq, CI.lower.b), color = "blue", size = 0.8, linetype = "dashed") +
            geom_line(aes(timeseq, CI.upper.b), color = "blue", size = 0.8, linetype = "dashed") +
            labs(title = "Plotting the beta coefficients",
                 x = "Time Sequence",
                 y = "Beta") +
            scale_x_continuous(breaks = seq(l, u, i))
          
          # Fourth plot: plotting tao coefficients (smoothed) across the time using ggplot
          plot4_t <- ggplot(data = final_results, aes(timeseq, tao_hat)) +
            geom_line(color = "red", size = 0.75) +
            geom_line(aes(timeseq, CI.lower.t), color = "blue", size = 0.8, linetype = "dashed") +
            geom_line(aes(timeseq, CI.upper.t), color = "blue", size = 0.8, linetype = "dashed") +
            labs(title = "Plotting the tao coefficients",
                 x = "Time Sequence",
                 y = "Tao") +
            scale_x_continuous(breaks = seq(l, u, i))
          
          # Fifth plot: plotting the mediation effects across time using ggplot
          plot5 <- ggplot(data = final_results, aes(timeseq, medEffect)) +
            geom_line(color = "red", size = 0.75) +
            labs(title = "Plotting the mediation effect (estimated via product method)",
                 x = "Time Sequence",
                 y = "Mediation Effect") +
            scale_x_continuous(breaks = seq(l, u, i))
          
          if(CI == "boot"){
            # Sixth plot: plotting the mediation effect with 95% CIs
            plot6 <- ggplot(data = final_results, aes(timeseq, medEffect)) +
              geom_line(size = 1, color = "red") +
              geom_line(aes(timeseq, CI.lower), color = "blue", size = 0.8, linetype = "dashed") +
              geom_line(aes(timeseq, CI.upper), color = "blue", size = 0.8, linetype = "dashed") +
              geom_line(aes(timeseq, 0)) +
              labs(title = "Mediation Effect (estimated via product method) with 95% CIs (computed with percentile bootstrap)",
                   x = "Time Sequence",
                   y = "Mediation Effect") + 
              theme(legend.position = "none") +
              scale_x_continuous(breaks = seq(l, u, i))
            
            # Seventh plot: plotting the mediation effect from bootstrap samples
            plot7 <- ggplot(data = IE_t, aes(t.seq.b, V2)) + 
              geom_line() + 
              labs(title = "Estimated Mediation Effect(s) from Bootstrapping",
                   x = "Time Sequence",
                   y = "Mediation Effect") +
              scale_x_continuous(breaks = seq(l, u, i))
            for (i in 2:ncol(IE_t)) {
              x <- data.frame(cbind(t.seq.b,IE_t[,i]))
              names(x)[2] <- "val"
              plot7 <- plot7 + geom_line(data = x, aes(t.seq.b, val))
            }
            
            plot_results <- list("plot1_a" = plot1_a,
                                 "plot2_g" = plot2_g,
                                 "plot3_b" = plot3_b,
                                 "plot4_t" = plot4_t,
                                 "MedEff" = plot5,
                                 "MedEff_CI" = plot6,
                                 "bootstrap" = plot7)
          }else{
            plot_results <- list("plot1_a" = plot1_a,
                                 "plot2_g" = plot2_g,
                                 "plot3_b" = plot3_b,
                                 "plot4_t" = plot4_t,
                                 "MedEff" = plot5)
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
                          "plot1_a" = plot1_a,
                          "plot2_g" = plot2_g,
                          "plot3_b" = plot3_b,
                          "plot4_t" = plot4_t,
                          "MedEff" = plot5,
                          "MedEff_CI" = plot6,
                          "bootstrap" = plot7)
        }
        else if(plot == TRUE & CI != "boot"){
          results <- list("Estimates" = final_results,
                          "plot1_a" = plot1_a,
                          "plot2_g" = plot2_g,
                          "plot3_b" = plot3_b,
                          "plot4_t" = plot4_t,
                          "MedEff" = plot5)  
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
