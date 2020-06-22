#' Main function for time varying mediation function for continuous Outcome and three treatment (Exposure) Groups
#' 
#' Part of the set of internal functions to estimate the time-varying mediation effect and bootstrap standard errors, involving three treatment groups and continuous outcome.
#' 
#' @param NRT1        a vector with treatment arm 1 values
#' @param NRT2        a vector with treatment arm 2 values
#' @param t.seq       a vector of time points for each observation
#' @param x           matrix of mediator values in wide format
#' @param y           matrix of outcome outcomes in wide format
#' @param t.est       time points to make the estimation             Default = t.seq.
#' 
#' @return \item{hat.alpha1}{estimated NRT1 effect on mediator}
#' @return \item{hat.alpha2}{estimated NRT2 effect on mediator}
#' @return \item{hat.beta3}{estimated mediation effect on outcome}
#' @return \item{hat.mediation1}{time varying mediation effect - NRT1 on outcome}
#' @return \item{hat.mediation2}{time varying mediation effect - NRT2 on outcome}
#' 
#' @export


tvmcurve_3trt<-function(NRT1, NRT2, t.seq, x, y, t.est)
{
  # Estimate time-varying mediation effect and bootstrap standard errors
  # for continuous outcome and 3 treatment arms
  #
  # Args:
  #   NRT1        -->   a vector with treatment arm 1 values
  #   NRT2        -->   a vector with treatment arm 2 values
  #   t.seq       -->   a vector of time points for each observation
  #   x           -->   matrix of mediator values in wide format
  #   y           -->   matrix of outcome outcomes in wide format
  #   t.est       -->   time points to make the estimation                              Default = t.seq
  #   plot        -->   TRUE or FALSE for plotting mediation effect                     Default = "FALSE"
  #   CI          -->   "none" or "boot" method of deriving confidence intervals.       Default = "boot"
  #   replicates  -->   Number of replicates for bootstrapping confidence intervals.    Default = 1000
  #   verbose     -->   TRUE or FALSE for printing results to screen.                   Default = "FALSE"
  #
  # Returns:
  #   hat.alpha.1       -->   estimated treatment effect on mediator
  #   hat.beta.2        -->   estimated mediation effect on outcome
  #   mediation.effect  -->   time varying mediation effect
  #   sd                -->   estimated standard deviation of the mediation effect (only for asymptotic estimation)
  #   CI.upper          -->   Upper confidence intervals
  #   CI.lower          -->   Lower confidence intervals
  #
  ##
  
  deltat=max(diff(t.seq))/2
  #temp.coeff stores all the estimated coefficient values at t-seq
  t.coeff=NULL
  for(l in 2:length(t.seq))
  {
    t.coeff=cbind(t.coeff, coeff(l, NRT1, NRT2, x, y)$coeff.est)
  }
  
  bw_alpha1 <- thumbBw(t.seq[-1], t.coeff[1,], deg=1, kernel=gaussK)
  bw_alpha2 <- thumbBw(t.seq[-1], t.coeff[2,], deg=1, kernel=gaussK)
  bw_beta3 <- thumbBw(t.seq[-1], t.coeff[5,], deg=1, kernel=gaussK)
  hat.alpha.1=locPolSmootherC(t.seq[-1], t.coeff[1,], t.est-deltat, bw_alpha1, deg=1, kernel=gaussK)$beta0
  hat.alpha.2=locPolSmootherC(t.seq[-1], t.coeff[2,], t.est-deltat, bw_alpha2, deg=1, kernel=gaussK)$beta0
  hat.beta.3=locPolSmootherC(t.seq[-1], t.coeff[5,], t.est, bw_beta3, deg=1, kernel=gaussK)$beta0
  
  list(hat.alpha1 = hat.alpha.1, hat.alpha2 = hat.alpha.2, hat.beta3 = hat.beta.3,
       hat.mediation1=hat.alpha.1*hat.beta.3, hat.mediation2=hat.alpha.2*hat.beta.3)
}