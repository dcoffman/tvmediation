#' Main function for time varying mediation function for continuous Outcome and three treatment arms (Exposure Groups)
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
#' @return \item{hat.alpha1}{estimated NRT1 (treatment arm 1) effect on mediator}
#' @return \item{hat.alpha2}{estimated NRT2 (treatment arm 2) effect on mediator}
#' @return \item{hat.beta1}{estimated NRT1 (treatment arm 1) effect on outcome}
#' @return \item{hat.beta2}{estimated NRT2 (treatment arm 2) effect on outcome}
#' @return \item{hat.beta3}{estimated mediator effect on outcome}
#' @return \item{hat.mediation1}{time varying mediation effect - NRT1 (treatment arm 1) on outcome}
#' @return \item{hat.mediation2}{time varying mediation effect - NRT2 (treatment arm 2) on outcome}
#' 
#' @export


tvmcurve_3trt<-function(NRT1, NRT2, t.seq, x, y, t.est)
{
  # Estimate time-varying mediation effect and bootstrap standard errors
  # for continuous outcome and 3 treatment arms
  #
  # Arguments:
  #   NRT1        -->   a vector with treatment arm 1 values
  #   NRT2        -->   a vector with treatment arm 2 values
  #   t.seq       -->   a vector of time points for each observation
  #   x           -->   matrix of mediator values in wide format
  #   y           -->   matrix of outcome outcomes in wide format
  #   t.est       -->   time points to make the estimation                  Default = t.seq
  #
  # Returns:
  #   hat.alpha1        -->   estimated NRT1 (treatment) effect on mediator
  #   hat.alpha2        -->   estimated NRT2 (treatment) effect on mediator
  #   hat.beta1         -->   estimated NRT1 (treatment) effect on outcome
  #   hat.beta2         -->   estimated NRT2 (treatment) effect on outcome
  #   hat.beta3         -->   estimated mediator effect on outcome
  #   hat.mediation1    -->   time varying mediation effect - NRT1 on outcome
  #   hat.mediation2    -->   time varying mediation effect - NRT2 on outcome
  #
  ##
  
  deltat=max(diff(t.seq))/2
  #temp.coeff stores all the estimated coefficient values at t-seq
  t.coeff=NULL
  for(l in 2:length(t.seq))
  {
    t.coeff=cbind(t.coeff, coeff(l, NRT1, NRT2, x, y)$coeff.est)
  }
  
  bw_alpha1 <- locpol::thumbBw(t.seq[-1], t.coeff[1,], deg=1, kernel= locpol::gaussK)
  bw_alpha2 <- locpol::thumbBw(t.seq[-1], t.coeff[2,], deg=1, kernel= locpol::gaussK)
  bw_beta1 <- locpol::thumbBw(t.seq[-1], t.coeff[3,], deg=1, kernel= locpol::gaussK)
  bw_beta2 <- locpol::thumbBw(t.seq[-1], t.coeff[4,], deg=1, kernel= locpol::gaussK)
  bw_beta3 <- locpol::thumbBw(t.seq[-1], t.coeff[5,], deg=1, kernel= locpol::gaussK)
  
  hat.alpha.1=locpol::locPolSmootherC(t.seq[-1], t.coeff[1,], t.est-deltat, bw_alpha1, deg=1, kernel= locpol::gaussK)$beta0
  hat.alpha.2=locpol::locPolSmootherC(t.seq[-1], t.coeff[2,], t.est-deltat, bw_alpha2, deg=1, kernel= locpol::gaussK)$beta0
  hat.beta.1=locpol::locPolSmootherC(t.seq[-1], t.coeff[3,], t.est, bw_beta1, deg=1, kernel= locpol::gaussK)$beta0
  hat.beta.2=locpol::locPolSmootherC(t.seq[-1], t.coeff[4,], t.est, bw_beta2, deg=1, kernel= locpol::gaussK)$beta0
  hat.beta.3=locpol::locPolSmootherC(t.seq[-1], t.coeff[5,], t.est, bw_beta3, deg=1, kernel= locpol::gaussK)$beta0
  
  list(hat.alpha1 = hat.alpha.1, hat.alpha2 = hat.alpha.2,
       hat.beta1 = hat.beta.1, hat.beta2 = hat.beta.2, hat.beta3 = hat.beta.3,
       hat.mediation1=hat.alpha.1*hat.beta.3, hat.mediation2=hat.alpha.2*hat.beta.3)
}