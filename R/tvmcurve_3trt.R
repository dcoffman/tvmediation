#' Main function for time varying mediation function for continuous Outcome and three treatment arms (Exposure Groups)
#' 
#' Part of the set of internal functions to estimate the time-varying mediation effect and 
#' bootstrap standard errors for three treatment groups and continuous outcome.
#' 
#' @param T1          a vector indicating assignment to treatment 1
#' @param T2          a vector indicating assignment to treatment 2
#' @param t.seq       a vector of time points for each observation
#' @param x           matrix of mediator values in wide format
#' @param y           matrix of outcome values in wide format
#' @param t.est       time points at which to make the estimation. Default = t.seq.
#' 
#' @return \item{hat.alpha1}{estimated Treatment 1 effect on mediator}
#' @return \item{hat.alpha2}{estimated Treatment 2 effect on mediator}
#' @return \item{hat.gamma1}{estimated Treatment 1 direct effect on outcome}
#' @return \item{hat.gamma2}{estimated Treatment 2 direct effect on outcome}
#' @return \item{hat.tau1}{estimated Treatment 1 total effect on outcome}
#' @return \item{hat.tau2}{estimated Treatment 2 total effect on outcome}
#' @return \item{hat.beta}{estimated mediator effect on outcome}
#' @return \item{hat.mediation1}{time varying mediation effect for Treatment 1 on outcome}
#' @return \item{hat.mediation2}{time varying mediation effect for Treatment 2 on outcome}
#' 
#' @export

tvmcurve_3trt<-function(T1, T2, t.seq, x, y, t.est)
{
  
  deltat <- max(diff(t.seq))/2
  #temp.coeff stores all the estimated coefficient values at t-seq
  t.coeff <- NULL
  for(l in 2:length(t.seq))
  {
    X.new.l <- cbind(T1, T2)
    X.new.l <- scale(X.new.l, center = TRUE, scale = FALSE)
    Y.new.l <- scale(y[l-1,], center=TRUE, scale=FALSE)
    
    nomissing.x <- complete.cases(X.new.l)
    nomissing.y <- complete.cases(Y.new.l)
    nomissing.index <- nomissing.x*nomissing.y
    
    X.new.l <- X.new.l[which(nomissing.index == 1),]
    Y.new.l <- Y.new.l[which(nomissing.index == 1)]
   
    coeff.est.t <- solve(t(X.new.l)%*%(X.new.l))%*%t(X.new.l)%*%(Y.new.l)
    
    coeff.all <- rbind(coeff(l, T1, T2, x, y)$coeff.est, coeff.est.t)
    t.coeff <- cbind(t.coeff, coeff.all)
  }
  
  bw_alpha1 <- locpol::thumbBw(t.seq[-1], t.coeff[1,], deg=1, kernel=locpol::gaussK)
  bw_alpha2 <- locpol::thumbBw(t.seq[-1], t.coeff[2,], deg=1, kernel=locpol::gaussK)
  bw_gamma1 <- locpol::thumbBw(t.seq[-1], t.coeff[3,], deg=1, kernel=locpol::gaussK)
  bw_gamma2 <- locpol::thumbBw(t.seq[-1], t.coeff[4,], deg=1, kernel=locpol::gaussK)
  bw_beta <- locpol::thumbBw(t.seq[-1], t.coeff[5,], deg=1, kernel=locpol::gaussK)
  bw_tau1 <- locpol::thumbBw(t.seq[-1], t.coeff[6,], deg=1, kernel=locpol::gaussK)
  bw_tau2 <- locpol::thumbBw(t.seq[-1], t.coeff[7,], deg=1, kernel=locpol::gaussK)
  
  hat.alpha.1 <- locpol::locPolSmootherC(t.seq[-1], t.coeff[1,], t.est-deltat, 
                                         bw_alpha1, deg=1, kernel=locpol::gaussK)$beta0
  hat.alpha.2 <- locpol::locPolSmootherC(t.seq[-1], t.coeff[2,], t.est-deltat, 
                                         bw_alpha2, deg=1, kernel=locpol::gaussK)$beta0
  hat.gamma.1 <- locpol::locPolSmootherC(t.seq[-1], t.coeff[3,], t.est, bw_gamma1, deg=1, 
                                         kernel=locpol::gaussK)$beta0
  hat.gamma.2 <- locpol::locPolSmootherC(t.seq[-1], t.coeff[4,], t.est, bw_gamma2, 
                                         deg=1, kernel=locpol::gaussK)$beta0
  hat.beta <- locpol::locPolSmootherC(t.seq[-1], t.coeff[5,], t.est, bw_beta, deg=1, 
                                      kernel=locpol::gaussK)$beta0
  hat.tau.1 <- locpol::locPolSmootherC(t.seq[-1], t.coeff[6,], t.est-deltat, bw_tau1, 
                                       deg=1, kernel=locpol::gaussK)$beta0
  hat.tau.2 <- locpol::locPolSmootherC(t.seq[-1], t.coeff[7,], t.est-deltat, bw_tau2, 
                                       deg=1, kernel=locpol::gaussK)$beta0
    
  list(hat.alpha1 = hat.alpha.1, hat.alpha2 = hat.alpha.2,
       hat.gamma1 = hat.gamma.1, hat.gamma2 = hat.gamma.2,
       hat.tau1 = hat.tau.1, hat.tau2 = hat.tau.2, hat.beta = hat.beta,
       hat.mediation1=hat.alpha.1*hat.beta, hat.mediation2=hat.alpha.2*hat.beta)
}
