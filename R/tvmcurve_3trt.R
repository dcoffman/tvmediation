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
#' @return \item{hat.gamma1}{estimated NRT1 (treatment arm 1) effect on outcome}
#' @return \item{hat.gamma2}{estimated NRT2 (treatment arm 2) effect on outcome}
#' @return \item{hat.tao1}{estimated NRT1 (treatment arm 1) effect on outcome, excluding adjustment for mediator}
#' @return \item{hat.tao2}{estimated NRT2 (treatment arm 2) effect on outcome, excluding adjustment for mediator}
#' @return \item{hat.beta}{estimated mediator effect on outcome}
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
  #   hat.gamma1        -->   estimated NRT1 (treatment) effect on outcome
  #   hat.gamma2        -->   estimated NRT2 (treatment) effect on outcome
  #   hat.tao1          -->   estimated NRT1 (treatment) effect on outcome, excluding adjustment for mediator
  #   hat.tao2          -->   estimated NRT2 (treatment) effect on outcome, excluding adjustment for mediator
  #   hat.beta          -->   estimated mediator effect on outcome
  #   hat.mediation1    -->   time varying mediation effect - NRT1 on outcome
  #   hat.mediation2    -->   time varying mediation effect - NRT2 on outcome
  #
  ##
  
  deltat=max(diff(t.seq))/2
  #temp.coeff stores all the estimated coefficient values at t-seq
  t.coeff=NULL
  for(l in 2:length(t.seq))
  {
    X.new.l<-cbind(NRT1, NRT2)
    X.new.l = scale(X.new.l, center = TRUE, scale = FALSE)
    # NRT1.new.l=scale(NRT1,center=TRUE, scale=FALSE)
    # NRT2.new.l=scale(NRT2,center=TRUE, scale=FALSE)
    Y.new.l=scale(y[l-1,],center=TRUE, scale=FALSE)
    
    nomissing.x <- complete.cases(X.new.l)
    # nomissing.NRT1<-complete.cases(NRT1.new.l)
    # nomissing.NRT2<-complete.cases(NRT2.new.l)
    nomissing.y<-complete.cases(Y.new.l)
    
    nomissing.index <- nomissing.x * nomissing.y
    # nomissing.index.NRT1 <- nomissing.NRT1 * nomissing.y
    # nomissing.index.NRT2 <- nomissing.NRT2 * nomissing.y
    
    X.new.l = X.new.l[which(nomissing.index == 1), ]
    Y.new.l = Y.new.l[which(nomissing.index == 1)]
    # NRT1.new.l = NRT1.new.l[which(nomissing.index.NRT1 == 1),]
    # NRT2.new.l = NRT2.new.l[which(nomissing.index.NRT2 == 1),]
    # Y.new.l.NRT1 = Y.new.l[which(nomissing.index.NRT1 == 1)]
    # Y.new.l.NRT2 = Y.new.l[which(nomissing.index.NRT2 == 1)]
    
    coeff.est.t <- solve(t(X.new.l)%*%(X.new.l))%*%t(X.new.l)%*%(Y.new.l)
    # coeff.est.t1 <- solve(t(NRT1.new.l)%*%(NRT1.new.l))%*%t(NRT1.new.l)%*%(Y.new.l.NRT1)
    # coeff.est.t2 <- solve(t(NRT2.new.l)%*%(NRT2.new.l))%*%t(NRT2.new.l)%*%(Y.new.l.NRT2)
    
    coeff.all <- rbind(coeff(l, NRT1, NRT2, x, y)$coeff.est, coeff.est.t)
    t.coeff=cbind(t.coeff, coeff.all)
  }
  
  bw_alpha1 <- locpol::thumbBw(t.seq[-1], t.coeff[1,], deg=1, kernel= locpol::gaussK)
  bw_alpha2 <- locpol::thumbBw(t.seq[-1], t.coeff[2,], deg=1, kernel= locpol::gaussK)
  bw_gamma1 <- locpol::thumbBw(t.seq[-1], t.coeff[3,], deg=1, kernel= locpol::gaussK)
  bw_gamma2 <- locpol::thumbBw(t.seq[-1], t.coeff[4,], deg=1, kernel= locpol::gaussK)
  bw_beta <- locpol::thumbBw(t.seq[-1], t.coeff[5,], deg=1, kernel= locpol::gaussK)
  bw_tao1 <- locpol::thumbBw(t.seq[-1], t.coeff[6,], deg=1, kernel= locpol::gaussK)
  bw_tao2 <- locpol::thumbBw(t.seq[-1], t.coeff[7,], deg=1, kernel= locpol::gaussK)
  
  
  hat.alpha.1=locpol::locPolSmootherC(t.seq[-1], t.coeff[1,], t.est-deltat, bw_alpha1, deg=1, kernel= locpol::gaussK)$beta0
  hat.alpha.2=locpol::locPolSmootherC(t.seq[-1], t.coeff[2,], t.est-deltat, bw_alpha2, deg=1, kernel= locpol::gaussK)$beta0
  hat.gamma.1=locpol::locPolSmootherC(t.seq[-1], t.coeff[3,], t.est, bw_gamma1, deg=1, kernel= locpol::gaussK)$beta0
  hat.gamma.2=locpol::locPolSmootherC(t.seq[-1], t.coeff[4,], t.est, bw_gamma2, deg=1, kernel= locpol::gaussK)$beta0
  hat.beta=locpol::locPolSmootherC(t.seq[-1], t.coeff[5,], t.est, bw_beta, deg=1, kernel= locpol::gaussK)$beta0
  hat.tao.1=locpol::locPolSmootherC(t.seq[-1], t.coeff[6,], t.est-deltat, bw_tao1, deg=1, kernel= locpol::gaussK)$beta0
  hat.tao.2=locpol::locPolSmootherC(t.seq[-1], t.coeff[7,], t.est-deltat, bw_tao2, deg=1, kernel= locpol::gaussK)$beta0
    
  list(hat.alpha1 = hat.alpha.1, hat.alpha2 = hat.alpha.2,
       hat.gamma1 = hat.gamma.1, hat.gamma2 = hat.gamma.2,
       hat.tao1 = hat.tao.1, hat.tao2 = hat.tao.2, hat.beta = hat.beta,
       hat.mediation1=hat.alpha.1*hat.beta, hat.mediation2=hat.alpha.2*hat.beta)
}