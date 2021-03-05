#' Function to compute local polynomial estimation using rule of thumb for bandwidth selection
#' 
#' Part of the set of internal functions called within the \code{tvma} function to assist 
#' in the estimation of the time varying mediation effect.
#' 
#' @param t.seq     a vector of time points at each observation
#' @param t.coeff   estimated coefficients
#' @param t.est     time points at which to make the estimation
#' @param deltat    a small constant which controls the time-lag of the effect of the 
#'                  mediator on the outcome, half the time between two time points
#' 
#' @return \item{bw_alpha}{a number computed via Fan and Gijbels' (1996) rule of thumb for 
#'                         bandwidth selector for alpha coefficient.}
#' @return \item{bw_gamma}{a number computed via Fan and Gijbels' (1996) rule of thumb for 
#'                         bandwidth selector for gamma coefficient.}
#' @return \item{bw_beta}{a number computed via Fan and Gijbels' (1996) rule of thumb for 
#'                        bandwidth selector for beta coefficient.}
#' @return \item{bw_beta}{a number computed via Fan and Gijbels' (1996) rule of thumb for 
#'                        bandwidth selector for tau coefficient.}
#' @return \item{hat.alpha}{estimated treatment effect on mediator}
#' @return \item{hat.gamma}{estimated treatment effect on outcome, adjusted for mediator}
#' @return \item{hat.beta}{estimated mediator effect on outcome}
#' @return \item{hat.tau}{estimated treatment effect on outcome, not adjusting for mediator}
#' @return \item{est.M}{estimated mediation effect, product of hat.alpha and hat.beta}
#'
#'

smoothest <- function(t.seq, t.coeff, t.est, deltat) {
  bw_alpha <- locpol::thumbBw(t.seq[-1], t.coeff[1, ], deg = 1, kernel = locpol::gaussK)
  bw_gamma <- locpol::thumbBw(t.seq[-1], t.coeff[2, ], deg = 1, kernel = locpol::gaussK)
  bw_beta  <- locpol::thumbBw(t.seq[-1], t.coeff[3, ], deg = 1, kernel = locpol::gaussK)
  bw_tau   <- locpol::thumbBw(t.seq[-1], t.coeff[4, ], deg = 1, kernel = locpol::gaussK)
  
  hat.alpha <- locpol::locPolSmootherC(t.seq[-1], t.coeff[1, ], t.est - deltat, bw_alpha,
                                        deg = 1, kernel = locpol::gaussK)$beta0
  hat.gamma <- locpol::locPolSmootherC(t.seq[-1], t.coeff[2, ], t.est, bw_gamma,
                                        deg = 1, kernel = locpol::gaussK)$beta0
  hat.beta <- locpol::locPolSmootherC(t.seq[-1], t.coeff[3, ], t.est, bw_beta, 
                                        deg = 1, kernel = locpol::gaussK)$beta0
  hat.tau <- locpol::locPolSmootherC(t.seq[-1], t.coeff[4, ], t.est, bw_tau, 
                                      deg = 1, kernel = locpol::gaussK)$beta0
  
  return(list(bw_alpha = bw_alpha, bw_gamma = bw_gamma, bw_beta = bw_beta, bw_tau = bw_tau,
              hat.alpha = hat.alpha, hat.gamma = hat.gamma, hat.beta = hat.beta, 
              hat.tau = hat.tau, est.M = hat.alpha*hat.beta))
}
