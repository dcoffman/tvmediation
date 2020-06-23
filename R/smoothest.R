#' Function to compute local polynomial estimation using rule of thumb for bandwidth selection
#' 
#' Part of the set of internal functions called within the \code{tvma} function to assist in the estimation of time varying mediation effect.
#' 
#' @param t.seq     a vector of time points at each observation
#' @param t.coeff   estimated coefficients
#' @param t.est     time points to make the estimation
#' @param deltat    a small constant which controls the time-lag of the effect of the mediator on the outcome.
#' 
#' @return \item{bw_alpha1}{a number computed via Fan and Gijbels(1996)'s Rule of thumb for bandwidth selector for alpha1 coefficient.}
#' @return \item{bw_beta2}{a number computed via Fan and Gijbels(1996)'s Rule of thumb for bandwidth selector for beta2 coefficient.}
#' @return \item{hat.alpha.1}{estimated treatment effect on mediator}
#' @return \item{hat.beta.2}{estimated mediator effect on outcome}
#' @return \item{est.M}{estimated medition effect, product of hat.alpha.1 and hat.beta.2}
#'
#' @export
#'

smoothest <- function(t.seq, t.coeff, t.est, deltat) {
  # Compute local polynomial estimation using rule of thumb for bandwidth selection
  #
  # Args:
  #   t.seq     -->   a vector of time points at each obs
  #   t.coeff   -->   estimated coefficients
  #   t.est     -->   time points to make the estimation
  #   deltat    -->   half the time between two time points
  #
  # Returns:
  #  smoothest  --> list containing ********
  #
  ##
  
  # Equations 4 & 5
  bw_alpha1 <- locpol::thumbBw(t.seq[-1], t.coeff[1, ], deg = 1, kernel = gaussK)
  bw_beta2 <- locpol::thumbBw(t.seq[-1], t.coeff[3, ], deg = 1, kernel = gaussK)
  hat.alpha.1 = locpol::locPolSmootherC(t.seq[-1], t.coeff[1, ], t.est - deltat,
                                bw_alpha1, deg = 1, kernel = gaussK)$beta0
  hat.beta.2 = locpol::locPolSmootherC(t.seq[-1], t.coeff[3, ], t.est, bw_beta2,
                               deg = 1, kernel = gaussK)$beta0
  
  return(list(bw_alpha1 = bw_alpha1, bw_beta2 = bw_beta2,
              hat.alpha.1 = hat.alpha.1, hat.beta.2 = hat.beta.2,
              est.M = hat.alpha.1*hat.beta.2))
}
