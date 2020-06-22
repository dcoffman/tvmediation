#' Function to compute new Mediator and Outcome using time T and T-1 mean centered on the individual
#' 
#' Part of the set of internal functions called within the \code{tvma} function to assist in the estimation of time varying mediation effect.
#' 
#' @param trt             numeric binary treatment schedule for each individual
#' @param mediator        (t.seq x N) matrix where N = number of observations. Column 1 is mediator at time T-1. Column 2 is mediator at time T.
#' @param outcome         (Nx1) matrix were N = number of observations. Column 1 is outcome at time T-1.
#' @param centerOutcome   outcomes must be centered estimating coefficients and covariance of error terms.
#' 
#' @return \item{newMO}{list containing new mediators, outcomes, and index of complete cases}
#' 
#' @export
#' 

newMediatorOutcome <- function(trt, M, Y) {
  # Compute new Mediator and Outcome using time T and T-1
  # mean centered on the individual
  #
  # Args:
  #   trt           -->   numeric binary treatment schedule for each individual.
  #   mediator      -->   (t.seq x N) matrix where N = number of obs.
  #                           Column 1 should be Mediator at time T-1
  #                           Column 2 should be Mediator at time T
  #   outcome       -->   (Nx1) matrix were N = number of obs.
  #                           Column 1 should be Outcome at time T-1
  #   centerOutcome -->   Outcomes must be centered estimating coefficients
  #                       & covariance of error terms.
  #
  # Returns:
  #   newMO         -->   list containing new Mediators, Outcomes, and Index of complete cases (nomissing)
  ##
  
  if (dim(M)[2] != 2) {
    stop("Argument Mediator must be Nx2 matrix: ",
         "\n\tMatrix is Nx", dim(M)[2], ".")
  }
  
  N <- length(trt)
  
  M.new <- cbind(c(trt, rep(0, N)),
                 c(rep(0, N), trt),
                 c(rep(0, N), M[,1]))
  Y.new <- c(M[,2], Y)
  
  M.new <- scale(M.new, center = TRUE, scale = FALSE)
  Y.new <- scale(Y.new, center = TRUE, scale = FALSE)
  
  nomissing.M <- complete.cases(M.new)
  nomissing.Y <- complete.cases(Y.new)
  nomissing.index <- nomissing.M * nomissing.Y
  M.new <- M.new[which(nomissing.index == 1), ]
  Y.new <- Y.new[which(nomissing.index == 1)]
  
  newMO <- list(M = M.new, Y = Y.new, nomissingIndex = nomissing.index)
  return(newMO)
  
}
