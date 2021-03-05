#' Function to estimate coefficients at time t
#' 
#' Part of the set of internal functions called within the \code{tvmcurve_3trt} function 
#' to assist in the estimation of time varying mediation effect.
#' 
#' @param j           a number indicating time point of observation
#' @param T1          a vector indicating assignment to treatment 1
#' @param T2          a vector indicating assignment to treatment 2
#' @param x           matrix of mediator values in wide format
#' @param y           matrix of outcome outcomes in wide format
#' 
#' @return \item{coeff.est}{estimated coefficients of the mediation model}
#' @return \item{nomissing.index}{index of complete cases}
#' 
#' 

coeff<-function(j, T1, T2, x, y)
{
  X.new.j <- cbind(c(T1, rep(0, length(T1))), c(T2, rep(0, length(T2))),
                 c(rep(0, length(T1)), T1), c(rep(0, length(T2)), T2),
                 c(rep(0, length(T1)), as.numeric(x[j-1,])))
  Y.new.j <- c(x[j,], y[j-1,])
  
  X.new.j <- scale(X.new.j, center=TRUE, scale=FALSE)
  Y.new.j <- scale(Y.new.j, center=TRUE, scale=FALSE)
  
  nomissing.x <- complete.cases(X.new.j)
  nomissing.y <- complete.cases(Y.new.j)
  
  nomissing.index <- nomissing.x * nomissing.y
  
  X.new.j <- X.new.j[which(nomissing.index==1),]
  Y.new.j <- Y.new.j[which(nomissing.index==1)]
  
  coeff.est <- solve(t(X.new.j)%*%(X.new.j))%*%t(X.new.j)%*%(Y.new.j)
  list(coeff.est=coeff.est, nomissing.index=nomissing.index)
}
